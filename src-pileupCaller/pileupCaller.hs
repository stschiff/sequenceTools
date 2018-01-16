{-# LANGUAGE OverloadedStrings #-}

import SeqTools.OrderedZip (orderedZip)
import SeqTools.VCF (liftParsingErrors)

import Control.Error (headErr)
import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad.Trans.Class (lift)
import Control.Monad.Trans.Reader (ReaderT, runReaderT, ask, asks)
import qualified Data.Attoparsec.Text as A
import Data.Char (isSpace, isDigit, toUpper)
import Data.List (sortBy, partition, sortOn, group, sort)
import Data.Maybe (catMaybes, fromMaybe)
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Data.Version (showVersion)
import Debug.Trace (trace)
import qualified Options.Applicative as OP
import Paths_sequenceTools (version)
import Pipes (Consumer, Pipe, yield, (>->), runEffect, Producer, Pipe, for, cat)
import Pipes.Attoparsec (parsed)
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, SafeT)
import Pipes.Safe.Prelude (withFile)
import Pipes.Text.Encoding (decodeUtf8)
import qualified Pipes.Text.IO as PT
import Prelude hiding (FilePath)
import System.IO (IOMode(..))
import System.Random (randomRIO, mkStdGen, setStdGen)
import Turtle hiding (tab, cat, stderr, err)

data ProgOpt = ProgOpt {
    optCallingModeString :: String,
    optSeed :: Maybe Int,
    optMinDepth :: Int,
    optMinSupport :: Int,
    optTransversionsOnly :: TransversionMode,
    optSnpFile :: Maybe FilePath,
    optOutChrom :: Maybe Text,
    optOutFormat :: OutFormat,
    optSampleNames :: Either [Text] FilePath,
    optSamplePopName :: Text,
    optEigenstratOutPrefix :: Maybe FilePath
}

data CallingMode =
        MajorityCalling Int |
        RandomCalling Int |
        RareCalling Int |
        RandomDiploidCalling
    deriving (Show, Read)
data TransversionMode =
    TransitionsMissing | SkipTransitions | AllSites deriving (Show, Read)
data OutFormat = EigenStrat | FreqSumFormat deriving (Show, Read)
data FreqSumRow = FreqSumRow T.Text Int Char Char [Int] deriving (Show)
data SnpEntry = SnpEntry T.Text Int Char Char deriving (Show)
                      -- Chrom  Pos Ref    Alt
data PileupRow = PileupRow T.Text Int Char [Text] deriving (Show)
type App = ReaderT ProgOpt (SafeT IO)
data Call = HaploidCall Char | DiploidCall Char Char | MissingCall

main :: IO ()
main = OP.execParser parser >>= runSafeT . runReaderT runWithOpts
  where
    parser = OP.info (pure (.) <*> versionInfoOpt <*> OP.helper <*> argParser)
                     (OP.progDesc "A program to perform genotype calling from \
                     \a pileup file")
    versionInfoOpt = OP.infoOption ("This is pileupCaller from sequenceTools \
        \Version " ++ showVersion version)
        (OP.long "version" <> OP.help "Print version and exit")

runWithOpts :: App ()
runWithOpts = do
    setRandomSeed
    let pileupProducer = parsed pileupParser PT.stdin >>= liftParsingErrors
    snpFile <- asks optSnpFile
    freqSumProducer <- case snpFile of
            Nothing -> simpleCalling pileupProducer
            Just fn -> snpListCalling fn pileupProducer
    outFormat <- asks optOutFormat
    case outFormat of
        FreqSumFormat -> printFreqSum freqSumProducer
        EigenStrat -> printEigenStrat freqSumProducer

setRandomSeed :: App ()
setRandomSeed = do
    seed <- asks optSeed
    case seed of
        Nothing -> return ()
        Just seed_ -> liftIO . setStdGen $ mkStdGen seed_

pileupParser :: A.Parser PileupRow
pileupParser = do
    chrom <- word
    _ <- A.space
    pos <- A.decimal
    _ <- A.space
    refA <- A.satisfy (A.inClass "ACTGNactgnM")
     -- for some reason, there is an M in the human reference at
     -- position 3:60830534 (both in hs37d5 and in hg19)
    _ <- A.space
    entries <- parsePileupPerSample chrom pos refA `A.sepBy1`
        A.satisfy A.isHorizontalSpace
    A.endOfLine
    let ret = PileupRow chrom pos refA entries
    --trace (show ret) $ return ret
    return ret
  where
    parsePileupPerSample chrom pos refA =
        processPileupEntry chrom pos refA <$> A.decimal <* A.space <*> word <*
            A.space <*> word

processPileupEntry :: Text -> Int -> Char -> Int -> Text -> Text -> Text
processPileupEntry chrom pos refA cov readBaseString _ =
    if cov == 0
    then ""
    else
        let returnString = T.pack $ go (T.unpack readBaseString)
        in  returnString
        -- in  if (T.length returnString /= cov && readBaseString /= "*")
        --     then
        --         trace ("Warning at " ++ show chrom ++ ", " ++
        --             show pos ++ ": readBaseString " ++
        --             show readBaseString ++ " does not match coverage number "
        --             ++ show cov) returnString
        --     else returnString
  where
    go (x:xs)
        | x `elem` (".," :: String) = refA : go xs
        | x `elem` ("ACTGNactgn" :: String) = toUpper x : go xs
        | x `elem` ("$*" :: String) = go xs
        | x == '^' = go (drop 1 xs)
        | x `elem` ("+-" :: String) =
            let [(num, rest)] = reads xs in go (drop num rest)
        | otherwise = error $ "cannot parse read base string: " ++ (x:xs)
    go [] = []

word :: A.Parser T.Text
word = A.takeTill isSpace

simpleCalling :: Producer PileupRow (SafeT IO) () ->
    App (Producer FreqSumRow (SafeT IO) ())
simpleCalling pileupProducer = do
    mode <- getCallingMode
    minDepth <- asks optMinDepth
    return $ for pileupProducer $ \pileupRow -> do
        let PileupRow chrom pos refA entryPerSample = pileupRow
        calls <- liftIO $ mapM (callGenotype mode minDepth refA)
            [T.unpack alleles | alleles <- entryPerSample]
        let altAlleles = findAlternativeAlleles refA calls
        when (length altAlleles == 1) $ do
            let altA = head altAlleles
                genotypes = map (callToGenotype refA altA) calls
            when (refA /= 'N' && any (>0) genotypes) $
                yield (FreqSumRow chrom pos refA altA genotypes)

getCallingMode :: App CallingMode
getCallingMode = do
    callingModeString <- asks optCallingModeString
    minSupport <- asks optMinSupport
    case callingModeString of
        "MajorityCalling" -> return $ MajorityCalling minSupport
        "RandomCalling" -> return $ RandomCalling minSupport
        "RareCalling" -> return $ RareCalling minSupport
        "RandomDiploidCalling" -> return RandomDiploidCalling
        _ -> error ("illegal calling mode " ++ callingModeString)

callGenotype :: CallingMode -> Int -> Char -> String -> IO Call
callGenotype mode minDepth refA alleles =
    if length alleles < minDepth then return MissingCall else
        case mode of
            MajorityCalling minSupport -> do
                let groupedAlleles = sortOn fst
                        [(length g, head g) | g <- group . sort $ alleles]
                    majorityCount = fst . head $ groupedAlleles
                    majorityAlleles =
                        [a | (n, a) <- groupedAlleles, n == majorityCount]
                a <- case majorityAlleles of
                    [a'] -> return a'
                    listA -> do
                        rn <- randomRIO (0, length listA - 1)
                        return (listA !! rn)
                let nrSupportingReads =
                        length . filter (==a) $ alleles
                if nrSupportingReads >= minSupport
                then
                    return $ HaploidCall a
                else
                    return MissingCall
            RandomCalling minSupport -> do
                res <- sampleWithoutReplacement alleles 1
                case res of
                    Nothing -> return MissingCall
                    Just [a] -> do
                        let nrSupportingReads =
                                length . filter (==a) $ alleles
                        if nrSupportingReads >= minSupport
                        then
                            return $ HaploidCall a
                        else
                            return MissingCall
            RareCalling minSupport -> do
                let groupedNonRefAlleles =
                        [(length g, head g) |
                         g <- group . sort . filter (/=refA) $ alleles,
                         length g >= minSupport]
                case groupedNonRefAlleles of
                    [] -> return $ DiploidCall refA refA
                    [(n, a)] -> return $ DiploidCall refA a
                    _ -> return MissingCall
            RandomDiploidCalling -> do
                res <- sampleWithoutReplacement alleles 2
                case res of
                    Nothing -> return MissingCall
                    Just [a1, a2] -> return $ DiploidCall a1 a2

sampleWithoutReplacement :: [a] -> Int -> IO (Maybe [a])
sampleWithoutReplacement = go []
  where
    go res xs 0 = return $ Just res
    go res xs n =
        if n > length xs
        then return Nothing
        else do
            rn <- randomRIO (0, length xs - 1)
            let a = xs !! rn
                xs' = let (ys, zs) = splitAt rn xs in ys ++ tail zs
            go (a:res) xs' (n - 1)

callToGenotype :: Char -> Char -> Call -> Int
callToGenotype refA altA call = case call of
    HaploidCall a | a == refA -> 0
                  | a == altA -> 1
                  | otherwise -> -1
    DiploidCall a1 a2 | (a1, a2) == (refA, refA) -> 0
                      | (a1, a2) == (refA, altA) -> 1
                      | (a1, a2) == (altA, refA) -> 1
                      | (a1, a2) == (altA, altA) -> 2
                      | otherwise                -> -1
    MissingCall -> -1

findAlternativeAlleles :: Char -> [Call] -> String
findAlternativeAlleles refA calls =
    let allAlleles = do
            c <- calls
            case c of
                HaploidCall a -> [a]
                DiploidCall a b -> [a, b]
                MissingCall -> []
        groupedNonRefAlleles = sortOn fst
            [(length g, head g) |
             g <- group . sort . filter (/=refA) $ allAlleles]
    in  if null groupedNonRefAlleles
        then []
        else
            let majorityCount = fst. head $ groupedNonRefAlleles
            in  [a | (n, a) <- groupedNonRefAlleles, n == majorityCount]

snpListCalling :: FilePath -> Producer PileupRow (SafeT IO) () ->
    App (Producer FreqSumRow (SafeT IO) ())
snpListCalling snpFileName pileupProducer = do
    sampleNameSpec <- asks optSampleNames
    sampleNames <- case sampleNameSpec of
        Left list -> return list
        Right fn -> T.lines <$> (liftIO . T.readFile . T.unpack . format fp) fn
    let snpTextProd = PT.readFile ((T.unpack . format fp) snpFileName)
        snpProd = parsed snpParser snpTextProd >>= liftParsingErrors
        jointProd = orderedZip cmp snpProd pileupProducer
    mode <- getCallingMode
    minDepth <- asks optMinDepth
    let ret = for jointProd $ \jointEntry ->
            case jointEntry of
                (Just (SnpEntry snpChrom snpPos snpRef snpAlt), Nothing) ->
                    yield $ FreqSumRow snpChrom snpPos snpRef snpAlt
                            (replicate (length sampleNames) (-1))
                (Just (SnpEntry snpChrom snpPos snpRef snpAlt),
                 Just pileupRow) -> do
                    let PileupRow chrom pos refA entryPerSample = pileupRow
                    calls <- liftIO $ mapM (callGenotype mode minDepth refA)
                        [T.unpack alleles | alleles <- entryPerSample]
                    let genotypes = map (callToGenotype snpRef snpAlt) calls
                    yield (FreqSumRow snpChrom snpPos snpRef snpAlt genotypes)
                _ -> return ()
    return (fst <$> ret)
  where
    cmp (SnpEntry snpChrom snpPos _ _) (PileupRow pChrom pPos _ _) =
        case snpChrom `chromNameCompare` pChrom of
            LT -> LT
            GT -> GT
            EQ -> snpPos `compare` pPos
    chromNameCompare c1 c2 =
        let (c1Nums, c1NonNums) = partition isDigit . T.unpack $ c1
            (c2Nums, c2NonNums) = partition isDigit . T.unpack $ c2
            c1Num = read c1Nums :: Int
            c2Num = read c2Nums :: Int
        in  case c1NonNums `compare` c2NonNums of
                LT -> LT
                GT -> GT
                EQ -> c1Num `compare` c2Num

snpParser :: A.Parser SnpEntry
snpParser = do
    _ <- many A.space
    _ <- word
    _ <- A.many1 A.space
    chrom <- word
    _ <- A.many1 A.space
    _ <- A.double
    _ <- A.many1 A.space
    pos <- A.decimal
    _ <- A.many1 A.space
    ref <- A.satisfy (A.inClass "ACTG")
    _ <- A.many1 A.space
    alt <- A.satisfy (A.inClass "ACTG")
    _ <- A.satisfy (\c -> c == '\r' || c == '\n')
    let ret = SnpEntry chrom pos ref alt
    return ret

printFreqSum :: Producer FreqSumRow (SafeT IO) () -> App ()
printFreqSum freqSumProducer = do
    outChrom <- asks optOutChrom
    transversionsOnly <- asks optTransversionsOnly
    sampleNameSpec <- asks optSampleNames
    callingMode <- getCallingMode
    sampleNames <- case sampleNameSpec of
        Left list -> return list
        Right fn -> T.lines <$> (liftIO . T.readFile . T.unpack . format fp) fn
    let nrHaplotypes = case callingMode of
            MajorityCalling _ -> 1
            RandomCalling _ -> 1
            RareCalling _ -> 2
            RandomDiploidCalling -> 2
    echo . unsafeTextToLine . format ("#CHROM\tPOS\tREF\tALT\t"%s) $
        T.intercalate "\t"
            [format (s%"("%d%")") n nrHaplotypes | n <- sampleNames]
    lift . runEffect $ freqSumProducer >->
        filterTransitions transversionsOnly >->
        P.map (showFreqSum outChrom) >-> printToStdOut
  where
    showFreqSum outChrom (FreqSumRow chrom pos ref alt calls) =
        let newChrom = fromMaybe chrom outChrom
        in  format (s%"\t"%d%"\t"%s%"\t"%s%"\t"%s) newChrom pos
                (T.singleton ref) (T.singleton alt) callsStr
      where
        callsStr = (T.intercalate "\t" . map (format d)) calls

printEigenStrat :: Producer FreqSumRow (SafeT IO) () -> App ()
printEigenStrat freqSumProducer = do
    outChrom <- asks optOutChrom
    transversionsOnly <- asks optTransversionsOnly
    sampleNameSpec <- asks optSampleNames
    callingMode <- getCallingMode
    let diploidizeCall = case callingMode of
            RandomCalling _ -> True
            MajorityCalling _ -> True
            RandomDiploidCalling -> False
            RareCalling _ -> False
    sampleNames <- case sampleNameSpec of
        Left list -> return list
        Right fn -> T.lines <$> (liftIO . T.readFile . T.unpack . format fp) fn
    eigenStratOutPrefix <- asks optEigenstratOutPrefix
    case eigenStratOutPrefix of
        Nothing -> liftIO . throwIO $
            AssertionFailed "need an eigenstratPrefix for EigenStratFormat"
        Just fn -> do
            let snpOut = fn <.> "snp.txt"
                indOut = fn <.> "ind.txt"
            p <- asks optSamplePopName
            lift . withFile (T.unpack . format fp $ indOut) WriteMode $
                \indOutHandle ->
                    sequence_ [liftIO $ T.hPutStrLn indOutHandle
                               (format (s%"\tU\t"%s) n p) | n <- sampleNames]
            lift . withFile (T.unpack . format fp $ snpOut) WriteMode $
                \snpOutHandle -> runEffect $
                    freqSumProducer >-> filterTransitions transversionsOnly >->
                        printEigenStratRow diploidizeCall outChrom snpOutHandle >-> printToStdOut

printToStdOut :: Consumer T.Text (SafeT IO) ()
printToStdOut = for cat (liftIO . T.putStrLn)

filterTransitions :: TransversionMode ->
    Pipe FreqSumRow FreqSumRow (SafeT IO) ()
filterTransitions transversionsMode =
    case transversionsMode of
        SkipTransitions ->
            P.filter (\(FreqSumRow _ _ ref alt _) -> isTransversion ref alt)
        TransitionsMissing ->
            P.map (\(FreqSumRow chrom pos ref alt calls) -> FreqSumRow chrom pos ref alt [-1 | c <- calls])
        AllSites -> cat
  where
    isTransversion ref alt = not $ isTransition ref alt
    isTransition ref alt = ((ref == 'A') && (alt == 'G')) ||
        ((ref == 'G') && (alt == 'A')) || ((ref == 'C') && (alt == 'T')) ||
        ((ref == 'T') && (alt == 'C'))

printEigenStratRow :: Bool -> Maybe Text -> Handle ->
    Pipe FreqSumRow Text (SafeT IO) r
printEigenStratRow diploidizeCall outChrom snpOutHandle =
    for cat $ \(FreqSumRow chrom pos ref alt calls) -> do
        let newChrom = fromMaybe chrom outChrom
        let n = format (s%"_"%d) newChrom pos
            snpLine = format (s%"\t"%s%"\t0\t"%d%"\t"%s%"\t"%s) n newChrom pos
                (T.singleton ref) (T.singleton alt)
        liftIO . T.hPutStrLn snpOutHandle $ snpLine
        yield . T.concat . map (format d . toEigenStratNum) $ calls
  where
    toEigenStratNum c =
        if diploidizeCall then
            case c of
                0 -> 2 :: Int
                1 -> 0
                -1 -> 9
                _ -> error "illegal call for pseudo-haploid Calling method"
        else
            case c of
                0 -> 2 :: Int
                1 -> 1
                2 -> 0
                -1 -> 9
                _ -> error ("unknown genotype " ++ show c)

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseCallingMode <*> parseSeed <*> parseMinDepth <*>
    parseMinSupport <*> parseTransversionsMode <*> parseSnpFile <*>
    parseOutChrom <*>
    parseFormat <*> parseSampleNames <*> parseSamplePopName <*>
    parseEigenstratOutPrefix
  where
    parseCallingMode = OP.strOption (OP.long "mode" <> OP.short 'm' <>
        OP.value "RandomCalling" <> OP.showDefault <> OP.metavar "<MODE>" <>
        OP.help "specify the mode of calling. \
        \MajorityCalling: Pick the allele supported by the \
        \most reads. If an equal numbers of Alleles fulfil this, pick one at \
        \random. This results in a haploid call.\
        \RandomCalling: Pick one read at random. This results in a haploid \
        \call;\
        \RareCalling (deprecated!): If at least n reads support the \
        \non-reference allele, call a heterozygote,  \
        \otherwise call homozygous reference, where n is set by --minSupport;\
        \RandomDiploid: Sample two random reads at random and represent the \
        \individual by the diploid genotype constructed from those two random \
        \picks. This will always assign missing data to positions where only \
        \one read is present, even if minDepth=1.\n")
    parseSeed = OP.option (Just <$> OP.auto) (OP.long "seed" <>
        OP.value Nothing <> OP.metavar "<RANDOM_SEED>" <>
        OP.help "random seed used for random calling. If not given, use \
        \system clock to seed the random number generator")
    parseMinDepth = OP.option OP.auto (OP.long "minDepth" <> OP.short 'd' <>
        OP.value 1 <> OP.showDefault <> OP.metavar "<DEPTH>" <>
        OP.help "specify the minimum depth for a call. For sites with fewer \
            \reads than this number, declare Missing")
    parseMinSupport = OP.option OP.auto (OP.long "minSupport" <>
        OP.value 1 <> OP.showDefault <> OP.metavar "<MIN_SUPPORT>" <>
        OP.help "specify the minimum number of supporting reads for the \
        \RandomCalling, MajorityCalling and RareCalling (deprecated) methods. For calling rare variants with either of those methods, you should set \
        \--minSupport 2 or higher.")
    parseTransversionsMode = OP.option OP.auto (OP.long "transversionsMode" <>
        OP.short 't' <> OP.value AllSites <> OP.metavar "MODE" <>
        OP.showDefault <> OP.help "Three \
        \options possible: SkipTransitions: skip transitions in the output; \
        \TransitionsMissing: output transition sites as missing data in all \
        \samples; AllSites: output all sites including transitions.")
    parseSnpFile = OP.option (Just . fromText . T.pack <$> OP.str)
        (OP.long "snpFile" <> OP.short 'f' <> OP.value Nothing <>
        OP.metavar "<FILE>" <> OP.help "specify an Eigenstrat SNP file for \
        \the positions and alleles to call. All \
        \positions in the SNP file will be output, adding missing data where \
        \necessary. Note that simpleBamCaller automatically checks whether \
        \alleles in the SNP file are flipped with respect to the human \
        \reference. But it assumes that the strand-orientation is the same.")
    parseOutChrom = OP.option (Just . T.pack <$> OP.str) (OP.long "outChrom" <>
        OP.metavar "<CHROM>" <> OP.help "specify the output chromosome name. \
        \This can be useful if the input chromosome name is something like \
        \'chr1' and you would like to merge with a dataset that has just \
        \'1'." <> OP.value Nothing)
    parseFormat = OP.option OP.auto (OP.metavar "<OUT_FORMAT>" <>
        OP.long "format" <> OP.short 'o' <> OP.value FreqSumFormat <>
        OP.showDefault <>
        OP.help "specify output format: EigenStrat or FreqSum")
    parseSampleNames = parseSampleNameList <|> parseSampleNameFile
    parseSampleNameList = OP.option (Left . T.splitOn "," . T.pack <$> OP.str)
        (OP.long "sampleNames" <> OP.metavar "NAME1,NAME2,..." <>
        OP.help "give the names of the samples as comma-separated list")
    parseSampleNameFile = OP.option (Right . fromText . T.pack <$> OP.str)
        (OP.long "sampleNameFile" <> OP.metavar "<FILE>" <>
        OP.help "give the names of the samples in a file with one name per \
        \line")
    parseSamplePopName = OP.option (T.pack <$> OP.str)
        (OP.long "samplePopName" <> OP.value "Unknown" <> OP.showDefault <>
        OP.metavar "POP" <>
        OP.help "specify the population name of the samples, which is included\
        \ in the output *.ind.txt file. This will be ignored if the output \
        \format is not Eigenstrat")
    parseEigenstratOutPrefix = OP.option (Just . fromText . T.pack <$> OP.str)
        (OP.long "eigenstratOutPrefix" <> OP.short 'e' <>
        OP.value Nothing <> OP.metavar "<FILE_PREFIX>" <>
        OP.help "specify the filenames for the EigenStrat SNP and IND \
        \file outputs: <FILE_PREFIX>.snp.txt and <FILE_PREFIX>.ind.txt \
        \Ignored if Output format is not Eigenstrat")

{-# LANGUAGE OverloadedStrings #-}

import SeqTools.OrderedZip (orderedZip)
import SequenceFormats.Utils (liftParsingErrors)
import SequenceFormats.Eigenstrat (readEigenstratSnpFile, EigenstratSnpEntry(..), GenoLine, 
    EigenstratIndEntry(..), Sex(..), writeEigenstrat, GenoEntry(..))
import SequenceFormats.FreqSum(FreqSumEntry(..), printFreqSumStdOut, FreqSumHeader(..))

import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad.Trans.Class (lift)
import Control.Monad.Trans.Reader (ReaderT, runReaderT, asks)
import qualified Data.Attoparsec.Text as A
import Data.Char (isSpace, isDigit, toUpper)
import Data.List (partition, sortOn, group, sort)
import Data.Maybe (fromMaybe, catMaybes)
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Data.Version (showVersion)
import Data.Vector (fromList)
-- import Debug.Trace (trace)
import qualified Options.Applicative as OP
import Paths_sequenceTools (version)
import Pipes (Pipe, yield, (>->), runEffect, Producer, Pipe, for, cat)
import Pipes.Attoparsec (parsed)
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, SafeT)
import qualified Pipes.Text.IO as PT
import Prelude hiding (FilePath)
import System.Random (randomRIO, mkStdGen, setStdGen)
import Turtle hiding (tab, cat, stderr, err, sort, sortOn, x, g)

data ProgOpt = ProgOpt {
    optCallingModeString :: String,
    optSeed :: Maybe Int,
    optMinDepth :: Int,
    optMinSupport :: Int,
    optDownSampling :: Bool,
    optTransitionsOnly :: TransitionsMode,
    optSnpFile :: Maybe FilePath,
    optOutChrom :: Maybe Text,
    optOutFormat :: OutFormat,
    optSampleNames :: Either [Text] FilePath,
    optSamplePopName :: Text,
    optEigenstratOutPrefix :: Maybe FilePath
}

data CallingMode = MajorityCalling Bool | RandomCalling | RareCalling Int |
        RandomDiploidCalling
    deriving (Show, Read)
data TransitionsMode =
    TransitionsMissing | SkipTransitions | AllSites deriving (Show, Read)
data OutFormat = EigenStrat | FreqSumFormat deriving (Show, Read)
data PileupRow = PileupRow T.Text Int Char [Text] deriving (Show)
type App = ReaderT ProgOpt (SafeT IO)
data Call = HaploidCall Char | DiploidCall Char Char | MissingCall

progInfo :: String
progInfo = "A program to perform genotype calling from a pileup file. Part of sequenceTools \
    \version " ++ showVersion version

main :: IO ()
main = OP.execParser parser >>= runSafeT . runReaderT runWithOpts
  where
    parser = OP.info (pure (.) <*> versionInfoOpt <*> OP.helper <*> argParser)
                     (OP.progDesc progInfo)
    versionInfoOpt = OP.infoOption (showVersion version)
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
        FreqSumFormat -> outputFreqSum freqSumProducer
        EigenStrat -> outputEigenStrat freqSumProducer

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
    entries <- parsePileupPerSample refA `A.sepBy1`
        A.satisfy A.isHorizontalSpace
    A.endOfLine
    let ret = PileupRow chrom pos refA entries
    --trace (show ret) $ return ret
    return ret
  where
    parsePileupPerSample refA =
        processPileupEntry refA <$> A.decimal <* A.space <*> word <*
            A.space <*> word

processPileupEntry :: Char -> Int -> Text -> Text -> Text
processPileupEntry refA cov readBaseString _ =
    if cov == 0 then "" else T.pack $ go (T.unpack readBaseString)
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
    App (Producer FreqSumEntry (SafeT IO) ())
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
            when (refA /= 'N' && (sum . catMaybes) genotypes > 0) $
                yield (FreqSumEntry chrom pos refA altA genotypes)

getCallingMode :: App CallingMode
getCallingMode = do
    callingModeString <- asks optCallingModeString
    downSampling <- asks optDownSampling
    minSupport <- asks optMinSupport
    case callingModeString of
        "MajorityCalling" -> return $ MajorityCalling downSampling
        "RandomCalling" -> return RandomCalling
        "RareCalling" -> return $ RareCalling minSupport
        "RandomDiploidCalling" -> return RandomDiploidCalling
        _ -> error ("illegal calling mode " ++ callingModeString)

callGenotype :: CallingMode -> Int -> Char -> String -> IO Call
callGenotype mode minDepth refA alleles =
    if length alleles < minDepth then return MissingCall else
        case mode of
            MajorityCalling downSampling -> do
                alleles' <-
                    if downSampling
                    then do
                        Just a <- sampleWithoutReplacement alleles minDepth
                        return a
                    else return alleles
                let groupedAlleles = sortOn fst
                        [(length g, head g) | g <- group . sort $ alleles']
                    majorityCount = fst . last $ groupedAlleles
                    majorityAlleles =
                        [a | (n, a) <- groupedAlleles, n == majorityCount]
                a <- case majorityAlleles of
                    [a'] -> return a'
                    listA -> do
                        rn <- randomRIO (0, length listA - 1)
                        return (listA !! rn)
                return $ HaploidCall a
            RandomCalling -> do
                res <- sampleWithoutReplacement alleles 1
                case res of
                    Nothing -> return MissingCall
                    Just [a] -> return $ HaploidCall a
                    _ -> error "should not happen"
            RareCalling minSupport -> do
                let groupedNonRefAlleles =
                        [(length g, head g) |
                         g <- group . sort . filter (/=refA) $ alleles,
                         length g >= minSupport]
                case groupedNonRefAlleles of
                    [] -> return $ DiploidCall refA refA
                    [(_, a)] -> return $ DiploidCall refA a
                    _ -> return MissingCall
            RandomDiploidCalling -> do
                res <- sampleWithoutReplacement alleles 2
                case res of
                    Nothing -> return MissingCall
                    Just [a1, a2] -> return $ DiploidCall a1 a2
                    _ -> error "should not happen"

sampleWithoutReplacement :: [a] -> Int -> IO (Maybe [a])
sampleWithoutReplacement = go []
  where
    go res _ 0 = return $ Just res
    go res xs n
        | n > length xs = return Nothing
        | n == length xs = return $ Just (xs ++ res)
        | otherwise = do
                rn <- randomRIO (0, length xs - 1)
                let a = xs !! rn
                    xs' = remove rn xs
                go (a:res) xs' (n - 1)
    remove i xs = let (ys, zs) = splitAt i xs in ys ++ tail zs

callToGenotype :: Char -> Char -> Call -> Maybe Int
callToGenotype refA altA call = case call of
    HaploidCall a | a == refA -> Just 0
                  | a == altA -> Just 1
                  | otherwise -> Nothing
    DiploidCall a1 a2 | (a1, a2) == (refA, refA) -> Just 0
                      | (a1, a2) == (refA, altA) -> Just 1
                      | (a1, a2) == (altA, refA) -> Just 1
                      | (a1, a2) == (altA, altA) -> Just 2
                      | otherwise                -> Nothing
    MissingCall -> Nothing

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
    App (Producer FreqSumEntry (SafeT IO) ())
snpListCalling snpFileName pileupProducer = do
    sampleNameSpec <- asks optSampleNames
    sampleNames <- case sampleNameSpec of
        Left list -> return list
        Right fn -> T.lines <$> (liftIO . T.readFile . T.unpack . format fp) fn
    let snpProd = readEigenstratSnpFile (T.unpack . format fp $ snpFileName)
        jointProd = orderedZip cmp snpProd pileupProducer
    mode <- getCallingMode
    minDepth <- asks optMinDepth
    let ret = for jointProd $ \jointEntry ->
            case jointEntry of
                (Just (EigenstratSnpEntry snpChrom snpPos snpRef snpAlt), Nothing) ->
                    yield $ FreqSumEntry snpChrom snpPos snpRef snpAlt
                            (replicate (length sampleNames) Nothing)
                (Just (EigenstratSnpEntry snpChrom snpPos snpRef snpAlt),
                 Just pileupRow) -> do
                    let PileupRow _ _ refA entryPerSample = pileupRow
                    calls <- liftIO $ mapM (callGenotype mode minDepth refA)
                        [T.unpack alleles | alleles <- entryPerSample]
                    let genotypes = map (callToGenotype snpRef snpAlt) calls
                    yield (FreqSumEntry snpChrom snpPos snpRef snpAlt genotypes)
                _ -> return ()
    return (fst <$> ret)
  where
    cmp (EigenstratSnpEntry snpChrom snpPos _ _) (PileupRow pChrom pPos _ _) =
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

outputFreqSum :: Producer FreqSumEntry (SafeT IO) () -> App ()
outputFreqSum freqSumProducer = do
    outChrom <- asks optOutChrom
    transitionsOnly <- asks optTransitionsOnly
    sampleNameSpec <- asks optSampleNames
    callingMode <- getCallingMode
    sampleNames <- case sampleNameSpec of
        Left list -> return list
        Right fn -> T.lines <$> (liftIO . T.readFile . T.unpack . format fp) fn
    let nrHaplotypes = case callingMode of
            MajorityCalling _ -> 1 :: Int
            RandomCalling -> 1
            RareCalling _ -> 2
            RandomDiploidCalling -> 2
    let header' = FreqSumHeader sampleNames [nrHaplotypes | _ <- sampleNames]
        outProd = freqSumProducer >-> filterTransitions transitionsOnly >->
            P.map (correctChrom outChrom)
    lift . runEffect $ outProd >-> printFreqSumStdOut header'
  where
    correctChrom outChrom =
        case outChrom of
            Just c -> \fse -> fse {fsChrom = c}
            Nothing -> id

outputEigenStrat :: Producer FreqSumEntry (SafeT IO) () -> App ()
outputEigenStrat freqSumProducer = do
    outChrom <- asks optOutChrom
    transitionsMode <- asks optTransitionsOnly
    sampleNameSpec <- asks optSampleNames
    callingMode <- getCallingMode
    let diploidizeCall = case callingMode of
            RandomCalling -> True
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
            let snpOut = T.unpack . format fp $ fn <.> "snp.txt"
                indOut = T.unpack . format fp $ fn <.> "ind.txt"
                genoOut = T.unpack . format fp $ fn <.> "geno.txt"
            p <- asks optSamplePopName
            let indEntries = [EigenstratIndEntry n Unknown p | n <- sampleNames]
            lift . runEffect $ freqSumProducer >-> filterTransitions transitionsMode >->
                        toEigenstrat diploidizeCall outChrom >->
                        writeEigenstrat genoOut snpOut indOut indEntries

filterTransitions :: TransitionsMode ->
    Pipe FreqSumEntry FreqSumEntry (SafeT IO) ()
filterTransitions transversionsMode =
    case transversionsMode of
        SkipTransitions ->
            P.filter (\(FreqSumEntry _ _ ref alt _) -> isTransversion ref alt)
        TransitionsMissing ->
            P.map (\(FreqSumEntry chrom pos ref alt calls) ->
                let calls' = if isTransversion ref alt then calls else [Nothing | _ <- calls]
                in  FreqSumEntry chrom pos ref alt calls')
        AllSites -> cat
  where
    isTransversion ref alt = not $ isTransition ref alt
    isTransition ref alt = ((ref == 'A') && (alt == 'G')) ||
        ((ref == 'G') && (alt == 'A')) || ((ref == 'C') && (alt == 'T')) ||
        ((ref == 'T') && (alt == 'C'))

toEigenstrat :: Bool -> Maybe Text -> Pipe FreqSumEntry (EigenstratSnpEntry, GenoLine) (SafeT IO) ()
toEigenstrat diploidizeCall outChrom = P.map toEigenstrat'
  where
    toEigenstrat' (FreqSumEntry chrom pos ref alt calls) =
        let newChrom = fromMaybe chrom outChrom
            snpEntry = EigenstratSnpEntry newChrom pos ref alt
            geno = fromList . map toGenoCall $ calls
        in  (snpEntry, geno)
    toGenoCall c =
        if diploidizeCall then
            case c of
                Just 0 -> HomRef
                Just 1 -> HomAlt
                Nothing -> Missing
                _ -> error "illegal call for pseudo-haploid Calling method"
        else
            case c of
                Just 0 -> HomRef
                Just 1 -> Het
                Just 2 -> HomAlt
                Nothing -> Missing
                _ -> error ("unknown genotype " ++ show c)

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseCallingMode <*> parseSeed <*> parseMinDepth <*>
    parseMinSupport <*> parseDownSampling <*> parseTransitionsMode <*> 
    parseSnpFile <*> parseOutChrom <*> parseFormat <*> parseSampleNames <*> 
    parseSamplePopName <*> parseEigenstratOutPrefix
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
        \RareCalling (deprecated) method. This option is ignored for other \
        \calling methods. For RareCalling, you should use \
        \--minSupport 2 or higher.")
    parseDownSampling = OP.switch (OP.long "withDownSampling" <> OP.help
        "When this switch is given, the MajorityCalling mode with downsample \
        \from the total number of reads a number of reads \
        \(without replacement) equal to the --minDepth given. This mitigates a \
        \subtle reference bias in the MajorityCalling model.")
    parseTransitionsMode = OP.option OP.auto (OP.long "transitionsMode" <>
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
        \necessary. Note that pileupCaller automatically checks whether \
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

{-# LANGUAGE OverloadedStrings #-}

import OrderedZip (orderedZip)

import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad (forM, void)
import Control.Monad.Random (evalRandIO)
import Control.Monad.Trans.Class (lift)
import qualified Data.Attoparsec.Text as A
import Data.Char (isSpace)
import Data.List (sortBy)
import Data.Maybe (catMaybes)
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Debug.Trace (trace)
import qualified Options.Applicative as OP
import Pipes (Pipe, await, yield, (>->), runEffect, Producer, Pipe, for, cat, next)
import Pipes.Attoparsec (parsed)
import Pipes.Cliff (CreateProcess(..), CmdSpec(..), pipeOutput, NonPipe(..))
import Pipes.Cliff.Core (defaultHandler)
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, SafeT, MonadSafe)
import Pipes.Text.Encoding (decodeUtf8)
import qualified Pipes.Text.IO as PT
import Prelude hiding (FilePath)
import System.Random (randomRIO, mkStdGen, setStdGen)
import System.Random.Shuffle (shuffleM)
import Turtle hiding (tab, cat)

data ProgOpt = ProgOpt {
    optCallingMode :: CallingMode,
    optSeed :: Maybe Int,
    optMinDepth :: Int,
    optFilter :: FilterMode,
    optSnpFile :: Maybe FilePath,
    optRegion :: Text,
    optReference :: FilePath,
    optOutFormat :: OutFormat,
    optBamFiles :: [FilePath]
}

data CallingMode = MajorityCalling | RandomCalling deriving (Show, Read)
data FilterMode = NoFilter | Transversions | TransversionsMissing deriving (Show, Read, Eq)
data OutFormat = EigenStrat | FreqSumFormat deriving (Show, Read)
data FreqSumRow = FreqSumRow Text Int Char Char [Int] deriving (Show)
data VCFentry = VCFentry Text Int [Char] [[Int]] deriving (Show) -- Chrom Pos Alleles Number_of_reads_per_individual
data SnpEntry = SnpEntry Text Int Char Char deriving (Show)-- Chrom Pos Ref Alt

main = OP.execParser parser >>= runWithOpts
  where
    parser = OP.info (OP.helper <*> argParser) (OP.progDesc "A program to perform simple genotype calling directly \ 
                                                             \from BAM")

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseCallingMode <*> parseSeed <*> parseMinDepth <*> parseFilter <*> parseSnpFile <*>
                        parseChrom <*> parseRef <*> parseFormat <*> OP.some parseBams
  where
    parseCallingMode = OP.option OP.auto (OP.long "mode" <> OP.short 'm' <> OP.value RandomCalling <> OP.showDefault <>
                                          OP.metavar "<MODE>" <>
                                          OP.help "specify the mode of calling: MajorityCalling or RandomCalling")
    parseSeed = OP.option (Just <$> OP.auto) (OP.long "seed" <> OP.value Nothing <> OP.metavar "<RANDOM_SEED>" <>
                                              OP.help "random seed used for random calling. If not given, use system \  
                                                       \clock to seed the random number generator")
    parseMinDepth = OP.option OP.auto (OP.long "minDepth" <> OP.short 'd' <> OP.value 1 <> OP.showDefault <>
                                       OP.metavar "<DEPTH>" <>
                                       OP.help "specify the minimum depth for a call")
    parseFilter = OP.option OP.auto (OP.long "tfilter" <> OP.short 't' <> OP.value NoFilter <> OP.showDefault <>
                            OP.help "filter transversions. Three options are available: NoFilter (call all sites), \ 
                                     \Tranversions (remove transition SNPs from the output), TransversionsMissing \
                                     \(like Tranversions, but only mark Transitions as missing, do not remove them). \
                                     \The second and third option are equivalent if no SNP file is provided. ")
    parseSnpFile = OP.option (Just . fromText . T.pack <$> OP.str)
                   (OP.long "snpFile" <> OP.short 'f' <> OP.value Nothing <> OP.metavar "<FILE>" <>
                    OP.help "specify a SNP file for the positions and alleles to call. All \
                             \positions in the SNP file will be output, adding missing data where necessary. \
                             \The file should have three columns (space- or tab-separated): Chromosome, position and \
                             \alleles, where the alleles should be one reference and one alternative allele \
                             \separated by a comma.")
    parseChrom = OP.option (T.pack <$> OP.str) (OP.long "chrom" <> OP.short 'c' <> OP.metavar "<CHROM>" <>
                                    OP.help "specify the region in the BAM file to call from. Can be just the \
                                             \chromosome, or a string of the form CHROM:START-END")
    parseRef = OP.option (fromText . T.pack <$> OP.str) (OP.long "reference" <> OP.short 'r' <> OP.metavar "<REF>" <>
                                  OP.help "the reference fasta file")
    parseBams = OP.argument (fromText . T.pack <$> OP.str) (OP.metavar "<BAM_FILE>" <> OP.help "input file")
    parseFormat = OP.option OP.auto (OP.metavar "<OUT_FORMAT>" <> OP.long "format" <> OP.short 'o' <>
                                     OP.value FreqSumFormat <> OP.showDefault <>
                                     OP.help "specify output format: EigenStrat or FreqSum")
    
runWithOpts :: ProgOpt -> IO ()
runWithOpts (ProgOpt mode seed minDepth filter_ snpFile region reference outFormat bamFiles) = do
    case seed of
        Nothing -> return ()
        Just seed_ -> setStdGen $ mkStdGen seed_
    freqSumProducer <- case snpFile of
        Nothing -> do
            let cmd = format ("samtools mpileup -q30 -Q30 -C50 -I -f "%fp%" -g -t DPR -r "%s%" "%s%" | bcftools \
                              \view -v snps -H") reference region bams
            vcfTextProd <- produceFromCommand cmd
            let vcfProd = parsed vcfParser vcfTextProd
            return $ vcfProd >-> processVcfSimple mode minDepth filter_
        Just fn -> do
            let cmd = format ("samtools mpileup -q30 -Q30 -C50 -I -f "%fp%" -g -t DPR -r "%s%" -l "%fp%" "%s%
                           " | bcftools view -H") reference region fn bams
            vcfTextProd <- produceFromCommand cmd
            let snpTextProd = PT.readFile ((T.unpack . format fp) fn)
            let vcfProd = parsed vcfParser vcfTextProd
            let chrom = head (T.splitOn ":" region)
            let snpProd = parsed snpParser snpTextProd >-> P.filter (\(SnpEntry c _ _ _) -> c == chrom)
            let jointProd = orderedZip cmp snpProd vcfProd
            return . fmap snd $ jointProd >-> processVcfWithSnpFile (length bamFiles) mode minDepth filter_
    let freqSumToText = case outFormat of
            FreqSumFormat -> printFreqSum
            EigenStrat -> printEigenStrat
    res <- runSafeT . runEffect $ freqSumProducer >-> P.map freqSumToText >-> printToStdOut
    case res of
        Left (e, rest) -> do
            err . format w $ e
            Right (c, _) <- runSafeT $ next rest
            err c
        Right () -> return ()
  where
    bams = (T.intercalate " " (map (format fp) bamFiles))
    cmp (SnpEntry _ snpPos _ _) (VCFentry _ vcfPos _ _) = snpPos `compare` vcfPos
    printToStdOut = for cat (lift . lift . T.putStrLn)

produceFromCommand :: Text -> IO (Producer Text (SafeT IO) ())
produceFromCommand cmd = do
    let createProcess = CreateProcess (ShellCommand (T.unpack cmd)) Nothing Nothing False False False defaultHandler
    (p, _) <- pipeOutput Inherit Inherit createProcess
    return . void . decodeUtf8 $ p

vcfParser :: A.Parser VCFentry
vcfParser = do
    chrom <- word
    tab
    pos <- A.decimal
    tab >> word >> tab
    ref <- A.satisfy (A.inClass "NACTG")
    tab
    alt <- (altAllele <|> xAllele) `A.sepBy1` (A.char ',')
    tab
    _ <- A.count 3 (word >> tab)
    _ <- A.string "PL:DPR" >> tab
    coverages <- coverage `A.sepBy1` tab
    _ <- A.satisfy (\c -> c == '\r' || c == '\n')
    -- trace (show (chrom, pos, ref, alt, coverages)) $ return ()
    let filteredAlt = catMaybes alt
    let filteredCoverages = [[c | (c, a) <- zip cov (Just ref:alt), a /= Nothing] | cov <- coverages]
    return $ VCFentry chrom pos (ref:filteredAlt) filteredCoverages
  where
    altAllele = Just <$> A.satisfy (A.inClass "ACTG")
    xAllele = A.string "<X>" >> return Nothing
    coverage = do
        _ <- A.decimal `A.sepBy1` (A.char ',')
        _ <- A.char ':'
        A.decimal `A.sepBy1` (A.char ',')

word :: A.Parser Text
word = T.pack <$> A.many1 (A.satisfy (not . isSpace))

tab :: A.Parser ()
tab = A.char '\t' >> return ()

processVcfSimple :: CallingMode -> Int -> FilterMode -> Pipe VCFentry FreqSumRow (SafeT IO) r
processVcfSimple mode minDepth filter_ = for cat $ \vcfEntry -> do
    let (VCFentry chrom pos alleles covNums) = vcfEntry
    (normalizedAlleles, normalizedCovNums) <- case alleles of
        [ref] -> lift . lift . throwIO $ AssertionFailed "should not happen, need at least one alternative allele"
        [ref, alt] -> return ([ref, alt], covNums)
        _ -> do
            let altNumPairs = [(alleles!!i, sum [c!!i | c <- covNums]) | i <- [1 .. (length alleles - 1)]]
            shuffledAltNumPairs <- lift . lift $ shuffle altNumPairs
            let (alt, _) = head . sortBy (\a b -> snd b `compare` snd a) $ shuffledAltNumPairs
            let altIndex = snd . head . filter ((==alt) . fst) $ zip alleles [0..]
            when (altIndex == 0) $ (lift . lift . throwIO) (AssertionFailed "should not happen, altIndex==0")
            return ([head alleles, alt], [[c !! 0, c !! altIndex] | c <- covNums])
    let [ref, alt] = normalizedAlleles
    when (ref /= 'N' && (not transversionsOnly || isTransversion ref alt)) $ do
        genotypes <- lift . lift $ mapM (callGenotype mode minDepth) normalizedCovNums
        when (any (>0) genotypes) $ yield (FreqSumRow chrom pos ref alt genotypes)
  where
    shuffle list = evalRandIO (shuffleM list)
    transversionsOnly = filter_ /= NoFilter
    
isTransversion :: Char -> Char -> Bool
isTransversion ref alt = not isTransition
  where
    isTransition = ((ref == 'A') && (alt == 'G')) || ((ref == 'G') && (alt == 'A')) ||
                   ((ref == 'C') && (alt == 'T')) || ((ref == 'T') && (alt == 'C'))

callGenotype :: CallingMode -> Int -> [Int] -> IO Int
callGenotype mode minDepth covs = do
    case covs of
        [numRef] -> return 0
        [numRef, numAlt] -> do
            if (numRef + numAlt < minDepth) then
                return (-1)
                else
                    case mode of
                        MajorityCalling -> case numRef `compare` numAlt of
                            LT -> return 2
                            GT -> return 0
                            EQ -> do
                                rn <- randomRIO (1, numRef + numAlt)
                                if rn <= numRef then return 0 else return 2
                        RandomCalling -> do
                                rn <- randomRIO (1, numRef + numAlt)
                                if rn <= numRef then return 0 else return 2
        _ -> return (-1)

processVcfWithSnpFile :: Int -> CallingMode -> Int -> FilterMode ->
                                  Pipe (Maybe SnpEntry, Maybe VCFentry) FreqSumRow (SafeT IO) r
processVcfWithSnpFile nrInds mode minDepth filter_ = for cat $ \jointEntry -> do
    -- trace (show jointEntry) (return ())
    case jointEntry of
        (Just (SnpEntry snpChrom snpPos snpRef snpAlt), Nothing) -> do
            yield $ FreqSumRow snpChrom snpPos snpRef snpAlt (replicate nrInds (-1))
        (Just (SnpEntry snpChrom snpPos snpRef snpAlt), Just (VCFentry vcfChrom vcfPos vcfAlleles vcfNums)) -> do
            let normalizedAlleleI = map snd . filter (\(a, _) -> a == snpRef || a == snpAlt) $ zip vcfAlleles [0..]
                normalizedVcfAlleles = map (vcfAlleles!!) normalizedAlleleI
                normalizedVcfNums = [map (v!!) normalizedAlleleI | v <- vcfNums]
            genotypes <- case normalizedVcfAlleles of
                    [ref] -> if ref == snpRef then return (replicate nrInds 0) else return (replicate nrInds 2)
                    [ref, alt] -> if [ref, alt] == [snpRef, snpAlt]
                        then lift . lift $ mapM (callGenotype mode minDepth) normalizedVcfNums
                        else lift . lift $ mapM (callGenotype mode minDepth) (reverse normalizedVcfNums)
                    _ -> lift . lift . throwIO $ AssertionFailed "should not happen, can only have two alleles after normalization"
            case filter_ of
                NoFilter -> yield (FreqSumRow snpChrom snpPos snpRef snpAlt genotypes)
                Transversions -> when (isTransversion snpRef snpAlt) $
                        yield (FreqSumRow snpChrom snpPos snpRef snpAlt genotypes)
                TransversionsMissing -> if isTransversion snpRef snpAlt
                        then yield (FreqSumRow snpChrom snpPos snpRef snpAlt genotypes)
                        else yield (FreqSumRow snpChrom snpPos snpRef snpAlt (replicate nrInds (-1)))
        _ -> return ()

snpParser :: A.Parser SnpEntry
snpParser = do
    chrom <- word
    tab
    pos <- A.decimal
    tab
    ref <- A.satisfy (A.inClass "ACTG")
    _ <- A.char ','
    alt <- A.satisfy (A.inClass "ACTG")
    _ <- A.satisfy (\c -> c == '\r' || c == '\n')
    let ret = SnpEntry chrom pos ref alt
    -- trace (show ret) $ return ()
    return $ SnpEntry chrom pos ref alt

printFreqSum :: FreqSumRow -> Text
printFreqSum (FreqSumRow chrom pos ref alt calls) =
    format (s%"\t"%d%"\t"%s%"\t"%s%"\t"%s) chrom pos (T.singleton ref) (T.singleton alt) callsStr
  where
    callsStr = (T.intercalate "\t" . map (format d)) calls

printEigenStrat :: FreqSumRow -> Text
printEigenStrat (FreqSumRow _ _ _ _ calls) = (T.concat . map (format d) . map toEigenStratNum) calls
  where
    toEigenStratNum c = case c of
        0 -> 2
        1 -> 1
        2 -> 0
        -1 -> 9
        _ -> error ("unknown genotype " ++ show c)

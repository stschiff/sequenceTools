{-# LANGUAGE OverloadedStrings #-}

import SeqTools.OrderedZip (orderedZip)
import SeqTools.VCF (readVCF, VCFheader(..), VCFentry(..))

import Control.Error (headErr)
import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad.Random (evalRandIO)
import Control.Monad.Trans.Class (lift)
import Control.Monad.Trans.Reader (ReaderT, runReaderT, asks)
import qualified Data.Attoparsec.Text as A
import Data.Char (isSpace)
import Data.List (sortBy)
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Debug.Trace (trace)
import qualified Options.Applicative as OP
import Pipes (Pipe, yield, (>->), runEffect, Producer, Pipe, for, cat)
import Pipes.Attoparsec (parsed)
import Pipes.Cliff (CreateProcess(..), CmdSpec(..), pipeOutput, NonPipe(..))
import Pipes.Cliff.Core (defaultHandler)
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, SafeT)
import Pipes.Safe.Prelude (withFile)
import Pipes.Text.Encoding (decodeUtf8)
import qualified Pipes.Text.IO as PT
import Prelude hiding (FilePath)
import System.IO (IOMode(..))
import System.Random (randomRIO, mkStdGen, setStdGen)
import System.Random.Shuffle (shuffleM)
import Turtle hiding (tab, cat, stderr)

data ProgOpt = ProgOpt {
    optCallingMode :: CallingMode,
    optSeed :: Maybe Int,
    optMinDepth :: Int,
    optTransversionsOnly :: Bool,
    optSnpFile :: Maybe FilePath,
    optRegion :: Text,
    optOutChrom :: Maybe Text,
    optReference :: FilePath,
    optOutFormat :: OutFormat,
    optEigenstratOutPrefix :: Maybe FilePath,
    optSamtoolsExe :: FilePath,
    optBcftoolsExe :: FilePath,
    optBamFiles :: [FilePath]
}

data CallingMode = MajorityCalling | RandomCalling | RareCalling deriving (Show, Read)
data OutFormat = EigenStrat | FreqSumFormat deriving (Show, Read)
data VCFentryReduced = VCFentryReduced T.Text Int [T.Text] [[Int]] deriving (Show)
data FreqSumRow = FreqSumRow T.Text Int T.Text T.Text [Int] deriving (Show)
data SnpEntry = SnpEntry T.Text Int T.Text T.Text deriving (Show)-- Chrom Pos Ref Alt

main :: IO ()
main = OP.execParser parser >>= runSafeT . runReaderT runWithOpts
  where
    parser = OP.info (OP.helper <*> argParser)
                     (OP.progDesc "A program to perform simple genotype calling directly from BAM")

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseCallingMode <*> parseSeed <*> parseMinDepth <*>
                        parseTransversionsOnly <*> parseSnpFile <*> parseChrom <*>
                        parseOutChrom <*> parseRef <*> parseFormat <*> parseEigenstratOutPrefix <*> 
                        parseSamtoolsExe <*> parseBcftoolsExe <*> OP.some parseBams
  where
    parseCallingMode = OP.option OP.auto (OP.long "mode" <> OP.short 'm' <>
                       OP.value RandomCalling <> OP.showDefault <> OP.metavar "<MODE>" <>
                       OP.help "specify the mode of calling: MajorityCalling, RandomCalling or \
                                \RareCalling. MajorityCalling: Pick the allele supported by the \   
                                \most reads. If equal numbers of Alleles fulfil this, pick one at \
                                \random. RandomCalling: Pick one read at random. RareCalling: \
                                \Require a number of reads equal to the minDepth supporting the \
                                \alternative allele to call a heterozygote. Otherwise call \
                                \homozygous reference or missing depending on depth. For \
                                \RareCalling you should use --minDepth 2.")
    parseSeed = OP.option (Just <$> OP.auto) (OP.long "seed" <> OP.value Nothing <>
                OP.metavar "<RANDOM_SEED>" <> OP.help "random seed used for random calling. If \
                \not given, use system clock to seed the random number generator")
    parseMinDepth = OP.option OP.auto (OP.long "minDepth" <> OP.short 'd' <> OP.value 1 <> 
                                       OP.showDefault <>
                                       OP.metavar "<DEPTH>" <>
                                       OP.help "specify the minimum depth for a call. This has a \
                                                \special meaning for RareCalling, see --mode")
    parseTransversionsOnly = OP.switch (OP.long "transversionsOnly" <> OP.short 't' <>
                    OP.help "Remove transition SNPs from the output)")
    parseSnpFile = OP.option (Just . fromText . T.pack <$> OP.str)
                   (OP.long "snpFile" <> OP.short 'f' <> OP.value Nothing <> OP.metavar "<FILE>" <>
                    OP.help "specify a SNP file for the positions and alleles to call. All \
                             \positions in the SNP file will be output, adding missing data where \
                             \necessary. The file should have three columns (space- or \
                             \tab-separated): Chromosome, position and \
                             \alleles, where the alleles should be one reference and one \
                             \alternative allele \
                             \separated by a comma. Note that this file needs to have a difference \
                             \format than \
                             \the snp files in Eigenstrat because of samtools specifications. Note \
                             \also that simpleBamCaller automatically checks whether alleles in \
                             \the SNP file are flipped with respect to the human reference. But \
                             \it assumes that the strand-orientation is the same.")
    parseChrom = OP.option (T.pack <$> OP.str) (OP.long "chrom" <> OP.short 'c' <>
                            OP.metavar "<CHROM>" <> OP.help "specify the region in the BAM file to \
                            \call from. Can be just the chromosome, or a string of the form \
                            \CHROM:START-END")
    parseOutChrom = OP.option (Just . T.pack <$> OP.str) (OP.long "outChrom" <>
                                    OP.metavar "<CHROM>" <>
                                   OP.help "specify the output chromosome name" <> OP.value Nothing)
    parseRef = OP.option (fromText . T.pack <$> OP.str) (OP.long "reference" <> OP.short 'r' <> 
                                   OP.metavar "<REF>" <>
                                  OP.help "the reference fasta file")
    parseBams = OP.argument (fromText . T.pack <$> OP.str) (OP.metavar "<BAM_FILE(S)>" <>
                                    OP.help "input file, give multiple files for multiple samples")
    parseFormat = OP.option OP.auto (OP.metavar "<OUT_FORMAT>" <> OP.long "format" <>
                                     OP.short 'o' <>
                                     OP.value FreqSumFormat <> OP.showDefault <>
                                     OP.help "specify output format: EigenStrat or FreqSum")
    parseEigenstratOutPrefix = OP.option (Just . fromText . T.pack <$> OP.str)
                                 (OP.long "eigenstratOutPrefix" <> 
                                    OP.short 'e' <>
                                 OP.value Nothing <> OP.metavar "<FILE_PREFIX>" <>
                                 OP.help "specify the filenames for the EigenStrat SNP and IND \
                                 \file outputs: <FILE_PREFIX>.snp.txt and <FILE_PREFIX>.ind.txt \
                                 \Ignored if Output format is not Eigenstrat")
    parseSamtoolsExe = OP.option (fromText . T.pack <$> OP.str) (OP.long "samtools" <>
                                OP.value "samtools" <> OP.showDefault <>
                                  OP.metavar "<SAMTOOLS_PATH>" <>
                                  OP.help "path to the samtools executable, version >= 1.2")
    parseBcftoolsExe = OP.option (fromText . T.pack <$> OP.str) (OP.long "bcftools" <>
                                  OP.value "bcftools" <> OP.showDefault <>
                                  OP.metavar "<BCFTOOLS_PATH>" <> OP.help "path to the \
                                  \bcftools executable, version >= 1.2")
    
runWithOpts :: ReaderT ProgOpt (SafeT IO) ()
runWithOpts = do
    seed <- asks optSeed
    region <- asks optRegion
    outChrom <- asks optOutChrom
    outFormat <- asks optOutFormat
    eigenStratOutPrefix <- asks optEigenstratOutPrefix
    transversionsOnly <- asks optTransversionsOnly
    case seed of
        Nothing -> return ()
        Just seed_ -> liftIO . setStdGen $ mkStdGen seed_
    let chrom = head (T.splitOn ":" region)
    let outputChrom = case outChrom of
            Nothing -> chrom
            Just label -> label
    (vcfHeader, freqSumProducer) <- runPileup
    -- liftIO $ print vcfHeader
    case outFormat of
        FreqSumFormat -> do
            let VCFheader _ n = vcfHeader
            echo $ format ("#CHROM\tPOS\tREF\tALT\t"%s)
                          (T.intercalate "\t" . map (format (s%"(2)")) $ n)
            lift . runEffect $ freqSumProducer >-> filterTransitions transversionsOnly >->
                               P.map (showFreqSum outputChrom) >-> printToStdOut
        EigenStrat -> case eigenStratOutPrefix of
            Nothing -> liftIO . throwIO $ AssertionFailed "need an eigenstratPrefix for \
                                                           \EigenStratFormat"
            Just fn -> do
                let snpOut = fn <.> "snp.txt"
                    indOut = fn <.> "ind.txt"
                lift . withFile (T.unpack . format fp $ indOut) WriteMode $ \indOutHandle -> do
                    let VCFheader _ sampleNames = vcfHeader
                    mapM_ (\n -> liftIO $ T.hPutStrLn indOutHandle (format (s%"\tU\tUnknown") n))
                          sampleNames
                lift . withFile (T.unpack . format fp $ snpOut) WriteMode $ \snpOutHandle -> do
                    runEffect $ freqSumProducer >-> filterTransitions transversionsOnly >->
                                printEigenStrat outputChrom snpOutHandle >-> printToStdOut
  where
    printToStdOut = for cat (liftIO . T.putStrLn)
    filterTransitions transversionsOnly =
        if transversionsOnly
        then P.filter (\(FreqSumRow _ _ ref alt _) -> isTransversion ref alt)
        else cat
    isTransversion ref alt = not $ isTransition ref alt
    isTransition ref alt = ((ref == "A") && (alt == "G")) || ((ref == "G") && (alt == "A")) ||
                           ((ref == "C") && (alt == "T")) || ((ref == "T") && (alt == "C"))
    
runPileup :: ReaderT ProgOpt (SafeT IO) (VCFheader, Producer FreqSumRow (SafeT IO) ())
runPileup = do
    snpFile <- asks optSnpFile
    case snpFile of
        Nothing -> runPileupSimple
        Just fn -> runPileupSnpFile fn

runPileupSimple :: ReaderT ProgOpt (SafeT IO) (VCFheader, Producer FreqSumRow (SafeT IO) ())
runPileupSimple = do
    samtools <- asks optSamtoolsExe
    bcftools <- asks optBcftoolsExe
    reference <- asks optReference
    region <- asks optRegion
    bamFiles <- asks optBamFiles
    mode <- asks optCallingMode
    minDepth <- asks optMinDepth
    let bams = (T.intercalate " " (map (format fp) bamFiles))
    let cmd = format (fp%" mpileup -q30 -Q30 -C50 -I -f "%fp%" -g -t DPR -r "%s%" "%s%
                      " | "%fp%" view -v snps") samtools reference region bams bcftools
    vcfTextProd <- liftIO $ produceFromCommand cmd
    (vcfHeader_, vcfProd) <- lift $ readVCF vcfTextProd
    let vcfProdPipe = vcfProd >-> processVcfSimple (length bamFiles) mode minDepth
    return (vcfHeader_, vcfProdPipe)

runPileupSnpFile :: FilePath ->
                    ReaderT ProgOpt (SafeT IO) (VCFheader, Producer FreqSumRow (SafeT IO) ())
runPileupSnpFile fn = do
    samtools <- asks optSamtoolsExe
    bcftools <- asks optBcftoolsExe
    reference <- asks optReference
    region <- asks optRegion
    bamFiles <- asks optBamFiles
    mode <- asks optCallingMode
    minDepth <- asks optMinDepth
    let bams = (T.intercalate " " (map (format fp) bamFiles))
    let chrom = head (T.splitOn ":" region)
    let cmd = format (fp%" mpileup -q30 -Q30 -C50 -I -f "%fp%" -g -t DPR -r "%s%" -l "%fp%
                      " "%s%" | "%fp%" view") samtools reference region fn bams bcftools
    vcfTextProd <- liftIO $ produceFromCommand cmd
    let snpTextProd = PT.readFile ((T.unpack . format fp) fn)
    (vcfHeader_, vcfProd) <- lift $ readVCF vcfTextProd
    let snpProd = parsed snpParser snpTextProd >-> P.filter (\(SnpEntry c _ _ _) -> c == chrom)
        jointProd = orderedZip cmp snpProd vcfProd
        jointProdPipe = jointProd >-> processVcfWithSnpFile (length bamFiles) mode minDepth
    return (vcfHeader_, fmap snd jointProdPipe)
  where
    cmp (SnpEntry _ snpPos _ _) vcfEntry = snpPos `compare` (vcfPos vcfEntry)
    
produceFromCommand :: Text -> IO (Producer Text (SafeT IO) ())
produceFromCommand cmd = do
    let createProcess = CreateProcess (ShellCommand (T.unpack cmd)) Nothing Nothing False False 
                                                        False False False False Nothing Nothing defaultHandler
    (p, _) <- pipeOutput Inherit Inherit createProcess
    return . void . decodeUtf8 $ p

processVcfSimple :: Int -> CallingMode -> Int -> Pipe VCFentry FreqSumRow (SafeT IO) r
processVcfSimple nrInds mode minDepth = for cat $ \vcfEntry -> do
    let Right (VCFentryReduced chrom pos alleles covNums) = makeReducedVCF vcfEntry
    -- trace (show vcfEntry) $ return ()
    -- trace (show (VCFentryReduced chrom pos alleles covNums)) $ return ()
    when (length covNums /= nrInds) $ (liftIO . throwIO) (AssertionFailed "inconsistent number \
            \of genotypes. Check that bam files have different readgroup sample names")
    (normalizedAlleles, normalizedCovNums) <- case alleles of
        [_] -> liftIO . throwIO $ AssertionFailed "should not happen, need at least one \
                            \alternative allele"
        [ref, alt] -> return ([ref, alt], covNums)
        _ -> do
            let altNumPairs =
                    [(alleles!!i, sum [c!!i | c <- covNums]) | i <- [1 .. (length alleles - 1)]]
            shuffledAltNumPairs <- liftIO $ shuffle altNumPairs
            let (alt, _) = head . sortBy (\a b -> snd b `compare` snd a) $ shuffledAltNumPairs
            let altIndex = snd . head . filter ((==alt) . fst) $ zip alleles [0..]
            when (altIndex == 0) $ (liftIO . throwIO)
                            (AssertionFailed "should not happen, altIndex==0")
            return ([head alleles, alt], [[c !! 0, c !! altIndex] | c <- covNums])
    let [ref, alt] = normalizedAlleles
    when (ref /= "N") $ do
        genotypes <- liftIO $ mapM (callGenotype mode minDepth) normalizedCovNums
        when (any (>0) genotypes) $ yield (FreqSumRow chrom pos ref alt genotypes)
  where
    shuffle list = evalRandIO (shuffleM list)

makeReducedVCF :: VCFentry -> Either String VCFentryReduced
makeReducedVCF (VCFentry chrom pos _ ref alt _ _ _ formatS genotypes) = do
    dprIndex <- fmap fst . headErr "Did not find DPR tag in format string" .
                        filter ((=="DPR") . snd) . zip [0..] $ formatS
    let covNums = map (getCorrectCovNums . (!!dprIndex)) genotypes
    return $ VCFentryReduced chrom pos normalizedAlleles covNums
  where
    getCorrectCovNums covNumStr =
        (\v -> map (v!!) normalizedAlleleIndices) . map (read . T.unpack) . T.splitOn "," $ 
        covNumStr
    normalizedAlleleIndexPairs = filter (\a -> snd a /= "<X>" && snd a /= "<*>") . zip [0..] $ 
                                 ref : alt
    normalizedAlleles = map snd normalizedAlleleIndexPairs
    normalizedAlleleIndices = map fst normalizedAlleleIndexPairs
    

callGenotype :: CallingMode -> Int -> [Int] -> IO Int
callGenotype mode minDepth covs = do
    if sum covs < minDepth then return (-1) else do
        case covs of
            [_] -> return 0
            [numRef, numAlt] -> do
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
                    RareCalling -> do
                            if numAlt >= 2 then return 1 else return 0
            _ -> throwIO (AssertionFailed "should not happen. CallGenotype called with more \
                                            \than two alleles")

processVcfWithSnpFile :: Int -> CallingMode -> Int -> Pipe (Maybe SnpEntry, Maybe VCFentry) 
                                                      FreqSumRow (SafeT IO) r
processVcfWithSnpFile nrInds mode minDepth = for cat $ \jointEntry -> do
    -- trace (show jointEntry) (return ())
    case jointEntry of
        (Just (SnpEntry snpChrom snpPos snpRef snpAlt), Nothing) -> do
            yield $ FreqSumRow snpChrom snpPos snpRef snpAlt (replicate nrInds (-1))
        (Just (SnpEntry snpChrom snpPos snpRef snpAlt), Just vcfEntry) -> do
            let Right (VCFentryReduced _ _ vcfAlleles vcfNums) = makeReducedVCF vcfEntry
            when (length vcfNums /= nrInds) $ (liftIO . throwIO) (AssertionFailed "inconsistent \ 
                            \number of genotypes. Check that bam files have different \
                            \readgroup sample names")
            let normalizedAlleleI =
                    map snd . filter (\(a, _) -> a == snpRef || a == snpAlt) $ zip vcfAlleles [0..]
                normalizedVcfAlleles = map (vcfAlleles!!) normalizedAlleleI
                normalizedVcfNums = [map (v!!) normalizedAlleleI | v <- vcfNums]
            -- trace (show (VCFentry vcfChrom vcfPos vcfAlleles vcfNums)) (return ())
            genotypes <- case normalizedVcfAlleles of
                    [] -> return (replicate nrInds (-1))
                    [ref] -> if ref == snpRef
                        then return [if sum c >= minDepth then 0 else (-1) | c <- normalizedVcfNums]
                        else return [if sum c >= minDepth then 2 else (-1) | c <- normalizedVcfNums]
                    [ref, alt] -> if [ref, alt] == [snpRef, snpAlt]
                        then liftIO $ mapM (callGenotype mode minDepth) normalizedVcfNums
                        else liftIO $ mapM (callGenotype mode minDepth)
                                        (map reverse normalizedVcfNums)
                    _ -> liftIO . throwIO $ AssertionFailed ("should not happen, can only have \
                                    \two alleles after normalization: " ++ show jointEntry)
            yield (FreqSumRow snpChrom snpPos snpRef snpAlt genotypes)
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
    let ret = SnpEntry chrom pos (T.singleton ref) (T.singleton alt)
    -- trace (show ret) $ return ()
    return ret
  where
    word = T.pack <$> A.many1 (A.satisfy (not . isSpace))
    tab = A.char '\t' >> return ()
    
showFreqSum :: Text -> FreqSumRow -> Text
showFreqSum outChrom (FreqSumRow _ pos ref alt calls) =
    format (s%"\t"%d%"\t"%s%"\t"%s%"\t"%s) outChrom pos ref alt callsStr
  where
    callsStr = (T.intercalate "\t" . map (format d)) calls

printEigenStrat :: Text -> Handle -> Pipe FreqSumRow Text (SafeT IO) r
printEigenStrat outChrom snpOutHandle = for cat $ \(FreqSumRow _ pos ref alt calls) -> do
    let n = format (s%"_"%d) outChrom pos
        snpLine = format (s%"\t"%s%"\t0\t"%d%"\t"%s%"\t"%s) n outChrom pos ref alt
    liftIO . T.hPutStrLn snpOutHandle $ snpLine
    yield . T.concat . map (format d . toEigenStratNum) $ calls
  where
    toEigenStratNum c = case c of
        0 -> 2 :: Int
        1 -> 1
        2 -> 0
        -1 -> 9
        _ -> error ("unknown genotype " ++ show c)

{-# LANGUAGE OverloadedStrings #-}

import Pipes.OrderedZip (orderedZip)
import SequenceFormats.Utils (liftParsingErrors, Chrom(..))
import SequenceFormats.Eigenstrat (readEigenstratSnpFile, EigenstratSnpEntry(..), GenoLine, 
    EigenstratIndEntry(..), Sex(..), writeEigenstrat, GenoEntry(..))
import SequenceFormats.FreqSum(FreqSumEntry(..), printFreqSumStdOut, FreqSumHeader(..))
import SequenceFormats.Pileup (PileupRow(..), readPileupFromStdIn)

import SequenceTools.Utils (sampleWithoutReplacement)

import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad.IO.Class (liftIO)
import Control.Monad.Trans.Class (lift)
import Control.Monad.Trans.Reader (ReaderT, runReaderT, asks)
import qualified Data.Attoparsec.ByteString.Char8 as A
import qualified Data.ByteString.Char8 as B
import Data.Char (isSpace, isDigit, toUpper)
import Data.List (partition, sortOn, group, sort)
import Data.Maybe (fromMaybe, catMaybes)
import Data.Version (showVersion)
import Data.Vector (fromList)
-- import Debug.Trace (trace)
import qualified Options.Applicative as OP
import Paths_sequenceTools (version)
import Pipes (Pipe, yield, (>->), runEffect, Producer, Pipe, for, cat)
import Pipes.Attoparsec (parsed)
import qualified Pipes.ByteString as PB
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, SafeT)
import System.Random (randomRIO, mkStdGen, setStdGen)

data ProgOpt = ProgOpt {
    optCallingMode :: CallingMode,
    optSeed :: Maybe Int,
    optMinDepth :: Int,
    optTransitionsMode :: TransitionsMode,
    optSnpFile :: Maybe FilePath,
    optOutFormat :: OutFormat,
    optSampleNames :: Either [String] FilePath,
    optSamplePopName :: String
}

data TransitionsMode = TransitionsMissing | SkipTransitions | AllSites
data OutFormat = EigenStratFormat String | FreqSumFormat
type App = ReaderT ProgOpt (SafeT IO)

main :: IO ()
main = OP.execParser parserInfo >>= runSafeT . runReaderT runWithOpts
  where
    parserInfo = OP.info (pure (.) <*> versionInfoOpt <*> OP.helper <*> argParser) (OP.progDesc progInfo)
    versionInfoOpt = OP.infoOption (showVersion version) (OP.long "version" <> OP.help "Print version and exit")

progInfo :: String
progInfo = "A program to perform genotype calling from a pileup file. Part of sequenceTools \
    \version " ++ showVersion version

runWithOpts :: App ()
runWithOpts = do
    setRandomSeed
    let pileupProducer = readPileupFromStdIn
    snpFile <- asks optSnpFile
    freqSumProducer <- pileupToFreqSum snpFile pileupProducer
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

pileupToFreqSum :: FilePath -> Producer PileupRow (SafeT IO) () -> App (Producer FreqSumEntry (SafeT IO) ())
pileupToFreqSum snpFileName pileupProducer = do
    sampleNameSpec <- asks optSampleNames
    sampleNames <- case sampleNameSpec of
        Left list -> return list
        Right fn -> lines <$> readFile fn
    let snpProd = readEigenstratSnpFile snpFileName
        jointProd = orderedZip cmp snpProd pileupProducer
    mode <- optCallingMode
    minDepth <- asks optMinDepth
    let ret = for jointProd $ \jointEntry ->
            case jointEntry of
                (Just (EigenstratSnpEntry snpChrom snpPos _ _ snpRef snpAlt), Nothing) ->
                    yield $ FreqSumEntry snpChrom snpPos snpRef snpAlt
                            (replicate (length sampleNames) Nothing)
                (Just (EigenstratSnpEntry snpChrom snpPos _ _ snpRef snpAlt),
                 Just pileupRow) -> do
                    let PileupRow _ _ refA entryPerSample = pileupRow
                    calls <- liftIO $ mapM (callGenotype mode minDepth refA) (map B.unpack entryPerSample)
                    let genotypes = map (callToGenotype snpRef snpAlt) calls
                    yield (FreqSumEntry snpChrom snpPos snpRef snpAlt genotypes)
                _ -> return ()
    return (fst <$> ret)
  where
    cmp (EigenstratSnpEntry snpChrom snpPos _ _ _ _) (PileupRow pChrom pPos _ _) =
        (snpChrom, snpPos) `compare` (pChrom, pPos)

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
            RandomDiploidCalling -> do
                res <- sampleWithoutReplacement alleles 2
                case res of
                    Nothing -> return MissingCall
                    Just [a1, a2] -> return $ DiploidCall a1 a2
                    _ -> error "should not happen"


outputFreqSum :: Producer FreqSumEntry (SafeT IO) () -> App ()
outputFreqSum freqSumProducer = do
    outChrom <- asks optOutChrom
    transitionsOnly <- asks optTransitionsOnly
    sampleNameSpec <- asks optSampleNames
    callingMode <- getCallingMode
    sampleNames <- case sampleNameSpec of
        Left list -> return list
        Right fn -> liftIO $ lines (readFile fn)
    let nrHaplotypes = case callingMode of
            MajorityCalling _ -> 1 :: Int
            RandomCalling -> 1
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
    sampleNames <- case sampleNameSpec of
        Left list -> return list
        Right fn -> lines <$> (liftIO . readFile) fn
    eigenStratOutPrefix <- asks optEigenstratOutPrefix
    case eigenStratOutPrefix of
        Nothing -> liftIO . throwIO $
            AssertionFailed "need an eigenstratPrefix for EigenStratFormat"
        Just fn -> do
            let snpOut = fn <> "snp.txt"
                indOut = fn <> "ind.txt"
                genoOut = fn <> "geno.txt"
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

toEigenstrat :: Bool -> Maybe Chrom ->
    Pipe FreqSumEntry (EigenstratSnpEntry, GenoLine) (SafeT IO) ()
toEigenstrat diploidizeCall outChrom = P.map toEigenstrat'
  where
    toEigenstrat' (FreqSumEntry chrom pos ref alt calls) =
        let newChrom = fromMaybe chrom outChrom
            snpId = format (s%"_"%d) (unChrom newChrom) pos
            snpEntry = EigenstratSnpEntry newChrom pos 0.0 snpId ref alt
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


-- ProgOpt {
--     optCallingMode :: CallingMode,
--     optSeed :: Maybe Int,
--     optMinDepth :: Int,
--     optTransitionsMode :: TransitionsMode,
--     optSnpFile :: Maybe FilePath,
--     optOutFormat :: OutFormat,
--     optSampleNames :: Either [String] FilePath,
--     optSamplePopName :: String
-- }

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseCallingMode <*> parseSeed <*> parseMinDepth <*>
    parseTransitionsMode <*> parseSnpFile <*> parseFormat <*> parseSampleNames <*> 
    parseSamplePopName
  where
    parseCallingMode = parseRandomCalling <|> parseMajorityCalling <|> parseRandomDiploidCalling
    parseRandomCalling = OP.flag' RandomCalling (OP.long "randomHaploid" <>
        OP.help "This method samples one read at random at each site, and uses the allele \
        \on that read as the one for the actual genotype. This results in a haploid \
        \call")
    parseRandomDiploidCalling = OP.flag' RandomDiploiCalling (OP.long "randomDiploid" <>
        OP.help "Sample two random reads at random and represent the \
        \individual by the diploid genotype constructed from those two random \
        \picks. This will always assign missing data to positions where only \
        \one read is present, even if minDepth=1.")
    parseMajorityCalling = MajorityCalling <$> (parseMajorityCallingFlag *> parseDownsamplingFlag)
    parseMajorityCallingFlag = OP.flag' True (OP.long "majorityCall" <> OP.help
        "Pick the allele supported by the \
        \most reads at a site. If an equal numbers of Alleles fulfil this, pick one at \
        \random. This results in a haploid call.")
    parseDownsamplingFlag = OP.switch (OP.long "downSampling" <> OP.help "When this switch is given, the MajorityCalling mode with downsample \
        \from the total number of reads a number of reads \
        \(without replacement) equal to the --minDepth given. This mitigates \
        \reference bias in the MajorityCalling model, which increases with higher coverage.")
    parseSeed = OP.option (Just <$> OP.auto) (OP.long "seed" <>
        OP.value Nothing <> OP.metavar "<RANDOM_SEED>" <>
        OP.help "random seed used for random calling. If not given, use \
        \system clock to seed the random number generator")
    parseMinDepth = OP.option OP.auto (OP.long "minDepth" <> OP.short 'd' <>
        OP.value 1 <> OP.showDefault <> OP.metavar "<DEPTH>" <>
        OP.help "specify the minimum depth for a call. For sites with fewer \
            \reads than this number, declare Missing")
    parseTransitionsMode = parseSkipTransitions <|> parseTransitionsMissing <|> pure AllSites
    parseSkipTransitions = OP.flag' (OP.long "skipTransitions" <> OP.help "skip transition SNPs entirely in the output, resulting in a dataset with fewer sites.\
        \ If neither option --skipTransitions nor --transitionsMissing is set, output all sites")
    parseTransitionsMissing = OP.flag' (OP.long "transitionsMissing" <> OP.help "mark transitions as missing in the output, but do output the sites. \
        \If neither option --skipTransitions nor --transitionsMissing is set, output all sites")
    parseSnpFile = OP.strOption (OP.long "snpFile" <> OP.short 'f' <> OP.value Nothing <>
        OP.metavar "<FILE>" <> OP.help "specify an Eigenstrat SNP file for \
        \the positions and alleles to call. All \
        \positions in the SNP file will be output, adding missing data where \
        \necessary. Note that pileupCaller automatically checks whether \
        \alleles in the SNP file are flipped with respect to the human \
        \reference. But it assumes that the strand-orientation is the same.")
    parseFormat = parseEigenstratPrefix <|> pure FreqSumFormat
    parseEigenstratPrefix = OP.strOption (OP.long "eigenstratOut" <> OP.short 'e' <>
        OP.metavar "<FILE_PREFIX>" <>
        OP.help "Set Eigenstrat as output format. Specify the filenames for the EigenStrat SNP and IND \
        \file outputs: <FILE_PREFIX>.snp.txt and <FILE_PREFIX>.ind.txt \
        \If not set, output will be FreqSum (Default)")
    parseSampleNames = parseSampleNameList <|> parseSampleNameFile
    parseSampleNameList = OP.option (Left . splitOn "," <$> OP.str)
        (OP.long "sampleNames" <> OP.metavar "NAME1,NAME2,..." <>
        OP.help "give the names of the samples as comma-separated list")
    parseSampleNameFile = OP.strOption (OP.long "sampleNameFile" <> OP.metavar "<FILE>" <>
        OP.help "give the names of the samples in a file with one name per \
        \line")
    parseSamplePopName = OP.strOption (OP.long "samplePopName" <> OP.value "Unknown" <> OP.showDefault <>
        OP.metavar "POP" <>
        OP.help "specify the population name of the samples, which is included\
        \ in the output *.ind.txt file. This will be ignored if the output \
        \format is not Eigenstrat")

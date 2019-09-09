{-# LANGUAGE OverloadedStrings, NoImplicitPrelude #-}

import SequenceFormats.Eigenstrat (readEigenstratSnpFile, EigenstratSnpEntry(..), 
    EigenstratIndEntry(..), Sex(..), writeEigenstrat)
import SequenceFormats.FreqSum(FreqSumEntry(..), printFreqSumStdOut, FreqSumHeader(..))
import SequenceFormats.Pileup (PileupRow(..), readPileupFromStdIn)
import SequenceTools.Utils (versionInfoOpt, versionInfoText)
import SequenceTools.PileupCaller (CallingMode(..), callGenotypeFromPileup, callToDosage,
    freqSumToEigenstrat, filterTransitions, TransitionsMode(..), cleanSSdamageAllSamples)

import qualified Data.ByteString.Char8 as B
import Data.List.Split (splitOn)
import qualified Data.Vector.Unboxed.Mutable as V
import qualified Options.Applicative as OP
import Pipes (yield, (>->), runEffect, Producer, for)
import Pipes.OrderedZip (orderedZip, orderCheckPipe)
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, SafeT)
import RIO
import System.IO (readFile, hPutStrLn, stderr)
import System.Random (mkStdGen, setStdGen)
import Text.Printf (printf)
import qualified Text.PrettyPrint.ANSI.Leijen as PP

data OutFormat = EigenstratFormat String String | FreqSumFormat

data ProgOpt = ProgOpt CallingMode Bool (Maybe Int) Int TransitionsMode FilePath OutFormat (Either [String] FilePath)
    --optCallingMode :: CallingMode,
    --optKeepInCongruentReads :: Bool,
    --optSeed :: Maybe Int,
    --optMinDepth :: Int,
    --optTransitionsMode :: TransitionsMode,
    --optSnpFile :: FilePath,
    --optOutFormat :: OutFormat,
    --optSampleNames :: Either [String] FilePath

data ReadStats = ReadStats {
    rsTotalSites :: IORef Int,
    rsNonMissingSites :: V.IOVector Int,
    rsRawReads :: V.IOVector Int,
    rsReadsCleanedSS :: V.IOVector Int,
    rsReadsCongruent :: V.IOVector Int
}

data Env = Env {
    envCallingMode :: CallingMode,
    envKeepInCongruentReads :: Bool,
    envMinDepth :: Int,
    envTransitionsMode :: TransitionsMode,
    envOutFormat :: OutFormat,
    envSnpFile :: FilePath,
    envSampleNames :: [String],
    envStats :: ReadStats
}

type App = ReaderT Env (SafeT IO)

main :: IO ()
main = do
    args <- OP.execParser parserInfo
    env <- initialiseEnvironment args
    runSafeT $ runReaderT runMain env

parserInfo :: OP.ParserInfo ProgOpt
parserInfo = OP.info (pure (.) <*> versionInfoOpt <*> OP.helper <*> argParser)
    (OP.progDescDoc (Just programHelpDoc))

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseCallingMode <*> parseKeepIncongruentReads <*> parseSeed <*> parseMinDepth <*>
    parseTransitionsMode <*> parseSnpFile <*> parseFormat <*> parseSampleNames
  where
    parseCallingMode = parseRandomCalling <|> parseMajorityCalling <|> parseRandomDiploidCalling
    parseRandomCalling = OP.flag' RandomCalling (OP.long "randomHaploid" <>
        OP.help "This method samples one read at random at each site, and uses the allele \
        \on that read as the one for the actual genotype. This results in a haploid \
        \call")
    parseRandomDiploidCalling = OP.flag' RandomDiploidCalling (OP.long "randomDiploid" <>
        OP.help "Sample two reads at random (without replacement) at each site and represent the \
        \individual by a diploid genotype constructed from those two random \
        \picks. This will always assign missing data to positions where only \
        \one read is present, even if minDepth=1. The main use case for this \
        \option is for estimating mean heterozygosity across sites.")
    parseMajorityCalling = MajorityCalling <$> (parseMajorityCallingFlag *> parseDownsamplingFlag)
    parseMajorityCallingFlag = OP.flag' True (OP.long "majorityCall" <> OP.help
        "Pick the allele supported by the \
        \most reads at a site. If an equal numbers of alleles fulfil this, pick one at \
        \random. This results in a haploid call. See --downSampling for best practices \
        \for calling rare variants")
    parseDownsamplingFlag = OP.switch (OP.long "downSampling" <> OP.help "When this switch is given, \
        \the MajorityCalling mode will downsample \
        \from the total number of reads a number of reads \
        \(without replacement) equal to the --minDepth given. This mitigates \
        \reference bias in the MajorityCalling model, which increases with higher coverage. \
        \The recommendation for rare-allele calling is --majorityCall --downsampling --minDepth 3")
    parseKeepIncongruentReads = OP.switch (OP.long "keepIncongruentReads" <> OP.help "By default, \
        \pileupCaller now removes reads with tri-allelic alleles that are neither of the two alleles specified in the SNP file. \
        \To keep those reads for sampling, set this flag. With this option given, if \
        \the sampled read has a tri-allelic allele that is neither of the two given alleles in the SNP file, a missing genotype is generated. \
        \IMPORTANT NOTE: The default behaviour has changed in pileupCaller version 1.4.0. If you want to emulate the previous \
        \behaviour, use this flag. I recommend now to NOT set this flag and use the new behaviour.")
    parseSeed = OP.option (Just <$> OP.auto) (OP.long "seed" <>
        OP.value Nothing <> OP.metavar "<RANDOM_SEED>" <>
        OP.help "random seed used for the random number generator. If not given, use \
        \system clock to seed the random number generator.")
    parseMinDepth = OP.option OP.auto (OP.long "minDepth" <> OP.short 'd' <>
        OP.value 1 <> OP.showDefault <> OP.metavar "<DEPTH>" <>
        OP.help "specify the minimum depth for a call. For sites with fewer \
            \reads than this number, declare Missing")
    parseTransitionsMode = parseSkipTransitions <|> parseTransitionsMissing <|> parseSingleStrandMode <|> pure AllSites
    parseSkipTransitions = OP.flag' SkipTransitions (OP.long "skipTransitions" <>
        OP.help "skip transition SNPs entirely in the output, resulting in a dataset with fewer sites.")
    parseTransitionsMissing = OP.flag' TransitionsMissing (OP.long "transitionsMissing" <>
        OP.help "mark transitions as missing in the output, but do output the sites.")
    parseSingleStrandMode = OP.flag' SingleStrandMode (OP.long "singleStrandMode" <>
        OP.help "At C/T polymorphisms, ignore reads aligning to the forward strand. \
        \At G/A polymorphisms, ignore reads aligning to the reverse strand. This should \
        \remove post-mortem damage in ancient DNA libraries prepared with the non-UDG single-stranded protocol.")
    parseSnpFile = OP.strOption (OP.long "snpFile" <> OP.short 'f' <>
        OP.metavar "<FILE>" <> OP.help "an Eigenstrat-formatted SNP list file for \
        \the positions and alleles to call. All \
        \positions in the SNP file will be output, adding missing data where \
        \there is no data. Note that pileupCaller automatically checks whether \
        \alleles in the SNP file are flipped with respect to the human \
        \reference, and in those cases flips the genotypes accordingly. \
        \But it assumes that the strand-orientation of the SNPs given in the SNP list is the one \
        \in the reference genome used in the BAM file underlying the pileup input. \
        \Note that both the SNP file and the incoming pileup data have to be \
        \ordered by chromosome and position, and this is checked. The chromosome order in humans is 1-22,X,Y,MT. \
        \Chromosome can generally begin with \"chr\". In case of non-human data with different chromosome \
        \names, you should convert all names to numbers. They will always considered to \
        \be numerically ordered, even beyond 22. Finally, I note that for internally, \
        \X is converted to 23, Y to 24 and MT to 90. This is the most widely used encoding in Eigenstrat \
        \databases for human data, so using a SNP file with that encoding will automatically be correctly aligned \
        \to pileup data with actual chromosome names X, Y and MT (or chrX, chrY and chrMT, respectively).")
    parseFormat = (EigenstratFormat <$> parseEigenstratPrefix <*> parseSamplePopName) <|> pure FreqSumFormat
    parseEigenstratPrefix = OP.strOption (OP.long "eigenstratOut" <> OP.short 'e' <>
        OP.metavar "<FILE_PREFIX>" <>
        OP.help "Set Eigenstrat as output format. Specify the filenames for the EigenStrat \
        \SNP and IND file outputs: <FILE_PREFIX>.snp.txt and <FILE_PREFIX>.ind.txt \
        \If not set, output will be FreqSum (Default). Note that freqSum format, described at \
        \https://rarecoal-docs.readthedocs.io/en/latest/rarecoal-tools.html#vcf2freqsum, \
        \is useful for testing your pipeline, since it's output to standard out")
    parseSampleNames = parseSampleNameList <|> parseSampleNameFile
    parseSampleNameList = OP.option (Left . splitOn "," <$> OP.str)
        (OP.long "sampleNames" <> OP.metavar "NAME1,NAME2,..." <>
        OP.help "give the names of the samples as comma-separated list (no spaces)")
    parseSampleNameFile = OP.option (Right <$> OP.str) (OP.long "sampleNameFile" <> OP.metavar "<FILE>" <>
        OP.help "give the names of the samples in a file with one name per \
        \line")
    parseSamplePopName = OP.strOption (OP.long "samplePopName" <> OP.value "Unknown" <> OP.showDefault <>
        OP.metavar "POP" <>
        OP.help "specify the population name of the samples, which is included\
        \ in the output *.ind.txt file in Eigenstrat output. This will be ignored if the output \
        \format is not Eigenstrat")

programHelpDoc :: PP.Doc
programHelpDoc =
    part1 PP.<$>
    PP.enclose PP.line PP.line (PP.indent 4 samtoolsExample) PP.<$>
    part2 PP.<$> PP.line PP.<$>
    PP.string versionInfoText
  where
    part1 = PP.fillSep . map PP.text . words $
        "PileupCaller is a tool to create genotype calls from bam files using read-sampling methods. \
        \To use this tool, you need to convert bam files into the mpileup-format, specified at \
        \http://www.htslib.org/doc/samtools.html (under \"mpileup\"). The recommended command line \
        \to create a multi-sample mpileup file to be processed with pileupCaller is"
    samtoolsExample = PP.hang 4 . PP.fillSep . map PP.text . words $
        "samtools mpileup -B -q30 -Q30 -l <BED_FILE> -R -f <FASTA_REFERENCE_FILE> \
        \Sample1.bam Sample2.bam Sample3.bam | pileupCaller ..."
    part2 = PP.fillSep . map PP.text . words $
        "You can lookup what these options do in the samtools documentation. \
        \Note that flag -B in samtools is very important to reduce reference \
        \bias in low coverage data."

initialiseEnvironment :: ProgOpt -> IO Env
initialiseEnvironment args = do
    let (ProgOpt callingMode keepInCongruentReads seed minDepth
            transitionsMode snpFile outFormat sampleNames) = args
    case seed of
        Nothing -> return ()
        Just seed_ -> liftIO . setStdGen $ mkStdGen seed_
    sampleNamesList <- case sampleNames of
        Left list -> return list
        Right fn -> lines <$> readFile fn
    let n = length sampleNamesList
    readStats <- ReadStats <$> newIORef 0 <*> makeVec n <*> makeVec n <*> makeVec n <*> makeVec n
    return $ Env callingMode keepInCongruentReads minDepth transitionsMode
        outFormat snpFile sampleNamesList readStats
  where
    makeVec n = do
        v <- V.new n 
        V.set v 0
        return v

runMain :: App ()
runMain = do
    let pileupProducer = readPileupFromStdIn
    snpFile <- asks envSnpFile
    freqSumProducer <- pileupToFreqSum snpFile pileupProducer
    outFormat <- asks envOutFormat
    case outFormat of
        FreqSumFormat -> outputFreqSum freqSumProducer
        EigenstratFormat outPrefix popName -> outputEigenStrat outPrefix popName freqSumProducer
    outputStats

pileupToFreqSum :: FilePath -> Producer PileupRow (SafeT IO) () ->
    App (Producer FreqSumEntry (SafeT IO) ())
pileupToFreqSum snpFileName pileupProducer = do
    nrSamples <- length <$> asks envSampleNames
    let snpProdOrderChecked =
            readEigenstratSnpFile snpFileName >-> orderCheckPipe cmpSnpPos
        pileupProdOrderChecked =
            pileupProducer >-> orderCheckPipe cmpPileupPos
        jointProd =
            orderedZip cmpSnpToPileupPos snpProdOrderChecked pileupProdOrderChecked
    mode <- asks envCallingMode
    keepInCongruentReads <- asks envKeepInCongruentReads
    transisitionsMode <- asks envTransitionsMode
    let singleStrandMode = (transisitionsMode == SingleStrandMode)
    minDepth <- asks envMinDepth
    readStats <- asks envStats
    let ret = Pipes.for jointProd $ \jointEntry ->
            case jointEntry of
                (Just esEntry, Nothing) -> do
                    let (EigenstratSnpEntry chr pos _ id_ ref alt) = esEntry
                        dosages = (replicate nrSamples Nothing) 
                    liftIO $ addOneSite readStats
                    yield $ FreqSumEntry chr pos (Just . B.unpack $ id_) ref alt dosages
                (Just esEntry, Just pRow) -> do
                    let (EigenstratSnpEntry chr pos _ id_ ref alt) = esEntry
                        (PileupRow _ _ _ rawPileupBasesPerSample rawStrandInfoPerSample) = pRow
                    let cleanBasesPerSample =
                            if   singleStrandMode
                            then cleanSSdamageAllSamples ref alt rawPileupBasesPerSample rawStrandInfoPerSample
                            else rawPileupBasesPerSample
                    let congruentBasesPerSample = 
                            if   keepInCongruentReads
                            then cleanBasesPerSample
                            else map (filter (\c -> c == ref || c == alt)) cleanBasesPerSample
                    liftIO $ addOneSite readStats
                    liftIO $ updateStatsAllSamples readStats (map length rawPileupBasesPerSample)
                        (map length cleanBasesPerSample) (map length congruentBasesPerSample)
                    calls <- liftIO $ mapM (callGenotypeFromPileup mode minDepth) congruentBasesPerSample
                    let genotypes = map (callToDosage ref alt) calls
                    yield (FreqSumEntry chr pos (Just . B.unpack $ id_) ref alt genotypes)
                _ -> return ()
    return (fst <$> ret)
  where
    cmpSnpPos :: EigenstratSnpEntry -> EigenstratSnpEntry -> Ordering
    cmpSnpPos es1 es2 = (snpChrom es1, snpPos es1) `compare` (snpChrom es2, snpPos es2)
    cmpPileupPos :: PileupRow -> PileupRow -> Ordering
    cmpPileupPos pr1 pr2 = (pileupChrom pr1, pileupPos pr1) `compare` (pileupChrom pr2, pileupPos pr2)
    cmpSnpToPileupPos :: EigenstratSnpEntry -> PileupRow -> Ordering
    cmpSnpToPileupPos es pr = (snpChrom es, snpPos es) `compare` (pileupChrom pr, pileupPos pr)

addOneSite :: ReadStats -> IO ()
addOneSite readStats = modifyIORef' (rsTotalSites readStats) (+1)

updateStatsAllSamples :: ReadStats -> [Int] -> [Int] -> [Int] -> IO ()
updateStatsAllSamples readStats rawBaseCounts damageCleanedBaseCounts congruencyCleanedBaseCounts = do
    sequence_ [V.modify (rsRawReads readStats) (+n) i | (i, n) <- zip [0..] rawBaseCounts]
    sequence_ [V.modify (rsReadsCleanedSS readStats) (+n) i | (i, n) <- zip [0..] damageCleanedBaseCounts]
    sequence_ [V.modify (rsReadsCongruent readStats) (+n) i | (i, n) <- zip [0..] congruencyCleanedBaseCounts]
    let nonMissingSites = [if n > 0 then 1 else 0 | n <- congruencyCleanedBaseCounts]
    sequence_ [V.modify (rsNonMissingSites readStats) (+n) i | (i, n) <- zip [0..] nonMissingSites]

outputFreqSum :: Producer FreqSumEntry (SafeT IO) () -> App ()
outputFreqSum freqSumProducer = do
    transitionsOnly <- asks envTransitionsMode
    callingMode <- asks envCallingMode
    sampleNames <- asks envSampleNames
    let nrHaplotypes = case callingMode of
            MajorityCalling _ -> 1 :: Int
            RandomCalling -> 1
            RandomDiploidCalling -> 2
    let header' = FreqSumHeader (map B.pack sampleNames) [nrHaplotypes | _ <- sampleNames]
        outProd = freqSumProducer >-> filterTransitions transitionsOnly
    lift . runEffect $ outProd >-> printFreqSumStdOut header'

outputEigenStrat :: FilePath -> String -> Producer FreqSumEntry (SafeT IO) () -> App ()
outputEigenStrat outPrefix popName freqSumProducer = do
    transitionsMode <- asks envTransitionsMode
    sampleNames <- asks envSampleNames
    callingMode <- asks envCallingMode
    let diploidizeCall = case callingMode of
            RandomCalling -> True
            MajorityCalling _ -> True
            RandomDiploidCalling -> False
    let snpOut = outPrefix <> "snp.txt"
        indOut = outPrefix <> "ind.txt"
        genoOut = outPrefix <> "geno.txt"
    let indEntries = [EigenstratIndEntry n Unknown popName | n <- sampleNames]
    lift . runEffect $ freqSumProducer >-> filterTransitions transitionsMode >->
                P.map (freqSumToEigenstrat diploidizeCall) >->
                writeEigenstrat genoOut snpOut indOut indEntries

outputStats :: App ()
outputStats = do
    ReadStats totalSites nonMissingSitesVec rawReadsVec damageCleanedReadsVec congruentReadsVec <- asks envStats
    sampleNames <- asks envSampleNames
    liftIO $ hPutStrLn stderr
        "# Summary Statistics per sample \n\
        \# SampleName: Name of the sample as given by the user \n\
        \# TotalSites: Total number of sites in the given Snp file (before transition filtering) \n\
        \# NonMissingCalls: Total number of sites output with a non-Missing call (before transition filtering) \n\
        \# avgRawReads: mean coverage of raw pileup input data across total sites (incl. missing sites) \n\
        \# avgDamageCleanedReads: mean coverage of pileup after single-stranded damage removal \n\
        \# avgSampledFrom: mean coverage of pileup after removing reads with tri-allelic alleles \n\
        \SampleName\tTotalSites\tNonMissingCalls\tavgRawReads\tavgDamageCleanedReads\tavgSampledFrom"
    forM_ (zip [0..] sampleNames) $ \(i, name) -> do
        totalS <- readIORef totalSites
        nonMissingSites <- V.read nonMissingSitesVec i
        rawReads <- V.read rawReadsVec i
        damageCleanedReads <- V.read damageCleanedReadsVec i
        congruentReads <- V.read congruentReadsVec i
        let avgRawReads = (fromIntegral rawReads / fromIntegral totalS) :: Double
            avgDamageCleanedReads = (fromIntegral damageCleanedReads / fromIntegral totalS) :: Double
            avgCongruentReads = (fromIntegral congruentReads / fromIntegral totalS) :: Double
        liftIO . hPutStrLn stderr $ printf "%s\t%d\t%d\t%g\t%g\t%g" name totalS nonMissingSites
            avgRawReads avgDamageCleanedReads avgCongruentReads

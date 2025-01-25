{-# LANGUAGE OverloadedStrings #-}

import           SequenceFormats.Eigenstrat   (EigenstratIndEntry (..),
                                               EigenstratSnpEntry (..),
                                               Sex (..), readEigenstratSnpFile,
                                               writeEigenstrat)
import           SequenceFormats.FreqSum      (FreqSumEntry (..),
                                               FreqSumHeader (..),
                                               printFreqSumStdOut)
import           SequenceFormats.Pileup       (PileupRow (..),
                                               readPileupFromStdIn)
import           SequenceFormats.Plink        (PlinkPopNameMode (..),
                                               eigenstratInd2PlinkFam,
                                               writePlink)
import           SequenceFormats.Utils        (Chrom (..),
                                               SeqFormatException (..))
import           SequenceFormats.VCF          (VCFentry (..), VCFheader(..), printVCFtoStdOut)
import           SequenceTools.PileupCaller   (CallingMode (..),
                                               TransitionsMode (..),
                                               callGenotypeFromPileup,
                                               callToDosage,
                                               cleanSSdamageAllSamples,
                                               filterTransitions,
                                               computeAlleleFreq)
import           SequenceTools.Utils          (UserInputException (..),
                                               freqSumToEigenstrat,
                                               versionInfoOpt, versionInfoText)

import           Control.Applicative          ((<|>))
import           Control.Exception            (catch, throwIO)
import           Control.Monad                (forM_, when)
import           Control.Monad.IO.Class       (liftIO)
import           Control.Monad.Trans.Class    (lift)
import           Control.Monad.Trans.Reader   (ReaderT, ask, asks, runReaderT)
import qualified Data.ByteString.Char8        as B
import           Data.IORef                   (IORef, modifyIORef', newIORef,
                                               readIORef)
import           Data.List                    (intercalate)
import           Data.List.Split              (splitOn)
import qualified Data.Vector.Unboxed.Mutable  as V
import           Data.Version                 (Version, showVersion)
import qualified Options.Applicative          as OP
import           Paths_sequenceTools          (version)
import           Pipes                        (Producer, for, runEffect, yield,
                                               (>->))
import           Pipes.OrderedZip             (orderCheckPipe, orderedZip)
import qualified Pipes.Prelude                as P
import           Pipes.Safe                   (SafeT, runSafeT)
import           System.Environment           (getArgs, getProgName)
import           System.IO                    (hPutStrLn, stderr)
import           System.Random                (mkStdGen, setStdGen)
import qualified Text.PrettyPrint.ANSI.Leijen as PP
import           Text.Printf                  (printf)

data OutFormat = EigenstratFormat FilePath Bool
               | PlinkFormat FilePath PlinkPopNameMode Bool
               | VCFformat
               | FreqSumFormat deriving (Show)

data ProgOpt = ProgOpt
    CallingMode                 --optCallingMode
    Bool                        --optKeepIncongruentReads
    (Maybe Int)                 --optSeed
    Int                         --optMinDepth
    TransitionsMode             --optTransitionsMode
    FilePath                    --optSnpFile
    OutFormat                   --optOutFormat
    (Either [String] FilePath)  --optSampleNames
    (Either String [String])    --optPopName

data ReadStats = ReadStats {
    rsTotalSites      :: IORef Int,
    rsNonMissingSites :: V.IOVector Int,
    rsRawReads        :: V.IOVector Int,
    rsReadsCleanedSS  :: V.IOVector Int,
    rsReadsCongruent  :: V.IOVector Int
}

data Env = Env {
    envCallingMode          :: CallingMode,
    envKeepInCongruentReads :: Bool,
    envMinDepth             :: Int,
    envTransitionsMode      :: TransitionsMode,
    envOutFormat            :: OutFormat,
    envSnpFile              :: FilePath,
    envSampleNames          :: [String],
    envPopName              :: Either String [String],
    envVersion              :: Version,
    envStats                :: ReadStats
}

instance Show Env where
    show e = show (
        envCallingMode e,
        envKeepInCongruentReads e,
        envMinDepth e,
        envTransitionsMode e,
        envOutFormat e,
        envSnpFile e,
        envSampleNames e,
        envPopName e,
        envVersion e)

type App = ReaderT Env (SafeT IO)

main :: IO ()
main = do
    args <- OP.execParser parserInfo
    env <- initialiseEnvironment args
    let handler = \(SeqFormatException msg) -> do
            throwIO $ SeqFormatException (take 200 msg)
    catch (runSafeT $ runReaderT runMain env) handler

parserInfo :: OP.ParserInfo ProgOpt
parserInfo = OP.info (pure (.) <*> versionInfoOpt <*> OP.helper <*> argParser)
    (OP.progDescDoc (Just programHelpDoc))

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseCallingMode
                    <*> parseKeepIncongruentReads
                    <*> parseSeed
                    <*> parseMinDepth
                    <*> parseTransitionsMode
                    <*> parseSnpFile
                    <*> parseFormat
                    <*> parseSampleNames
                    <*> parsePopName
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
        OP.help "[THIS IS CURRENTLY AN EXPERIMENTAL FEATURE]. At C/T polymorphisms, ignore reads aligning to the forward strand. \
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
    parseFormat = (EigenstratFormat <$> parseEigenstratPrefix <*> parseZipOut) <|>
        (PlinkFormat <$> parsePlinkPrefix <*> parsePlinkPopMode <*> parseZipOut) <|>
        parseVCFformat <|> pure FreqSumFormat
    parseZipOut = OP.switch (OP.long "--zip" <> OP.short 'z' <> OP.help "GZip the output genotype files. Filenames will be appended with '.gz'.")
    parseEigenstratPrefix = OP.strOption (OP.long "eigenstratOut" <> OP.short 'e' <>
        OP.metavar "<FILE_PREFIX>" <>
        OP.help "Set Eigenstrat as output format. Specify the filenames for the EigenStrat \
        \SNP, IND and GENO file outputs: <FILE_PREFIX>.snp, <FILE_PREFIX>.ind and <FILE_PREFIX>.geno. \
        \If not set, output will be FreqSum (Default). Note that freqSum format, described at \
        \https://rarecoal-docs.readthedocs.io/en/latest/rarecoal-tools.html#vcf2freqsum, \
        \is useful for testing your pipeline, since it's output to standard out")
    parsePlinkPrefix = OP.strOption (OP.long "plinkOut" <> OP.short 'p' <>
        OP.metavar "<FILE_PREFIX>" <>
        OP.help "Set Plink as output format. Specify the filenames for the Plink \
        \BIM, FAM and BED file outputs: <FILE_PREFIX>.bim, <FILE_PREFIX>.fam and <FILE_PREFIX>.bed. \
        \If not set, output will be FreqSum (Default). Note that freqSum format, described at \
        \https://rarecoal-docs.readthedocs.io/en/latest/rarecoal-tools.html#vcf2freqsum, \
        \is useful for testing your pipeline, since it's output to standard out")
    parsePlinkPopMode = parsePlinkPopPhenotype <|> parsePlinkPopBoth <|> pure PlinkPopNameAsFamily
    parsePlinkPopPhenotype = OP.flag' PlinkPopNameAsPhenotype (OP.long "popNameAsPhenotype" <> OP.help "Only valid for Plink Output: \
        \Write the population name into the last column of the fam file, as a Phenotype according to the Plink Spec. \
        \By default, the population name is specified as the first column only (family name in the Plink spec)")
    parsePlinkPopBoth = OP.flag' PlinkPopNameAsBoth (OP.long "popNameAsBoth" <> OP.help "Only valid for Plink Output: \
        \Write the population name into both the first and last column of the fam file, so both as Family-ID and as a \
        \Phenotype according to the Plink Spec. By default, the population name is specified only as the first column (family name in the Plink spec)")
    parseVCFformat = OP.flag' VCFformat (OP.long "vcf" <> OP.help "output VCF format to stdout")
    parseSampleNames = parseSampleNameList <|> parseSampleNameFile
    parseSampleNameList = OP.option (Left . splitOn "," <$> OP.str)
        (OP.long "sampleNames" <> OP.metavar "NAME1,NAME2,..." <>
        OP.help "give the names of the samples as comma-separated list (no spaces)")
    parseSampleNameFile = OP.option (Right <$> OP.str) (OP.long "sampleNameFile" <> OP.metavar "<FILE>" <>
        OP.help "give the names of the samples in a file with one name per \
        \line")
    parsePopName = OP.option (processPopNames . splitOn "," <$> OP.str) (OP.long "samplePopName" <> OP.value (Left "Unknown") <> OP.showDefault <>
        OP.metavar "POP(s)" <>
        OP.help "specify the population name(s) of the samples, which are included\
        \ in the output *.ind.txt file in Eigenstrat output. This will be ignored if the output \
        \format is not Eigenstrat. If a single name is given, it is applied to all samples, if multiple are given, their number must match the \
        \the number of samples")
    processPopNames names = if length names == 1 then Left (head names) else Right names

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
            transitionsMode snpFile outFormat sampleNames popName) = args
    case seed of
        Nothing    -> return ()
        Just seed_ -> liftIO . setStdGen $ mkStdGen seed_
    sampleNamesList <- case sampleNames of
        Left list -> return list
        Right fn  -> lines <$> readFile fn
    let n = length sampleNamesList
    readStats <- ReadStats <$> newIORef 0 <*> makeVec n <*> makeVec n <*> makeVec n <*> makeVec n
    return $ Env callingMode keepInCongruentReads minDepth transitionsMode
        outFormat snpFile sampleNamesList popName version readStats
  where
    makeVec n = do
        v <- V.new n
        V.set v 0
        return v

runMain :: App ()
runMain = do
    let pileupProducer = readPileupFromStdIn
    outFormat <- asks envOutFormat
    popNameSpec <- asks envPopName
    n <- length <$> asks envSampleNames
    let popNames = case popNameSpec of
            Left singlePopName -> replicate n singlePopName
            Right p -> if length p /= n then error "number of specified populations must equal sample size" else p
    case outFormat of
        FreqSumFormat -> do
            freqSumProducer <- pileupToFreqSum pileupProducer
            outputFreqSum (freqSumProducer >-> P.map snd)
        EigenstratFormat outPrefix zipOut -> do
            freqSumProducer <- pileupToFreqSum pileupProducer
            outputEigenStratOrPlink outPrefix zipOut popNames Nothing (freqSumProducer >-> P.map snd)
        PlinkFormat outPrefix popNameMode zipOut -> do
            freqSumProducer <- pileupToFreqSum pileupProducer
            outputEigenStratOrPlink outPrefix zipOut popNames (Just popNameMode) (freqSumProducer >-> P.map snd)
        VCFformat -> do
            freqSumProducer <- pileupToFreqSum pileupProducer
            outputVCF popNames freqSumProducer
    outputStats

pileupToFreqSum :: Producer PileupRow (SafeT IO) () -> App (Producer (Maybe PileupRow, FreqSumEntry) (SafeT IO) ())
pileupToFreqSum pileupProducer = do
    snpFileName <- asks envSnpFile
    nrSamples <- length <$> asks envSampleNames
    let snpProdOrderChecked =
            readEigenstratSnpFile snpFileName >-> orderCheckPipe cmpSnpPos
        pileupProdOrderChecked =
            pileupProducer >-> orderCheckPipe cmpPileupPos
        jointProd =
            orderedZip cmpSnpToPileupPos snpProdOrderChecked pileupProdOrderChecked
    mode <- asks envCallingMode
    keepInCongruentReads <- asks envKeepInCongruentReads
    transitionsMode <- asks envTransitionsMode
    let singleStrandMode = (transitionsMode == SingleStrandMode)
    minDepth <- asks envMinDepth
    readStats <- asks envStats
    let ret = Pipes.for jointProd $ \jointEntry ->
            case jointEntry of
                (Just esEntry, Nothing) -> do
                    let (EigenstratSnpEntry chr pos gpos id_ ref alt) = esEntry
                        dosages = (replicate nrSamples Nothing)
                    liftIO $ addOneSite readStats
                    yield $ (Nothing, FreqSumEntry chr pos (Just id_) (Just gpos) ref alt dosages)
                (Just esEntry, Just pRow) -> do
                    let (EigenstratSnpEntry chr pos gpos id_ ref alt) = esEntry
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
                    yield $ (Just pRow, FreqSumEntry chr pos (Just id_) (Just gpos) ref alt genotypes)
                _ -> return ()
    return $ (fst <$> ret) >-> filterTransitions transitionsMode
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
    let nSamples = V.length (rsRawReads readStats)
    when (length rawBaseCounts /= nSamples) . throwIO . UserInputException $
        "number of individuals specified (" ++ show nSamples ++
        ") differs from number of individuals in the pileup input (" ++
        show (length rawBaseCounts) ++ ")"
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
            MajorityCalling _    -> 1 :: Int
            RandomCalling        -> 1
            RandomDiploidCalling -> 2
    let header' = FreqSumHeader sampleNames [nrHaplotypes | _ <- sampleNames]
    lift . runEffect $ freqSumProducer >-> printFreqSumStdOut header'

outputEigenStratOrPlink :: FilePath -> Bool -> [String] -> Maybe PlinkPopNameMode -> Producer FreqSumEntry (SafeT IO) () -> App ()
outputEigenStratOrPlink outPrefix zipOut popNames maybePlinkPopMode freqSumProducer = do
    transitionsMode <- asks envTransitionsMode
    sampleNames <- asks envSampleNames
    let (snpOut, indOut, genoOut) = case (maybePlinkPopMode, zipOut) of
            (Just _, False)  -> (outPrefix <> ".bim", outPrefix <> ".fam", outPrefix <> ".bed")
            (Just _, True)  -> (outPrefix <> ".bim.gz", outPrefix <> ".fam", outPrefix <> ".bed.gz")
            (Nothing, False) -> (outPrefix <> ".snp", outPrefix <> ".ind", outPrefix <> ".geno")
            (Nothing, True) -> (outPrefix <> ".snp.gz", outPrefix <> ".ind", outPrefix <> ".geno.gz")
    let indEntries = [EigenstratIndEntry n Unknown p | (n, p) <- zip sampleNames popNames]
    let writeFunc = case maybePlinkPopMode of
            Nothing -> (\g s i -> writeEigenstrat g s i indEntries)
            Just popMode ->
                let famEntries = map (eigenstratInd2PlinkFam popMode) indEntries
                in  (\g s i -> writePlink g s i famEntries)
    lift . runEffect $ freqSumProducer >-> P.map freqSumToEigenstrat >-> writeFunc genoOut snpOut indOut

outputVCF :: [String] -> Producer (Maybe PileupRow, FreqSumEntry) (SafeT IO) () -> App ()
outputVCF popNames freqSumProd = do
    sampleNames <- map B.pack <$> asks envSampleNames
    version <- asks envVersion
    env <- ask
    prog_name <- liftIO $ getProgName
    prog_args <- liftIO $ getArgs
    let command_line = prog_name ++ " " ++ intercalate " " prog_args
    let metaInfoLines = map B.pack [
            "##fileformat=VCFv4.2",
            "##source=pileupCaller_v" ++ showVersion version,
            "##command_line=" ++ command_line,
            "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
            "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">",
            "##FILTER=<ID=s10,Description=\"Less than 10% of samples have data\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">",
            "##FORMAT=<ID=DP2,Number=2,Type=Integer,Description=\"Nr of Reads supporting each of the two alleles\">"]
        vcfh = VCFheader metaInfoLines sampleNames
    lift . runEffect $ freqSumProd >-> P.map (\(mpr, fse) -> createVcfEntry mpr fse) >-> printVCFtoStdOut vcfh

createVcfEntry :: Maybe PileupRow -> FreqSumEntry -> VCFentry
createVcfEntry maybePileupRow (FreqSumEntry chrom@(Chrom c) pos maybeSnpId maybeGeneticPos ref alt calls) =
    VCFentry chrom pos maybeSnpId (B.pack [ref]) [B.pack [alt]] Nothing (Just filterString) infoFields genotypeInfos
  where
    nrMissing = length . filter (==Nothing) $ calls
    nrSamples = length calls
    filterString =
        if nrMissing * 10 < nrSamples then "s10;s50"
        else if nrMissing * 2 < nrSamples then "s50"
        else "PASS"
    totalDepth = case maybePileupRow of
        Nothing -> 0
        Just pr -> sum . map length . pileupBases $ pr
    nrSamplesWithData = nrSamples - nrMissing
    alleleFreq = computeAlleleFreq calls
    infoFields = [
        B.pack $ "NS=" ++ show nrSamplesWithData,
        B.pack $ "DP=" ++ show totalDepth,
        B.pack $ "AF=" ++ show alleleFreq]
    formatField = case maybePileupRow of
        Nothing -> ["GT"]
        Just pr -> ["GT", "DP", "DP2"]
    genotypeFields = do -- list monad over samples
        i <- [0 .. (nrSamples - 1)]
        let c = calls !! i
        let gt = case c of
                Nothing     -> "."
                Just (0, 1) -> "0"
                Just (1, 1) -> "1"
                Just (0, 2) -> "0/0"
                Just (1, 2) -> "0/1"
                Just (2, 2) -> "1/1"
                _           -> error "should never happen"
        case maybePileupRow of
            Nothing -> return [gt]
            Just pr -> do
                let bases = pileupBases pr !! i
                let dp = length bases
                let dpRef = length . filter (==ref) $ bases
                let dpAlt = length . filter (==alt) $ bases
                return [gt, B.pack $ show dp, B.pack $ show dpRef ++ "," ++ show dpAlt]
    genotypeInfos = Just (formatField, genotypeFields)


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
        totalS <- liftIO $ readIORef totalSites
        nonMissingSites <- V.read nonMissingSitesVec i
        rawReads <- V.read rawReadsVec i
        damageCleanedReads <- V.read damageCleanedReadsVec i
        congruentReads <- V.read congruentReadsVec i
        let avgRawReads = (fromIntegral rawReads / fromIntegral totalS) :: Double
            avgDamageCleanedReads = (fromIntegral damageCleanedReads / fromIntegral totalS) :: Double
            avgCongruentReads = (fromIntegral congruentReads / fromIntegral totalS) :: Double
        liftIO . hPutStrLn stderr $ printf "%s\t%d\t%d\t%g\t%g\t%g" name totalS nonMissingSites
            avgRawReads avgDamageCleanedReads avgCongruentReads

{-# LANGUAGE OverloadedStrings #-}

import OrderedZip (orderedZip)

import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad.Random (evalRandIO)
import Control.Monad.Trans.Class (lift)
import Control.Monad.Trans.Reader (ReaderT, runReaderT, asks)
import Control.Monad.Trans.State.Strict (runStateT)
import qualified Data.Attoparsec.Text as A
import Data.Char (isSpace)
import Data.List (sortBy, intercalate)
import Data.Maybe (catMaybes)
import qualified Data.Text as T
import qualified Data.Text.IO as T
-- import Debug.Trace (trace)
import qualified Options.Applicative as OP
import Pipes (Pipe, yield, (>->), runEffect, Producer, Pipe, for, cat, next)
import Pipes.Attoparsec (parsed, parse)
import Pipes.Cliff (CreateProcess(..), CmdSpec(..), pipeOutput, NonPipe(..))
import Pipes.Cliff.Core (defaultHandler)
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, SafeT)
import Pipes.Safe.Prelude (withFile)
import Pipes.Text.Encoding (decodeUtf8)
import qualified Pipes.Text.IO as PT
import Prelude hiding (FilePath)
import System.IO (IOMode(..), hPutStrLn)
import System.Random (randomRIO, mkStdGen, setStdGen)
import System.Random.Shuffle (shuffleM)
import Turtle hiding (tab, cat, stderr)

data ProgOpt = ProgOpt {
    optTransversionsOnly :: Bool,
    optFillHomRef :: Bool,
    optSnpFile :: FilePath,
    optChrom :: Text,
    optOutChrom :: Maybe Text,
    optOutPrefix :: FilePath,
    optBcftoolsExe :: FilePath,
    optVcfFile :: FilePath
}

data VCFentry = VCFentry Text Int [Char] [[Int]] deriving (Show)
data VCFheader = VCFheader [Text] [String] deriving (Show)-- simple comment lines, sample names
data SnpEntry = SnpEntry Text Int Char Char deriving (Show)-- Chrom Pos Ref Alt

main :: IO ()
main = OP.execParser parser >>= runSafeT . runReaderT runWithOpts
  where
    parser = OP.info (OP.helper <*> argParser)
                     (OP.progDesc "A program to convert a VCF file to Eigenstrat")

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseTransversionsOnly <*> parseFillHomRef <*> parseSnpFile <*>  
                        parseChrom <*> parseOutChrom <*> parseOutPrefix <*> 
                        parseBcftoolsExe <*> parseVcfFileName
  where
    parseTransversionsOnly = OP.switch (OP.long "transversionsOnly" <> OP.short 't' <>
                    OP.help "Remove transition SNPs from the output)")
    parseFillHomRef = OP.switch (OP.long "fillHomRef" <> OP.short 'h' <>
                    OP.help "At missing sites in the VCF, output hom-ref genotypes instead of \
                             \missing data")
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
                             \the snp files in Eigenstrat because of samtools specifications.")
    parseChrom = OP.option (T.pack <$> OP.str) (OP.long "chrom" <> OP.short 'c' <>
                            OP.metavar "<CHROM>" <> OP.help "specify the region in the BAM file to \
                            \call from. Can be just the chromosome, or a string of the form \
                            \CHROM:START-END")
    parseOutChrom = OP.option (Just . T.pack <$> OP.str) (OP.long "outChrom" <>
                                    OP.metavar "<CHROM>" <>
                                   OP.help "specify the output chromosome name" <> OP.value Nothing)
    parseVcfFileName = OP.argument (fromText . T.pack <$> OP.str) (OP.metavar "<VCF_FILE>" <>
                                    OP.help "input VCF file")
    parseOutPrefix = OP.option (Just . fromText . T.pack <$> OP.str)
                                 (OP.long "outPrefix" <> OP.short 'e' <>
                                 OP.value Nothing <> OP.metavar "<FILE_PREFIX>" <>
                                 OP.help "specify the filenames for the EigenStrat SNP and IND \
                                 \file outputs: <FILE_PREFIX>.snp.txt and <FILE_PREFIX>.ind.txt \
                                 \Ignored if Output format is not Eigenstrat")
    parseBcftoolsExe = OP.option (fromText . T.pack <$> OP.str) (OP.long "bcftools" <>
                                  OP.value "bcftools" <> OP.showDefault <>
                                  OP.metavar "<BCFTOOLS_PATH>" <> OP.help "path to the \
                                  \bcftools executable, version >= 1.2")
    
-- runWithOpts :: ReaderT ProgOpt (SafeT IO) ()
-- runWithOpts = do
--     seed <- asks optSeed
--     region <- asks optRegion
--     outChrom <- asks optOutChrom
--     outFormat <- asks optOutFormat
--     eigenStratOutPrefix <- asks optEigenstratOutPrefix
--     transversionsOnly <- asks optTransversionsOnly
--     case seed of
--         Nothing -> return ()
--         Just seed_ -> liftIO . setStdGen $ mkStdGen seed_
--     let chrom = head (T.splitOn ":" region)
--     let outputChrom = case outChrom of
--             Nothing -> chrom
--             Just label -> label
--     (vcfHeader, freqSumProducer) <- runPileup
--     case outFormat of
--         FreqSumFormat -> do
--             let VCFheader _ n = vcfHeader
--             liftIO . putStrLn $ "#CHROM\tPOS\tREF\tALT\t" ++ (intercalate "\t" . map (++"(2)") $ n)
--             lift . runEffect $ freqSumProducer >-> filterTransitions transversionsOnly >->
--                                P.map (showFreqSum outputChrom) >-> printToStdOut
--         EigenStrat -> case eigenStratOutPrefix of
--             Nothing -> liftIO . throwIO $ AssertionFailed "need an eigenstratPrefix for \
--                                                            \EigenStratFormat"
--             Just fn -> do
--                 let snpOut = fn <.> "snp.txt"
--                     indOut = fn <.> "ind.txt"
--                 lift . withFile (T.unpack . format fp $ indOut) WriteMode $ \indOutHandle -> do
--                     let VCFheader _ sampleNames = vcfHeader
--                     mapM_ (\n -> liftIO $ hPutStrLn indOutHandle (n ++ "\tU\tUnknown")) sampleNames
--                 lift . withFile (T.unpack . format fp $ snpOut) WriteMode $ \snpOutHandle -> do
--                     runEffect $ freqSumProducer >-> filterTransitions transversionsOnly >->
--                                 printEigenStrat outputChrom snpOutHandle >-> printToStdOut
--   where
--     printToStdOut = for cat (liftIO . T.putStrLn)
--     filterTransitions transversionsOnly =
--         if transversionsOnly
--         then P.filter (\(FreqSumRow _ _ ref alt _) -> isTransversion ref alt)
--         else cat
--     isTransversion ref alt = not $ isTransition ref alt
--     isTransition ref alt = ((ref == 'A') && (alt == 'G')) || ((ref == 'G') && (alt == 'A')) ||
--                            ((ref == 'C') && (alt == 'T')) || ((ref == 'T') && (alt == 'C'))
--
-- runPileup :: ReaderT ProgOpt (SafeT IO) (VCFheader, Producer FreqSumRow (SafeT IO) ())
-- runPileup = do
--     snpFile <- asks optSnpFile
--     case snpFile of
--         Nothing -> runPileupSimple
--         Just fn -> runPileupSnpFile fn
--
-- runPileupSimple :: ReaderT ProgOpt (SafeT IO) (VCFheader, Producer FreqSumRow (SafeT IO) ())
-- runPileupSimple = do
--     samtools <- asks optSamtoolsExe
--     bcftools <- asks optBcftoolsExe
--     reference <- asks optReference
--     region <- asks optRegion
--     bamFiles <- asks optBamFiles
--     mode <- asks optCallingMode
--     minDepth <- asks optMinDepth
--     let bams = (T.intercalate " " (map (format fp) bamFiles))
--     let cmd = format (fp%" mpileup -q30 -Q30 -C50 -I -f "%fp%" -g -t DPR -r "%s%" "%s%
--                       " | "%fp%" view -v snps") samtools reference region bams bcftools
--     vcfTextProd <- liftIO $ produceFromCommand cmd
--     (vcfHeader_, vcfProd) <- lift $ parseVCF vcfTextProd
--     let vcfProdPipe = vcfProd >-> processVcfSimple (length bamFiles) mode minDepth
--     return (vcfHeader_, vcfProdPipe)
--
-- runPileupSnpFile :: FilePath ->
--                     ReaderT ProgOpt (SafeT IO) (VCFheader, Producer FreqSumRow (SafeT IO) ())
-- runPileupSnpFile fn = do
--     samtools <- asks optSamtoolsExe
--     bcftools <- asks optBcftoolsExe
--     reference <- asks optReference
--     region <- asks optRegion
--     bamFiles <- asks optBamFiles
--     mode <- asks optCallingMode
--     minDepth <- asks optMinDepth
--     let bams = (T.intercalate " " (map (format fp) bamFiles))
--     let chrom = head (T.splitOn ":" region)
--     let cmd = format (fp%" mpileup -q30 -Q30 -C50 -I -f "%fp%" -g -t DPR -r "%s%" -l "%fp%
--                       " "%s%" | "%fp%" view") samtools reference region fn bams bcftools
--     vcfTextProd <- liftIO $ produceFromCommand cmd
--     let snpTextProd = PT.readFile ((T.unpack . format fp) fn)
--     (vcfHeader_, vcfProd) <- lift $ parseVCF vcfTextProd
--     let snpProd =
--             parsed snpParser snpTextProd >-> P.filter (\(SnpEntry c _ _ _) -> c == chrom)
--     let jointProd = orderedZip cmp snpProd vcfProd
--         jointProdPipe = jointProd >-> processVcfWithSnpFile (length bamFiles) mode minDepth
--     return (vcfHeader_, fmap snd jointProdPipe)
--   where
--     cmp (SnpEntry _ snpPos _ _) (VCFentry _ vcfPos _ _) = snpPos `compare` vcfPos
--
--
-- produceFromCommand :: Text -> IO (Producer Text (SafeT IO) ())
-- produceFromCommand cmd = do
--     let createProcess = CreateProcess (ShellCommand (T.unpack cmd)) Nothing Nothing False False
--                                                         False defaultHandler
--     (p, _) <- pipeOutput Inherit Inherit createProcess
--     return . void . decodeUtf8 $ p
--
-- parseVCF :: Producer Text (SafeT IO) () -> SafeT IO (VCFheader, Producer VCFentry (SafeT IO) ())
-- parseVCF prod = do
--     (res, rest) <- runStateT (parse vcfHeaderParser) prod
--     header <- case res of
--         Nothing -> liftIO . throwIO $ AssertionFailed "vcf header not readible. VCF file empty?"
--         Just (Left e_) -> do
--             Right (chunk, _) <- next rest
--             let msg = show e_ ++ T.unpack chunk
--             liftIO . throwIO $ AssertionFailed ("VCF header parsing error: " ++ msg)
--         Just (Right h) -> return h
--     return (header, parsed vcfParser rest >>= liftErrors)
--   where
--     liftErrors res = case res of
--         Left (e_, prod_) -> do
--             Right (chunk, _) <- lift $ next prod_
--             let msg = show e_ ++ T.unpack chunk
--             lift . liftIO . throwIO $ AssertionFailed msg
--         Right () -> return ()
--
-- vcfHeaderParser :: A.Parser VCFheader
-- vcfHeaderParser = VCFheader <$> A.many' doubleCommentLine <*> singleCommentLine
--   where
--     doubleCommentLine = do
--         c1 <- A.string "##"
--         s_ <- A.takeTill A.isEndOfLine <* A.endOfLine
--         return $ T.append c1 s_
--     singleCommentLine = do
--         void $ A.char '#'
--         s_ <- A.takeTill A.isEndOfLine <* A.endOfLine
--         let fields = T.splitOn "\t" s_
--         return . drop 9 . map T.unpack $ fields
--
-- vcfParser :: A.Parser VCFentry
-- vcfParser = do
--     chrom <- word
--     tab
--     pos <- A.decimal
--     tab >> word >> tab
--     ref <- A.satisfy (A.inClass "NACTG")
--     tab
--     alt <- (altAllele <|> xAllele) `A.sepBy1` (A.char ',')
--     tab
--     _ <- A.count 3 (word >> tab)
--     _ <- A.string "PL:DPR" >> tab
--     coverages <- coverage `A.sepBy1` tab
--     _ <- A.satisfy (\c -> c == '\r' || c == '\n')
--     -- trace (show (chrom, pos, ref, alt, coverages)) $ return ()
--     let filteredAlt = catMaybes alt
--     let filteredCoverages =
--             [[c | (c, a) <- zip cov (Just ref:alt), a /= Nothing] | cov <- coverages]
--     return $ VCFentry chrom pos (ref:filteredAlt) filteredCoverages
--   where
--     altAllele = Just <$> A.satisfy (A.inClass "ACTG")
--     xAllele = (A.string "<X>" <|> A.string "<*>") >> return Nothing
--     coverage = do
--         _ <- A.decimal `A.sepBy1` (A.char ',') :: A.Parser [Int]
--         _ <- A.char ':'
--         A.decimal `A.sepBy1` (A.char ',')
--
-- word :: A.Parser Text
-- word = T.pack <$> A.many1 (A.satisfy (not . isSpace))
--
-- tab :: A.Parser ()
-- tab = A.char '\t' >> return ()
--
-- processVcfSimple :: Int -> CallingMode -> Int -> Pipe VCFentry FreqSumRow (SafeT IO) r
-- processVcfSimple nrInds mode minDepth = for cat $ \vcfEntry -> do
--     let (VCFentry chrom pos alleles covNums) = vcfEntry
--     when (length covNums /= nrInds) $ (liftIO . throwIO) (AssertionFailed "inconsistent number \
--             \of genotypes. Check that bam files have different readgroup sample names")
--     (normalizedAlleles, normalizedCovNums) <- case alleles of
--         [_] -> liftIO . throwIO $ AssertionFailed "should not happen, need at least one \
--                             \alternative allele"
--         [ref, alt] -> return ([ref, alt], covNums)
--         _ -> do
--             let altNumPairs =
--                     [(alleles!!i, sum [c!!i | c <- covNums]) | i <- [1 .. (length alleles - 1)]]
--             shuffledAltNumPairs <- liftIO $ shuffle altNumPairs
--             let (alt, _) = head . sortBy (\a b -> snd b `compare` snd a) $ shuffledAltNumPairs
--             let altIndex = snd . head . filter ((==alt) . fst) $ zip alleles [0..]
--             when (altIndex == 0) $ (liftIO . throwIO)
--                             (AssertionFailed "should not happen, altIndex==0")
--             return ([head alleles, alt], [[c !! 0, c !! altIndex] | c <- covNums])
--     let [ref, alt] = normalizedAlleles
--     when (ref /= 'N') $ do
--         genotypes <- liftIO $ mapM (callGenotype mode minDepth) normalizedCovNums
--         when (any (>0) genotypes) $ yield (FreqSumRow chrom pos ref alt genotypes)
--   where
--     shuffle list = evalRandIO (shuffleM list)
--
-- callGenotype :: CallingMode -> Int -> [Int] -> IO Int
-- callGenotype mode minDepth covs = do
--     if sum covs < minDepth then return (-1) else do
--         case covs of
--             [_] -> return 0
--             [numRef, numAlt] -> do
--                 case mode of
--                     MajorityCalling -> case numRef `compare` numAlt of
--                         LT -> return 2
--                         GT -> return 0
--                         EQ -> do
--                             rn <- randomRIO (1, numRef + numAlt)
--                             if rn <= numRef then return 0 else return 2
--                     RandomCalling -> do
--                             rn <- randomRIO (1, numRef + numAlt)
--                             if rn <= numRef then return 0 else return 2
--                     RareCalling -> do
--                             if numAlt >= 2 then return 1 else return 0
--             _ -> throwIO (AssertionFailed "should not happen. CallGenotype called with more \
--                                             \than two alleles")
--
-- processVcfWithSnpFile :: Int -> CallingMode -> Int -> Pipe (Maybe SnpEntry, Maybe VCFentry)
--                                                       FreqSumRow (SafeT IO) r
-- processVcfWithSnpFile nrInds mode minDepth = for cat $ \jointEntry -> do
--     -- trace (show jointEntry) (return ())
--     case jointEntry of
--         (Just (SnpEntry snpChrom snpPos snpRef snpAlt), Nothing) -> do
--             yield $ FreqSumRow snpChrom snpPos snpRef snpAlt (replicate nrInds (-1))
--         (Just (SnpEntry snpChrom snpPos snpRef snpAlt),
--          Just (VCFentry _ _ vcfAlleles vcfNums)) -> do
--             when (length vcfNums /= nrInds) $ (liftIO . throwIO) (AssertionFailed "inconsistent \
--                             \number of genotypes. Check that bam files have different \
--                             \readgroup sample names")
--             let normalizedAlleleI =
--                     map snd . filter (\(a, _) -> a == snpRef || a == snpAlt) $ zip vcfAlleles [0..]
--                 normalizedVcfAlleles = map (vcfAlleles!!) normalizedAlleleI
--                 normalizedVcfNums = [map (v!!) normalizedAlleleI | v <- vcfNums]
--             -- trace (show (VCFentry vcfChrom vcfPos vcfAlleles vcfNums)) (return ())
--             genotypes <- case normalizedVcfAlleles of
--                     [] -> return (replicate nrInds (-1))
--                     [ref] -> if ref == snpRef
--                         then return [if sum c >= minDepth then 0 else (-1) | c <- normalizedVcfNums]
--                         else return [if sum c >= minDepth then 2 else (-1) | c <- normalizedVcfNums]
--                     [ref, alt] -> if [ref, alt] == [snpRef, snpAlt]
--                         then liftIO $ mapM (callGenotype mode minDepth) normalizedVcfNums
--                         else liftIO $ mapM (callGenotype mode minDepth)
--                                         (map reverse normalizedVcfNums)
--                     _ -> liftIO . throwIO $ AssertionFailed ("should not happen, can only have \
--                                     \two alleles after normalization: " ++ show jointEntry)
--             yield (FreqSumRow snpChrom snpPos snpRef snpAlt genotypes)
--         _ -> return ()
--
-- snpParser :: A.Parser SnpEntry
-- snpParser = do
--     chrom <- word
--     tab
--     pos <- A.decimal
--     tab
--     ref <- A.satisfy (A.inClass "ACTG")
--     _ <- A.char ','
--     alt <- A.satisfy (A.inClass "ACTG")
--     _ <- A.satisfy (\c -> c == '\r' || c == '\n')
--     let ret = SnpEntry chrom pos ref alt
--     -- trace (show ret) $ return ()
--     return ret
--
-- showFreqSum :: Text -> FreqSumRow -> Text
-- showFreqSum outChrom (FreqSumRow _ pos ref alt calls) =
--     format (s%"\t"%d%"\t"%s%"\t"%s%"\t"%s) outChrom pos (T.singleton ref) (T.singleton alt) callsStr
--   where
--     callsStr = (T.intercalate "\t" . map (format d)) calls
--
-- printEigenStrat :: Text -> Handle -> Pipe FreqSumRow Text (SafeT IO) r
-- printEigenStrat outChrom snpOutHandle = for cat $ \(FreqSumRow _ pos ref alt calls) -> do
--     let n = format (s%"_"%d) outChrom pos
--         snpLine = format (s%"\t"%s%"\t0\t"%d%"\t"%s%"\t"%s) n outChrom pos (T.singleton ref)
--                             (T.singleton alt)
--     liftIO . T.hPutStrLn snpOutHandle $ snpLine
--     yield . T.concat . map (format d . toEigenStratNum) $ calls
--   where
--     toEigenStratNum c = case c of
--         0 -> 2 :: Int
--         1 -> 1
--         2 -> 0
--         -1 -> 9
--         _ -> error ("unknown genotype " ++ show c)

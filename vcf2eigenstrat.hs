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
main = readOptions >>= runMain

readOptions :: IO ProgOpt
readOptions = OP.execParser parserInfo
  where
    parserInfo = OP.info (OP.helper <*> argParser)
                         (OP.progDesc "A program to convert a VCF file to Eigenstrat")

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseTransversionsOnly <*> parseSnpFile <*>  
                        parseChrom <*> parseOutChrom <*> parseOutPrefix <*> 
                        parseBcftoolsExe <*> parseVcfFileName
  where
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
                             \the snp files in Eigenstrat because of samtools specifications.")
    parseChrom = OP.option (T.pack <$> OP.str) (OP.long "chrom" <> OP.short 'c' <>
                            OP.metavar "<CHROM>" <> OP.help "specify the chromosome in the VCF \
                            \file to call from")
    parseOutChrom = OP.option (Just . T.pack <$> OP.str) (OP.long "outChrom" <>
                                    OP.metavar "<CHROM>" <>
                                   OP.help "specify the output chromosome name" <> OP.value Nothing)
    parseVcfFileName = OP.argument (fromText . T.pack <$> OP.str) (OP.metavar "<VCF_FILE>" <>
                                    OP.help "input VCF file")
    parseOutPrefix = OP.option (fromText . T.pack <$> OP.str)
                                 (OP.long "outPrefix" <> OP.short 'e' <>
                                 OP.metavar "<FILE_PREFIX>" <>
                                 OP.help "specify the filenames for the EigenStrat SNP and IND \
                                 \file outputs: <FILE_PREFIX>.snp.txt and <FILE_PREFIX>.ind.txt")
    parseBcftoolsExe = OP.option (fromText . T.pack <$> OP.str) (OP.long "bcftools" <>
                                  OP.value "bcftools" <> OP.showDefault <>
                                  OP.metavar "<BCFTOOLS_PATH>" <> OP.help "path to the \
                                  \bcftools executable, version >= 1.2")
    
runMain :: ProgOpt -> IO ()
runMain = runSafeT . runReaderT $ do --ReaderT ProgOpt (SafeT IO) ()
    chrom <- asks optChrom
    outChrom <- asks optOutChrom
    outPrefix <- asks optOutPrefix
    transversionsOnly <- asks optTransversionsOnly
    let outputChrom = case outChrom of
            Nothing -> chrom
            Just label -> label
    (vcfHeader, freqSumProducer) <- readVCF
    let snpOut = outPrefix <.> "snp.txt"
        indOut = outPrefix <.> "ind.txt"
    lift . withFile (T.unpack . format fp $ indOut) WriteMode $ \indOutHandle -> do
        let VCFheader _ sampleNames = vcfHeader
        mapM_ (\n -> liftIO $ hPutStrLn indOutHandle (n ++ "\tU\tUnknown")) sampleNames
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
    isTransition ref alt = ((ref == 'A') && (alt == 'G')) || ((ref == 'G') && (alt == 'A')) ||
                           ((ref == 'C') && (alt == 'T')) || ((ref == 'T') && (alt == 'C'))

readVCF :: ReaderT ProgOpt (SafeT IO) (VCFheader, Producer FreqSumRow (SafeT IO) ())
readVCF = do
    bcftools <- asks optBcftoolsExe
    vcfFN <- asks optVcfFile
    snpFN <- asks optSnpFile
    chrom <- asks optChrom
    let cmd = format (fp%" view -v snps -M2 -r "%s%" "%fp) bcftools chrom vcfFN
    vcfTextProd <- liftIO $ produceFromCommand cmd
    let snpTextProd = PT.readFile ((T.unpack . format fp) snpFN)
    (vcfHeader_, vcfProd) <- lift $ parseVCF vcfTextProd
    let snpProd = parsed snpParser snpTextProd >-> P.filter (\(SnpEntry c _ _ _) -> c == chrom)
    let jointProd = orderedZip cmp snpProd vcfProd >-> processVcfWithSnpFile
    return (vcfHeader_, fmap snd jointProd)
  where
    cmp (SnpEntry _ snpPos _ _) (VCFentry _ vcfPos _ _) = snpPos `compare` vcfPos

produceFromCommand :: Text -> IO (Producer Text (SafeT IO) ())
produceFromCommand cmd = do
    let createProcess = CreateProcess (ShellCommand (T.unpack cmd)) Nothing Nothing False False
                                                        False defaultHandler
    (p, _) <- pipeOutput Inherit Inherit createProcess
    return . void . decodeUtf8 $ p

parseVCF :: Producer Text (SafeT IO) () -> SafeT IO (VCFheader, Producer VCFentry (SafeT IO) ())
parseVCF prod = do
    (res, rest) <- runStateT (parse vcfHeaderParser) prod
    header <- case res of
        Nothing -> liftIO . throwIO $ AssertionFailed "vcf header not readible. VCF file empty?"
        Just (Left e_) -> do
            Right (chunk, _) <- next rest
            let msg = show e_ ++ T.unpack chunk
            liftIO . throwIO $ AssertionFailed ("VCF header parsing error: " ++ msg)
        Just (Right h) -> return h
    return (header, parsed vcfParser rest >>= liftErrors)
  where
    liftErrors res = case res of
        Left (e_, prod_) -> do
            Right (chunk, _) <- lift $ next prod_
            let msg = show e_ ++ T.unpack chunk
            lift . liftIO . throwIO $ AssertionFailed msg
        Right () -> return ()

vcfHeaderParser :: A.Parser VCFheader
vcfHeaderParser = VCFheader <$> A.many' doubleCommentLine <*> singleCommentLine
  where
    doubleCommentLine = do
        c1 <- A.string "##"
        s_ <- A.takeTill A.isEndOfLine <* A.endOfLine
        return $ T.append c1 s_
    singleCommentLine = do
        void $ A.char '#'
        s_ <- A.takeTill A.isEndOfLine <* A.endOfLine
        let fields = T.splitOn "\t" s_
        return . drop 9 . map T.unpack $ fields

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
    let filteredCoverages =
            [[c | (c, a) <- zip cov (Just ref:alt), a /= Nothing] | cov <- coverages]
    return $ VCFentry chrom pos (ref:filteredAlt) filteredCoverages
  where
    altAllele = Just <$> A.satisfy (A.inClass "ACTG")
    xAllele = (A.string "<X>" <|> A.string "<*>") >> return Nothing
    coverage = do
        _ <- A.decimal `A.sepBy1` (A.char ',') :: A.Parser [Int]
        _ <- A.char ':'
        A.decimal `A.sepBy1` (A.char ',')

word :: A.Parser Text
word = T.pack <$> A.many1 (A.satisfy (not . isSpace))

tab :: A.Parser ()
tab = A.char '\t' >> return ()

processVcfWithSnpFile :: Pipe (Maybe SnpEntry, Maybe VCFentry) FreqSumRow (SafeT IO) r
processVcfWithSnpFile = for cat $ \jointEntry -> do
    -- trace (show jointEntry) (return ())
    case jointEntry of
        (Just snpEntry, Just vcfEntry) -> do
            let SnpEntry snpChrom snpPos snpRef snpAlt = snpEntry
                VCFentry _ _ vcfAlleles vcfNums = vcfEntry  
                nrInds = length vcfNums
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
    let ret = SnpEntry chrom pos ref alt
    -- trace (show ret) $ return ()
    return ret

printEigenStrat :: Text -> Handle -> Pipe FreqSumRow Text (SafeT IO) r
printEigenStrat outChrom snpOutHandle = for cat $ \(FreqSumRow _ pos ref alt calls) -> do
    let n = format (s%"_"%d) outChrom pos
        snpLine = format (s%"\t"%s%"\t0\t"%d%"\t"%s%"\t"%s) n outChrom pos (T.singleton ref)
                            (T.singleton alt)
    liftIO . T.hPutStrLn snpOutHandle $ snpLine
    yield . T.concat . map (format d . toEigenStratNum) $ calls
  where
    toEigenStratNum c = case c of
        0 -> 2 :: Int
        1 -> 1
        2 -> 0
        -1 -> 9
        _ -> error ("unknown genotype " ++ show c)

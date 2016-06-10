{-# LANGUAGE OverloadedStrings #-}

import OrderedZip (orderedZip)

import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Applicative ((<|>))
import Control.Monad (forM_, void, when)
import Control.Monad.IO.Class (liftIO, MonadIO)
import Control.Monad.Trans.Class (lift)
import Control.Monad.Trans.State.Strict (runStateT)
import qualified Data.Attoparsec.Text as A
import Data.Char (isSpace)
import Data.Monoid ((<>))
import qualified Data.Text as T
import qualified Data.Text.IO as T
-- import Debug.Trace (trace)
import qualified Options.Applicative as OP
import Pipes (Pipe, yield, (>->), runEffect, Producer, Pipe, for, cat, next)
import Pipes.Attoparsec (parsed, parse, ParsingError)
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, SafeT, MonadSafe)
import qualified Pipes.Safe.Prelude as S
import qualified Pipes.Text.IO as PT
import System.IO (IOMode(..), withFile, Handle)
import Turtle.Format (format, d, s, (%))

data ProgOpt = ProgOpt {
    optSnpPosFile :: Maybe FilePath,
    optFillHomRef :: Bool,
    optOutPrefix :: FilePath,
    optChrom :: String,
    optOutChrom :: Maybe String,
    optTransversionsOnly :: Bool
}

data VCFheader = VCFheader [T.Text] [T.Text] deriving (Show)-- simple comment lines, sample names
data VCFentry = VCFentry T.Text Int Char Char [Int] deriving (Show) -- chrom, pos, ref, alt, dosages
data SnpEntry = SnpEntry T.Text Int Char Char deriving (Show)-- Chrom Pos Ref Alt

main :: IO ()
main = readOptions >>= runMain

readOptions :: IO ProgOpt
readOptions = OP.execParser parserInfo
  where
    parserInfo = OP.info (OP.helper <*> argParser)
                         (OP.progDesc "A program to convert a VCF file (stdin) to Eigenstrat")

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseSnpPosFile <*> parseFillHomRef <*> parseOutPrefix <*> parseChrom <*> 
                        parseOutChrom <*> parseTransversionsOnly
  where
    parseSnpPosFile = OP.option (Just <$> OP.str)
                   (OP.long "snpFile" <> OP.short 'f' <> OP.value Nothing <> OP.metavar "<FILE>" <>
                    OP.help "specify an Eigenstrat SNP file with the positions and alleles of a \
                             \reference set. \
                             \All  positions in the SNP file will be output, adding missing data \
                             \or hom-ref where necessary.")
    parseFillHomRef = OP.switch (OP.long "fillHomRef" <> OP.short 'r' <>
                             OP.help "Declare missing sites in VCF as Hom-Ref instead of missing. \ 
                                      \This is useful if you have high coverage data but your \
                                      \VCF file only contains non-ref sites")
    parseOutPrefix = OP.strOption (OP.long "outPrefix" <> OP.short 'e' <>
                                  OP.metavar "<FILE_PREFIX>" <>
                                  OP.help "specify the filenames for the EigenStrat SNP and IND \
                                  \file outputs: <FILE_PREFIX>.snp.txt and <FILE_PREFIX>.ind.txt")
    parseChrom = OP.strOption (OP.long "chrom" <> OP.short 'c' <>
                            OP.metavar "<CHROM>" <> OP.help "specify the chromosome in the VCF \
                            \file to \
                            \call from. This is important if a SNP file has been given.")
    parseOutChrom = OP.option (Just <$> OP.str) (OP.long "outChrom" <>
                                    OP.metavar "<CHROM>" <>
                                   OP.help "specify the output chromosome name" <> OP.value Nothing)
    parseTransversionsOnly = OP.switch (OP.long "transversionsOnly" <> OP.short 't' <>
                             OP.help "Remove transition SNPs from the output")
    
runMain :: ProgOpt -> IO ()
runMain (ProgOpt snpPosFile fillHomRef outPrefix chrom maybeOutChrom transversionsOnly) =
    runSafeT $ do
        (vcfHeader, vcfBody) <- readVCF
        let snpOut = outPrefix ++ ".snp.txt"
            indOut = outPrefix ++ ".ind.txt"
            VCFheader _ sampleNames = vcfHeader
            nrInds = length sampleNames
        S.withFile indOut WriteMode $ \indOutHandle -> do
            forM_ sampleNames $ \n -> do
                liftIO . T.hPutStrLn indOutHandle . T.intercalate "\t" $ [n, "U", "Unknown"]
        let vcfProducer = case snpPosFile of
                Just fn -> runJointly vcfBody nrInds chrom fn fillHomRef
                Nothing -> runSimple vcfBody chrom
        let outChrom = case maybeOutChrom of
                Just c -> c
                Nothing -> chrom
        let eigenStratPipe = S.withFile snpOut WriteMode (printEigenStrat outChrom)
        runEffect $
            vcfProducer >-> filterTransitions transversionsOnly >-> eigenStratPipe >-> printToStdOut
  where
    printToStdOut = for cat (liftIO . T.putStrLn)
    filterTransitions transversionsOnly =
        if transversionsOnly 
        then P.filter (\(VCFentry _ _ ref alt _) -> isTransversion ref alt)
        else cat
    isTransversion ref alt = not $ isTransition ref alt
    isTransition ref alt = ((ref == 'A') && (alt == 'G')) || ((ref == 'G') && (alt == 'A')) ||
                           ((ref == 'C') && (alt == 'T')) || ((ref == 'T') && (alt == 'C'))

readVCF :: (MonadIO m) => m (VCFheader, Producer VCFentry m ())
readVCF = do
    (res, rest) <- runStateT (parse vcfHeaderParser) PT.stdin
    header <- case res of
        Nothing -> liftIO . throwIO $ AssertionFailed "vcf header not readible. VCF file empty?"
        Just (Left e_) -> do
            Right (chunk, _) <- next rest
            let msg = show e_ ++ T.unpack chunk
            liftIO . throwIO $ AssertionFailed ("VCF header parsing error: " ++ msg)
        Just (Right h) -> return h
    return (header, parsed vcfEntryParser rest >>= liftErrors)

liftErrors :: (MonadIO m) => Either (ParsingError, Producer T.Text m r) () -> Producer a m ()
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
        return . drop 9 $ fields

vcfEntryParser :: A.Parser VCFentry
vcfEntryParser = do
    chrom <- word
    void tab
    pos <- A.decimal
    void tab
    void word
    void tab
    ref <- allele
    void tab
    alt <- allele
    void tab
    void $ A.count 4 (word >> tab)
    genotypes <- genotype `A.sepBy1` tab
    void A.endOfLine
    return $ VCFentry chrom pos ref alt genotypes
  where
    allele = A.satisfy (\c -> c `elem` actg)
    actg :: String
    actg = "ACTG"

tab :: A.Parser Char
tab = A.char '\t'

word :: A.Parser T.Text
word = A.takeTill isSpace

genotype :: A.Parser Int
genotype = do
    gen1 <- (A.char '0' <|> A.char '1' <|> A.char '.')
    void (A.char '/' <|> A.char '|')
    gen2 <- (A.char '0' <|> A.char '1' <|> A.char '.')
    _ <- A.takeTill (\a -> (a `elem` ['\r', '\t', '\n', ' ']))
    if   gen1 == '.' || gen2 == '.'
    then return (-1) 
    else return . length . filter (=='1') $ [gen1, gen2]

runJointly :: (MonadIO m, MonadSafe m) => Producer VCFentry m r -> Int -> String -> FilePath -> 
                                          Bool -> Producer VCFentry m r
runJointly vcfBody nrInds chrom snpPosFile fillHomRef =
    let snpProd = parsed snpParser (PT.readFile snpPosFile) >->
                  P.filter (\(SnpEntry c _ _ _) -> c == T.pack chrom) >>= liftErrors
        jointProd = snd <$> orderedZip cmp snpProd vcfBody
    in  jointProd >-> processVcfWithSnpFile nrInds fillHomRef
  where
    cmp (SnpEntry _ snpPos _ _) (VCFentry _ vcfPos _ _ _) = snpPos `compare` vcfPos

processVcfWithSnpFile :: (MonadIO m) => Int -> Bool ->
                         Pipe (Maybe SnpEntry, Maybe VCFentry) VCFentry m r
processVcfWithSnpFile nrInds fillHomRef = for cat $ \jointEntry -> do
    case jointEntry of
        (Just (SnpEntry snpChrom snpPos snpRef snpAlt), Nothing) -> do
            let dosages = if fillHomRef then replicate nrInds 0 else replicate nrInds (-1)
            yield $ VCFentry snpChrom snpPos snpRef snpAlt dosages
        (Just (SnpEntry snpChrom snpPos snpRef snpAlt),
         Just (VCFentry vcfChrom _ vcfRef vcfAlt vcfNums)) -> do
            when (length vcfNums /= nrInds) $ (liftIO . throwIO) (AssertionFailed "inconsistent \ 
                            \number of genotypes. Check that bam files have different \
                            \readgroup sample names")
            when (snpChrom /= vcfChrom) $ (liftIO . throwIO) (AssertionFailed "wrong chromosome \
                                            \name in VCF")
            let normalizedDosages =
                    if (vcfRef, vcfAlt) == (snpRef, snpAlt)
                    then vcfNums
                    else
                        if (vcfRef, vcfAlt) == (snpAlt, snpRef)
                        then map flipDosages vcfNums
                        else
                            replicate nrInds (-1)
            yield (VCFentry snpChrom snpPos snpRef snpAlt normalizedDosages)
        _ -> return ()
  where
    flipDosages d = case d of
        0 -> 2
        1 -> 1
        2 -> 0
        -1 -> -1

snpParser :: A.Parser SnpEntry
snpParser = do
    _ <- word
    tab
    chrom <- word
    tab
    _ <- word
    pos <- A.decimal
    tab
    ref <- A.satisfy (A.inClass "ACTG")
    tab
    alt <- A.satisfy (A.inClass "ACTG")
    _ <- A.satisfy (\c -> c == '\r' || c == '\n')
    let ret = SnpEntry chrom pos ref alt
    return ret

runSimple :: (MonadIO m) => Producer VCFentry m r -> String -> Producer VCFentry m r
runSimple vcfBody chrom = for vcfBody $ \v@(VCFentry vcfChrom pos ref alt dosages) -> do
    when (vcfChrom /= T.pack chrom) $ (liftIO . throwIO) (AssertionFailed "wrong chromosome in VCF")
    yield v

printEigenStrat :: (MonadIO m) => String -> Handle -> Pipe VCFentry T.Text m r
printEigenStrat outChrom snpOutHandle = for cat $ \(VCFentry chrom pos ref alt dosages) -> do
    let n = format (s%"_"%d) (T.pack outChrom) pos
        snpLine = format (s%"\t"%s%"\t0\t"%d%"\t"%s%"\t"%s) n (T.pack outChrom) pos
                         (T.singleton ref) (T.singleton alt)
    liftIO . T.hPutStrLn snpOutHandle $ snpLine
    yield . T.concat . map (format d . toEigenStratNum) $ dosages
  where
    toEigenStratNum :: Int -> Int
    toEigenStratNum c = case c of
        0 -> 2
        1 -> 1
        2 -> 0
        -1 -> 9
        _ -> error ("unknown dosage " ++ show c)

{-# LANGUAGE OverloadedStrings #-}

import SeqTools.OrderedZip (orderedZip)
import SeqTools.Fasta (loadFastaChrom)
import SeqTools.VCF (readVCF, VCFheader(..), VCFentry(..), SimpleVCFentry(..),
                     isBiallelicSnp, isTransversionSnp, liftParsingErrors, getDosages,
                     makeSimpleVCFentry)

import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad (forM_, void, when)
import Control.Monad.IO.Class (liftIO, MonadIO)
import qualified Data.Attoparsec.Text as A
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.ByteString.Char8 as B
import Data.Char (isSpace)
import Data.Monoid ((<>))
import qualified Data.Text as T
import qualified Data.Text.IO as T
-- import Debug.Trace (trace)
import qualified Options.Applicative as OP
import Pipes (Pipe, yield, (>->), runEffect, Producer, Pipe, for, cat)
import Pipes.Attoparsec (parsed)
import qualified Pipes.ByteString as PB
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, MonadSafe)
import qualified Pipes.Safe.Prelude as S
import qualified Pipes.Text.IO as PT
import System.IO (IOMode(..), Handle, withFile)
import Turtle.Format (format, d, s, (%))
import Turtle.Prelude (err)

data ProgOpt = ProgOpt (Maybe FilePath) (Maybe FilePath) FilePath String (Maybe String) Bool

data SnpEntry = SnpEntry T.Text Int Char Char deriving (Show)-- Chrom Pos Ref Alt

-- main :: IO ()
-- main = do
--     prod <- loadFastaChrom "/data/schiffels/ReferenceGenome/hs37d5.fa" "10"
--     let loop = for (prod >-> PB.take 200) $ \b -> liftIO $ print b
--     runEffect loop
    
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
    parseFillHomRef = OP.option (Just <$> OP.str) (OP.long "fillHomRef" <> OP.value Nothing <>
                             OP.short 'r' <>
                             OP.help "Input a reference sequence (uncompressed fasta format) to \
                                      \use to declare missing sites in the VCF as Hom-Ref instead of \
                                      \missing. This is useful if your VCF only contains non-ref \
                                      \sites. This option only makes sense if you use a SNP file.")
    parseOutPrefix = OP.strOption (OP.long "outPrefix" <> OP.short 'e' <>
                                  OP.metavar "<FILE_PREFIX>" <>
                                  OP.help "specify the filenames for the EigenStrat SNP and IND \
                                  \file outputs: <FILE_PREFIX>.snp.txt and <FILE_PREFIX>.ind.txt")
    parseChrom = OP.strOption (OP.long "chrom" <> OP.short 'c' <>
                            OP.metavar "<CHROM>" <> OP.help "specify the chromosome in the VCF \
                            \file to \
                            \call from. This is important if a SNP file has been given.")
    parseOutChrom = OP.option (Just <$> OP.str) (OP.long "outChrom" <> OP.value Nothing <>
                                    OP.metavar "<CHROM>" <>
                                   OP.help "specify the output chromosome name" <> OP.value Nothing)
    parseTransversionsOnly = OP.switch (OP.long "transversionsOnly" <> OP.short 't' <>
                             OP.help "Remove transition SNPs from the output")

runMain :: ProgOpt -> IO ()
runMain (ProgOpt snpPosFile fillHomRef outPrefix chrom maybeOutChrom transversionsOnly) =
    runSafeT $ do
        (vcfHeader, vcfBody) <- readVCF PT.stdin
        let snpOut = outPrefix ++ ".snp.txt"
            indOut = outPrefix ++ ".ind.txt"
            VCFheader _ sampleNames = vcfHeader
            nrInds = length sampleNames
        S.withFile indOut WriteMode $ \indOutHandle -> do
            forM_ sampleNames $ \n -> do
                liftIO . T.hPutStrLn indOutHandle . T.intercalate "\t" $ [n, "U", "Unknown"]
        let vcfBodyBiAllelic = vcfBody >-> P.filter (\e -> isBiallelicSnp (vcfRef e) (vcfAlt e))
        vcfProducer <- case snpPosFile of
                Just fn -> do
                    refSeq <- case fillHomRef of
                            Just fp -> do
                                S.withFile fp ReadMode $ \fh -> do
                                    bs <- liftIO $ loadFastaChrom fh chrom >>= PB.toLazyM
                                    return $ Just (BL.toStrict bs)
                            Nothing -> return Nothing
                    return $ runJointly vcfBodyBiAllelic nrInds chrom fn refSeq
                Nothing -> return $ runSimple vcfBodyBiAllelic chrom
        let outChrom = case maybeOutChrom of
                Just c -> c
                Nothing -> chrom
        let eigenStratPipe = S.withFile snpOut WriteMode (printEigenStrat outChrom)
        runEffect $ vcfProducer >-> filterTransitions >-> eigenStratPipe >-> printToStdOut
  where
    printToStdOut = for cat (liftIO . T.putStrLn)
    filterTransitions = if transversionsOnly
                        then P.filter (\e -> isTransversionSnp (sVCFref e) (sVCFalt e))
                        else cat

runJointly :: (MonadIO m, MonadSafe m) => Producer VCFentry m r -> Int -> String -> FilePath ->
                                          Maybe B.ByteString -> Producer SimpleVCFentry m r
runJointly vcfBody nrInds chrom snpPosFile refSeq =
    let snpProd = parsed snpParser (PT.readFile snpPosFile) >->
                  P.filter (\(SnpEntry c _ _ _) -> c == T.pack chrom) >>= liftParsingErrors
        jointProd = snd <$> orderedZip cmp snpProd vcfBody
    in  jointProd >-> processVcfWithSnpFile nrInds refSeq
  where
    cmp (SnpEntry _ snpPos _ _) vcfEntry = snpPos `compare` (vcfPos vcfEntry)

processVcfWithSnpFile :: (MonadIO m) => Int -> Maybe B.ByteString ->
                         Pipe (Maybe SnpEntry, Maybe VCFentry) SimpleVCFentry m r
processVcfWithSnpFile nrInds refSeq = for cat $ \jointEntry -> do
    case jointEntry of
        (Just (SnpEntry snpChrom snpPos snpRef snpAlt), Nothing) -> do
            let dosages = case refSeq of
                              Just seq_ -> let nuc = seq_ `B.index` (fromIntegral snpPos - 1)
                                           in  if nuc == snpRef
                                               then replicate nrInds (Just 0)
                                               else
                                                   if nuc == snpAlt
                                                   then replicate nrInds (Just 2)
                                                   else replicate nrInds (Nothing)
                              Nothing -> replicate nrInds Nothing
            yield $ SimpleVCFentry snpChrom snpPos (T.singleton snpRef) [T.singleton snpAlt]
                                   dosages
        (Just (SnpEntry snpChrom snpPos snpRef snpAlt), Just vcfEntry) -> do
            dosages <- case getDosages vcfEntry of
                Right dos -> return dos
                Left err -> liftIO . throwIO $ AssertionFailed err
            when (length dosages /= nrInds) $ (liftIO . throwIO) (AssertionFailed "inconsistent \
                            \number of genotypes. Check that bam files have different \
                            \readgroup sample names")
            when (snpChrom /= vcfChrom vcfEntry) $ do
                liftIO . throwIO $ AssertionFailed "wrong chromosome name in VCF"
            let normalizedDosages =
                    case vcfAlt vcfEntry of
                        [alt] -> if (vcfRef vcfEntry, alt) ==
                                        (T.singleton snpRef, T.singleton snpAlt)
                                 then dosages
                                 else
                                     if (vcfRef vcfEntry, alt) ==
                                             (T.singleton snpAlt, T.singleton snpRef)
                                     then map flipDosages dosages
                                     else replicate nrInds Nothing
                        _ -> replicate nrInds Nothing
            yield $ SimpleVCFentry snpChrom snpPos (T.singleton snpRef) [T.singleton snpAlt]
                                   normalizedDosages
        _ -> return ()
  where
    flipDosages dos = case dos of
        Just 0 -> Just 2
        Just 1 -> Just 1
        Just 2 -> Just 0
        _ -> Nothing

snpParser :: A.Parser SnpEntry
snpParser = do
    A.skipMany A.space
    void word
    A.skipMany1 A.space
    chrom <- word
    A.skipMany1 A.space
    void word
    A.skipMany1 A.space
    pos <- A.decimal
    A.skipMany1 A.space
    ref <- A.satisfy (A.inClass "ACTGX")
    A.skipMany1 A.space
    alt <- A.satisfy (A.inClass "ACTGX")
    void A.endOfLine
    let ret = SnpEntry chrom pos ref alt
    return ret
  where
    word = A.takeTill isSpace

runSimple :: (MonadIO m) => Producer VCFentry m r -> String -> Producer SimpleVCFentry m r
runSimple vcfBody chrom = for vcfBody $ \e -> do
    when (vcfChrom e /= T.pack chrom) $ (liftIO . throwIO) (AssertionFailed "wrong chromosome in VCF")
    case makeSimpleVCFentry e of
        Right e' -> do
            liftIO $ print e'
            yield e'
        Left err -> (liftIO . throwIO) (AssertionFailed err)

printEigenStrat :: (MonadIO m) => String -> Handle -> Pipe SimpleVCFentry T.Text m r
printEigenStrat outChrom snpOutHandle = for cat $ \(SimpleVCFentry _ pos ref alt dosages) -> do
    let n = format (s%"_"%d) (T.pack outChrom) pos
        snpLine = format (s%"\t"%s%"\t0\t"%d%"\t"%s%"\t"%s) n (T.pack outChrom) pos ref (head alt)
    liftIO . T.hPutStrLn snpOutHandle $ snpLine
    yield . T.concat . map (format d . toEigenStratNum) $ dosages
  where
    toEigenStratNum :: Maybe Int -> Int
    toEigenStratNum c = case c of
        Just 0 -> 2
        Just 1 -> 1
        Just 2 -> 0
        Nothing -> 9
        _ -> error ("unknown dosage " ++ show c)

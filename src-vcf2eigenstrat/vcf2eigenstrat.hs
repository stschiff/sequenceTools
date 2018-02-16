{-# LANGUAGE OverloadedStrings #-}

import SeqTools.OrderedZip (orderedZip)
import SequenceFormats.Fasta (loadFastaChrom)
import SequenceFormats.VCF (readVCF, VCFheader(..), VCFentry(..), SimpleVCFentry(..),
                     isBiallelicSnp, isTransversionSnp, getDosages,
                     makeSimpleVCFentry)
import SequenceFormats.Eigenstrat (EigenstratSnpEntry(..), readEigenstratSnpFile, writeEigenstrat, 
    GenoLine, GenoEntry(..), Sex(..), EigenstratIndEntry(..))

import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad (when)
import Control.Monad.IO.Class (liftIO, MonadIO)
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.ByteString.Char8 as B
import Data.Monoid ((<>))
import qualified Data.Text as T
-- import Debug.Trace (trace)
import Data.Vector (fromList)
import Data.Version (showVersion)
import qualified Options.Applicative as OP
import Paths_sequenceTools (version)
import Pipes (Pipe, yield, (>->), runEffect, Producer, Pipe, for, cat)
import qualified Pipes.ByteString as PB
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, MonadSafe)
import qualified Pipes.Safe.Prelude as S
import qualified Pipes.Text.IO as PT
import System.IO (IOMode(..))

data ProgOpt = ProgOpt (Maybe FilePath) (Maybe FilePath) FilePath T.Text (Maybe T.Text) Bool

main :: IO ()
main = readOptions >>= runMain

readOptions :: IO ProgOpt
readOptions = OP.execParser parserInfo
  where
    parserInfo = OP.info
        (pure (.) <*> versionInfoOpt <*> OP.helper <*> argParser)
        (OP.progDesc "A program to convert a VCF file (stdin) to Eigenstrat")
    versionInfoOpt = OP.infoOption ("This is vcf2eigenstrat from sequenceTools \
        \Version " ++ showVersion version)
        (OP.long "version" <> OP.help "Print version and exit")


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
    parseChrom = OP.option (T.pack <$> OP.str) (OP.long "chrom" <> OP.short 'c' <>
                            OP.metavar "<CHROM>" <> OP.help "specify the chromosome in the VCF \
                            \file to \
                            \call from. This is important if a SNP file has been given.")
    parseOutChrom = OP.option (Just . T.pack <$> OP.str) (OP.long "outChrom" <> OP.value Nothing <>
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
            genoOut = outPrefix ++ ".geno.txt"
            VCFheader _ sampleNames = vcfHeader
            nrInds = length sampleNames
            indEntries = [EigenstratIndEntry n Unknown "Unknown" | n <- sampleNames]
        let vcfBodyBiAllelic = vcfBody >-> P.filter (\e -> isBiallelicSnp (vcfRef e) (vcfAlt e))
        vcfProducer <- case snpPosFile of
                Just fn -> do
                    refSeq <- case fillHomRef of
                            Just fp -> do
                                S.withFile fp ReadMode $ \fh -> do
                                    bs <- liftIO $ loadFastaChrom fh (T.unpack chrom) >>= PB.toLazyM
                                    return $ Just (BL.toStrict bs)
                            Nothing -> return Nothing
                    return $ runJointly vcfBodyBiAllelic nrInds chrom fn refSeq
                Nothing -> return $ runSimple vcfBodyBiAllelic (T.unpack chrom)
        let outChrom = case maybeOutChrom of
                Just c -> c
                Nothing -> chrom
        runEffect $ vcfProducer >-> filterTransitions >-> eigenStratPipe outChrom >->
            writeEigenstrat genoOut snpOut indOut indEntries
  where
    filterTransitions = if transversionsOnly
                        then P.filter (\e -> isTransversionSnp (sVCFref e) (sVCFalt e))
                        else cat

runJointly :: (MonadIO m, MonadSafe m) => Producer VCFentry m r -> Int -> T.Text -> FilePath ->
                                          Maybe B.ByteString -> Producer SimpleVCFentry m r
runJointly vcfBody nrInds chrom snpPosFile refSeq =
    let snpProd = readEigenstratSnpFile snpPosFile >->
            P.filter (\(EigenstratSnpEntry c _ _ _) -> c == chrom)
        jointProd = snd <$> orderedZip cmp snpProd vcfBody
    in  jointProd >-> processVcfWithSnpFile nrInds refSeq
  where
    cmp (EigenstratSnpEntry _ snpPos _ _) vcfEntry = snpPos `compare` (vcfPos vcfEntry)

processVcfWithSnpFile :: (MonadIO m) => Int -> Maybe B.ByteString ->
                         Pipe (Maybe EigenstratSnpEntry, Maybe VCFentry) SimpleVCFentry m r
processVcfWithSnpFile nrInds refSeq = for cat $ \jointEntry -> do
    case jointEntry of
        (Just (EigenstratSnpEntry snpChrom snpPos snpRef snpAlt), Nothing) -> do
            let dosages = case refSeq of
                              Just seq_ -> let nuc = seq_ `B.index` (fromIntegral snpPos - 1)
                                           in  if nuc == snpRef
                                               then replicate nrInds (Just 0)
                                               else
                                                   if nuc == snpAlt
                                                   then replicate nrInds (Just 2)
                                                   else replicate nrInds (Nothing)
                              Nothing -> replicate nrInds Nothing
            yield $ SimpleVCFentry snpChrom snpPos (T.singleton snpRef)
                [T.singleton snpAlt] dosages
        (Just (EigenstratSnpEntry snpChrom snpPos snpRef snpAlt), Just vcfEntry) -> do
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
            yield $ SimpleVCFentry snpChrom snpPos (T.singleton snpRef)
                [T.singleton snpAlt] normalizedDosages
        _ -> return ()
  where
    flipDosages dos = case dos of
        Just 0 -> Just 2
        Just 1 -> Just 1
        Just 2 -> Just 0
        _ -> Nothing

runSimple :: (MonadIO m) => Producer VCFentry m r -> String -> Producer SimpleVCFentry m r
runSimple vcfBody chrom = for vcfBody $ \e -> do
    when (vcfChrom e /= T.pack chrom) $ (liftIO . throwIO) (AssertionFailed "wrong chromosome in VCF")
    case makeSimpleVCFentry e of
        Right e' -> do
            liftIO $ print e'
            yield e'
        Left err -> (liftIO . throwIO) (AssertionFailed err)

eigenStratPipe :: (MonadIO m) => T.Text -> Pipe SimpleVCFentry (EigenstratSnpEntry, GenoLine) m r
eigenStratPipe outChrom = P.map vcfToEigenstrat
  where
    vcfToEigenstrat (SimpleVCFentry _ pos ref alt dosages) =
        let snpEntry = EigenstratSnpEntry outChrom pos (T.head ref) (T.head . head $ alt)
            genoLine = fromList [dosageToCall d | d <- dosages]
        in  (snpEntry, genoLine)
    dosageToCall d = case d of
        Just 0 -> HomRef
        Just 1 -> Het
        Just 2 -> HomAlt
        Nothing -> Missing
        _ -> error ("unknown dosage " ++ show d)

{-# LANGUAGE OverloadedStrings #-}

import Pipes.OrderedZip (orderedZip)
import SequenceFormats.VCF (readVCFfromStdIn, VCFheader(..), VCFentry(..),
                     isBiallelicSnp, getDosages, vcfToFreqSumEntry)
import SequenceFormats.Eigenstrat (EigenstratSnpEntry(..), readEigenstratSnpFile, writeEigenstrat, 
    GenoLine, GenoEntry(..), Sex(..), EigenstratIndEntry(..))
import SequenceFormats.FreqSum (FreqSumEntry(..), freqSumEntryToText)
import SequenceFormats.Utils (Chrom(..))

import SequenceTools.Utils (versionInfoText)

import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad (when)
import Control.Monad.IO.Class (liftIO, MonadIO)
import qualified Data.ByteString.Char8 as B
import Data.Monoid ((<>))
-- import Debug.Trace (trace)
import Data.Vector (fromList)
import Data.Version (showVersion)
import qualified Options.Applicative as OP
import Paths_sequenceTools (version)
import Pipes (Pipe, yield, (>->), runEffect, Producer, Pipe, for, cat)
import qualified Pipes.Prelude as P
import Pipes.Safe (runSafeT, MonadSafe)

-- snpPosFile outPrefix
data ProgOpt = ProgOpt (Maybe FilePath) FilePath

main :: IO ()
main = readOptions >>= runMain

readOptions :: IO ProgOpt
readOptions = OP.execParser parserInfo
  where
    parserInfo = OP.info (pure (.) <*> versionInfoOpt <*> OP.helper <*> argParser)
        (OP.progDesc ("A program to convert a VCF file (stdin) to Eigenstrat. Part of \
            \sequenceTools version " ++ versionInfoText))
    versionInfoOpt = OP.infoOption (showVersion version)
        (OP.long "version" <> OP.help "Print version and exit")


argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseSnpPosFile <*> parseOutPrefix
  where
    parseSnpPosFile = OP.option (Just <$> OP.str)
                   (OP.long "snpFile" <> OP.short 'f' <> OP.value Nothing <> OP.metavar "<FILE>" <>
                    OP.help "specify an Eigenstrat SNP file with the positions and alleles of a \
                             \reference set. \
                             \All  positions in the SNP file will be output, adding missing data \
                             \or hom-ref where necessary.")
    parseOutPrefix = OP.strOption (OP.long "outPrefix" <> OP.short 'e' <>
                                  OP.metavar "<FILE_PREFIX>" <>
                                  OP.help "specify the filenames for the EigenStrat SNP and IND \
                                  \file outputs: <FILE_PREFIX>.snp.txt and <FILE_PREFIX>.ind.txt")

runMain :: ProgOpt -> IO ()
runMain (ProgOpt maybeSnpPosFile outPrefix) =
    runSafeT $ do
        (vcfHeader, vcfBody) <- readVCFfromStdIn
        let snpOut = outPrefix ++ ".snp.txt"
            indOut = outPrefix ++ ".ind.txt"
            genoOut = outPrefix ++ ".geno.txt"
            VCFheader _ sampleNames = vcfHeader
            nrInds = length sampleNames
            indEntries = [EigenstratIndEntry n Unknown "Unknown" | n <- sampleNames]
        let vcfBodyBiAllelic = vcfBody >-> P.filter (\e -> isBiallelicSnp (vcfRef e) (vcfAlt e))
        vcfProducer <- case maybeSnpPosFile of
                Just snpPosFile ->
                    return $ runJointly vcfBodyBiAllelic nrInds snpPosFile
                Nothing -> return $ runSimple vcfBodyBiAllelic
        runEffect $ vcfProducer >-> eigenStratPipe >-> writeEigenstrat genoOut snpOut indOut indEntries

runJointly :: (MonadIO m, MonadSafe m) => Producer VCFentry m r -> Int -> FilePath -> Producer FreqSumEntry m r
runJointly vcfBody nrInds snpPosFile =
    let snpProd = readEigenstratSnpFile snpPosFile
        jointProd = snd <$> orderedZip cmp snpProd vcfBody
    in  jointProd >-> processVcfWithSnpFile nrInds
  where
    cmp (EigenstratSnpEntry snpChrom' snpPos' _ _ _ _) vcfEntry = (snpChrom', snpPos') `compare` (vcfChrom vcfEntry, vcfPos vcfEntry)

processVcfWithSnpFile :: (MonadIO m) => Int -> Pipe (Maybe EigenstratSnpEntry, Maybe VCFentry) FreqSumEntry m r
processVcfWithSnpFile nrInds = for cat $ \jointEntry -> do
    case jointEntry of
        (Just (EigenstratSnpEntry snpChrom' snpPos' _ _ snpRef' snpAlt'), Nothing) -> do
            let dosages = replicate nrInds Nothing
            yield $ FreqSumEntry snpChrom' snpPos' snpRef' snpAlt' dosages
        (Just (EigenstratSnpEntry snpChrom' snpPos' _ _ snpRef' snpAlt'), Just vcfEntry) -> do
            dosages <- case getDosages vcfEntry of
                Right dos -> return dos
                Left err -> liftIO . throwIO $ AssertionFailed err
            when (length dosages /= nrInds) $ (liftIO . throwIO) (AssertionFailed "inconsistent number of genotypes.")
            let normalizedDosages =
                    case vcfAlt vcfEntry of
                        [alt] -> if (vcfRef vcfEntry, alt) == (B.singleton snpRef', B.singleton snpAlt')
                                 then dosages
                                 else
                                     if (vcfRef vcfEntry, alt) == (B.singleton snpAlt', B.singleton snpRef')
                                     then map flipDosages dosages
                                     else replicate nrInds Nothing
                        _ -> replicate nrInds Nothing
            yield $ FreqSumEntry snpChrom' snpPos' snpRef' snpAlt' normalizedDosages
        _ -> return ()
  where
    flipDosages dos = case dos of
        Just 0 -> Just 2
        Just 1 -> Just 1
        Just 2 -> Just 0
        _ -> Nothing

runSimple :: (MonadIO m) => Producer VCFentry m r -> Producer FreqSumEntry m r
runSimple vcfBody = for vcfBody $ \e -> do
    case vcfToFreqSumEntry e of
        Right e' -> do
            liftIO . B.putStr . freqSumEntryToText $ e'
            yield e'
        Left err -> (liftIO . throwIO) (AssertionFailed err)

eigenStratPipe :: (MonadIO m) => Pipe FreqSumEntry (EigenstratSnpEntry, GenoLine) m r
eigenStratPipe = P.map vcfToEigenstrat
  where
    vcfToEigenstrat (FreqSumEntry chrom pos ref alt dosages) =
        let snpId' = B.pack (unChrom chrom <> show pos)
            snpEntry = EigenstratSnpEntry chrom pos 0.0 snpId' ref alt
            genoLine = fromList [dosageToCall d | d <- dosages]
        in  (snpEntry, genoLine)
    dosageToCall d = case d of
        Just 0 -> HomRef
        Just 1 -> Het
        Just 2 -> HomAlt
        Nothing -> Missing
        _ -> error ("unknown dosage " ++ show d)

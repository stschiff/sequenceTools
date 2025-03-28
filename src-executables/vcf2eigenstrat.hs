{-# LANGUAGE OverloadedStrings #-}

import           Pipes.OrderedZip           (orderedZip)
import           SequenceFormats.Eigenstrat (EigenstratIndEntry (..),
                                             EigenstratSnpEntry (..), Sex (..),
                                             readEigenstratSnpFile,
                                             writeEigenstrat)
import           SequenceFormats.FreqSum    (FreqSumEntry (..))
import           SequenceFormats.VCF        (VCFentry (..), VCFheader (..),
                                             getDosages, isBiallelicSnp,
                                             readVCFfromStdIn,
                                             vcfToFreqSumEntry)

import           SequenceTools.Utils        (freqSumToEigenstrat,
                                             versionInfoOpt, versionInfoText)

import           Control.Exception.Base     (AssertionFailed (..), throwIO)
import           Control.Monad              (when)
import           Control.Monad.IO.Class     (MonadIO, liftIO)
import qualified Data.ByteString.Char8      as B
-- import Debug.Trace (trace)
import qualified Options.Applicative        as OP
import           Pipes                      (Pipe, Producer, cat, for,
                                             runEffect, yield, (>->))
import qualified Pipes.Prelude              as P
import           Pipes.Safe                 (MonadSafe, runSafeT)

-- snpPosFile outPrefix
data ProgOpt = ProgOpt (Maybe FilePath) FilePath

main :: IO ()
main = readOptions >>= runMain

readOptions :: IO ProgOpt
readOptions = OP.execParser parserInfo
  where
    parserInfo = OP.info (pure (.) <*> versionInfoOpt <*> OP.helper <*> argParser)
        (OP.progDesc ("A program to convert a VCF file (stdin) to Eigenstrat. " ++ versionInfoText))

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseSnpPosFile <*> parseOutPrefix
  where
    parseSnpPosFile = OP.option (Just <$> OP.str)
                   (OP.long "snpFile" <> OP.short 'f' <> OP.value Nothing <> OP.metavar "<FILE>" <>
                    OP.help "specify an Eigenstrat SNP file with the positions and alleles of a \
                             \reference set. If this option is given, only positions that are both in the SNP file \
                             \and in the VCF will be output. Without this option, all sites in the VCF will be output")
    parseOutPrefix = OP.strOption (OP.long "outPrefix" <> OP.short 'e' <>
                                  OP.metavar "<FILE_PREFIX>" <>
                                  OP.help "specify the prefix for the EigenStrat output files. Output files will then be  \
                                  \<FILE_PREFIX>.geno, <FILE_PREFIX>.snp and <FILE_PREFIX>.ind")

runMain :: ProgOpt -> IO ()
runMain (ProgOpt maybeSnpPosFile outPrefix) =
    runSafeT $ do
        (vcfHeader, vcfBody) <- readVCFfromStdIn
        let snpOut = outPrefix ++ ".snp"
            indOut = outPrefix ++ ".ind"
            genoOut = outPrefix ++ ".geno"
            VCFheader _ sampleNames = vcfHeader
            nrInds = length sampleNames
            indEntries = [EigenstratIndEntry n Unknown "Unknown" | n <- sampleNames]
        let vcfBodyBiAllelic = vcfBody >-> P.filter (\e -> isBiallelicSnp (vcfRef e) (vcfAlt e))
        vcfProducer <- case maybeSnpPosFile of
                Just snpPosFile ->
                    return $ runJointly vcfBodyBiAllelic nrInds snpPosFile
                Nothing -> return $ runSimple vcfBodyBiAllelic
        runEffect $ vcfProducer >-> P.map freqSumToEigenstrat >-> writeEigenstrat genoOut snpOut indOut indEntries

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
        (Just (EigenstratSnpEntry snpChrom' snpPos' gpos snpId' snpRef' snpAlt'), Nothing) -> do
            let dosages = replicate nrInds Nothing
            yield $ FreqSumEntry snpChrom' snpPos' (Just snpId') (Just gpos) snpRef' snpAlt' dosages
        (Just (EigenstratSnpEntry snpChrom' snpPos' gpos snpId' snpRef' snpAlt'), Just vcfEntry) -> do
            dosages <- liftIO $ getDosages vcfEntry
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
            yield $ FreqSumEntry snpChrom' snpPos' (Just snpId') (Just gpos) snpRef' snpAlt' normalizedDosages
        _ -> return ()
  where
    flipDosages dos = case dos of
        Just (0, 2) -> Just (2, 2)
        Just (1, 2) -> Just (1, 2)
        Just (2, 2) -> Just (0, 2)
        Just (0, 1) -> Just (1, 1)
        Just (1, 1) -> Just (0, 1)
        _           -> Nothing

runSimple :: (MonadIO m) => Producer VCFentry m r -> Producer FreqSumEntry m r
runSimple vcfBody = for vcfBody $ \e -> do
    fs <- liftIO $ vcfToFreqSumEntry e
    yield fs

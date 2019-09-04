{-# LANGUAGE OverloadedStrings #-}

import SequenceTools.Utils (versionInfoOpt, versionInfoText)

import Control.Monad.IO.Class (MonadIO, liftIO)
import qualified Data.ByteString.Char8 as B
import Data.Maybe (isJust, fromJust)
import SequenceFormats.VCF (readVCFfromFile, VCFentry(..), isBiallelicSnp)
import SequenceFormats.Eigenstrat (EigenstratSnpEntry(..), readBimFile, writeEigenstratSnp)
import SequenceFormats.Utils (Chrom(..))
import qualified Options.Applicative as OP
import Pipes (runEffect, (>->), Pipe, for, cat, yield)
import qualified Pipes.Prelude as P
import Pipes.OrderedZip (orderedZip)
import Pipes.Safe (runSafeT)
import System.IO (stdout, stderr, hPutStrLn)

type ProgOpt = (FilePath, FilePath)

main :: IO ()
main = OP.execParser optionSpec >>= runWithOptions

optionSpec :: OP.ParserInfo ProgOpt
optionSpec = OP.info (pure (.) <*> versionInfoOpt <*> OP.helper <*> optParser) (
    OP.fullDesc <>
    OP.progDesc ("Program to flip all alleles in a BIM file with a refrence VCF file \
        \into the correct REF and ALT order and the reference strand" <> versionInfoText)) 

optParser :: OP.Parser ProgOpt
optParser = (,) <$> OP.strOption (OP.long "BIM-file" <> OP.metavar "FILE") <*> OP.strOption (OP.long "VCF-file" <> OP.metavar "FILE")

runWithOptions :: ProgOpt -> IO ()
runWithOptions (bimFileName, vcfFileName) = do
    runSafeT $ do
        (_, vcfProd) <- readVCFfromFile vcfFileName
        let vcfProdFiltered = vcfProd >-> P.filter isValidSnp
            bimProd = readBimFile bimFileName
            mergedProd = orderedZip comp bimProd vcfProdFiltered >> return ()
        _ <- runEffect $ mergedProd >-> processJointEntries >-> writeEigenstratSnp stdout
        return ()
  where
    comp (EigenstratSnpEntry bimChrom bimPos _ _ _ _) vcfEntry = 
        compare (bimChrom, bimPos) (vcfChrom vcfEntry, vcfPos vcfEntry)
    isValidSnp vcf = isBiallelicSnp (vcfRef vcf) (vcfAlt vcf)

processJointEntries :: (MonadIO m) =>
    Pipe (Maybe EigenstratSnpEntry, Maybe VCFentry) EigenstratSnpEntry m ()
processJointEntries = for cat (\(mes,mvcf) -> do
    case (mes, mvcf) of
        (Just es, Just vcf) ->
            if isJust (vcfId vcf) && fromJust (vcfId vcf) /= snpId es
            then
                liftIO $ hPutStrLn stderr ("SKIP_ID_MISMATCH: " <> B.unpack (snpId es) <> " <> " <> B.unpack (fromJust (vcfId vcf)))
            else do
                if snpRef es == B.head (vcfRef vcf) && snpAlt es == B.head (head (vcfAlt vcf))
                then yield es
                else do
                    liftIO $ hPutStrLn stderr ("WARN_ALLELE_CHANGE at " <> B.unpack (snpId es) <> " (" <> unChrom (snpChrom es) <> ":" <> show (snpPos es) <>
                            "): (" <> [snpRef es] <> "," <> [snpAlt es] <> ") -> (" <> B.unpack (vcfRef vcf) <> "," <> (B.unpack . head) (vcfAlt vcf) <> ")") 
                    yield es {snpRef = B.head (vcfRef vcf), snpAlt = B.head (head (vcfAlt vcf))}
        (Just es, Nothing) ->
            liftIO $ hPutStrLn stderr ("SKIP_MISSING: Did not find position " <> B.unpack (snpId es) <> " (" <> unChrom (snpChrom es) <> ":" <> show (snpPos es) <> ") in VCF file")
        _ -> return ())
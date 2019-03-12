#!/usr/bin/env stack
-- stack runghc

{-# LANGUAGE OverloadedStrings #-}

import Data.Maybe (isJust, fromJust)
import Data.Ord (Ordering(..))
import qualified Data.Text as T
import Data.Version (showVersion)
import Filesystem.Path.CurrentOS (encodeString)
import Turtle hiding (stdout, cat)
import SequenceFormats.VCF (readVCFfromFile, VCFentry(..), isBiallelicSnp)
import SequenceFormats.Eigenstrat (EigenstratSnpEntry(..), readBimFile, writeEigenstratSnp)
import SequenceFormats.Utils (Chrom(..))
import Paths_sequenceTools (version)
import Pipes (runEffect, (>->), await, Consumer, Producer, Pipe, for, cat, yield)
import qualified Pipes.Prelude as P
import Pipes.OrderedZip (orderedZip)
import Pipes.Safe (runSafeT, MonadSafe)
import Pipes.Safe.Prelude (withFile)
import Pipes.Text.Encoding (decodeAscii)
import Prelude hiding (FilePath, withFile, hGetContents)
import System.IO (IOMode(..), stdout, stderr)
import qualified Text.PrettyPrint.ANSI.Leijen as PT

optParser :: Parser (FilePath, FilePath)
optParser = (,) <$> argPath "BIM-file" "BIM-file" <*> argPath "VCF-file" "VCF-file"

main = do
    (bimFileName, vcfFileName) <- options helpMessage optParser
    runSafeT $ do
        (vcfHeader, vcfProd) <- readVCFfromFile (encodeString vcfFileName)
        let vcfProdFiltered = vcfProd >-> P.filter isValidSnp
            bimProd = readBimFile (encodeString bimFileName)
            mergedProd = orderedZip comp bimProd vcfProd >> return ()
        _ <- runEffect $ mergedProd >-> processJointEntries >-> writeEigenstratSnp stdout
        return ()
  where
    helpMessage = Description $ PT.text ("program to normalise a BIM file with respect to a VCF file. Part of sequenceTools version " ++ showVersion version)
    comp (EigenstratSnpEntry bimChrom bimPos _ _ _ _) vcfEntry = 
        compare (bimChrom, bimPos) (vcfChrom vcfEntry, vcfPos vcfEntry)
    isValidSnp vcf = isBiallelicSnp (vcfRef vcf) (vcfAlt vcf)

processJointEntries :: (MonadIO m) =>
    Pipe (Maybe EigenstratSnpEntry, Maybe VCFentry) EigenstratSnpEntry m ()
processJointEntries = for cat (\(mes,mvcf) -> do
    case (mes, mvcf) of
        (Just es, Just vcf) ->
            if isJust (vcfId vcf) && fromJust (vcfId vcf) /= snpId es
            then liftIO . err . unsafeTextToLine $
                format ("SKIP_ID_MISMATCH: "%s%" <> "%s) (snpId es) (fromJust $ vcfId vcf)
            else do
                if snpRef es == T.head (vcfRef vcf) && snpAlt es == T.head (head (vcfAlt vcf))
                then yield es
                else do
                    liftIO . err . unsafeTextToLine $
                        format ("WARN_ALLELE_CHANGE at "%s%" ("%s%":"%d%
                            "): ("%s%","%s%") -> ("%s%","%s%")") (snpId es) 
                            (unChrom $ snpChrom es) (snpPos es) (T.singleton $ snpRef es) 
                            (T.singleton $ snpAlt es) (vcfRef vcf)
                            (head $ vcfAlt vcf)
                    yield es {snpRef = T.head (vcfRef vcf), snpAlt = T.head (head (vcfAlt vcf))}
        (Just es, Nothing) ->
            liftIO . err . unsafeTextToLine $
                format ("SKIP_MISSING: Did not find position "%s%" ("%s%":"%d%
                    ") in VCF file") (snpId es) (unChrom $ snpChrom es) (snpPos es)
        _ -> return ())
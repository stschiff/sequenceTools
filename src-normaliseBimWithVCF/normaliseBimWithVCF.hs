#!/usr/bin/env stack
-- stack runghc

{-# LANGUAGE OverloadedStrings #-}

import Data.Version (showVersion)
import Filesystem.Path.CurrentOS (encodeString)
import Turtle
import SequenceFormats.VCF (readVCFfromFile, VCFentry(..), isBiallelicSnp)
import SequenceFormats.Eigenstrat (EigenstratSnpEntry(..), readBimFile)
import Paths_sequenceTools (version)
import Pipes (runEffect, (>->), await, Consumer, Producer)
import qualified Pipes.Prelude as P
import Pipes.OrderedZip (orderedZip)
import Pipes.Safe (runSafeT, MonadSafe)
import Pipes.Safe.Prelude (withFile)
import Pipes.Text.Encoding (decodeAscii)
import Prelude hiding (FilePath, withFile, hGetContents)
import System.IO (IOMode(..))
import qualified Text.PrettyPrint.ANSI.Leijen as PT

optParser :: Parser (FilePath, FilePath)
optParser = (,) <$> argPath "BIM-file" "BIM-file" <*> argPath "VCF-file" "VCF-file"

main = do
    (bimFileName, vcfFileName) <- options helpMessage optParser
    runSafeT $ do
        (vcfHeader, vcfProd) <- readVCFfromFile (encodeString vcfFileName)
        let vcfProdFiltered = vcfProd >-> P.filter isValidSnp
            bimProd = readBimFile (encodeString bimFileName)
            mergedProd = orderedZip comp bimProd vcfProd
        _ <- runEffect $ mergedProd >-> P.print
        return ()
  where
    helpMessage = Description $ PT.text ("program to normalise a BIM file with respect to a VCF file. Part of sequenceTools version " ++ showVersion version)
    comp (EigenstratSnpEntry bimChrom bimPos _ _ _ _) vcfEntry =
        compare (vcfChrom vcfEntry, vcfPos vcfEntry) (bimChrom, bimPos)
    isValidSnp vcf = isBiallelicSnp (vcfRef vcf) (vcfAlt vcf)


-- outConsumer :: (Monad m) => Consumer (Maybe EigenstratSnpEntry, Maybe VCFentry) m ((),())
-- outConsumer =  do
--     (eigenstratEntry, vcf) <- await
--     liftIO $ print (eigenstratEntry, vcf)
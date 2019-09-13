{-# LANGUAGE OverloadedStrings #-}
module SequenceTools.Utils (versionInfoOpt, versionInfoText, sampleWithoutReplacement,
    freqSumToEigenstrat, dosageToEigenstratGeno) where 

import SequenceFormats.FreqSum (FreqSumEntry(..))
import SequenceFormats.Eigenstrat (EigenstratSnpEntry(..), GenoLine, GenoEntry(..))
import SequenceFormats.Utils (Chrom(..))

import qualified Data.ByteString.Char8 as B
import Data.Vector (fromList)
import Data.Version (showVersion)
import qualified Options.Applicative as OP
import Paths_sequenceTools (version)
import System.Random (randomRIO)

versionInfoOpt :: OP.Parser (a -> a)
versionInfoOpt = OP.infoOption (showVersion version) (OP.long "version" <> OP.help "Print version and exit")

versionInfoText :: String
versionInfoText = "This tool is part of sequenceTools version " ++ showVersion version

sampleWithoutReplacement :: [a] -> Int -> IO (Maybe [a])
sampleWithoutReplacement = go []
  where
    go res _ 0 = return $ Just res
    go res xs n
        | n > length xs = return Nothing
        | n == length xs = return $ Just (xs ++ res)
        | otherwise = do
                rn <- randomRIO (0, length xs - 1)
                let a = xs !! rn
                    xs' = remove rn xs
                go (a:res) xs' (n - 1)
    remove i xs = let (ys, zs) = splitAt i xs in ys ++ tail zs

-- |convert a freqSum entry to an eigenstrat SNP entry
freqSumToEigenstrat :: Bool -> FreqSumEntry -> (EigenstratSnpEntry, GenoLine)
freqSumToEigenstrat diploidizeCall (FreqSumEntry chrom@(Chrom c) pos maybeSnpId maybeGeneticPos ref alt calls) =
    let snpId_ = case maybeSnpId of 
            Just id_ -> id_
            Nothing -> c <> "_" <> B.pack (show pos)
        geneticPos = case maybeGeneticPos of
            Just p -> p
            Nothing -> 0.0
        snpEntry = EigenstratSnpEntry chrom pos geneticPos snpId_ ref alt
        geno = fromList . map (dosageToEigenstratGeno diploidizeCall) $ calls
    in  (snpEntry, geno)

-- |convert a Dosage to an eigenstrat-encoded genotype
dosageToEigenstratGeno :: Bool -> Maybe Int -> GenoEntry
dosageToEigenstratGeno diploidizeCall c =
    if diploidizeCall then
        case c of
            Just 0 -> HomRef
            Just 1 -> HomAlt
            Nothing -> Missing
            _ -> error "illegal call for pseudo-haploid Calling method"
    else
        case c of
            Just 0 -> HomRef
            Just 1 -> Het
            Just 2 -> HomAlt
            Nothing -> Missing
            _ -> error ("unknown genotype " ++ show c)

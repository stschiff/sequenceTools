{-# LANGUAGE OverloadedStrings #-}

import qualified SequenceFormats.Eigenstrat as E
import SequenceFormats.FreqSum (readFreqSumFile, readFreqSumStdIn, FreqSumHeader(..), 
    FreqSumEntry(..))
import SequenceFormats.Eigenstrat (readEigenstrat)

import Data.Text (Text, pack, unpack)
import Data.Version (showVersion)
import Paths_sequenceTools (version)
import Pipes (for, Producer, runEffect, (>->), yield)
import qualified Pipes.Prelude as P
import Pipes.Safe (MonadSafe, runSafeT)
import Prelude hiding (FilePath)
import Turtle hiding (stdin, x)

data ProgOpt = ProgOpt InputOption Bool

data InputOption = FreqsumInput (Maybe FilePath) | EigenstratInput FilePath FilePath FilePath

data InputEntry = InputEntry Text [Genotype] deriving (Show)
data Genotype = HomRef | HomAlt | Het | Missing deriving (Show)

main :: IO ()
main = options descString argParser >>= runWithOpts
  where
    descString = "A program to evaluate per-chromosome and total statistics of genotyping \
            \data, read either as Eigenstrat or FreqSum"

argParser :: Parser ProgOpt
argParser = ProgOpt <$> parseInputOpt <*> parseVersionOpt
  where
    parseInputOpt = parseFreqsumInput <|> parseEigenstratInput
    parseVersionOpt = switch "version" 'v' "Print the version and exit"

parseFreqsumInput :: Parser InputOption
parseFreqsumInput =
    process <$> optText "freqsum" 'f' "a freqsum file to read as input. Use - to read from stdin"
  where
    process p =
        if p == "-"
        then FreqsumInput Nothing
        else FreqsumInput . Just . fromText $ p

parseEigenstratInput :: Parser InputOption
parseEigenstratInput = EigenstratInput <$> parseGenoFile <*> parseSnpFile <*> parseIndFile
  where
    parseGenoFile = optPath "eigenstratGeno" 'g' "Eigenstrat Geno File"
    parseSnpFile = optPath "eigenstratSnp" 's' "Eigenstrat Snp File"
    parseIndFile = optPath "eigenstratInd" 'i' "Eigenstrat Ind File"

runWithOpts :: ProgOpt -> IO ()
runWithOpts (ProgOpt inputOpt optVersion) = do
    if optVersion
    then do
        let v = pack . showVersion $ version
        echo . repr $ format ("This is genoStats from sequenceTools version "%s) v
    else runSafeT $ do
        (names, entryProducer) <- case inputOpt of 
            FreqsumInput fsFile -> runWithFreqSum fsFile
            EigenstratInput genoFile snpFile indFile -> runWithEigenstrat genoFile snpFile indFile
        liftIO . print $ names
        runEffect $ entryProducer >-> P.print
        
        -- collectStats names entryProducer >>= reportResults

runWithFreqSum :: (MonadSafe m) => Maybe FilePath -> m ([Text], Producer InputEntry m ())
runWithFreqSum fsFile = do
    (FreqSumHeader names nrHaps, fsProd) <- case fsFile of
        Nothing -> readFreqSumStdIn
        Just fn -> readFreqSumFile (unpack $ format fp fn)
    let prod = for fsProd $ \(FreqSumEntry chrom _ _ _ counts) -> do
            let genotypes = do
                    (count', nrHap) <- zip counts nrHaps
                    case count' of
                        0 -> return HomRef
                        x | x == nrHap -> return HomAlt
                          | x > 0 && x < nrHap -> return Het
                          | x < 0 -> return Missing
                        _ -> error "should not happen"
            yield $ InputEntry chrom genotypes
    return (names, prod)

runWithEigenstrat :: (MonadSafe m) =>
    FilePath -> FilePath -> FilePath -> m ([Text], Producer InputEntry m ())
runWithEigenstrat genoFile snpFile indFile = do
    let genoFile' = unpack $ format fp genoFile
        snpFile' = unpack $ format fp snpFile
        indFile' = unpack $ format fp indFile
    (indEntries, genoStream) <- readEigenstrat genoFile' snpFile' indFile'
    let names = [name | E.EigenstratIndEntry name _ _ <- indEntries]
    let prod = for genoStream $ \(E.EigenstratSnpEntry chrom _ _ _, genoLine) -> do
            let genotypes = do
                    geno <- genoLine
                    case geno of
                        E.HomRef -> return HomRef
                        E.HomAlt -> return HomAlt
                        E.Het -> return Het
                        E.Missing -> return Missing
            yield $ InputEntry chrom genotypes
    return (names, prod)

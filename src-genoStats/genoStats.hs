{-# LANGUAGE OverloadedStrings #-}

import qualified SequenceFormats.Eigenstrat as E
import SequenceFormats.FreqSum (parseFreqSum, FreqSumHeader(..), FreqSumEntry(..), 
    streamEigenstratGeno, streamEigenstratSnp)

import Control.Error (Script, runScript, scriptIO)
import Data.Text (Text, pack, unpack)
import Data.Version (showVersion)
import Paths_sequenceTools (version)
import Pipes (for, Producer, runEffect, (>->), yield)
import qualified Pipes.Prelude as P
import Pipes.Text.IO (fromHandle)
import Prelude hiding (FilePath)
import System.IO (IOMode(..), openFile, stdin)
import Turtle hiding (stdin)

data ProgOpt = ProgOpt {
    optInput :: InputOption,
    optVersion :: Bool
}

data InputOption = FreqsumInput (Maybe FilePath) | EigenstratInput FilePath FilePath FilePath

data InputEntry = InputEntry Int [Genotype]
data Genotype = HomRef | HomAlt | Het | Missing

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
parseFreqsumInput = do
    p <- optText "freqsum" 'f' "a freqsum file to read as input. Use - to read from stdin"
    if p == "-"
    then return $ FreqsumInput Nothing
    else return . FreqsumInput . Just . fromText $ p

parseEigenstratInput :: Parser InputOption
parseEigenstratInput = EigenstratInput <$> parseSnpFile <*> parseIndFile <*> parseGenoFile
  where
    parseSnpFile = optPath "eigenstratSnp" 's' "Eigenstrat Snp File"
    parseIndFile = optPath "eigenstratInd" 'i' "Eigenstrat Ind File"
    parseGenoFile = optPath "eigenstratGeno" 'g' "Eigenstrat Geno File"

runWithOpts :: ProgOpt -> IO ()
runWithOpts (ProgOpt inputOpt optVersion) = do
    if optVersion
    then do
        let v = pack . showVersion $ version
        echo . repr $ format ("This is genoStats from sequenceTools version "%s) v
    else do
        (names, entryProducer) <- case inputOpt of 
            FreqsumInput fsFile -> runWithFreqSum fsFile
            EigenstratInput genoFile snpFile indFile -> runWithEigenstrat genoFile snpFile indFile
        print names
        runEffect $ entryProducer >-> P.print
        
        -- collectStats names entryProducer >>= reportResults

runWithFreqSum :: Maybe FilePath -> IO ([Text], Producer InputEntry IO ())
runWithFreqSum fsFile = do
    handle <- case fsFile of
        Nothing -> return stdin
        Just f -> openFile (repr f) ReadMode
    (FreqSumHeader namesStr nrHaps, fsProd) <- parseFreqSum handle
    let names = map pack namesStr
    let prod = for fsProd $ \(FreqSumEntry chrom _ _ _ counts) -> do
            let genotypes = do
                    (count, nrHap) <- zip counts nrHaps
                    case count of
                        0 -> return HomRef
                        x | x == nrHap -> return HomAlt
                          | x > 0 && x < nrHap -> return Het
                          | x < 0 -> return Missing
            yield $ InputEntry chrom genotypes
    return (names, prod)

runWithEigenstrat :: FilePath -> FilePath -> FilePath -> IO ([Text], Producer InputEntry IO ())
runWithEigenstrat genoFile snpFile indFile = do
    indEntries <- E.readEigenstratInd (repr indFile)
    let names = [name | E.EigenstratIndEntry name _ _ <- indEntries]
    genoHandle <- openFile genoFile ReadMode
    snpHandle <- openFile snpFile ReadMode
    let genoProd = streamEigenstratGeno genoHandle
        snpProd = streamEigenstratSnp snpHandle
        zippedProducer = P.zip snpProd genoProd
    let prod = for zippedProducer $ \(E.EigenstratSnpEntry chrom _ _ _, genoLine) -> do
            let genotypes = do
                    geno <- genoLine
                    case geno of
                        E.HomRef -> return HomRef
                        E.HomAlt -> return HomAlt
                        E.Het -> return Het
                        E.Missing -> return Missing
            yield $ InputEntry chrom genotypes
    return (names, prod)

        
        
    
    
    
{-# LANGUAGE OverloadedStrings #-}

import qualified SequenceFormats.Eigenstrat as E
import SequenceFormats.FreqSum (readFreqSumFile, readFreqSumStdIn, FreqSumHeader(..), 
    FreqSumEntry(..))
import SequenceFormats.Eigenstrat (readEigenstrat)

import Control.Foldl (purely, Fold)
import Control.Monad (forM_)
import Data.Text (Text, pack, unpack)
import Data.Text.IO (putStrLn)
import Data.Version (showVersion)
import qualified Data.Vector as V
import Lens.Family2 (view)
import Paths_sequenceTools (version)
import Pipes (for, Producer, runEffect, (>->), yield, Consumer, cat)
import Pipes.Group (groupsBy, folds)
import Pipes.Safe (MonadSafe, runSafeT)
import Prelude hiding (FilePath, putStrLn)
import Turtle hiding (stdin, x, view, cat)

data ProgOpt = ProgOpt InputOption Bool

data InputOption = FreqsumInput (Maybe FilePath) | EigenstratInput FilePath FilePath FilePath

data InputEntry = InputEntry Text (V.Vector Genotype) deriving (Show)
data Genotype = HomRef | HomAlt | Het | Missing deriving (Show)

type StatsReportAllSamples = [StatsReport]

data StatsReport = StatsReport {
    srNrSitesMissing :: Int,
    srNrSitesHomRef :: Int,
    srNrSitesHomAlt :: Int,
    srNrSitesHet :: Int
} deriving (Show)

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
        err . unsafeTextToLine $ format ("processing samples: "%w) names
        runEffect $ runStats names entryProducer >-> reportStats names

runWithFreqSum :: (MonadSafe m) => Maybe FilePath -> m ([Text], Producer InputEntry m ())
runWithFreqSum fsFile = do
    (FreqSumHeader names nrHaps, fsProd) <- case fsFile of
        Nothing -> readFreqSumStdIn
        Just fn -> readFreqSumFile (unpack $ format fp fn)
    let prod = for fsProd $ \(FreqSumEntry chrom _ _ _ counts) -> do
            let genotypes = V.fromList $ do
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
            let genotypes = V.fromList $ do
                    geno <- genoLine
                    case geno of
                        E.HomRef -> return HomRef
                        E.HomAlt -> return HomAlt
                        E.Het -> return Het
                        E.Missing -> return Missing
            yield $ InputEntry chrom genotypes
    return (names, prod)

runStats :: (MonadIO m) => [Text] -> Producer InputEntry m () ->
    Producer (Text, StatsReportAllSamples) m ()
runStats names entryProducer =
    let groupedProd = view (groupsBy (\(InputEntry c1 _) (InputEntry c2 _) -> c1 == c2))        
            entryProducer
    in  purely folds (runStatsPerChrom (length names)) groupedProd

runStatsPerChrom :: Int -> Fold InputEntry (Text, StatsReportAllSamples)
runStatsPerChrom nrSamples = (,) <$> getChrom <*>
    traverse runStatsPerChromPerSample [0..(nrSamples - 1)]    
    
getChrom :: Fold InputEntry Text
getChrom = Fold (\_ (InputEntry c _) -> c) "" id

runStatsPerChromPerSample :: Int -> Fold InputEntry StatsReport
runStatsPerChromPerSample i = Fold step initial extract
  where
    step :: StatsReport -> InputEntry -> StatsReport
    step accum@(StatsReport miss homr homa het) (InputEntry _ line) = case (line V.! i) of
        Missing -> accum {srNrSitesMissing = miss + 1}
        HomRef -> accum {srNrSitesHomRef = homr + 1}
        HomAlt -> accum {srNrSitesHomAlt = homa + 1}
        Het -> accum {srNrSitesHet = het + 1}
    initial :: StatsReport
    initial = StatsReport 0 0 0 0
    extract :: StatsReport -> StatsReport
    extract = id
    
reportStats :: (MonadIO m) => [Text] -> Consumer (Text, StatsReportAllSamples) m ()
reportStats names = do
    liftIO . putStrLn $ format ("Chrom\tSample\tMissing\tHomRef\tHomAlt\tHet")
    for cat $ \(chrom, reports) -> do
        forM_ (zip names reports) $ \(n, StatsReport mis ref alt het) -> do
            liftIO . putStrLn $
                format (s%"\t"%s%"\t"%d%"\t"%d%"\t"%d%"\t"%d) chrom n mis ref alt het

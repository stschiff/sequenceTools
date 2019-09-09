{-# LANGUAGE OverloadedStrings #-}

import SequenceFormats.FreqSum (readFreqSumFile, readFreqSumStdIn, FreqSumHeader(..), 
    FreqSumEntry(..))
import SequenceFormats.Eigenstrat (readEigenstrat, GenoEntry(..), GenoLine, EigenstratSnpEntry(..), EigenstratIndEntry(..))
import SequenceFormats.Utils (Chrom(..))
import SequenceTools.Utils (versionInfoOpt, versionInfoText)

import Control.Applicative ((<|>))
import Control.Foldl (purely, Fold(..))
import Control.Monad (forM_)
import Control.Monad.IO.Class (MonadIO, liftIO)
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector as V
import Lens.Family2 (view)
import qualified Options.Applicative as OP
import Pipes (for, Producer, (>->), yield, Consumer, cat)
import Pipes.Group (groupsBy, folds)
import Pipes.Safe (MonadSafe, runSafeT)
import qualified Pipes.Prelude as P
import System.IO (hPutStrLn, stderr)
import Text.Printf (printf)

data ProgOpt = ProgOpt InputOption

data InputOption = FreqsumInput (Maybe FilePath) | EigenstratInput FilePath FilePath FilePath

data InputEntry = InputEntry Chrom GenoLine deriving (Show)

type StatsReportAllSamples = [StatsReport]

data StatsReport = StatsReport {
    srNrSitesMissing :: Int,
    srNrSitesHomRef :: Int,
    srNrSitesHomAlt :: Int,
    srNrSitesHet :: Int
} deriving (Show)

main :: IO ()
main = OP.execParser optionSpec >>= runWithOpts

optionSpec :: OP.ParserInfo ProgOpt
optionSpec = OP.info (pure (.) <*> versionInfoOpt <*> OP.helper <*> argParser) (
    OP.progDesc ("A program \
    \to evaluate per-chromosome and total statistics \
    \of genotyping data, read either as Eigenstrat (by specifying options -g, -s and -i) or FreqSum (default, or by specifying option -f). " ++ versionInfoText))

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> (parseFreqsumInput <|> parseEigenstratInput)

parseFreqsumInput :: OP.Parser InputOption
parseFreqsumInput =
    process <$> OP.strOption (OP.long "freqsum" <> OP.short 'f' <> OP.help "a freqsum file to read as input. Use - to read from stdin (the default)" <>
        OP.value "-" <> OP.showDefault <> OP.metavar "FILEPATH")
  where
    process p =
        if p == "-"
        then FreqsumInput Nothing
        else FreqsumInput (Just p)

parseEigenstratInput :: OP.Parser InputOption
parseEigenstratInput = EigenstratInput <$> parseGenoFile <*> parseSnpFile <*> parseIndFile
  where
    parseGenoFile = OP.strOption (OP.long "eigenstratGeno" <> OP.short 'g' <> OP.help "Eigenstrat Geno File" <> OP.metavar "FILEPATH")
    parseSnpFile =  OP.strOption (OP.long "eigenstratSnp" <> OP.short 's' <> OP.help "Eigenstrat Snp File" <> OP.metavar "FILEPATH")
    parseIndFile =  OP.strOption (OP.long "eigenstratInd" <> OP.short 'i' <> OP.help "Eigenstrat Ind File" <> OP.metavar "FILEPATH")

runWithOpts :: ProgOpt -> IO ()
runWithOpts (ProgOpt inputOpt) = runSafeT $ do
    (names, entryProducer) <- case inputOpt of 
        FreqsumInput fsFile -> runWithFreqSum fsFile
        EigenstratInput genoFile snpFile indFile -> runWithEigenstrat genoFile snpFile indFile
    liftIO $ hPutStrLn stderr ("processing samples: " <> show names)
    let p = runStats names entryProducer >-> P.tee (reportStats names)
    totalReport <- purely P.fold (accumulateAllChromStats names) p
    printReports names totalReport

runWithFreqSum :: (MonadSafe m) => Maybe FilePath -> m ([String], Producer InputEntry m ())
runWithFreqSum fsFile = do
    (FreqSumHeader names nrHaps, fsProd) <- case fsFile of
        Nothing -> readFreqSumStdIn
        Just fn -> readFreqSumFile fn
    let prod = for fsProd $ \(FreqSumEntry chrom _ _ _ _ counts) -> do
            let genotypes = V.fromList $ do
                    (count', nrHap) <- zip counts nrHaps
                    case count' of
                        Just 0 -> return HomRef
                        Just x | x == nrHap -> return HomAlt
                               | x > 0 && x < nrHap -> return Het
                        Nothing -> return Missing
                        _ -> error "should not happen"
            yield $ InputEntry chrom genotypes
    return (map B.unpack names, prod)

runWithEigenstrat :: (MonadSafe m) =>
    FilePath -> FilePath -> FilePath -> m ([String], Producer InputEntry m ())
runWithEigenstrat genoFile snpFile indFile = do
    (indEntries, genoStream) <- readEigenstrat genoFile snpFile indFile
    let names = [name | EigenstratIndEntry name _ _ <- indEntries]
    let prod = genoStream >-> P.map (\(EigenstratSnpEntry chrom _ _ _ _ _, genoLine) -> InputEntry chrom genoLine)
    return (names, prod)

runStats :: (MonadIO m) => [String] -> Producer InputEntry m () ->
    Producer (Chrom, StatsReportAllSamples) m ()
runStats names entryProducer =
    let groupedProd = view (groupsBy (\(InputEntry c1 _) (InputEntry c2 _) -> c1 == c2))        
            entryProducer
    in  purely folds (runStatsPerChrom (length names)) groupedProd

runStatsPerChrom :: Int -> Fold InputEntry (Chrom, StatsReportAllSamples)
runStatsPerChrom nrSamples = (,) <$> getChrom <*>
    traverse runStatsPerChromPerSample [0..(nrSamples - 1)]    
    
getChrom :: Fold InputEntry Chrom
getChrom = Fold (\_ (InputEntry c _) -> c) (Chrom "") id

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
    
reportStats :: (MonadIO m) => [String] -> Consumer (Chrom, StatsReportAllSamples) m ()
reportStats names = do
    liftIO . putStrLn $ "Chrom\tSample\tMissing\tHomRef\tHomAlt\tHet"
    for cat $ \(chrom, reports) -> printReports names (chrom, reports)

printReports :: (MonadIO m) => [String] -> (Chrom, StatsReportAllSamples) -> m ()
printReports names (chrom, reports) =
    forM_ (zip names reports) $ \(n, StatsReport mis ref alt het) -> do
        let total = mis + ref + alt + het
            misPerc = round $ (fromIntegral mis / fromIntegral total) * (100.0 :: Double) :: Int
        liftIO . putStrLn $ printf "%s\t%s\t%d(%d%%)\t%d\t%d\t%d" (unChrom chrom) n mis misPerc ref alt het

accumulateAllChromStats :: [String] ->
    Fold (Chrom, StatsReportAllSamples) (Chrom, StatsReportAllSamples)
accumulateAllChromStats names = Fold step initial extract
  where
    step :: StatsReportAllSamples -> (Chrom, StatsReportAllSamples) -> StatsReportAllSamples
    step sumReports (_, newReports) = do
        (StatsReport smiss shomr shoma shet, StatsReport miss homr homa het) <-
                zip sumReports newReports
        return $ StatsReport (smiss + miss) (shomr + homr) (shoma + homa) (shet + het)
    initial :: StatsReportAllSamples
    initial = [StatsReport 0 0 0 0 | _ <- names]
    extract :: StatsReportAllSamples -> (Chrom, StatsReportAllSamples)
    extract r = (Chrom "Total", r)

        

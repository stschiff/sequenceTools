{-# LANGUAGE OverloadedStrings #-}

import Control.Monad (forM)
import qualified Data.Text as T
import Data.Text.Read (decimal)
import qualified Options.Applicative as OP
import Prelude hiding (FilePath)
import System.Random (randomRIO, mkStdGen, setStdGen)
import Turtle hiding (decimal)

data ProgOpt = ProgOpt {
    optCallingMode :: CallingMode,
    optSeed :: Maybe Int,
    optMinDepth :: Int,
    optSnpFile :: Maybe FilePath,
    optRegion :: Text,
    optReference :: FilePath,
    optBamFiles :: [FilePath],
    optOutFormat :: OutFormat
}

data CallingMode = MajorityCalling | RandomCalling deriving (Show, Read)
data OutFormat = EigenStrat | FreqSumFormat deriving (Show, Read)

data FreqSumRow = FreqSumRow Text Int Char Char [Maybe Int] deriving (Show)

main = OP.execParser parser >>= runWithOpts
  where
    parser = OP.info (OP.helper <*> argParser) (OP.progDesc "A program to perform simple genotype calling directly from BAM")

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseCallingMode <*> parseSeed <*> parseMinDepth <*> parseSnpFile <*> parseChrom <*>
                        parseRef <*> OP.some parseBams <*> parseFormat
  where
    parseCallingMode = OP.option OP.auto (OP.long "mode" <> OP.short 'm' <> OP.value RandomCalling <> OP.showDefault <>
                                          OP.metavar "<MODE>" <>
                                          OP.help "specify the mode of calling: MajorityCalling or RandomCalling")
    parseSeed = OP.option (Just <$> OP.auto) (OP.long "seed" <> OP.value Nothing <> OP.metavar "<RANDOM_SEED>" <>
                                              OP.help "random seed used for random calling. If not given, use system clock to seed the random number generator")
    parseMinDepth = OP.option OP.auto (OP.long "minDepth" <> OP.short 'd' <> OP.value 1 <> OP.showDefault <>
                                       OP.metavar "<DEPTH>" <>
                                       OP.help "specify the minimum depth for a call")
    parseSnpFile = OP.option (Just . fromText <$> OP.auto)
                   (OP.long "snpFile" <> OP.short 'f' <> OP.value Nothing <> OP.metavar "<FILE>" <>
                    OP.help "specify a SNP file (EigenStrat format) for the positions and alleles to call")
    parseChrom = OP.option OP.auto (OP.long "chrom" <> OP.short 'c' <> OP.metavar "<CHROM>" <>
                                    OP.help "specify the region in the BAM file to call from. Can be just the chromosome, or a string of the form CHROM:START-END")
    parseRef = OP.option (fromText <$> OP.auto) (OP.long "reference" <> OP.short 'r' <> OP.metavar "<REF>" <>
                                  OP.help "the reference fasta file")
    parseBams = OP.argument (fromText <$> OP.auto) (OP.metavar "<BAM_FILE>" <> OP.help "input file")
    parseFormat = OP.option OP.auto (OP.metavar "<OUT_FORMAT>" <> OP.long "format" <> OP.short 'o' <>
                                     OP.value FreqSumFormat <> OP.showDefault <>
                                     OP.help "specify output format: EigenStrat or FreqSum. EigenStrat requires a SnpFile given via --snpFile")
    
runWithOpts :: ProgOpt -> IO ()
runWithOpts opts = do
    case optSeed opts of
        Nothing -> return ()
        Just seed -> setStdGen $ mkStdGen seed
    let args = concat [["mpileup", "-q30", "-Q30", "-C50", "-f", format fp (optReference opts), "-g", "-t", 
                        "DPR", "-r", optRegion opts], map (format fp) (optBamFiles opts)]
        mpileupStream = inproc "samtools" args empty
        rawVcfStream = inproc "bcftools" ["view", "-m3", "-M3", "-v", "snps"] mpileupStream
        freqSumStream = processLines (optCallingMode opts) (optMinDepth opts) (optSnpFile opts) rawVcfStream
    view freqSumStream
    
processLines :: CallingMode -> Int -> Maybe FilePath -> Shell Text -> Shell FreqSumRow
processLines mode minDepth maybeSnpFile rawStream = do
    line <- rawStream
    let (chrom:posStr:_:refStr:altStr:_:_:_:_:genStrings) = T.words line
    pos <- liftIO $ case decimal posStr of
            Left e -> error e
            Right (p, _) -> return p
    genotypes <- liftIO $ mapM (callFromGenString mode minDepth) genStrings
    return $ FreqSumRow chrom pos (T.head refStr) (T.head altStr) genotypes

callFromGenString :: CallingMode -> Int -> Text -> IO (Maybe Int)
callFromGenString mode minDepth str = do
    let [_, dprStr] = T.splitOn ":" str
        numStrings = T.splitOn "," dprStr
    (numRef:numAlt:_) <- forM numStrings $ \numStr -> do
        case decimal numStr of
            Left e -> error e
            Right (n, _) -> return n
    if (numRef + numAlt < minDepth) then
        return Nothing
        else
            case mode of
                MajorityCalling -> case numRef `compare` numAlt of
                    LT -> return $ Just 2
                    GT -> return $ Just 0
                    EQ -> do
                        rn <- randomRIO (1, numRef + numAlt)
                        if rn <= numRef then return (Just 0) else return (Just 2)
                RandomCalling -> do
                        rn <- randomRIO (1, numRef + numAlt)
                        if rn <= numRef then return (Just 0) else return (Just 2)
    
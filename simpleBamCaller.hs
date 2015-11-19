{-# LANGUAGE OverloadedStrings #-}

import Control.Error (runScript, Script)
import qualified Data.Text as T
import qualified Options.Applicative as OP
import Prelude hiding (FilePath)
import Turtle

data ProgOpt = ProgOpt {
    optCallingMode :: CallingMode,
    optMinDepth :: Int,
    optSnpFile :: Maybe FilePath,
    optRegion :: Text,
    optReference :: FilePath,
    optBamFile :: FilePath,
    optOutFormat :: Format
}

data CallingMode = MajorityCalling | RandomCalling deriving (Show, Read)
data Format = EigenStrat | FreqSumFormat deriving (Show, Read)

data FreqSumRow = FreqSumRow Text Int Char Char [Maybe Int] deriving (Show)

main = OP.execParser parser >>= runWithOpts
  where
    parser = OP.info (OP.helper <*> argParser) (OP.progDesc "A program to perform simple genotype calling directly from BAM")

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseCallingMode <*> parseMinDepth <*> parseSnpFile <*> parseChrom <*> parseRef <*> parseBam <*>
                        parseFormat <*> parseSnpOut
  where
    parseCallingMode = OP.option OP.auto (OP.long "mode" <> OP.short 'm' <> OP.value RandomCalling <> OP.showDefault <>
                                          OP.metavar "<MODE>" <>
                                          OP.help "specify the mode of calling: MajorityCalling or RandomCalling")
    parseMinDepth = OP.option OP.auto (OP.long "minDepth" <> OP.short 'd' <> OP.value 1 <> OP.showDefault <>
                                       OP.metavar "<DEPTH>"
                                       OP.help "specify the minimum depth for a call")
    parseSnpFile = OP.option (Just . fromText <$> OP.auto)
                   (OP.long "snpFile" <> OP.short 'f' <> OP.value Nothing <> OP.metavar "<FILE>" <>
                    OP.help "specify a SNP file (EigenStrat format) for the positions and alleles to call")
    parseChrom = OP.option OP.auto (OP.long "chrom" <> OP.short 'c' <> OP.metavar "<CHROM>" <>
                                    OP.help "specify the region in the BAM file to call from. Can be just the chromosome, or a string of the form CHROM:START-END")
    parseRef = OP.option OP.auto (OP.long "reference" <> OP.short 'r' <> OP.metavar "<REF>" <>
                                  OP.help "the reference fasta file")
    parseBam = OP.argument OP.auto (OP.metavar "<BAM_FILE>" <> OP.help "input file")
    parseFormat = OP.argument OP.auto (OP.metavar "<OUT_FORMAT>" <> OP.long "format" <> OP.short 'o' <>
                                       OP.value FreqSumFormat <> OP.showDefault <>
                                       OP.help "specify output format: EigenStrat or FreqSum. EigenStrat requires a SnpFile given via --snpFile")
    
runWithOpts :: ProgOpt -> IO ()
runWithOpts opts = runScript $ do
    let rawVcfStream = inproc "samtools" ["mpileup", "-q", "30", "-Q", "30", "-C", "50", "-f",
                                          format fp (optReference opts), "-g", "-t", "DPR", "-r", optRegion opts,
                                          format fp (optBamFile opts)] empty
    let freqSumStream = processLines (optCallingMode opts) (optMinDepth opts) (optSnpFile opts) rawVcfStream
    view freqSumStream
    
processLines :: CallingMode -> Int -> Maybe FilePath -> Shell Text -> Shell FreqSumRow
processLines mode minDepth maybeSnpFile rawStream = do
    line <- rawStream
    let fields = T.words line
    
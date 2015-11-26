{-# LANGUAGE OverloadedStrings #-}

import OrderedZip (orderedZip)

import Control.Monad (forM)
import qualified Data.Text as T
import qualified Options.Applicative as OP
import Prelude hiding (FilePath)
import System.Random (randomRIO, mkStdGen, setStdGen)
import Turtle

data ProgOpt = ProgOpt {
    optCallingMode :: CallingMode,
    optSeed :: Maybe Int,
    optMinDepth :: Int,
    optTransversionOnly :: Bool,
    optSnpFile :: Maybe FilePath,
    optRegion :: Text,
    optReference :: FilePath,
    optOutFormat :: OutFormat,
    optBamFiles :: [FilePath]
}

data CallingMode = MajorityCalling | RandomCalling deriving (Show, Read)
data OutFormat = EigenStrat | FreqSumFormat deriving (Show, Read)

data FreqSumRow = FreqSumRow Text Int Char Char [Int] deriving (Show)

data VCFentry = VCFentry Text Int [Char] [[Int]] deriving (Show) -- Chrom Pos Alleles Number_of_reads_per_individual

data SnpEntry = SnpEntry Text Int Char Char deriving (Show)-- Chrom Pos Ref Alt

main = OP.execParser parser >>= runWithOpts
  where
    parser = OP.info (OP.helper <*> argParser) (OP.progDesc "A program to perform simple genotype calling directly \ 
                                                             \from BAM")

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseCallingMode <*> parseSeed <*> parseMinDepth <*> parseTrans <*> parseSnpFile <*>
                        parseChrom <*> parseRef <*> parseFormat <*> OP.some parseBams
  where
    parseCallingMode = OP.option OP.auto (OP.long "mode" <> OP.short 'm' <> OP.value RandomCalling <> OP.showDefault <>
                                          OP.metavar "<MODE>" <>
                                          OP.help "specify the mode of calling: MajorityCalling or RandomCalling")
    parseSeed = OP.option (Just <$> OP.auto) (OP.long "seed" <> OP.value Nothing <> OP.metavar "<RANDOM_SEED>" <>
                                              OP.help "random seed used for random calling. If not given, use system \  
                                                       \clock to seed the random number generator")
    parseMinDepth = OP.option OP.auto (OP.long "minDepth" <> OP.short 'd' <> OP.value 1 <> OP.showDefault <>
                                       OP.metavar "<DEPTH>" <>
                                       OP.help "specify the minimum depth for a call")
    parseTrans = OP.switch (OP.long "--transversionsOnly" <> OP.short 't' <>
                            OP.help "restrict calling to transversion SNPs. Replace Transition SNP calls with Missing \                                      \Data")
    parseSnpFile = OP.option (Just . fromText . T.pack <$> OP.str)
                   (OP.long "snpFile" <> OP.short 'f' <> OP.value Nothing <> OP.metavar "<FILE>" <>
                    OP.help "specify a SNP file for the positions and alleles to call. All \
                             \positions in the SNP file will be output, adding missing data where necessary. \
                             \The file should have three columns (space- or tab-separated): Chromosome, position and \
                             \alleles, where the alleles should be one reference and one alternative allele \
                             \separated by a comma.")
    parseChrom = OP.option (T.pack <$> OP.str) (OP.long "chrom" <> OP.short 'c' <> OP.metavar "<CHROM>" <>
                                    OP.help "specify the region in the BAM file to call from. Can be just the \
                                             \chromosome, or a string of the form CHROM:START-END")
    parseRef = OP.option (fromText . T.pack <$> OP.str) (OP.long "reference" <> OP.short 'r' <> OP.metavar "<REF>" <>
                                  OP.help "the reference fasta file")
    parseBams = OP.argument (fromText . T.pack <$> OP.str) (OP.metavar "<BAM_FILE>" <> OP.help "input file")
    parseFormat = OP.option OP.auto (OP.metavar "<OUT_FORMAT>" <> OP.long "format" <> OP.short 'o' <>
                                     OP.value FreqSumFormat <> OP.showDefault <>
                                     OP.help "specify output format: EigenStrat or FreqSum. EigenStrat requires a \     
                                              \SnpFile given via --snpFile")
    
runWithOpts :: ProgOpt -> IO ()
runWithOpts (ProgOpt mode seed minDepth transversionsOnly snpFile region reference outFormat bamFiles) = do
    case seed of
        Nothing -> return ()
        Just seed_ -> setStdGen $ mkStdGen seed_
    let vcfStream = parseVCF (inshell cmd empty)
        chrom = head (T.splitOn ":" region)
        freqSumStream = case snpFile of
            Nothing -> processLinesSimple mode minDepth transversionsOnly vcfStream
            Just fn -> processLinesWithSnpFile fn nrInds chrom mode minDepth transversionsOnly vcfStream
    view freqSumStream
  where
    cmd = case snpFile of
        Nothing -> format ("samtools mpileup -q30 -Q30 -C50 -I -f "%fp%" -g -t DPR -r "%s%" "%s%
                           " | bcftools view -m3 -M3 -v snps -H") reference region bams
        Just fn -> format ("samtools mpileup -q30 -Q30 -C50 -I -f "%fp%" -g -t DPR -r "%s%" -l "%fp%" "%s%
                           " | bcftools view -m2 -M3 -H") reference region fn bams
    nrInds = length bamFiles
    bams = (T.intercalate " " (map (format fp) bamFiles))

parseVCF :: Shell Text -> Shell VCFentry
parseVCF vcfTextStream = do
    line <- vcfTextStream
    let matchResult = match vcfPattern line
    case matchResult of
        [vcf] -> return vcf
        _ -> error $ "Could not parse VCF line " ++ show line

vcfPattern :: Pattern VCFentry
vcfPattern = do
    chrom <- word
    skip tab
    pos <- decimal
    skip tab >> skip word
    skip tab
    ref <- oneOf "NACTG"
    skip tab
    alt <- altAlleles <|> (text "<X>" >> return [])
    skip tab
    skip $ count 3 (word >> tab)
    skip (text "PL:DPR") >> skip tab
    coverages <- coverage `sepBy1` tab
    return $ VCFentry chrom pos (ref:alt) coverages
  where
    coverage = do
        skip word
        skip (char ':')
        init <$> decimal `sepBy1` (char ',')
    altAlleles = do
        a <- (oneOf "ACTG") `sepBy` (char ',')
        skip (text ",<X>")
        return a

word :: Pattern Text
word = plus $ noneOf "\r\t\n "

processLinesSimple :: CallingMode -> Int -> Bool -> Shell VCFentry -> Shell FreqSumRow
processLinesSimple mode minDepth transversionsOnly vcfTextStream = do
    e@(VCFentry chrom pos [ref, alt] covNums) <- vcfTextStream
    False <- return $ ref == 'N'
    when transversionsOnly $ do
        True <- return $ isTransversion ref alt 
        return ()
    genotypes <- liftIO $ mapM (callGenotype mode minDepth) covNums
    return $ FreqSumRow chrom pos ref alt genotypes

isTransversion :: Char -> Char -> Bool
isTransversion ref alt = not isTransition
  where
    isTransition = ((ref == 'A') && (alt == 'G')) || ((ref == 'G') && (alt == 'A')) ||
                   ((ref == 'C') && (alt == 'T')) || ((ref == 'T') && (alt == 'C'))

callGenotype :: CallingMode -> Int -> [Int] -> IO Int
callGenotype mode minDepth covs = do
    case covs of
        [numRef] -> return 0
        [numRef, numAlt] -> do
            if (numRef + numAlt < minDepth) then
                return (-1)
                else
                    case mode of
                        MajorityCalling -> case numRef `compare` numAlt of
                            LT -> return 2
                            GT -> return 0
                            EQ -> do
                                rn <- randomRIO (1, numRef + numAlt)
                                if rn <= numRef then return 0 else return 2
                        RandomCalling -> do
                                rn <- randomRIO (1, numRef + numAlt)
                                if rn <= numRef then return 0 else return 2
        _ -> return (-1)

processLinesWithSnpFile :: FilePath -> Int -> Text -> CallingMode -> Int -> Bool -> Shell VCFentry -> Shell FreqSumRow
processLinesWithSnpFile fn nrInds chrom mode minDepth transversionsOnly vcfStream = do
    x@(Just (SnpEntry snpChrom snpPos snpRef snpAlt), maybeVCF) <- orderedZip cmp snpStream vcfStream
    err (format w x)
    case maybeVCF of
        Nothing -> return $ FreqSumRow snpChrom snpPos snpRef snpAlt (replicate nrInds (-1))
        Just (VCFentry vcfChrom vcfPos vcfAlleles vcfNums) -> do
            genotypes <- if transversionsOnly && (not (isTransversion snpRef snpAlt))
                then return (replicate nrInds (-1))
                else do
                    case vcfAlleles of
                            [ref] -> if ref == snpRef
                                then return (replicate nrInds 0)
                                else return (replicate nrInds (-1))
                            [ref, alt] -> if [ref, alt] == [snpRef, snpAlt]
                                then liftIO $ mapM (callGenotype mode minDepth) vcfNums
                                else return (replicate nrInds (-1))
                            _ -> return (replicate nrInds (-1))
            return $ FreqSumRow snpChrom snpPos snpRef snpAlt genotypes
  where
    cmp (SnpEntry _ snpPos _ _) (VCFentry _ vcfPos _ _) = snpPos `compare` vcfPos
    snpStream = do
        line <- input fn
        let matchResult = match snpPattern line
        case matchResult of
            [e@(SnpEntry snpChrom _ _ _)] -> do
                True <- return (snpChrom == chrom)
                return e

snpPattern :: Pattern SnpEntry
snpPattern = do
    chrom <- word
    skip tab
    pos <- decimal
    skip tab
    ref <- oneOf "ACTG"
    skip (char ',')
    alt <- oneOf "ACTG"
    return $ SnpEntry chrom pos ref alt


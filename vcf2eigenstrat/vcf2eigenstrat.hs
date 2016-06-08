{-# LANGUAGE OverloadedStrings #-}

import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Applicative ((<|>))
import Control.Monad (forM_, void)
import Control.Monad.IO.Class (liftIO)
import Control.Monad.Trans.State.Strict (runStateT)
import qualified Data.Attoparsec.Text as A
import Data.Char (isSpace)
import Data.Monoid ((<>))
import qualified Data.Text as T
import qualified Data.Text.IO as T
-- import Debug.Trace (trace)
import qualified Options.Applicative as OP
import Pipes (Pipe, yield, (>->), runEffect, Producer, Pipe, for, cat, next)
import Pipes.Attoparsec (parsed, parse)
import qualified Pipes.Text.IO as PT
import System.IO (IOMode(..), withFile, Handle)
import Turtle.Format (format, d, s, (%))

data ProgOpt = ProgOpt FilePath (Maybe FilePath)-- outPrefix SnpFile

data VCFentry = VCFentry T.Text Int Char Char [Int] deriving (Show) -- chrom, pos, ref, alt, dosage
data VCFheader = VCFheader [T.Text] [T.Text] deriving (Show)-- simple comment lines, sample names

main :: IO ()
main = readOptions >>= runMain

readOptions :: IO ProgOpt
readOptions = OP.execParser parserInfo
  where
    parserInfo = OP.info (OP.helper <*> argParser)
                         (OP.progDesc "A program to convert a VCF file to Eigenstrat")

argParser :: OP.Parser ProgOpt
argParser = ProgOpt <$> parseOutPrefix
  where
    parseOutPrefix = OP.strOption (OP.long "outPrefix" <> OP.short 'e' <>
                                  OP.metavar "<FILE_PREFIX>" <>
                                  OP.help "specify the filenames for the EigenStrat SNP and IND \
                                  \file outputs: <FILE_PREFIX>.snp.txt and <FILE_PREFIX>.ind.txt")
    
runMain :: ProgOpt -> IO ()
runMain (ProgOpt outPrefix) = do
    (vcfHeader, vcfEntryProducer) <- readVCF
    let snpOut = outPrefix ++ ".snp.txt"
        indOut = outPrefix ++ ".ind.txt"
    liftIO . withFile indOut WriteMode $ \indOutHandle -> do
        let VCFheader _ sampleNames = vcfHeader
        forM_ sampleNames $ \n -> do
            T.hPutStrLn indOutHandle . T.intercalate "\t" $ [n, "U", "Unknown"]
    liftIO . withFile snpOut WriteMode $ \snpOutHandle -> do
        runEffect $ vcfEntryProducer >-> printEigenStrat snpOutHandle >-> printToStdOut
  where
    printToStdOut = for cat (liftIO . T.putStrLn)

readVCF :: IO (VCFheader, Producer VCFentry IO ())
readVCF = do
    (res, rest) <- runStateT (parse vcfHeaderParser) PT.stdin
    header <- case res of
        Nothing -> liftIO . throwIO $ AssertionFailed "vcf header not readible. VCF file empty?"
        Just (Left e_) -> do
            Right (chunk, _) <- next rest
            let msg = show e_ ++ T.unpack chunk
            liftIO . throwIO $ AssertionFailed ("VCF header parsing error: " ++ msg)
        Just (Right h) -> return h
    return (header, parsed vcfEntryParser rest >>= liftErrors)
  where
    liftErrors res = case res of
        Left (e_, prod_) -> do
            Right (chunk, _) <- liftIO $ next prod_
            let msg = show e_ ++ T.unpack chunk
            liftIO . throwIO $ AssertionFailed msg
        Right () -> return ()

vcfHeaderParser :: A.Parser VCFheader
vcfHeaderParser = VCFheader <$> A.many' doubleCommentLine <*> singleCommentLine
  where
    doubleCommentLine = do
        c1 <- A.string "##"
        s_ <- A.takeTill A.isEndOfLine <* A.endOfLine
        return $ T.append c1 s_
    singleCommentLine = do
        void $ A.char '#'
        s_ <- A.takeTill A.isEndOfLine <* A.endOfLine
        let fields = T.splitOn "\t" s_
        return . drop 9 $ fields

vcfEntryParser :: A.Parser VCFentry
vcfEntryParser = do
    chrom <- word
    void tab
    pos <- A.decimal
    void tab
    void word
    void tab
    ref <- allele
    void tab
    alt <- allele
    void tab
    void $ A.count 4 (word >> tab)
    genotypes <- genotype `A.sepBy1` tab
    void A.endOfLine
    return $ VCFentry chrom pos ref alt genotypes
  where
    word = A.takeTill isSpace
    allele = A.satisfy (\c -> c `elem` actg)
    actg :: String
    actg = "ACTG"
    tab = A.char '\t'

genotype :: A.Parser Int
genotype = do
    gen1 <- (A.char '0' <|> A.char '1' <|> A.char '.')
    void (A.char '/' <|> A.char '|')
    gen2 <- (A.char '0' <|> A.char '1' <|> A.char '.')
    _ <- A.takeTill (\a -> (a `elem` ['\r', '\t', '\n', ' ']))
    if   gen1 == '.' || gen2 == '.'
    then return (-1) 
    else return . length . filter (=='1') $ [gen1, gen2]

printEigenStrat :: Handle -> Pipe VCFentry T.Text IO r
printEigenStrat snpOutHandle = for cat $ \(VCFentry chrom pos ref alt dosages) -> do
    let n = format (s%"_"%d) chrom pos
        snpLine = format (s%"\t"%s%"\t0\t"%d%"\t"%s%"\t"%s) n chrom pos (T.singleton ref)
                            (T.singleton alt)
    liftIO . T.hPutStrLn snpOutHandle $ snpLine
    yield . T.concat . map (format d . toEigenStratNum) $ dosages
  where
    toEigenStratNum :: Int -> Int
    toEigenStratNum c = case c of
        0 -> 2
        1 -> 1
        2 -> 0
        -1 -> 9
        _ -> error ("unknown dosage " ++ show c)

{-# LANGUAGE OverloadedStrings #-}

module SeqTools.VCF (VCFheader(..),
                     VCFentry(..),
                     SimpleVCFentry(..),
                     readVCF,
                     vcfHeaderParser,
                     vcfEntryParser,
                     getGenotypes,
                     getDosages,
                     isTransversionSnp,
                     makeSimpleVCFentry,
                     isBiallelicSnp,
                     liftParsingErrors) where

import Control.Applicative ((<|>), empty)
import Control.Error (headErr)
import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad (void)
import Control.Monad.IO.Class (MonadIO, liftIO)
import Control.Monad.Trans.Class (lift)
import Control.Monad.Trans.State.Strict (runStateT)
import qualified Data.Attoparsec.Text as A
import Data.Char (isSpace)
import Data.Text (Text, count, unpack, append, splitOn)
import Pipes (Producer, next)
import Pipes.Attoparsec (parse, parsed, ParsingError(..))

data VCFheader = VCFheader {
    vcfHeaderComments :: [Text],
    vcfSampleNames :: [Text]
} deriving (Show)

data VCFentry = VCFentry {
    vcfChrom :: Text,
    vcfPos :: Int,
    vcfId :: Maybe Text,
    vcfRef :: Text,
    vcfAlt :: [Text],
    vcfQual :: Double,
    vcfFilter :: Maybe Text,
    vcfInfo :: [Text],
    vcfFormatString :: [Text],
    vcfGenotypeInfo :: [[Text]]
} deriving (Show)

data SimpleVCFentry = SimpleVCFentry {
    sVCFchrom :: Text,
    sVCFpos :: Int,
    sVCFref :: Text,
    sVCFalt :: [Text],
    sVCFdosages :: [Maybe Int]
} deriving (Show)

readVCF :: (MonadIO m) => Producer Text m () -> m (VCFheader, Producer VCFentry m ())
readVCF prod = do
    (res, rest) <- runStateT (parse vcfHeaderParser) prod
    header <- case res of
        Nothing -> liftIO . throwIO $ AssertionFailed "vcf header not readible. VCF file empty?"
        Just (Left e_) -> do
            Right (chunk, _) <- next rest
            let msg = show e_ ++ unpack chunk
            liftIO . throwIO $ AssertionFailed ("VCF header parsing error: " ++ msg)
        Just (Right h) -> return h
    return (header, parsed vcfEntryParser rest >>= liftParsingErrors)

liftParsingErrors :: (MonadIO m) => Either (ParsingError, Producer Text m r) () -> Producer a m ()
liftParsingErrors res = case res of
    Left (e_, prod_) -> do
        Right (chunk, _) <- lift $ next prod_
        let msg = show e_ ++ "\n" ++ unpack chunk
        lift . liftIO . throwIO $ AssertionFailed msg
    Right () -> return ()

vcfHeaderParser :: A.Parser VCFheader
vcfHeaderParser = VCFheader <$> A.many1' doubleCommentLine <*> singleCommentLine
  where
    doubleCommentLine = do
        c1 <- A.string "##"
        s_ <- A.takeWhile1 (not . A.isEndOfLine)
        A.endOfLine
        return $ append c1 s_
    singleCommentLine = do
        void $ A.char '#'
        s_ <- A.takeWhile1 (not . A.isEndOfLine)
        A.endOfLine
        let fields = splitOn "\t" s_
        return . drop 9 $ fields

vcfEntryParser :: A.Parser VCFentry
vcfEntryParser = VCFentry <$> word <* sp <*> A.decimal <* sp <*> parseId <* sp <*> word <* sp <*>
                              parseAlternativeAlleles <* sp <*> A.double <* sp <*>
                              parseFilter <* sp <*> parseInfoFields <* sp <*>
                              parseFormatStrings <* sp <*> parseGenotypeInfos <* A.endOfLine
  where
    word = A.takeTill A.isHorizontalSpace
    sp = A.satisfy A.isHorizontalSpace
    parseId = parseDot <|> (Just <$> word)
    parseDot = A.char '.' *> empty
    parseAlternativeAlleles = parseDot <|> (parseAllele `A.sepBy1` A.char ',')
    parseAllele = A.takeTill (\c -> c == ',' || A.isHorizontalSpace c)
    parseFilter = parseDot <|> (Just <$> word)
    parseInfoFields = parseDot <|> (parseInfoField `A.sepBy1` A.char ';')
    parseInfoField = A.takeTill (\c -> c == ';' || A.isHorizontalSpace c)
    parseFormatStrings = A.takeTill (\c -> c == ':' || A.isHorizontalSpace c) `A.sepBy1` A.char ':'
    parseGenotypeInfos = parseGenotype `A.sepBy1` (A.satisfy A.isHorizontalSpace)
    parseGenotype = parseGenoField `A.sepBy1` (A.char ':')
    parseGenoField = A.takeTill (\c -> c == ':' || isSpace c) 

isBiallelicSnp :: Text -> [Text] -> Bool
isBiallelicSnp ref alt = validRef && validAlt
  where
    validRef = (ref `elem` ["A", "C", "G", "T"])
    validAlt = case alt of
        [alt'] -> alt' `elem` ["A", "C", "G", "T"]
        _ -> False

isTransversionSnp :: Text -> [Text] -> Bool
isTransversionSnp ref alt =
    case alt of
        [alt'] -> isBiallelicSnp ref alt && (not $ isTransition ref alt')
        _ -> False
  where
    isTransition r a = ((r == "A") && (a == "G")) || ((r == "G") && (a == "A")) ||
                       ((r == "C") && (a == "T")) || ((r == "T") && (a == "C"))

getGenotypes :: VCFentry -> Either String [Text]
getGenotypes vcfEntry = do
    gtIndex <- fmap fst . headErr "GT format field not found" . filter ((=="GT") . snd) .
               zip [0..] . vcfFormatString $ vcfEntry
    return $ map (!!gtIndex) (vcfGenotypeInfo vcfEntry)

getDosages :: VCFentry -> Either String [Maybe Int]
getDosages vcfEntry = do
    genotypes <- getGenotypes vcfEntry
    let dosages = do
            gen <- genotypes
            case gen of
                "." -> return Nothing
                genString -> return . Just $ count "1" genString
    return dosages

makeSimpleVCFentry :: VCFentry -> Either String SimpleVCFentry
makeSimpleVCFentry vcfEntry = do
    dosages <- getDosages vcfEntry
    return $ SimpleVCFentry (vcfChrom vcfEntry) (vcfPos vcfEntry) (vcfRef vcfEntry)
                            (vcfAlt vcfEntry) dosages
    
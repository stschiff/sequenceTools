{-# LANGUAGE OverloadedStrings #-}

module SeqTools.Fasta (readFastaEntry) where

import Control.Monad.IO.Class (liftIO, MonadIO)
import Control.Monad.Trans.State.Strict (runStateT)
import qualified Data.Attoparsec.ByteString.Char8 as A
import qualified Data.ByteString.Char8 as B
import Data.Char (isAlphaNum)
import Pipes.Attoparsec (parse, parsed, ParsingError(..))
import Pipes.ByteString (fromHandle)
import Pipes.Safe (MonadSafe)
import System.IO (withFile, IOMode(..))
import Turtle.Prelude (err)

-- loadFasta :: (MonadIO m) => FilePath -> String -> m B.ByteString
-- loadFasta refFilePath chrom = do
--     err $ format ("loading reference sequence from "%s) (pack refFilePath)
--     withFile refFilePath ReadMode $ \fh -> do
--         let prod = fromHandle fh
--         (chrom, seqProd) <- readFasta

readFastaEntry :: (MonadIO m) => Producer B.ByteString m () -> m (String, Producer B.ByteString m (Producer B.ByString m ()))
readFastaEntry prod = do
    (res, rest) <- runStateT (parse fastaHeaderLineParser) prod
    header <- case res of
        Nothing -> liftIO . throwIO $ AssertionFailed "Fasta header not readible"
        Just (Left e_) -> do
            Right (chunk, _) <- next rest
            let msg = show e_ ++ unpack chunk
            liftIO . throwIO $ AssertionFailed ("Fasta header parsing error: " ++ msg)
        Just (Right h) -> return h
    return (header, parsed fastaBodyParser rest >>= liftParsingErrors)

fastaHeaderLineParser :: A.Parser String
fastaHeaderLineParser = do
    _ <- A.char '>'
    chrom <- A.takeWhile isAlphaNum
    A.skipSpace
    A.skipWhile (not . A.isEndOfLine)
    A.endOfLine
    return . B.unpack $ chrom

fastaBodyParser :: A.Parser B.ByteString
fastaBodyParser = A.takeTill (=='>')

liftParsingErrors :: (MonadIO m) => Either (ParsingError, Producer Text m r) r' -> Producer a m r'
liftParsingErrors res = case res of
    Left (e_, prod_) -> do
        Right (chunk, _) <- lift $ next prod_
        let msg = show e_ ++ "\n" ++ unpack chunk
        lift . liftIO . throwIO $ AssertionFailed msg
    Right r -> return r

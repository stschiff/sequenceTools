{-# LANGUAGE OverloadedStrings #-}

module SeqTools.Fasta (readNextFastaEntry, loadFastaChrom) where

import Control.Exception.Base (throwIO, AssertionFailed(..))
import Control.Monad (void)
import Control.Monad.IO.Class (liftIO, MonadIO)
import Control.Monad.Trans.State.Strict (runStateT)
import qualified Data.Attoparsec.ByteString.Char8 as A
import qualified Data.ByteString.Char8 as B
import Data.Char (isAlphaNum)
import Data.Text (pack)
import Lens.Family2 (view)
import Pipes (Producer, next, (>->), runEffect)
import Pipes.Attoparsec (parse)
import qualified Pipes.ByteString as P
import Pipes.Prelude (drain)
import System.IO (Handle)
import Turtle.Format (format, (%), s)
import Turtle.Prelude (err)

loadFastaChrom :: Handle -> String -> IO (Producer B.ByteString IO ())
loadFastaChrom refFileHandle chrom = do
    let prod = P.fromHandle refFileHandle
    go prod
  where
    go prod = do
        (chrom_, prod') <- readNextFastaEntry prod
        err $ format ("found chromosome "%s) (pack chrom_)
        if chrom_ == chrom
        then return (void prod')
        else do
            newProd <- runEffect $ prod' >-> drain
            go newProd

readNextFastaEntry :: (MonadIO m) => Producer B.ByteString m () ->
                      m (String, Producer B.ByteString m (Producer B.ByteString m ()))
readNextFastaEntry prod = do
    (res, rest) <- runStateT (parse fastaHeaderLineParser) prod
    header <- case res of
        Nothing -> liftIO . throwIO $ AssertionFailed "Could not find chromosome. Fasta file exhausted."
        Just (Left e_) -> do
            Right (chunk, _) <- next rest
            let msg = show e_ ++ B.unpack chunk
            liftIO . throwIO $ AssertionFailed ("Fasta header parsing error: " ++ msg)
        Just (Right h) -> return h
    return (header, view (P.break (==62)) rest >-> P.filter (\c -> c /= 10 && c /= 13))
-- '>' == 62, '\n' == 10, \r == 13

fastaHeaderLineParser :: A.Parser String
fastaHeaderLineParser = do
    _ <- A.char '>'
    chrom <- A.takeWhile isAlphaNum
    A.skipSpace
    A.skipWhile (\c -> c /= '\n' && c /= '\r')
    A.endOfLine
    return . B.unpack $ chrom


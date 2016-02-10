{-# LANGUAGE OverloadedStrings #-}

import Lib (searchReadEndDistance, mergeReads)

import Codec.Compression.GZip (compress)
import Control.Monad (void, when)
import Control.Monad.IO.Class (MonadIO)
import Control.Monad.Managed (runManaged, managed)
import Control.Monad.Trans.Class (lift)
import qualified Data.Attoparsec.ByteString.Char8 as A
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as B8
import qualified Data.ByteString.Lazy.Char8 as BL8
import Data.IORef (IORef, newIORef, readIORef, modifyIORef)
import Options.Applicative (execParser, info, fullDesc, progDesc, auto, strOption, value, option, long, short, metavar, help, showDefault, (<>))
import qualified Pipes.ByteString as PB
import Pipes (runEffect, (>->), for, liftIO, Effect, next)
import Pipes.Attoparsec (parsed)
import Pipes.GZip (decompress)
import qualified Pipes.Prelude as P
import System.IO (withFile, IOMode(..), stderr, hPutStrLn, Handle)

data MyOptions = MyOptions FilePath Int Double Int (Maybe Int) FilePath FilePath
data FastqEntry = FastqEntry B8.ByteString B8.ByteString B8.ByteString deriving (Show)

main :: IO ()
main = execParser (info parser (fullDesc <> progDesc "merge paired sequencing reads and clip any adapter sequence left or right of the merged fragment.")) >>= runWithOptions
  where
    parser = MyOptions <$> strOption (long "outPrefix" <> short 'o' <> metavar "<OUT_PREFIX>" <> help "Prefix for the output filename. Four \
                                      \files are output: <prefix>_1.fastq contains the forward reads that could not be merged. \
                                      \<prefix>_2.fastq contains the backward reads that could not be merged. \
                                      \<prefix>_merged.fastq contains the merged reads")
                       <*> option auto (long "minLength" <> short 'l' <> value 30 <> showDefault <> metavar "<MIN_LENGTH>" <>
                                        help "minimum length for merged reads. If shorter, do not merge")
                       <*> option auto (long "mismatchRate" <> short 'e' <> value 0.05 <> showDefault <> metavar "<MISMATCH_ERRORRATE>" <>
                                        help "fraction of base pairs allowed to differ between both reads in order to merge.")
                       <*> option auto (long "minOverlap" <> short 'm' <> value 20 <> showDefault <> metavar "<MIN_OVERLAP>" <>
                                        help "minimum overlap length for merged reads. If shorter, do not merge")
                       <*> option (Just <$> auto) (long "maxNr" <> value Nothing <> metavar "nr of reads after which to stop. For debugging purpose only")
                       <*> strOption (long "in1" <> short '1' <> metavar "<READ1_FASTQ>" <> help "fastq file containing read 1. Can be gzipped (with .gz file ending)")
                       <*> strOption (long "in2" <> short '2' <> metavar "<READ2_FASTQ>" <> help "fastq file containing read 2. Can be gzipped (with .gz file ending)")       

runWithOptions :: MyOptions -> IO ()
runWithOptions (MyOptions outPrefix minLength mismatchRate minOverlapSize maxNr read1FN read2FN) = runManaged $ do
    read1H <- managed $ withFile read1FN ReadMode
    read2H <- managed $ withFile read2FN ReadMode
    outR1H <- managed $ withFile (outPrefix ++ "_1.fastq.gz") WriteMode
    outR2H <- managed $ withFile (outPrefix ++ "_2.fastq.gz") WriteMode
    outMergedH <- managed $ withFile (outPrefix ++ "_merged.fastq.gz") WriteMode
    let read1Stream = parsed fastqParser . decompress $ PB.fromHandle read1H
        read2Stream = parsed fastqParser . decompress $ PB.fromHandle read2H
        combinedStream = P.zip read1Stream read2Stream
    statsRef <- liftIO . newIORef $ (0, 0)
    res <- runEffect $ for combinedStream (processFastq statsRef minLength mismatchRate minOverlapSize maxNr outR1H outR2H outMergedH)
    case res of
        Left (err, restProd) -> do
            liftIO $ hPutStrLn stderr ("Parsing error: " ++ show err)
            Right (chunk, _) <- next restProd
            liftIO $ hPutStrLn stderr (B8.unpack chunk)
        Right _ -> do
            (nrReads, nrMerged) <- liftIO . readIORef $ statsRef
            liftIO . putStrLn $ "Total Reads processed: " ++ show nrReads
            liftIO . putStrLn $ "Merged Reads: " ++ show nrMerged

fastqParser :: A.Parser FastqEntry
fastqParser = FastqEntry <$> line <* A.endOfLine <*> seq_ <* A.endOfLine <* A.char '+' <* A.endOfLine <*> line <* A.endOfLine 
  where
    line = A.takeTill (\c -> c == '\n' || c == '\r')
    seq_ = A.takeWhile1 (\c -> c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N')
    
processFastq :: (MonadIO m) => IORef (Int, Int) -> Int -> Double -> Int -> Maybe Int -> Handle -> Handle -> Handle -> (FastqEntry, FastqEntry) -> Effect m ()
processFastq statsRef minLength mismatchRate minOverlapSize maxNr outR1H outR2H outMergedH fqEntries = do
    liftIO $ modifyIORef statsRef (\(r, m) -> (r + 1, m))
    (nrReads, nrMerged) <- liftIO $ readIORef statsRef
    let stop = case maxNr of
            Just n -> nrReads > n
            Nothing -> False
    if stop then
        return ()
    else do
        when (nrReads `mod` 10000 == 0) $ (liftIO . hPutStrLn stderr) ("processing read " ++ show nrReads)
        let (FastqEntry header1 seq1 qual1, FastqEntry header2 seq2 qual2) = fqEntries
            readEndDistance = searchReadEndDistance seq1 seq2 mismatchRate minOverlapSize
        case readEndDistance of
            Just dist -> do
                let (mergedSeq, mergedQual) = mergeReads seq1 qual1 seq2 qual2 dist
                if (B8.length mergedSeq >= minLength) then do
                    liftIO $ modifyIORef statsRef (\(r, m) -> (r, m + 1))
                    liftIO . writeFQ outMergedH $ FastqEntry header1 mergedSeq mergedQual
                else do
                    liftIO . writeFQ outR1H $ FastqEntry header1 seq1 qual1
                    liftIO . writeFQ outR2H $ FastqEntry header2 seq2 qual2
            Nothing -> do
                liftIO . writeFQ outR1H $ FastqEntry header1 seq1 qual1
                liftIO . writeFQ outR2H $ FastqEntry header2 seq2 qual2
                    
writeFQ :: Handle -> FastqEntry -> IO ()
writeFQ outH (FastqEntry header seq_ qual) =
    mapM_ (BL8.hPutStr outH . compress . BL8.fromStrict . flip B8.snoc '\n') [header, seq_, "+", qual]

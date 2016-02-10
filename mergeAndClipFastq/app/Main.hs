import Lib

import Control.Monad (void)
import qualified Data.Attoparsec.Text as A
import Data.Text (Text)
import Options.Applicative (execParser, info, fullDesc, progDesc, strArgument, metavar, help, (<>))
import qualified Pipes.ByteString as PB
import Pipes.GZip (decompress)
import Pipes (runEffect, (>->), for, liftIO)
import Pipes.Attoparsec (parsed)
import qualified Pipes.Prelude as P
import Pipes.Text.Encoding (decodeUtf8)
import System.IO (withFile, IOMode(..))

data MyOptions = MyOptions FilePath FilePath

main :: IO ()
main = execParser (info parser (fullDesc <> progDesc "merge paired sequencing reads and clip adapters and barcodes.")) >>= runWithOptions
  where
    parser = MyOptions <$> strArgument (metavar "<READ1_FASTQ>" <> help "fastq file containing read 1. Can be gzipped (with .gz file ending)")
                       <*> strArgument (metavar "<READ2_FASTQ>" <> help "fastq file containing read 2. Can be gzipped (with .gz file ending)")

runWithOptions :: MyOptions -> IO ()
runWithOptions (MyOptions read1FN read2FN) = do
    withFile read1FN ReadMode $ \read1H -> do
        withFile read2FN ReadMode $ \read2H -> do
            let read1Stream = parsed fastqParser . decodeUtf8 . decompress $ PB.fromHandle read1H
                read2Stream = parsed fastqParser . decodeUtf8 . decompress $ PB.fromHandle read2H
                combinedStream = P.zip read1Stream read2Stream
            runEffect $ void combinedStream >-> P.take 10 >-> P.print

data FastqEntry = FastqEntry Text Text Text deriving (Show)

fastqParser :: A.Parser FastqEntry
fastqParser = FastqEntry <$> line <* A.endOfLine <*> seq_ <* A.endOfLine <* A.char '+' <* A.endOfLine <*> line <* A.endOfLine 
  where
    line = A.takeTill A.isEndOfLine
    seq_ = A.takeWhile1 (\c -> c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N')
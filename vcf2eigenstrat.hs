import Pipes (Producer, Pipe, yield, await, (>->), runEffect, lift)
import qualified System.IO as IO
import qualified Pipes.Prelude as P
import System.Environment (getArgs)
import Data.List.Split (splitOn)
import Data.List (intercalate)
import Control.Monad (liftM)

main = do
    snpFileName <- liftM head getArgs
    snpFileH <- IO.openFile snpFileName IO.WriteMode
    runEffect $ input >-> mainLoop snpFileH >-> P.stdoutLn
    IO.hClose snpFileH
    
input :: Producer String IO ()
input = P.stdinLn >-> P.filter (\l -> head l /= '#')

mainLoop :: IO.Handle -> Pipe String String IO ()
mainLoop snpFileH = do
    vcfString <- await
    let vcfFields = splitOn "\t" vcfString
    lift $ writeSnpLine snpFileH vcfFields
    let genFields = drop 9 vcfFields
    yield $ map getGenotype genFields
    mainLoop snpFileH

writeSnpLine :: IO.Handle -> [String] -> IO ()
writeSnpLine snpFileH fields = do
    let (chrom:pos:_:ref:alt:_) = fields
    let snpName = chrom ++ "_" ++ pos
    IO.hPutStrLn snpFileH (intercalate "\t" [snpName, chrom, "0.0", pos, ref, alt])

getGenotype :: String -> Char
getGenotype genField = case take 3 genField of
    ['0', _, '0'] -> '0'
    ['1', _, '0'] -> '1'
    ['0', _, '1'] -> '1'
    ['1', _, '1'] -> '2'
    _     -> '9'

    
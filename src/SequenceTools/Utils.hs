module SequenceTools.Utils (versionInfoOpt, versionInfoText, sampleWithoutReplacement) where 

import Paths_sequenceTools (version)
import Data.Version (showVersion)
import qualified Options.Applicative as OP
import System.Random (randomRIO)

versionInfoOpt :: OP.Parser (a -> a)
versionInfoOpt = OP.infoOption (showVersion version) (OP.long "version" <> OP.help "Print version and exit")

versionInfoText :: String
versionInfoText = "This tool is part of sequenceTools version " ++ showVersion version

sampleWithoutReplacement :: [a] -> Int -> IO (Maybe [a])
sampleWithoutReplacement = go []
  where
    go res _ 0 = return $ Just res
    go res xs n
        | n > length xs = return Nothing
        | n == length xs = return $ Just (xs ++ res)
        | otherwise = do
                rn <- randomRIO (0, length xs - 1)
                let a = xs !! rn
                    xs' = remove rn xs
                go (a:res) xs' (n - 1)
    remove i xs = let (ys, zs) = splitAt i xs in ys ++ tail zs

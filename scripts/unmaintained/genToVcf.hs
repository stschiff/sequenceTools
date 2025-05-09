import           Control.Monad      (liftM)
import           Data.List          (intercalate)
import           Data.List.Split    (chunksOf, splitOn)
import           System.Environment (getArgs)

convertLine :: String -> String -> String
convertLine chr line =
    let (_:_:pos:ref:alt:gens) = splitOn " " line
        genFields = map makeGenField $ chunksOf 3 (map read gens)
    in  intercalate "\t" $ [chr, pos, ".", ref, alt, "100", ".", ".", "GT:GP"] ++ genFields

makeGenField :: [Double] -> String
makeGenField [p1, p2, p3] =
    if p1 > p2 && p1 > p3 then "0/0:" ++ probStr else
        if p2 > p1 && p2 > p3 then "0/1:" ++ probStr else "1/1:" ++ probStr
    where probStr = intercalate "," $ map show [p1, p2, p3]

main = do
    chr <- liftM head getArgs
    interact $ unlines . map (convertLine chr) . lines

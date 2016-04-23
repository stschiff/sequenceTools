#!/usr/bin/env stack
-- stack runghc --package turtle --package text

{-# LANGUAGE OverloadedStrings #-}

import qualified Data.Text as T
import qualified Data.Text.IO as T
import Turtle
import qualified System.IO as IO

main :: IO ()
main = do
    snpFileName <- options "script to convert a vcf to eigenstrat format"
                   (argPath "snpFileName" "name of the SNP file")
    IO.withFile (T.unpack . format fp $ snpFileName) IO.WriteMode $ \snpFileH ->
        foldIO (grep (prefix (notChar '#')) stdin) (mainFold snpFileH)

mainFold :: IO.Handle -> FoldM IO Text ()
mainFold snpFileH = FoldM (step snpFileH) (return ()) (const (return ()))

step :: IO.Handle -> () -> Text -> IO ()
step snpFileH _ line = do
    let vcfFields = T.splitOn "\t" line
    writeSnpLine snpFileH vcfFields
    let genFields = drop 9 vcfFields
    putStrLn $ map getGenotype genFields

writeSnpLine :: IO.Handle -> [Text] -> IO ()
writeSnpLine snpFileH fields = do
    let (chrom:pos:_:ref:alt:_) = fields
        snpName = format (s%"_"%s) chrom pos
    T.hPutStrLn snpFileH (T.intercalate "\t" [snpName, chrom, "0", pos, ref, alt])

getGenotype :: Text -> Char
getGenotype genField = case T.unpack (T.take 3 genField) of
    ['0', _, '0'] -> '2'
    ['1', _, '0'] -> '1'
    ['0', _, '1'] -> '1'
    ['1', _, '1'] -> '0'
    _     -> '9'

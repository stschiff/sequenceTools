#!/usr/bin/env stack
-- stack runghc --package turtle

{-# LANGUAGE OverloadedStrings #-}

import Control.Applicative (optional)
import Prelude hiding (FilePath)
import Turtle

data Options = Options {
    optGeno :: FilePath,
    optSnp :: FilePath,
    optInd :: FilePath,
    optPopList :: FilePath,
    optLow :: Maybe Int,
    optHigh :: Maybe Int
}

main = do
    args <- options "Admixtools qpDstat wrapper" parser
    runManaged $ do
        paramFile <- mktempfile "." "qpDstat_wrapper"
        let content = return (format ("genotypename:\t"%fp) (optGeno args)) <|>
                      return (format ("snpname:\t"%fp) (optSnp args)) <|>
                      return (format ("indivname:\t"%fp) (optInd args)) <|>
                      return (format ("popfilename:\t"%fp) (optPopList args))
        output paramFile content
        let execParams = ["-p", format fp paramFile] ++
                         maybe [] (\low -> ["-l", format d low]) (optLow args) ++
                         maybe [] (\high -> ["-h", format d high]) (optHigh args)
        ec <- proc "qpDstat" execParams empty
        case ec of
            ExitSuccess -> return ()
            ExitFailure n -> err $ format ("qpDstat failed with exit code "%d) n

parser :: Parser Options
parser = Options <$> optPath "geno" 'g' "Genotype File"
                 <*> optPath "snp" 's' "Snp File"
                 <*> optPath "ind" 'i' "Ind File"
                 <*> optPath "popList" 'p' "give a list with all population triples"
                 <*> optional (optInt "lower" 'l' "analyse population quadruples from this line in the popList")
                 <*> optional (optInt "upper" 'u' "analyse population quadruples up to this line in the popList")

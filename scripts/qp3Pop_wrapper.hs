#!/usr/bin/env stack
-- stack --resolver lts-6.4 --install-ghc runghc --package turtle 

{-# LANGUAGE OverloadedStrings #-}

import Control.Applicative (optional)
import Prelude hiding (FilePath)
import Turtle

data Options = Options {
    optGeno :: FilePath,
    optSnp :: FilePath,
    optInd :: FilePath,
    optPopList :: FilePath
}

main = do
    args <- options "Admixtools qp3Pop wrapper" parser
    runManaged $ do
        paramFile <- mktempfile "/tmp" "qp3Pop_wrapper"
        let content = return (format ("genotypename:\t"%fp) (optGeno args)) <|>
                      return (format ("snpname:\t"%fp) (optSnp args)) <|>
                      return (format ("indivname:\t"%fp) (optInd args)) <|>
                      return (format ("popfilename:\t"%fp) (optPopList args))
        output paramFile content
        ec <- proc "qp3Pop" ["-p", format fp paramFile] empty
        case ec of
            ExitSuccess -> return ()
            ExitFailure n -> err $ format ("qp3Pop failed with exit code "%d) n

parser :: Parser Options
parser = Options <$> optPath "geno" 'g' "Genotype File"
                 <*> optPath "snp" 's' "Snp File"
                 <*> optPath "ind" 'i' "Ind File"
                 <*> optPath "popList" 'p' "give a list with all population triples"

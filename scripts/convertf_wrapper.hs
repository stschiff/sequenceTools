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
    optFormat :: FilePath,
    optPoplist :: FilePath,
    optOutGeno :: FilePath,
    optOutSnp :: FilePath,
    optOutInd :: FilePath
}

main = do
    args <- options "Eigensoft convertf wrapper" parser
    runManaged $ do
        paramFile <- mktempfile "." "convert_wrapper"
        let content = return (format ("genotypename:\t"%fp) (optGeno args)) <|>
                      return (format ("snpname:\t"%fp) (optSnp args)) <|>
                      return (format ("indivname:\t"%fp) (optInd args)) <|>
                      return (format ("outputformat:\t"%fp) (optFormat args)) <|>
                      return (format ("poplistname:\t"%fp) (optPoplist args)) <|>
                      return (format ("genotypeoutname:\t"%fp) (optOutGeno args)) <|>
                      return (format ("snpoutname:\t"%fp) (optOutSnp args)) <|>
                      return (format ("indivoutname:\t"%fp) (optOutInd args))
        output paramFile content
        ec <- proc "convertf" ["-p", format fp paramFile] empty
        case ec of
            ExitSuccess -> return ()
            ExitFailure n -> err $ format ("convertf failed with exit code "%d) n

parser :: Parser Options
parser = Options <$> optPath "geno" 'g' "Genotype File"
                 <*> optPath "snp" 's' "Snp File"
                 <*> optPath "ind" 'i' "Ind File"
                 <*> optPath "outFormat" 'f' "output format"
                 <*> optPath "popList" 'p' "population list"
                 <*> optPath "outGeno" 'G' "Output Genotype File"
                 <*> optPath "outSnp" 'S' "Output Snp File"
                 <*> optPath "outInd" 'I' "Output Ind File"

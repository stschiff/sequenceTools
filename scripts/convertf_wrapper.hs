#!/usr/bin/env stack
-- stack script --resolver lts-14.1 --package turtle
{-# LANGUAGE OverloadedStrings #-}

import           Control.Applicative (optional)
import           Prelude             hiding (FilePath)
import           Turtle

data Options = Options {
    optGeno    :: FilePath,
    optSnp     :: FilePath,
    optInd     :: FilePath,
    optFormat  :: FilePath,
    optPoplist :: Maybe FilePath,
    optOutGeno :: FilePath,
    optOutSnp  :: FilePath,
    optOutInd  :: FilePath
}

main = do
    args <- options "Eigensoft convertf wrapper" parser
    runManaged $ do
        paramFile <- mktempfile "." "convert_wrapper"
        let popListRow = case optPoplist args of
                Just popList -> return . unsafeTextToLine $ format ("poplistname:\t"%fp) popList
                Nothing -> empty
        let content = [(format ("genotypename:\t"%fp) (optGeno args)),
                       (format ("snpname:\t"%fp) (optSnp args)),
                       (format ("indivname:\t"%fp) (optInd args)),
                       (format ("outputformat:\t"%fp) (optFormat args)),
                       (format ("genotypeoutname:\t"%fp) (optOutGeno args)),
                       (format ("snpoutname:\t"%fp) (optOutSnp args)),
                       (format ("indivoutname:\t"%fp) (optOutInd args))]
        output paramFile $ select (map unsafeTextToLine content) <|> popListRow
        ec <- proc "convertf" ["-p", format fp paramFile] empty
        case ec of
            ExitSuccess -> return ()
            ExitFailure n -> err . unsafeTextToLine $ format ("convertf failed with exit code "%d) n

parser :: Parser Options
parser = Options <$> optPath "geno" 'g' "Genotype File"
                 <*> optPath "snp" 's' "Snp File"
                 <*> optPath "ind" 'i' "Ind File"
                 <*> optPath "outFormat" 'f' "output format"
                 <*> optional (optPath "popList" 'p' "population list")
                 <*> optPath "outGeno" 'G' "Output Genotype File"
                 <*> optPath "outSnp" 'S' "Output Snp File"
                 <*> optPath "outInd" 'I' "Output Ind File"

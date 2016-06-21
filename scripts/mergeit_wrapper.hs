#!/usr/bin/env stack
-- stack runghc --package turtle

{-# LANGUAGE OverloadedStrings #-}

import Control.Applicative (optional)
import Prelude hiding (FilePath)
import Turtle

data Options = Options {
    optGeno1 :: FilePath,
    optSnp1 :: FilePath,
    optInd1 :: FilePath,
    optGeno2 :: FilePath,
    optSnp2 :: FilePath,
    optInd2 :: FilePath,
    optOutPrefix :: FilePath,
    optOutFormat :: Maybe OutFormat
}

data OutFormat = ANCESTRYMAP | EIGENSTRAT | PED | PACKEDPED | PACKEDANCESTRYMAP
                 deriving (Show, Read)

main = do
    args <- options "Eigensoft mergeit wrapper" parser
    runManaged $ do
        paramFile <- mktempfile "/tmp" "mergeit_wrapper"
        let content = return (format ("geno1:\t"%fp) (optGeno1 args)) <|>
                      return (format ("snp1:\t"%fp) (optSnp1 args)) <|>
                      return (format ("ind1:\t"%fp) (optInd1 args)) <|>
                      return (format ("geno2:\t"%fp) (optGeno2 args)) <|>
                      return (format ("snp2:\t"%fp) (optSnp2 args)) <|>
                      return (format ("ind2:\t"%fp) (optInd2 args)) <|>
                      return (format ("genooutfilename:\t"%fp%".geno") (optOutPrefix args)) <|>
                      return (format ("snpoutfilename:\t"%fp%".snp") (optOutPrefix args)) <|>
                      return (format ("indoutfilename:\t"%fp%".ind") (optOutPrefix args))
        let content' = case optOutFormat args of
                Just outFormat -> content <|> return (format ("outputformat:\t"%w) outFormat)
                Nothing -> content
        output paramFile content'
        ec <- proc "mergeit" ["-p", format fp paramFile] empty
        case ec of
            ExitSuccess -> return ()
            ExitFailure n -> err $ format ("mergeit failed with exit code "%d) n

parser :: Parser Options
parser = Options <$> optPath "geno1" 'g' "First Genotype File"
                 <*> optPath "snp1" 's' "First Snp File"
                 <*> optPath "ind1" 'i' "First Ind File"
                 <*> optPath "geno2" 'G' "Second Genotype File"
                 <*> optPath "snp2" 'S' "Second Snp File"
                 <*> optPath "ind2" 'I' "Second Ind File"
                 <*> optPath "outPrefix" 'o' "Output prefix for *.geno, *.snp and *.ind \
                                              \output files"
                 <*> optional (optRead "outFormat" 'f' "Output format. One of ANCESTRYMAP, \
                               \EIGENSTRAT, PED, PACKEDPED, PACKEDANCESTRYMAP. Default is \
                               \PACKEDANCESTRYMAP")

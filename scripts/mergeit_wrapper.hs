#!/usr/bin/env stack
-- stack script --resolver lts-14.1 --package turtle
{-# LANGUAGE OverloadedStrings #-}

import           Control.Applicative (optional)
import           Prelude             hiding (FilePath)
import           Turtle

data Options = Options {
    optGeno1       :: FilePath,
    optSnp1        :: FilePath,
    optInd1        :: FilePath,
    optGeno2       :: FilePath,
    optSnp2        :: FilePath,
    optInd2        :: FilePath,
    optOutPrefix   :: FilePath,
    optAllowDups   :: Bool,
    optStrandCheck :: Bool,
    optOutFormat   :: Maybe OutFormat
}

data OutFormat = ANCESTRYMAP | EIGENSTRAT | PED | PACKEDPED | PACKEDANCESTRYMAP
                 deriving (Show, Read)

main = do
    args <- options "Eigensoft mergeit wrapper" parser
    runManaged $ do
        paramFile <- mktempfile "." "mergeit_wrapper"
        let content = [(format ("geno1:\t"%fp) (optGeno1 args)),
                       (format ("snp1:\t"%fp) (optSnp1 args)),
                       (format ("ind1:\t"%fp) (optInd1 args)),
                       (format ("geno2:\t"%fp) (optGeno2 args)),
                       (format ("snp2:\t"%fp) (optSnp2 args)),
                       (format ("ind2:\t"%fp) (optInd2 args)),
                       (format ("allowdups:\t"%s)
                               (if optAllowDups args then "YES" else "NO")),
                       (format ("strandcheck:\t"%s)
                               (if optStrandCheck args then "YES" else "NO")),
                       (format ("genooutfilename:\t"%fp%".geno") (optOutPrefix args)),
                       (format ("snpoutfilename:\t"%fp%".snp") (optOutPrefix args)),
                       (format ("indoutfilename:\t"%fp%".ind") (optOutPrefix args))]
        let outputFormatLine = case optOutFormat args of
                Just outFormat -> return . unsafeTextToLine $ format ("outputformat:\t"%w) outFormat
                Nothing -> empty
        output paramFile $ select (map unsafeTextToLine content) <|> outputFormatLine
        ec <- proc "mergeit" ["-p", format fp paramFile] empty
        case ec of
            ExitSuccess -> return ()
            ExitFailure n -> err . unsafeTextToLine $ format ("mergeit failed with exit code "%d) n

parser :: Parser Options
parser = Options <$> optPath "geno1" 'g' "First Genotype File"
                 <*> optPath "snp1" 's' "First Snp File"
                 <*> optPath "ind1" 'i' "First Ind File"
                 <*> optPath "geno2" 'G' "Second Genotype File"
                 <*> optPath "snp2" 'S' "Second Snp File"
                 <*> optPath "ind2" 'I' "Second Ind File"
                 <*> optPath "outPrefix" 'o' "Output prefix for *.geno, *.snp and *.ind \
                                              \output files"
                 <*> switch "allowDups" 'd' "Allow duplicates, leading for any duplicate \
                                             \individual in the second data set to be ignored"
                 <*> switch "strandcheck" 'c' "Check for strand misalignment. Warning: If set, \
                                               \removes all A/T and C/G SNPs"
                 <*> optional (optRead "outFormat" 'f' "Output format. One of ANCESTRYMAP, \
                               \EIGENSTRAT, PED, PACKEDPED, PACKEDANCESTRYMAP. Default is \
                               \PACKEDANCESTRYMAP")

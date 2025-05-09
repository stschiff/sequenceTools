#!/usr/bin/env stack
-- stack script --resolver lts-14.1 --package turtle
{-# LANGUAGE OverloadedStrings #-}

import           Control.Applicative (optional)
import           Prelude             hiding (FilePath)
import           Turtle

data Options = Options {
    optGeno       :: FilePath,
    optSnp        :: FilePath,
    optInd        :: FilePath,
    optOutPrefix  :: FilePath,
    optLSQproject :: Bool,
    optPopList    :: Maybe FilePath
}

main = do
    args <- options "Eigensoft smartpca wrapper" parser
    runManaged $ do
        paramFile <- mktempfile "." "smartpca_wrapper"
        let content        = [(format ("genotypename:\t"%fp) (optGeno args)),
                              (format ("snpname:\t"%fp) (optSnp args)),
                              (format ("indivname:\t"%fp) (optInd args)),
                              (format ("evecoutname:\t"%fp%".evec.txt") (optOutPrefix args)),
                              (format ("evaloutname:\t"%fp%".eval.txt") (optOutPrefix args))]
            lsqProjectLine = if (optLSQproject args) then return "lsqproject:\tYES" else empty
            popListLine    = case optPopList args of
                                Just popList -> return . unsafeTextToLine $ format ("poplistname:\t"%fp) popList
                                Nothing -> empty
        output paramFile $ select (map unsafeTextToLine content) <|> lsqProjectLine <|> popListLine
        ec <- proc "smartpca" ["-p", format fp paramFile] empty
        case ec of
            ExitSuccess -> return ()
            ExitFailure n -> err . unsafeTextToLine $ format ("mergeit failed with exit code "%d) n

parser :: Parser Options
parser = Options <$> optPath "geno" 'g' "Genotype File"
                 <*> optPath "snp" 's' "Snp File"
                 <*> optPath "ind" 'i' "Ind File"
                 <*> optPath "outPrefix" 'o' "Output prefix for *.evec.txt and *.eval.txt output \
                                              \files"
                 <*> switch "lsqProject" 'l' "set lsqproject option to YES"
                 <*> optional (optPath "popList" 'p' "give poplist file to restrict PCA to \
                                                      \populations listed.")

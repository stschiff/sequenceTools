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
    optOutPrefix :: FilePath,
    optLSQproject :: Bool,
    optPopList :: Maybe FilePath
}

main = do
    args <- options "Eigensoft mergeit wrapper" parser
    runManaged $ do
        paramFile <- mktempfile "." "mergeit_wrapper"
        let content = return (format ("genotypename:\t"%fp) (optGeno args)) <|>
                      return (format ("snpname:\t"%fp) (optSnp args)) <|>
                      return (format ("indivname:\t"%fp) (optInd args)) <|>
                      return (format ("evecoutname:\t"%fp%".evec.txt") (optOutPrefix args)) <|>
                      return (format ("evaloutname:\t"%fp%".eval.txt") (optOutPrefix args)) <|>
                      (if (optLSQproject args) then "lsqproject:\tYES" else "") <|>
                      case optPopList args of
                          Just popList -> return (format ("poplistname:\t"%fp) popList)
                          Nothing -> ""
        output paramFile content
        ec <- proc "smartpca" ["-p", format fp paramFile] empty
        case ec of
            ExitSuccess -> return ()
            ExitFailure n -> err $ format ("mergeit failed with exit code "%d) n

parser :: Parser Options
parser = Options <$> optPath "geno" 'g' "Genotype File"
                 <*> optPath "snp" 's' "Snp File"
                 <*> optPath "ind" 'i' "Ind File"
                 <*> optPath "outPrefix" 'o' "Output prefix for *.evec.txt and *.eval.txt output \  
                                              \files"
                 <*> switch "lsqProject" 'l' "set lsqproject option to YES"
                 <*> optional (optPath "popList" 'p' "give poplist file to restrict PCA to \
                                                      \populations listed.")

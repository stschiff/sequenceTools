#!/usr/bin/env stack
-- stack --resolver lts-6.4 --install-ghc runghc --package turtle --package foldl --package text

{-# LANGUAGE OverloadedStrings #-}

import Control.Foldl (list)
import Control.Monad (forM_)
import Data.List (sortOn, groupBy)
import Data.Maybe (fromJust)
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Prelude hiding (FilePath)
import Turtle

data Options = Options FilePath FilePath FilePath Int

main = do
    opts <- options "prepare Admixture data for DataGraph" optParser
    printData opts
  where
    optParser = Options <$> optPath "admixtureFile" 'a' "Input Admixture file (*.Q)"
                        <*> optPath "indFile" 'i' "Input *.ind file. This should be formatted in \
                                    \EigenStrat, which includes the individual name in the \
                                     \first, and the population name in the third column."
                        <*> optPath "popGroupFile" 'p' "A file with two columns. The first is the \
                              \population, the second population group, e.g. a continental group. \
                              \Note that the order of populations in this file also determines \
                              \the order in the output."
                        <*> optInt "blankLines" 'b' "Number of blank lines between populations"
                           
printData (Options admixtureF popF popGroupF blankLines) = do
    popGroupDat <- readPopGroupDat popGroupF
    admixtureDat <- fold (readAdmixtureDat popGroupDat admixtureF popF) list
    let (_, _, _, vals) = head admixtureDat
        k = length vals
        sortIndices = [(p, i) | ((p, _), i) <- zip popGroupDat [0..]]
        sortedDat = sortOn (\(_, p, _, _) -> fromJust $ lookup p sortIndices) admixtureDat
        legendedDat = putLegend sortedDat
    echo . T.intercalate "\t" $ ["Sample", "Pop", "PopGroup", "Label"] ++
                                [format ("Q"%d) i | i <- [1..k]]
    forM_ legendedDat $ \group -> do
        forM_ group $ \(sample, pop, popGroup, legend, vals) -> do
            echo . T.intercalate "\t" $ [sample, pop, popGroup, legend] ++ map (format g) vals
        replicateM_ blankLines (echo "")

readPopGroupDat :: FilePath -> IO [(Text, Text)]
readPopGroupDat popGroupF = do
    l <- fold (input popGroupF) list
    return [(p, pG) | [p, pG] <- map (cut (some space) . T.strip) l]

readAdmixtureDat :: [(Text, Text)] -> FilePath -> FilePath -> Shell (Text, Text, Text, [Double])
readAdmixtureDat popGroupDat admixtureF popF = do
    (admixtureL, indL) <- paste (input admixtureF) (input popF)
    let vals = map (read . T.unpack) . cut (some space) $ admixtureL
        [sample, _, pop] = cut (some space) . T.strip $ indL
    Just popGroup <- return $ pop `lookup` popGroupDat
    return (sample, pop, popGroup, vals)

putLegend :: [(Text, Text, Text, [Double])] -> [[(Text, Text, Text, Text, [Double])]]
putLegend admixtureDat = do
    group <- groups
    let l = length group
        (_, pop, _, _) = head group
        labels = [if i == l `div` 2 then pop else "" | i <- [0..(l - 1)]]
    return [(s, p, pg, l, v) | ((s, p, pg, v), l) <- zip group labels]
  where
    groups = groupBy (\(_, pop1, _, _) (_, pop2, _, _) -> pop1 == pop2) admixtureDat
        
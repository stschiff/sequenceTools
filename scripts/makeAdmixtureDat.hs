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

data Options = Options FilePath FilePath FilePath Int Bool

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
                        <*> switch "clusterSORT" 'c' "Sort by cluster, ignore order given in \
                                                      \popGroups file"
                           
printData (Options admixtureF popF popGroupF blankLines clusterSort) = do
    popGroupDat <- readPopGroupDat popGroupF
    admixtureDat <- fold (readAdmixtureDat popGroupDat admixtureF popF) list
    let (_, _, _, vals) = head admixtureDat
        k = length vals
        sortedDat = if clusterSort
                    then sortByCluster admixtureDat
                    else sortByPopGroupFile admixtureDat popGroupDat
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

sortByPopGroupFile :: [(Text, Text, Text, [Double])] -> [(Text, Text)]
                   -> [(Text, Text, Text, [Double])]
sortByPopGroupFile admixtureDat popGroupDat = 
    sortOn (\(_, p, _, _) -> fromJust $ lookup p sortIndices) admixtureDat
  where
    sortIndices = [(p, i) | ((p, _), i) <- zip popGroupDat [0..]]

sortByCluster :: [(Text, Text, Text, [Double])] -> [(Text, Text, Text, [Double])]
sortByCluster admixtureDat =
    let groups = groupBy (\(_, p1, _, _) (_, p2, _, _) -> p1 == p2) . sortOn (\(_, p, _, _) -> p) $ 
                 admixtureDat
        groupClusterWeights = [getColumnAverage . map (\(_, _, _, vals) -> vals) $ g | g <- groups]
        internallySortedGroups = zipWith sortInternally groupClusterWeights groups
    in  concat $ sortExternally groupClusterWeights internallySortedGroups
  where
    sortInternally weights group = sortOn (\(_, _, _, vals) -> negate (vals !! maxCluster)) group
      where
        maxCluster = last . map fst . sortOn snd . zip [0..] $ weights
    sortExternally weightMatrix groups = map snd . sortOn (maxIndex . fst) . sortOn (negate . maximum . fst) . zip weightMatrix $ groups
      where
        maxIndex = last . map fst . sortOn snd . zip [0..]
    
    
getColumnAverage :: [[Double]] -> [Double]
getColumnAverage mat = [(sum . map (!!i) $ mat) / fromIntegral n | i <- [0 .. (k - 1)]]
  where
    k = length . head $ mat
    n = length mat

putLegend :: [(Text, Text, Text, [Double])] -> [[(Text, Text, Text, Text, [Double])]]
putLegend admixtureDat = do
    group <- groups
    let l = length group
        (_, pop, _, _) = head group
        labels = [if i == l `div` 2 then pop else "" | i <- [0..(l - 1)]]
    return [(s, p, pg, l, v) | ((s, p, pg, v), l) <- zip group labels]
  where
    groups = groupBy (\(_, pop1, _, _) (_, pop2, _, _) -> pop1 == pop2) admixtureDat
        
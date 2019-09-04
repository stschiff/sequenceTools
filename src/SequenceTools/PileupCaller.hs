module SequenceTools.PileupCaller (callToDosage, Call(..), callGenotypeFromPileup,
    callMajorityAllele, findMajorityAlleles, callRandomAllele,
    callRandomDiploid, dosageToEigenstratGeno, freqSumToEigenstrat, CallingMode(..),
    TransitionsMode(..), filterTransitions) where

import SequenceFormats.Utils (Chrom(..))
import SequenceFormats.FreqSum (FreqSumEntry(..))
import SequenceFormats.Eigenstrat (EigenstratSnpEntry(..), GenoEntry(..), GenoLine)
import SequenceTools.Utils (sampleWithoutReplacement)

import qualified Data.ByteString.Char8 as B
import Data.List (sortOn, group, sort)
import Data.Vector (fromList)
import Pipes (Pipe, cat)
import qualified Pipes.Prelude as P

-- |A datatype to represent a single genotype call
data Call = HaploidCall Char | DiploidCall Char Char | MissingCall deriving (Show, Eq)

-- |A datatype to specify the calling mode
data CallingMode = MajorityCalling Bool | RandomCalling | RandomDiploidCalling

data TransitionsMode = TransitionsMissing | SkipTransitions | AllSites

-- |a function to turn a call into the dosage of non-reference alleles
callToDosage :: Char -> Char -> Call -> Maybe Int
callToDosage refA altA call = case call of
    HaploidCall a | a == refA -> Just 0
                  | a == altA -> Just 1
                  | otherwise -> Nothing
    DiploidCall a1 a2 | (a1, a2) == (refA, refA) -> Just 0
                      | (a1, a2) == (refA, altA) -> Just 1
                      | (a1, a2) == (altA, refA) -> Just 1
                      | (a1, a2) == (altA, altA) -> Just 2
                      | otherwise                -> Nothing
    MissingCall -> Nothing

-- |Make a call from alleles 
callGenotypeFromPileup :: CallingMode -> Int -> String -> IO Call
callGenotypeFromPileup mode minDepth alleles =
    if length alleles < minDepth then return MissingCall else
        case mode of
            MajorityCalling withDownsampling -> callMajorityAllele withDownsampling minDepth alleles
            RandomCalling -> callRandomAllele alleles
            RandomDiploidCalling -> callRandomDiploid alleles

-- |Sample the majority allele, or one of the majority alleles
callMajorityAllele :: Bool -> Int-> String -> IO Call
callMajorityAllele withDownsampling minDepth alleles = do
    maybeAlleles <- if withDownsampling
        then sampleWithoutReplacement alleles minDepth
        else return (Just alleles)
    case maybeAlleles of
        Nothing -> return MissingCall
        Just alleles' -> do
            a <- case findMajorityAlleles alleles' of
                [] -> error "should not happen"
                [a'] -> return a'
                listA -> do
                    r <- sampleWithoutReplacement listA 1
                    case r of
                        Just [r'] -> return r'
                        _ -> error "should not happen"
            return $ HaploidCall a

-- |Find the majority allele(s)
findMajorityAlleles :: String -> [Char]
findMajorityAlleles alleles =
    let groupedAlleles = sortOn fst [(length g, head g) | g <- group . sort $ alleles]
        majorityCount = fst . last $ groupedAlleles
    in  [a | (n, a) <- groupedAlleles, n == majorityCount]

-- |call a random allele
callRandomAllele :: String -> IO Call
callRandomAllele alleles = do
    res <- sampleWithoutReplacement alleles 1
    case res of
        Nothing -> return MissingCall
        Just [a] -> return $ HaploidCall a
        _ -> error "should not happen"

-- |call two random alleles
callRandomDiploid :: String -> IO Call
callRandomDiploid alleles = do
    res <- sampleWithoutReplacement alleles 2
    case res of
        Nothing -> return MissingCall
        Just [a1, a2] -> return $ DiploidCall a1 a2
        _ -> error "should not happen"

-- |convert a freqSum entry to an eigenstrat SNP entry
freqSumToEigenstrat :: Bool -> FreqSumEntry -> (EigenstratSnpEntry, GenoLine)
freqSumToEigenstrat diploidizeCall (FreqSumEntry chrom@(Chrom c) pos maybeSnpId ref alt calls) =
    let snpId_ = case maybeSnpId of 
            Just id_ -> B.pack id_
            Nothing -> B.pack $ c <> "_" <> show pos
        snpEntry = EigenstratSnpEntry chrom pos 0.0 snpId_ ref alt
        geno = fromList . map (dosageToEigenstratGeno diploidizeCall) $ calls
    in  (snpEntry, geno)

-- |convert a Dosage to an eigenstrat-encoded genotype
dosageToEigenstratGeno :: Bool -> Maybe Int -> GenoEntry
dosageToEigenstratGeno diploidizeCall c =
    if diploidizeCall then
        case c of
            Just 0 -> HomRef
            Just 1 -> HomAlt
            Nothing -> Missing
            _ -> error "illegal call for pseudo-haploid Calling method"
    else
        case c of
            Just 0 -> HomRef
            Just 1 -> Het
            Just 2 -> HomAlt
            Nothing -> Missing
            _ -> error ("unknown genotype " ++ show c)

filterTransitions :: (Monad m) => TransitionsMode -> Pipe FreqSumEntry FreqSumEntry m ()
filterTransitions transversionsMode =
    case transversionsMode of
        SkipTransitions ->
            P.filter (\(FreqSumEntry _ _ _ ref alt _) -> isTransversion ref alt)
        TransitionsMissing ->
            P.map (\(FreqSumEntry chrom pos id_ ref alt calls) ->
                let calls' = if isTransversion ref alt then calls else [Nothing | _ <- calls]
                in  FreqSumEntry chrom pos id_ ref alt calls')
        AllSites -> cat
  where
    isTransversion ref alt = not $ isTransition ref alt
    isTransition ref alt =
        ((ref == 'A') && (alt == 'G')) ||
        ((ref == 'G') && (alt == 'A')) ||
        ((ref == 'C') && (alt == 'T')) ||
        ((ref == 'T') && (alt == 'C'))

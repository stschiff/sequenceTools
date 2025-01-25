{-# LANGUAGE OverloadedStrings #-}
module SequenceTools.PileupCallerSpec (spec) where

import           SequenceTools.PileupCaller (Call (..), CallingMode (..),
                                             TransitionsMode (..),
                                             callGenotypeFromPileup,
                                             callMajorityAllele,
                                             callRandomAllele,
                                             callRandomDiploid, callToDosage,
                                             cleanSSdamageAllSamples,
                                             filterTransitions,
                                             findMajorityAlleles)
import           SequenceTools.Utils        (dosageToEigenstratGeno,
                                             freqSumToEigenstrat)

import           Control.Monad              (forM_, replicateM)
import qualified Data.ByteString.Char8      as B
import           Data.List                  (sort)
import           Data.Vector                (fromList)
import           Pipes                      (each, (>->))
import qualified Pipes.Prelude              as P
import           SequenceFormats.Eigenstrat (EigenstratSnpEntry (..),
                                             GenoEntry (..))
import           SequenceFormats.FreqSum    (FreqSumEntry (..))
import           SequenceFormats.Pileup     (Strand (..))
import           SequenceFormats.Utils      (Chrom (..))
import           Test.Hspec

spec :: Spec
spec = do
    testCallToDosage
    testCallRandomDiploid
    testCallRandomAllele
    testCallMajorityAllele
    testFindMajorityAlleles
    testCallGenotypeFromPileup
    testDosageToEigenstratGeno
    testFreqSumToEigenstrat
    testFilterTransitions
    testCleanSSdamageAllSamples

testCallToDosage :: Spec
testCallToDosage = describe "callToDosage" $ do
    it "should return Nothing for missing call" $
        callToDosage 'A' 'C' MissingCall `shouldBe` Nothing
    it "should return Nothing for haploid non-congruent call" $
        callToDosage 'A' 'C' (HaploidCall 'G') `shouldBe` Nothing
    it "should return 0 for haploid ref call" $
        callToDosage 'A' 'C' (HaploidCall 'A') `shouldBe` Just (0, 1)
    it "should return 1 for haploid alt call" $
        callToDosage 'A' 'C' (HaploidCall 'C') `shouldBe` Just (1, 1)
    it "should return Nothing for diploid non-congruent call" $
        callToDosage 'A' 'C' (DiploidCall 'A' 'G') `shouldBe` Nothing
    it "should return 0 for diploid hom-ref call" $
        callToDosage 'A' 'C' (DiploidCall 'A' 'A') `shouldBe` Just (0, 2)
    it "should return 1 for diploid het call" $ do
        callToDosage 'A' 'C' (DiploidCall 'A' 'C') `shouldBe` Just (1, 2)
        callToDosage 'A' 'C' (DiploidCall 'C' 'A') `shouldBe` Just (1, 2)
    it "should return 2 for diploid hom-alt call" $
        callToDosage 'A' 'C' (DiploidCall 'C' 'C') `shouldBe` Just (2, 2)


testCallGenotypeFromPileup :: Spec
testCallGenotypeFromPileup = describe "callGenotypeFromPileup" $ do
    it "should return missing if pileup below minDepth" $
        callGenotypeFromPileup RandomCalling 3 "A" `shouldReturn` MissingCall
    it "should not return missing if pileup above minDepth" $
        callGenotypeFromPileup RandomCalling 3 "AACCC" `shouldNotReturn` MissingCall


testCallMajorityAllele :: Spec
testCallMajorityAllele = describe "callMajorityAllele" $ do
    it "should call A from AAA" $
        callMajorityAllele False 1 "AAA" `shouldReturn` HaploidCall 'A'
    it "should call A from AAAAA with ds" $
        callMajorityAllele True 3 "AAAAA" `shouldReturn` HaploidCall 'A'
    it "should call Missing from AA with ds 3" $
        callMajorityAllele True 3 "AA" `shouldReturn` MissingCall
    it "should call A from AAC" $
        callMajorityAllele False 1 "AAC" `shouldReturn` HaploidCall 'A'
    it "should call C from ACC" $
        callMajorityAllele False 1 "ACC" `shouldReturn` HaploidCall 'C'
    it "should call 50/50 from AACC" $ do
        r <- replicateM 1000 (callMajorityAllele False 1 "AACC")
        let c = [rr | rr <- r, rr == HaploidCall 'A']
        length c `shouldSatisfy` (\c' -> c' >= 418 && c' <= 582) --p < 1e-7

testFindMajorityAlleles :: Spec
testFindMajorityAlleles = describe "findMajorityAllele" $ do
    it "should return A for AAC" $
        findMajorityAlleles "AAC" `shouldBe` "A"
    it "should return C for ACC" $
        findMajorityAlleles "ACC" `shouldBe` "C"
    it "should return AC for AACC" $
        findMajorityAlleles "AACC" `shouldBe` "AC"

testCallRandomAllele :: Spec
testCallRandomAllele = describe "callRandomAllele" $ do
    it "should return A for AAA" $
        callRandomAllele "AAA" `shouldReturn` HaploidCall 'A'
    it "should return C for C" $
        callRandomAllele "C" `shouldReturn` HaploidCall 'C'
    it "should return A,C or G for ACG roughly with 30% each" $ do
        r <- replicateM 1000 (callRandomAllele "ACG")
        forM_ ['A', 'C', 'G'] $ \nuc -> do
            let n = length . filter (==HaploidCall nuc) $ r
            n `shouldSatisfy` (\nn -> nn >= 257 && nn <= 412) --p < 1e-7

testCallRandomDiploid :: Spec
testCallRandomDiploid = describe "callRandomDiploid" $ do
    it "should return Missing for A" $
        callRandomDiploid "A" `shouldReturn` MissingCall
    it "should return AC for AC" $ do
        DiploidCall a1 a2 <- callRandomDiploid "AC"
        sort [a1, a2] `shouldBe` "AC"
    it "should return 50% het for AACC" $ do
        r <- replicateM 1000 (callRandomDiploid "AACC")
        let n = length ['A' | DiploidCall a1 a2 <- r, a1 /= a2]
        n `shouldSatisfy` (\nn -> nn >= 588 && nn < 743) --p < 1e-7

testDosageToEigenstratGeno :: Spec
testDosageToEigenstratGeno = describe "dosageToEigenstratGeno" $ do
    it "should give Hom-Ref for 0 pseudo-haploid" $
        dosageToEigenstratGeno (Just (0, 1)) `shouldBe` HomRef
    it "should give Hom-Alt for 1 pseudo-haploid" $
        dosageToEigenstratGeno (Just (1, 1)) `shouldBe` HomAlt
    it "should give Missing for Nothing pseudo-haploid" $
        dosageToEigenstratGeno Nothing `shouldBe` Missing
    it "should give Hom-Ref for 0 diploid" $
        dosageToEigenstratGeno (Just (0, 2)) `shouldBe` HomRef
    it "should give Het for 1 diploid" $
        dosageToEigenstratGeno (Just (1, 2)) `shouldBe` Het
    it "should give Hom-Alt for 2 diploid" $
        dosageToEigenstratGeno (Just (2, 2)) `shouldBe` HomAlt
    it "should give Missing for Nothing diploid" $
        dosageToEigenstratGeno Nothing `shouldBe` Missing

testFreqSumToEigenstrat :: Spec
testFreqSumToEigenstrat = describe "freqSumtoEigenstrat" $ do
    let fs = FreqSumEntry (Chrom "1") 1000 Nothing Nothing 'A' 'C' [Just (0, 1), Just (1, 1), Just (1, 1), Nothing, Just (0, 1)]
    let es = EigenstratSnpEntry (Chrom "1") 1000 0.0 (B.pack "1_1000") 'A' 'C'
        genoLine = fromList [HomRef, HomAlt, HomAlt, Missing, HomRef]
    it "should convert a freqSum example correctly to eigenstrat" $
        freqSumToEigenstrat fs `shouldBe` (es, genoLine)
    it "should convert a freqSum example with rsId correctly to eigenstrat" $
        freqSumToEigenstrat (fs {fsSnpId = Just "rs123"}) `shouldBe` (es {snpId = "rs123"}, genoLine)


mockFreqSumData :: [FreqSumEntry]
mockFreqSumData = [
    FreqSumEntry (Chrom "1") 1000 (Just "rs1") Nothing 'A' 'C' [Just (1, 2), Just (2, 2), Nothing, Just (0, 2), Just (0, 2)],
    FreqSumEntry (Chrom "1") 2000 (Just "rs2") Nothing 'C' 'T' [Just (1, 2), Just (2, 2), Nothing, Just (0, 2), Just (0, 2)],
    FreqSumEntry (Chrom "1") 3000 (Just "rs3") Nothing 'A' 'G' [Just (1, 2), Just (2, 2), Nothing, Just (0, 2), Just (0, 2)],
    FreqSumEntry (Chrom "2") 1000 (Just "rs4") Nothing 'A' 'G' [Just (1, 2), Just (2, 2), Nothing, Just (0, 2), Just (0, 2)],
    FreqSumEntry (Chrom "2") 2000 (Just "rs5") Nothing 'T' 'A' [Just (1, 2), Just (2, 2), Nothing, Just (0, 2), Just (0, 2)],
    FreqSumEntry (Chrom "2") 3000 (Just "rs6") Nothing 'T' 'C' [Just (1, 2), Just (2, 2), Nothing, Just (0, 2), Just (0, 2)]]

testFilterTransitions :: Spec
testFilterTransitions = describe "filterTransitions" $ do
    it "should remove transitions with SkipTransitions" $ do
        let r = P.toList $ each mockFreqSumData >-> filterTransitions SkipTransitions
        r `shouldBe` [
            FreqSumEntry (Chrom "1") 1000 (Just "rs1") Nothing 'A' 'C'
                [Just (1, 2), Just (2, 2), Nothing, Just (0, 2), Just (0, 2)],
            FreqSumEntry (Chrom "2") 2000 (Just "rs5") Nothing 'T' 'A'
                [Just (1, 2), Just (2, 2), Nothing, Just (0, 2), Just (0, 2)]]
    it "should mark transitions as missing with TransitionsMissing" $ do
        let r = P.toList $ each mockFreqSumData >-> filterTransitions TransitionsMissing
        r `shouldBe` [
            FreqSumEntry (Chrom "1") 1000 (Just "rs1") Nothing 'A' 'C'
                [Just (1, 2), Just (2, 2), Nothing, Just (0, 2), Just (0, 2)],
            FreqSumEntry (Chrom "1") 2000 (Just "rs2") Nothing 'C' 'T' [Nothing, Nothing, Nothing, Nothing, Nothing],
            FreqSumEntry (Chrom "1") 3000 (Just "rs3") Nothing 'A' 'G' [Nothing, Nothing, Nothing, Nothing, Nothing],
            FreqSumEntry (Chrom "2") 1000 (Just "rs4") Nothing 'A' 'G' [Nothing, Nothing, Nothing, Nothing, Nothing],
            FreqSumEntry (Chrom "2") 2000 (Just "rs5") Nothing 'T' 'A'
                [Just (1, 2), Just (2, 2), Nothing, Just (0, 2), Just (0, 2)],
            FreqSumEntry (Chrom "2") 3000 (Just "rs6") Nothing 'T' 'C' [Nothing, Nothing, Nothing, Nothing, Nothing]]
    it "should output all sites with AllSites" $ do
        let r = P.toList $ each mockFreqSumData >-> filterTransitions AllSites
        r `shouldBe` mockFreqSumData
    it "should output all sites with SingleStrandMode" $ do
        let r = P.toList $ each mockFreqSumData >-> filterTransitions SingleStrandMode
        r `shouldBe` mockFreqSumData

testCleanSSdamageAllSamples :: Spec
testCleanSSdamageAllSamples = describe "cleanSSdamageAllSamples" $ do
    let bases = ["AACATG", "AACATT", "AACTTG"]
        strands = [[f, r, r, f, r, r], [r, f, r, f, f, r], [f, f, r, f, f, r]]
    it "should not remove anything if not C/T or G/A SNP" $
        cleanSSdamageAllSamples 'C' 'A' bases strands `shouldBe` bases
    it "should remove forward reads from C/T SNPs" $ do
        cleanSSdamageAllSamples 'C' 'T' bases strands `shouldBe` ["ACTG", "ACT", "CG"]
        cleanSSdamageAllSamples 'T' 'C' bases strands `shouldBe` ["ACTG", "ACT", "CG"]
    it "should remove reverse reads from G/A SNPs" $ do
        cleanSSdamageAllSamples 'A' 'G' bases strands `shouldBe` ["AA", "AAT", "AATT"]
        cleanSSdamageAllSamples 'G' 'A' bases strands `shouldBe` ["AA", "AAT", "AATT"]
  where
    f = ForwardStrand
    r = ReverseStrand

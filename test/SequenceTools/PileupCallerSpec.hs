module SequenceTools.PileupCallerSpec (spec) where

import SequenceFormats.FreqSum (FreqSumEntry(..))
import SequenceFormats.Eigenstrat (GenoEntry(..), EigenstratSnpEntry(..))
import SequenceFormats.Utils (Chrom(..))
import SequenceTools.PileupCaller (callToDosage, Call(..), callGenotypeFromPileup,
    callMajorityAllele, findMajorityAlleles, callRandomAllele,
    callRandomDiploid, callToEigenstratGeno, freqSumToEigenstrat, CallingMode(..))

import Control.Monad (replicateM, forM_)
import qualified Data.ByteString.Char8 as B
import Data.List (sort)
import Data.Vector (fromList)
import Test.Hspec

spec :: Spec
spec = do
    testCallToDosage
    testCallRandomDiploid
    testCallRandomAllele
    testCallMajorityAllele
    testFindMajorityAlleles
    testCallGenotypeFromPileup
    testCallToEigenstratGeno
    testFreqSumToEigenstrat

testCallToDosage :: Spec
testCallToDosage = describe "callToDosage" $ do
    it "should return Nothing for missing call" $
        callToDosage 'A' 'C' MissingCall `shouldBe` Nothing
    it "should return Nothing for haploid non-congruent call" $
        callToDosage 'A' 'C' (HaploidCall 'G') `shouldBe` Nothing
    it "should return 0 for haploid ref call" $
        callToDosage 'A' 'C' (HaploidCall 'A') `shouldBe` Just 0
    it "should return 1 for haploid alt call" $
        callToDosage 'A' 'C' (HaploidCall 'C') `shouldBe` Just 1        
    it "should return Nothing for diploid non-congruent call" $
        callToDosage 'A' 'C' (DiploidCall 'A' 'G') `shouldBe` Nothing
    it "should return 0 for diploid hom-ref call" $
        callToDosage 'A' 'C' (DiploidCall 'A' 'A') `shouldBe` Just 0
    it "should return 1 for diploid het call" $ do
        callToDosage 'A' 'C' (DiploidCall 'A' 'C') `shouldBe` Just 1
        callToDosage 'A' 'C' (DiploidCall 'C' 'A') `shouldBe` Just 1
    it "should return 2 for diploid hom-alt call" $
        callToDosage 'A' 'C' (DiploidCall 'C' 'C') `shouldBe` Just 2


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

testCallToEigenstratGeno :: Spec
testCallToEigenstratGeno = describe "callToEigenstratGeno" $ do
    it "should give Hom-Ref for 0 pseudo-haploid" $
        callToEigenstratGeno True (Just 0) `shouldBe` HomRef
    it "should give Hom-Alt for 1 pseudo-haploid" $
        callToEigenstratGeno True (Just 1) `shouldBe` HomAlt
    it "should give Missing for Nothing pseudo-haploid" $
        callToEigenstratGeno True Nothing `shouldBe` Missing
    it "should give Hom-Ref for 0 diploid" $
        callToEigenstratGeno False (Just 0) `shouldBe` HomRef
    it "should give Het for 1 diploid" $
        callToEigenstratGeno False (Just 1) `shouldBe` Het
    it "should give Hom-Alt for 2 diploid" $
        callToEigenstratGeno False (Just 2) `shouldBe` HomAlt
    it "should give Missing for Nothing diploid" $
        callToEigenstratGeno False Nothing `shouldBe` Missing

testFreqSumToEigenstrat :: Spec
testFreqSumToEigenstrat = describe "freqSumtoEigenstrat" $ do
    let fs = FreqSumEntry (Chrom "1") 1000 'A' 'C' [Just 0, Just 1, Just 1, Nothing, Just 0]
    let es = EigenstratSnpEntry (Chrom "1") 1000 0.0 (B.pack "1_1000") 'A' 'C'
        genoLine = fromList [HomRef, HomAlt, HomAlt, Missing, HomRef]
    it "should convert a freqSum example correctly to eigenstrat" $
        freqSumToEigenstrat True fs `shouldBe` (es, genoLine)
    


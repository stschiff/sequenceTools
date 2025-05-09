module SequenceTools.UtilsSpec (spec) where

import           SequenceTools.Utils (sampleWithoutReplacement)

import           Control.Monad       (replicateM_)
import           Data.List           (nub, sort, union)
import           Test.Hspec

spec :: Spec
spec = testSampleWithoutReplacement

testSampleWithoutReplacement :: Spec
testSampleWithoutReplacement = describe "sampleWithoutReplacement" $ do
    it "should return Nothing if sample 1 from empty list" $
        sampleWithoutReplacement ([] :: [Char]) 1 `shouldReturn` Nothing
    it "should return one item if sample 1 from 1" $
        sampleWithoutReplacement ['A'] 1 `shouldReturn` Just ['A']
    it "should return an empty list of 0 sampled from 0" $
        sampleWithoutReplacement ([] :: [Char]) 0 `shouldReturn` Just []
    it "should return an empty list if 0 sampled from 1" $
        sampleWithoutReplacement ['A'] 0 `shouldReturn` Just []
    it "should return Nothing if 2 sampled from 1" $
        sampleWithoutReplacement ['A'] 2 `shouldReturn` Nothing
    it "should return 2 items if 2 sampled from 2" $ do
        r <- sampleWithoutReplacement ['A', 'C'] 2
        fmap sort r `shouldBe` Just (sort ['A', 'C'])
    it "should return a non-duplicate subset of ABCDEFGHIJ if 4 are sampled" $ do
        replicateM_ 10 $ do
            Just r <- sampleWithoutReplacement "ABCDEFGHIJ" 4
            length (nub r) `shouldBe` 4
            length (union r "ABCDEFGHIJ") `shouldBe` 10


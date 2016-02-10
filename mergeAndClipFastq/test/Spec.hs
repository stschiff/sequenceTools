{-# LANGUAGE OverloadedStrings #-}

import Lib (searchReadEndDistance, mergeReads, revComp, qualToChar, charToQual)

import qualified Data.ByteString.Char8 as B
import Test.Tasty (defaultMain, TestTree, testGroup)
import Test.Tasty.HUnit (Assertion, assertEqual, testCase)

main :: IO ()
main = defaultMain tests

tests :: TestTree
tests = testGroup "Tests" [ testCase "testing revComp" assertCorrectRevComp,
                            testCase "testing searchEndDistance" assertCorrectReadEndDistance,
                            testCase "testing mergeReads" assertCorrectMergeReads]

assertCorrectRevComp :: Assertion
assertCorrectRevComp = do
    assertEqual "revComp" "AAGGT" (revComp "ACCTT")
    assertEqual "revComp" "AANGT" (revComp "ACNTT")

assertCorrectReadEndDistance :: Assertion
assertCorrectReadEndDistance = do
    assertEqual "searchReadEndDist" (Just 6) (searchReadEndDistance "ACCTGCCTGC" "GCTGGTAATC" 0 4)
    assertEqual "searchReadEndDist" (Just 8) (searchReadEndDistance "ACCTGCCTGC" "AGCCAGGTCC" 0 4)
    assertEqual "searchReadEndDist" (Just 10) (searchReadEndDistance "ACCTGCCTGC" "GCAGGCTGGT" 0 8)
    assertEqual "searchReadEndDist" (Just 13) (searchReadEndDistance "ACCTGCCTGC" "TTTGCAGGCT" 0 4)
    assertEqual "searchReadEndDist" Nothing (searchReadEndDistance "ACCTGCCTGC" "TTGGTTGGCC" 0 4)

assertCorrectMergeReads :: Assertion
assertCorrectMergeReads = do
    let seq1 = "ACCTGCCTGC"
        seq2 = "GCTGGTAATC"
        qual1 = replicate 10 40
        qual2 = replicate 10 30
        (mergedSeq, mergedQual) = mergeReads seq1 (qualToText qual1) seq2 (qualToText qual2) 6
    assertEqual "mergedSeq" "ACCTGC" mergedSeq
    assertEqual "mergedQual" [40, 40, 30, 40, 40, 40] (textToQual mergedQual)
  where
    qualToText = B.pack . map qualToChar
    textToQual = map charToQual . B.unpack
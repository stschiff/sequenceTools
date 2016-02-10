module Lib
    ( searchReadEndDistance
    , mergeReads
    ) where

import qualified Data.Text as T

searchReadEndDistance :: T.Text -> T.Text -> Double -> Int -> Maybe Int
searchReadEndDistance seq1 seq2 mismatchRate minOverlapSize = undefined

mergeReads :: T.Text -> T.Text -> T.Text -> T.Text -> Int -> (T.Text, T.Text)
mergeReads seq1 qual1 seq2 qual2 dist = undefined

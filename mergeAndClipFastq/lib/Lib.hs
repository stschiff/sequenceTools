{-# LANGUAGE OverloadedStrings #-}

module Lib
    ( searchReadEndDistance
    , mergeReads
    , revComp
    , qualToChar
    , charToQual
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Char (ord, chr)

searchReadEndDistance :: B.ByteString -> B.ByteString -> Double -> Int -> Maybe Int
searchReadEndDistance seq1 seq2 mismatchRate minOverlapSize = undefined

mergeReads :: B.ByteString -> B.ByteString -> B.ByteString -> B.ByteString -> Int -> (B.ByteString, B.ByteString)
mergeReads seq1 qual1 seq2 qual2 dist = undefined

revComp :: B.ByteString -> B.ByteString
revComp = undefined

qualToChar :: Int -> Char
qualToChar = chr . (+33)

charToQual :: Char -> Int
charToQual = (-) 33 . ord

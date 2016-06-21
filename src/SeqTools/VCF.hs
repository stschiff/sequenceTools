module SeqTools.VCF (VCFheader(..),
                     VCFentry(..),
                     readVCF,
                     vcfHeaderParser,
                     vcfEntryParser,
                     getGenotypes,
                     getDosages,
                     isBiallelicSnp) where

{-# LANGUAGE OverloadedStrings #-}

import Data.Text (Text, count)

data VCFheader = VCFheader {
    vcfHeaderComments :: [Text],
    vcfSampleNames :: [Text]
}

data VCFentry = VCFentry {
    vcfChrom :: Text,
    vcfPos :: Int,
    vcfId :: Text,
    vcfRef :: Text,
    vcfAlt :: [Text],
    vcfQual :: Double,
    vcfFilter :: Text,
    vcfInfo :: [(Text, Text)],
    vcfFormatString :: [Text],
    vcfGenotypeInfo :: [[Text]]
}

readVCF :: (MonadIO m) => m (VCFheader, Producer VCFentry m ())
readVCF = do
    (res, rest) <- runStateT (parse vcfHeaderParser) PT.stdin
    header <- case res of
        Nothing -> liftIO . throwIO $ AssertionFailed "vcf header not readible. VCF file empty?"
        Just (Left e_) -> do
            Right (chunk, _) <- next rest
            let msg = show e_ ++ T.unpack chunk
            liftIO . throwIO $ AssertionFailed ("VCF header parsing error: " ++ msg)
        Just (Right h) -> return h
    return (header, parsed vcfEntryParser rest >>= liftErrors)

liftErrors :: (MonadIO m) => Either (ParsingError, Producer T.Text m r) () -> Producer a m ()
liftErrors res = case res of
    Left (e_, prod_) -> do
        Right (chunk, _) <- lift $ next prod_
        let msg = show e_ ++ T.unpack chunk
        lift . liftIO . throwIO $ AssertionFailed msg
    Right () -> return ()

vcfHeaderParser :: A.Parser VCFheader
vcfHeaderParser = VCFheader <$> A.many' doubleCommentLine <*> singleCommentLine
  where
    doubleCommentLine = do
        c1 <- A.string "##"
        s_ <- A.takeTill A.isEndOfLine <* A.endOfLine
        return $ T.append c1 s_
    singleCommentLine = do
        void $ A.char '#'
        s_ <- A.takeTill A.isEndOfLine <* A.endOfLine
        let fields = T.splitOn "\t" s_
        return . drop 9 $ fields

vcfEntryParser :: A.Parser VCFentry
vcfEntryParser = VCFentry <$> word <* A.space <*> A.decimal <* A.space <*> word <* A.space <*> 
                              alternativeAlleles <* A.space <*> A.double <* A.space <*> word <*
                              A.space <*> infoFields <* A.space <*> formatStrings <* A.space <*>
                              genotypeInfos <* A.endOfLine
  where
    word = A.takeTill isSpace
    alternativeAlleles = word `A.sepBy1` A.char ','
    infoFields = infoField `A.sepBy1` A.char ';'
    infoField = (,) <$> word <* A.char '=' <*> word
    formatStrings = word `A.sepBy` A.char ':'
    genotypeInfos = genotype `A.sepBy1` A.space
    genotype = word `A.sepBy1` ':'

isBiallelicSnp :: VCFentry -> Bool
isBiallelicSnp vcfEntry = validRef && validAlt
  where
    validRef = (vcfRef vcfEntry `elem` "ACTG")
    validAlt = length (vcfAlt vcfEntry) == 1 && head (vcfAlt vcfEntry) `elem` "ACTG"

getGenotypes :: VCFentry -> Either String [Text]
getGenotypes vcfEntry = do
    gtIndex <- fst <$> tryHead "GT format field not found" . filter ((=="GT") . snd) . zip [0..] . 
                       vcfFormatString $ vcfEntry
    return $ map (!!gtIndex) (vcfGenotypeInfo vcfEntry)

getDosages :: VCFentry -> Either String [Int]
getDosages vcfEntry = do
    genotypes <- getGenotypes vcfEntry
    return $ map (count "1") genotypes


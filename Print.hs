-- Print module.
-- By G.W. Schwartz

-- Collects all functions pertaining to the printing of data structures to
-- a string for saving to a file.

module Print where

-- Built in
import Data.List
import Data.Char
import qualified Data.Map as M
import Data.Ord

-- | Local
import Types
import FastaDiversity
import Diversity
import CompareDiversityMutationCount

-- Return the results of the mutation or stable counts as a string
printMutStabCounts :: Bool -> MutationMap AminoAcid -> String
printMutStabCounts mutBool mutationMap = header ++ body
  where
    header           = "position,count,count_weight,hydrophobicity\n"
    body             = unlines                          .
                       map mapLine                      .
                       M.toAscList                      $
                       mutationMap
    mapLine (x, xs) = show x                                            ++
                      ","                                               ++
                      (show . length . mutList $ xs)                    ++
                      ","                                               ++
                      (show . length . filterMutStab (\_ -> True) $ xs) ++
                      ","                                               ++
                      (classifyPosition . nub . getImportantAA . mutList $ xs)
    mutList          = map toUpper . mutUniquer mutBool
    mutUniquer True  = map snd . filterMutStab isMutation
    mutUniquer False = map fst . filterMutStab (not . isMutation)

-- Return the results of the different types of mutations or stable  as a
-- string. Can basically just add "nub" with "snd" to mapLine. Awwwwwww yeah.
-- The power of functional programming. :P
printMutStabTypeCounts :: Bool -> Double -> MutationMap AminoAcid -> String
printMutStabTypeCounts mutBool order mutationMap = header ++ body
  where
    header           = "position,count,count_weight,hydrophobicity\n"
    body             = unlines                          .
                       map mapLine                      .
                       M.toAscList                      $
                       mutationMap
    mapLine (x, xs) = show x                                            ++
                      ","                                               ++
                      (show . diversity order . mutList $ xs)           ++
                      ","                                               ++
                      (show . length . filterMutStab (\_ -> True) $ xs) ++
                      ","                                               ++
                      (classifyPosition . nub . getImportantAA . mutList $ xs)
    mutList          = map toUpper . mutUniquer mutBool
    mutUniquer True  = map snd . filterMutStab isMutation
    mutUniquer False = map fst . filterMutStab (not . isMutation)

printMutStabAAUse :: Bool -> MutationMap AminoAcid -> String
printMutStabAAUse mutBool mutationMap = header ++ body
  where
    header           = "position,aa_use,count,count_weight\n"
    body             = concat              .
                       map mapPositionLine .
                       M.toAscList         $
                       mutationMap
    mapPositionLine (x, xs) = unlines                 .
                              map (mapAALine x count) .
                              groupLengthSort         .
                              getImportantAA          $
                              mutList
      where
        count      = length mutList
        mutList    = map toUpper . mutUniquer mutBool $ xs
    mapAALine x weight xs = show x       ++
                            ","          ++
                            showAAUse xs ++
                            ","          ++
                            show weight
    showAAUse []           = "0, 0"
    showAAUse total@(x:xs) = x : ("," ++ (show . length $ total))
    mutUniquer True        = map snd . filterMutStab isMutation
    mutUniquer False       = map fst . filterMutStab (not . isMutation)

-- Return the results of the sample rarefaction percents as a string
printRarefaction :: MutationMap AminoAcid -> String
printRarefaction mutationMap = header ++ body
  where
    header           = "position,percent_above\n"
    body             = unlines                          .
                       map mapLine                      .
                       M.toAscList                      $
                       mutationMap
    mapLine (x, xs) = show x ++
                      ","    ++
                      (show . percent $ xs)
    percent         = rarefactionViable . rarefactionCurve . mutList
    mutList         = filterMutStab taut
    taut _          = True

-- Return the ChangedAAMap as a string for saving to a file, either with
-- positions or diversities specified by the input
printChangedAAMap :: DivPos -> ChangedAAMap -> String
printChangedAAMap divPos changedAAMap = header divPos ++ body
  where
    header Diversity      = "diversity,germline,clone,germline_codon,\
                            \clone_codon,germline_before,germline_after,\
                            \clone_before,clone_after,mutations,\
                            \mutation_position,size,type\n"
    header Position       = "position,germline,clone,germline_codon,\
                            \clone_codon,germline_before,germline_after,\
                            \clone_before,clone_after,mutations,\
                            \mutation_position,size,type\n"
    body                  = intercalate "\n"  .
                            map mapDiv        .
                            M.toAscList       .
                            M.map countGroups $
                            changedAAMap
    mapDiv (div, xs)      = intercalate "\n" . map (\ys -> mapLine div ys) $ xs
    mapLine div xs        = show div ++ "," ++ lineFormat xs
    lineFormat all@(x:xs) = intercalate "," [ germlineAA x : []
                                            , cloneAA x : []
                                            , germlineCodon x
                                            , cloneCodon x
                                            , germlineBefore x
                                            , germlineAfter x
                                            , cloneBefore x
                                            , cloneAfter x
                                            , (show . numMutations $ x)
                                            , mutPositions x
                                            , (show . length $ all)
                                            , mutType x ]
    countGroups           = groupBy groupFormat .
                            sortBy (comparing sortFormatCodon)
    groupFormat x y       = sortFormatCodon x == sortFormatCodon y
    mutType x             = if germlineCodon x == cloneCodon x
                                then
                                    if germlineAA x == cloneAA x
                                        then "silent"
                                        else "replacement"
                                else "na"

-- Return the SepAAMap as a string for saving to a file, either with
-- positions or diversities specified by the input
printSepAAMap :: DivPos -> SepAAMap -> String
printSepAAMap divPos sepAAMap = header divPos ++ body
  where
    header Diversity      = "diversity,seq_aa,seq_codon,\
                            \seq_before,seq_after,\
                            \size\n"
    header Position       = "position,seq_aa,seq_codon,\
                            \seq_before,seq_after,\
                            \size\n"
    body                  = intercalate "\n"
                          . map mapDiv
                          . M.toAscList
                          . M.map countGroups
                          $ sepAAMap
    mapDiv (div, xs)      = intercalate "\n" . map (\ys -> mapLine div ys) $ xs
    mapLine div xs        = show div ++ "," ++ lineFormat xs
    lineFormat all@(x:xs) = intercalate "," [ seqAA x : []
                                            , seqCodon x
                                            , seqBefore x
                                            , seqAfter x
                                            , (show . length $ all) ]
    countGroups           = groupBy groupFormat
                          . sortBy (comparing sortFormatSepCodon)
    groupFormat x y       = sortFormatSepCodon x == sortFormatSepCodon y

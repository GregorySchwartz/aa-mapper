-- CompareDiversityMutationCount module.
-- By G.W. Schwartz
--
-- Called by Main to pull together the functions pertaining to Main's
-- primary goal.

module CompareDiversityMutationCount where

-- Built in
import Data.List
import Data.Char
import qualified Data.Map as M
import Data.Ord
import Control.Applicative
import System.IO
import System.Environment

-- Cabal
import qualified Data.List.Split as Split

-- Local
import Diversity
import Translation

-- Algebraic
data GeneticUnit = AminoAcid | Codon

-- Basic
type ID                = Int
type Clone a           = Sequence a
type Germline a        = Sequence a
type Position          = Int
type Diversity         = Int
type Size              = Int

-- Advanced
type Mutation a    = (a, a)
type CloneMap a    = M.Map (ID, Germline a) [Clone a]
type MutationMap a = M.Map Position [Mutation a]
type CloneMutMap a = M.Map (ID, Germline a) (MutationMap a)
type DiversityMap  = M.Map Position Diversity
type ChangedAAMap  = M.Map Diversity [[(AminoAcid, AminoAcid, Size)]]


-- Takes a DW2 fasta file string and returns a CloneMap in order to
-- generate the basic building block for the mutation counting.
-- Note: Several repeating germlines, so they need a unique identifier (an
-- integer in this case).
generateCloneMap :: String -> CloneMap AminoAcid
generateCloneMap = M.fromList . getSequences
  where
    getSequences                   = map filterAssocList . assocList
    filterAssocList ((x1, x2), xs) = ((x1, last . filterHeaders $ x2),
                                      filterHeaders xs)
    filterHeaders                  = filter (\x -> head x /= '>')
    assocList                      = map assocMap . germlineSplit
    assocMap (x, y)                = ((x, germline y), clones y)
    germlineSplit                  = zip [0..]      .
                                     filter (/= "") .
                                     Split.splitOn ">>"
    germline                       = take 2 . lines
    clones                         = drop 2 . lines

-- Takes a DW2 fasta file string and returns a CloneMap in order to
-- generate the basic building block for the mutation counting.
-- Note: Several repeating germlines, so they need a unique identifier (an
-- integer in this case). This version is for Codons.
generateCodonCloneMap :: String -> CloneMap Codon
generateCodonCloneMap = M.fromList . map getCodonSplit . getSequences
  where
    getCodonSplit ((x1, x2), xs)   = ((x1, fullCodon . Split.chunksOf 3 $ x2),
                                      map (fullCodon . Split.chunksOf 3) xs)
    fullCodon                      = filter ((==3) . length)
    getSequences                   = map filterAssocList . assocList
    filterAssocList ((x1, x2), xs) = ((x1, last . filterHeaders $ x2),
                                      filterHeaders xs)
    filterHeaders                  = filter (\x -> head x /= '>')
    assocList                      = map assocMap . germlineSplit
    assocMap (x, y)                = ((x, germline y), clones y)
    germlineSplit                  = zip [0..]      .
                                     filter (/= "") .
                                     Split.splitOn ">>"
    germline                       = take 2 . lines
    clones                         = drop 2 . lines

-- Checks if a pair is actually a mutation
isMutation :: Mutation AminoAcid -> Bool
isMutation (x, y)
    | x == y    = False
    | otherwise = True

isCodonMutation :: Mutation Codon -> Bool
isCodonMutation (x, y)
    | codon2aa x == codon2aa y = False
    | otherwise                = True

-- Checks if a pair is actually a mutation. Can only be done with codons.
isSilentMutation :: Mutation Codon -> Bool
isSilentMutation (x, y)
    | x /= y && codon2aa x == codon2aa y = True
    | otherwise = False

-- Return the mutations if they exist between the germline and a clone
countMutations :: Germline a ->
                  Clone a    ->
                  [(Position, Mutation a)]
countMutations germline clone = mutation germline clone
  where
    mutation x y         = zip [1..] . zip x $ y

-- Takes a list of mutations and returns the mutations (or stable) that are
-- actual mutations (the mutation is not just the same character, which is
-- stable) with no gaps.
filterMutStab :: (Mutation AminoAcid -> Bool) ->
                 [Mutation AminoAcid]         ->
                 [Mutation AminoAcid]
filterMutStab isWhat = filter filterRules
  where
    filterRules x    = isWhat x            &&
                       not (inTuple '-' x) &&
                       not (inTuple '.' x) &&
                       not (inTuple '~' x)
    inTuple c (x, y) = if x == c || y == c then True else False
filterCodonMutStab :: (Mutation Codon -> Bool) ->
                      [Mutation Codon]         ->
                      [Mutation Codon]
filterCodonMutStab isWhat = filter filterRules
  where
    filterRules x    = isWhat x            &&
                       not (inTuple '-' x) &&
                       not (inTuple '.' x) &&
                       not (inTuple '~' x) &&
                       not (inTuple 'N' x)
    inTuple c (x, y)     = if c `elem` x || c `elem` y then True else False

-- Filters the cloneMap to remove clones with 30 mutations or greater
filterCloneMap :: CloneMap AminoAcid -> CloneMap AminoAcid
filterCloneMap = M.mapWithKey filterMutated
  where
    filterMutated k xs  = filter (not . isHighlyMutated (snd k)) xs
    isHighlyMutated k x = 30 <= (length . realMutations k $ x)
    realMutations k x   = filterMutStab isMutation .
                          map snd                      .
                          countMutations k             $
                          x

-- Filters the cloneMap to remove clones with 30 mutations or greater, this
-- is the nucleotide version. There is too much adhoc going on here,
-- I don't like the state of this file.
filterCodonCloneMap :: CloneMap Codon -> CloneMap Codon
filterCodonCloneMap = M.mapWithKey filterMutated
  where
    filterMutated k xs  = filter (not . isHighlyMutated (snd k)) xs
    isHighlyMutated k x = 30 <= (length . realMutations k $ x)
    realMutations k x   = filterCodonMutStab isCodonMutation .
                          map snd                      .
                          countMutations k             $
                          x

-- Join together mutation lists
joinMutations :: (Eq a) => [[(Position, Mutation a)]] -> MutationMap a
joinMutations = M.map nub . groupedMutations
  where
    groupedMutations = M.fromListWith (++) . map (\(x, y) -> (x, [y])) . concat

-- Generate a CloneMutMap which will then be printed to save files.
generateCloneMutMap :: (Eq a) => CloneMap a -> CloneMutMap a
generateCloneMutMap = M.mapWithKey gatherMutations
  where
    gatherMutations k xs = joinMutations . map (countMutations (snd k)) $ xs

-- Generate a ChangedAAMap which contains all of the aminoacids a certain
-- amino acid at a certain diversity goes to.
generateChangedAAMap :: [Int]             ->
                        DiversityMap      ->
                        MutationMap Codon ->
                        ChangedAAMap
generateChangedAAMap viablePos germDivMap = aaDivMap              .
                                            M.filter (not . null) .
                                            diversityMap          .
                                            filterNonviablePos    .
                                            realMutMap
  where
    aaDivMap                = M.map (\xs -> group . sort . map numMut $ xs)
    realMutMap              = M.map (filterCodonMutStab isCodonMutation)
    numMut (x, y)           = (codon2aa x, codon2aa y, hamming x y)
    filterNonviablePos      = M.filterWithKey (\k _ -> elem k viablePos)
    diversityMap            = M.fromListWith (++) .
                              map keysToDiversity .
                              M.toAscList
    keysToDiversity (x, xs) = (getDiversity x, xs)
    getDiversity x          = extractMaybe . M.lookup x $ germDivMap
    extractMaybe (Just x)   = x
    extractMaybe Nothing    = error "Clone position not in germline positions"

-- Returns a list of stuff ordered by how many of a type of stuff there
-- are, like [1,1,2,2,3,3,3,3,4,4,4,4,4] -> [[4,4,4,4,4], [3,3,3,3], [2,2],
-- [1,1]]
groupLengthSort :: (Ord a) => [a] -> [[a]]
groupLengthSort = reverse . sortBy (comparing length) . group . sort

-- Classify a list of amino acids as a certain kind of hydrophobicity
classifyAA :: Char -> String
classifyAA aa
    | aa `elem` "IVLFCMW" = "Hydrophobic"
    | aa `elem` "AGTSYPH" = "Neutral"
    | aa `elem` "NDQEKR"  = "Hydrophilic"
    | aa == '*'           = "Stop"
    | otherwise           = error ("Amino acid not found: " ++ [aa])

-- CLassify a position based on the numbers of types of hydrophobicity
-- I realize that the indeterminate should be grouped under otherwise, but
-- I just want to make obvious the rules here
classifyPosition :: [Char] -> String
classifyPosition aaList
    | dominate sortedList      = head . head $ sortedList
    | ignore sortedList        = head . head $ sortedList
    | weakFirst sortedList     = "Weak " ++ (head . head $ sortedList)
    | weakSecond sortedList    = "Weak " ++ (head . last $ sortedList)
    | indeterminate sortedList = "Indeterminate"
    | otherwise                = "Indeterminate"
  where
    dominate xs      = length xs == 1
    ignore xs        = length xs == 2 &&
                       length (head xs) >= 3 &&
                       length (last xs) == 1
    weakFirst xs     = length xs == 2 &&
                       length (last xs) > 1 &&
                       (head . last $ xs) == "Neutral"
    weakSecond xs    = length xs == 2 &&
                       length (last xs) > 1 &&
                       (head . head $ xs) == "Neutral"
    indeterminate xs = length xs == 3
    sortedList       = groupLengthSort    .
                       filter (/= "Stop") .
                       map classifyAA $ aaList

-- Return the Diversity Order 1 most important amino acids
getImportantAA :: (Ord a) => Sequence a -> Sequence a
getImportantAA mutList = concat . takeWhile biggerLength $ sortedMutList
  where
    biggerLength x = length x >= weight
    weight         = if shannonDiv == 0
                    then length . head $ sortedMutList
                    else length $ sortedMutList !! (shannonDiv - 1)
    shannonDiv     = round . diversity 1 $ mutList
    sortedMutList  = groupLengthSort mutList

-- Return the results of the mutation or stable counts as a string
printMutStabCounts :: Bool -> MutationMap AminoAcid -> String
printMutStabCounts mutBool mutationMap = header ++ body
  where
    header           = "position,count,count_weight,hydrophobicity\n"
    body             = unlines                          .
                       map mapLine                      .
                       M.toAscList                      $
                       mutationMap
    mapLine (x, xs) = show x                                    ++
                      ","                                       ++
                      (show . length . mutList $ xs)            ++
                      ","                                       ++
                      (show . length . mutList $ xs)            ++
                      ","                                       ++
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
    mapLine (x, xs) = show x                                  ++
                      ","                                     ++
                      (show . diversity order . mutList $ xs) ++
                      ","                                     ++
                      (show . length . mutList $ xs)          ++
                      ","                                     ++
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

printChangedAAMap :: ChangedAAMap -> String
printChangedAAMap changedAAMap = header ++ body
  where
    header                = "diversity,germline,clone,mutations,size\n"
    body                  = intercalate "\n" .
                            map mapDiv       .
                            M.toAscList      $
                            changedAAMap
    mapDiv (div, xs)      = intercalate "\n" . map (\ys -> mapLine div ys) $ xs
    mapLine div xs        = show div ++ "," ++ lineFormat xs
    lineFormat all@(x:xs) = [fst' x]          ++
                            ","               ++
                            [snd' x]          ++
                            ","               ++
                            (show . thd' $ x) ++
                            ","               ++
                            (show . length $ all)
    fst' (x, _, _) = x
    snd' (_, x, _) = x
    thd' (_, _, x) = x

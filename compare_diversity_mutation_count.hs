-- Compare Diversity with Mutation Counts
-- By G.W. Schwartz

-- Takes a DW2 fasta file and generates the mutation counts of the clones
-- from each listed germline in the file per position in a dataframe type
-- format. For use with comparing the counts with the diversities generated
-- already.

import Data.List
import Data.Char
import qualified Data.Map as M
import Data.Ord
import Control.Applicative
import System.IO
import System.Environment

import qualified Data.List.Split as Split

-- Basic
type ID                = Int
type Sequence          = String
type Clone             = Sequence
type Germline          = Sequence
type Position          = Int

-- Advanced
type Mutation          = (Char, Char)
type CloneMap    = M.Map (ID, Germline) [Clone]
type MutationMap = M.Map Position [Mutation]
type CloneMutMap = M.Map (ID, Germline) MutationMap

-- Returns the diversity of a list of things
diversity :: (Ord b) => Double -> [b] -> Double
diversity order sample
    | length sample == 0 = 0
    | order == 1         = exp . h $ speciesList
    | otherwise          = (sum . map ((** order) . p_i) $ speciesList) ** pow
  where
    pow          = 1 / (1 - order)
    h            = negate . sum . map (\x -> (p_i x) * (log (p_i x)))
    p_i x        = ((fromIntegral . length $ x) :: Double) /
                   ((fromIntegral . length $ sample) :: Double)
    speciesList  = group . sort $ sample

-- Calculates the binary coefficient
choose :: (Integral a) => a -> a -> a
choose n 0 = 1
choose 0 k = 0
choose n k = choose (n - 1) (k - 1) * n `div` k

-- Returns the rarefaction curve for each position in a list
rarefactionCurve :: (Eq a, Ord a) => [a] -> [Double]
rarefactionCurve xs = map rarefact [1..n_total]
  where
    rarefact n
        | n == 0       = 0
        | n == 1       = 1
        | n == n_total = k
        | otherwise    = k - ((1 / (fromIntegral (choose n_total n))) * inner n)
    inner n = fromIntegral                              .
              sum                                       .
              map (\g -> choose (n_total - length g) n) $
              grouped
    n_total = length xs 
    k       = genericLength grouped
    grouped = group . sort $ xs

-- Calculates the percent of the curve that is above 95% of height of the curve
rarefactionViable :: [Double] -> Double
rarefactionViable xs = (genericLength valid / genericLength xs) * 100
  where
    valid = dropWhile (< (0.95 * last xs)) xs

-- Takes a DW2 fasta file string and returns a CloneMap in order to
-- generate the basic building block for the mutation counting.
-- Note: Several repeating germlines, so they need a unique identifier (an
-- integer in this case).
generateCloneMap :: String -> CloneMap
generateCloneMap = M.fromList . getSequences
  where
    getSequences            = map filterAssocList . assocList
    filterAssocList ((x1, x2), xs) = ((x1, last . filterHeaders $ x2),
                                      filterHeaders xs)
    filterHeaders           = filter (\x -> head x /= '>')
    assocList               = map assocMap . germlineSplit
    assocMap (x, y)         = ((x, germline y), clones y)
    germlineSplit           = zip [0..] . filter (/= "") . Split.splitOn ">>"
    germline                = take 2 . lines
    clones                  = drop 2 . lines

-- Checks if a pair is actually a mutation
isMutation :: Mutation -> Bool
isMutation (x, y)
    | x == y    = False
    | otherwise = True

-- Return the mutations if they exist between the germline and a clone
countMutations :: Germline ->
                  Clone    ->
                  [(Position, Mutation)]
countMutations germline clone = mutation germline clone
  where
    mutation x y         = zip [1..] . zip x $ y

-- Takes a list of mutations and returns the mutations (or stable) that are
-- actual mutations (the mutation is not just the same character, which is
-- stable) with no gaps.
filterMutStab :: (Mutation -> Bool) ->
                 [Mutation]         ->
                 [Mutation]
filterMutStab isWhat = filter filterRules
  where
    filterRules x    = isWhat x            &&
                       not (inTuple '-' x) &&
                       not (inTuple '.' x) &&
                       not (inTuple '~' x)
    inTuple c (x, y) = if x == c || y == c then True else False

-- Filters the cloneMap to remove clones with 30 mutations or greater
filterCloneMap :: CloneMap -> CloneMap
filterCloneMap = M.mapWithKey filterMutated
  where
    filterMutated k xs  = filter (not . isHighlyMutated (snd k)) xs
    isHighlyMutated k x = 30 <= (length . realMutations k $ x)
    realMutations k x   = filterMutStab isMutation .
                          map snd                  .
                          countMutations k         $
                          x

-- Join together mutation lists
joinMutations :: [[(Position, Mutation)]] -> MutationMap
joinMutations = M.map nub . groupedMutations
  where
    groupedMutations = M.fromListWith (++) . map (\(x, y) -> (x, [y])) . concat

-- Generate a CloneMutMap which will then be printed to save files
generateCloneMutMap :: CloneMap -> CloneMutMap
generateCloneMutMap = M.mapWithKey gatherMutations
  where
    gatherMutations k xs = joinMutations . map (countMutations (snd k)) $ xs

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
getImportantAA :: [Char] -> [Char]
getImportantAA mutList = concat . takeWhile biggerLength $ sortedMutList
  where
    biggerLength x = length x >= weight
    weight         = if shannonDiv == 0
                    then length . head $ sortedMutList
                    else length $ sortedMutList !! (shannonDiv - 1)
    shannonDiv     = round . diversity 1 $ mutList
    sortedMutList  = groupLengthSort mutList


-- Return the results of the mutation or stable counts as a string
printMutStabCounts :: Bool -> MutationMap -> String
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
printMutStabTypeCounts :: Bool -> Double -> MutationMap -> String
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

printMutStabAAUse :: Bool -> MutationMap -> String
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
printRarefaction :: MutationMap -> String
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

main = do
    (orderContents:
     source:
     saveMutCounts:
     saveStabCounts:
     saveMutDiversityCounts:
     saveStabDiversityCounts:
     saveMutAAUse:
     saveStabAAUse:
     saveRarefaction:[]) <- getArgs

    contents <- readFile source
    let order = (read orderContents) :: Double

    let unfilteredCloneMap  = generateCloneMap contents
    let cloneMap            = filterCloneMap unfilteredCloneMap
    let cloneMutMap         = generateCloneMutMap cloneMap
    let combinedCloneMutMap = M.unionsWith (++) .
                              map snd           .
                              M.toAscList       $
                              cloneMutMap

    writeFile saveMutCounts $ printMutStabCounts True combinedCloneMutMap
    writeFile saveStabCounts $ printMutStabCounts False combinedCloneMutMap
    writeFile saveMutDiversityCounts $ printMutStabTypeCounts True order combinedCloneMutMap
    writeFile saveStabDiversityCounts $ printMutStabTypeCounts False order combinedCloneMutMap
    writeFile saveMutAAUse $ printMutStabAAUse True combinedCloneMutMap
    writeFile saveStabAAUse $ printMutStabAAUse False combinedCloneMutMap
    writeFile saveRarefaction $ printRarefaction combinedCloneMutMap

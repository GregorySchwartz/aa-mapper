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
import Types
import Diversity
import Translation
import FastaDiversity

-- Same as takeWhile, but also cuts off at a certain length and ignores
-- unsatisfied predicates
takeWhileN :: (a -> Bool) -> Int -> [a] -> [a]
takeWhileN _ _ [] = []
takeWhileN p n (x:xs)
    | p x && n > 0 = x : takeWhileN p (n - 1) xs
    | not (p x) && n > 0 = takeWhileN p n xs
    | otherwise     = []

-- Takes a fasta file string and removes newlines in the sequences to make
-- this compatible with the fasta parser. The lineCompress function should
-- get rid of any extra newlines messing with the code.
joinSeq :: String -> String
joinSeq = lineCompress
        . tail
        . concat
        . map newEntry
        . filter (/= "")
        . Split.splitOn ">>"
  where
    newEntry x             = if elem '>' x then cloneEntry x else germEntry x
    germEntry x            = newGerm x
    cloneEntry x           = newGerm (germline x)
                          ++ concat (map newClone . filter (/= "") . clone $ x)
    newGerm x
        | seq x /= ""      = "\n>>" ++ (header x) ++ "\n" ++ (seq x)
        | otherwise        = ""
    newClone x
        | seq x /= ""      = "\n>" ++ (header x) ++ "\n" ++ (seq x)
        | otherwise        = ""
    germline               = head . Split.splitOn ">"
    clone                  = tail . Split.splitOn ">"
    header                 = head . lines
    seq                    = concat . tail . lines
    lineCompress []        = []
    lineCompress ('\n':xs) = '\n' : (lineCompress $ dropWhile (== '\n') xs)
    lineCompress (x:xs)    = x : (lineCompress xs)

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
    germlineSplit                  = zip [0..]
                                   . filter (\x -> elem '>' x)  -- Only clones
                                   . filter (/= "")
                                   . Split.splitOn ">>"
    germline                       = take 2 . lines
    clones                         = drop 2 . lines

-- | Takes in a list of Codons and finds the Septamers according to the
-- position of the codons. For use with finding mutability for the amino
-- acid maps.
septize :: [Codon] -> [Septamer]
septize xs
    | length xs < 3 = []
    | otherwise     = ( Septamer { before = ""
                                 , after  = takeWhileN (/= '-') 2
                                          . concat
                                          . drop 1
                                          $ xs
                                 , codon  = head xs }
                      : septizeR [] xs )
  where
    septizeR ls (x:y:z:xs)
        | length y == 3 = Septamer { before = reverse
                                            . takeWhileN (/= '-') 2
                                            . concat
                                            $ (x:ls)
                                   , after  = takeWhileN (/= '-') 2
                                            . concat
                                            $ z
                                            : xs
                                   , codon  = y }
                        : septizeR (x:ls) (y:z:xs)
        | otherwise     = septizeR (x:ls) (y:z:xs)
    septizeR ls (x:y:[])
        | length y == 3 = Septamer { before = reverse
                                            . takeWhileN (/= '-') 2
                                            . concat
                                            $ x:ls
                                   , after  = ""
                                   , codon  = y }
                        : []
        | otherwise     = []

-- Takes a DW2 fasta file string and returns a CloneMap in order to
-- generate the basic building block for the mutation counting.
-- Note: Several repeating germlines, so they need a unique identifier (an
-- integer in this case). This version is for Codons.
generateCodonCloneMap :: String -> CloneMap Septamer
generateCodonCloneMap = M.fromList . map getCodonSplit . getSequences
  where
    getCodonSplit ((x1, x2), xs)   = ((x1, septize . Split.chunksOf 3 $ x2),
                                      map (septize . Split.chunksOf 3) xs)
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

isCodonMutation :: Mutation Septamer -> Bool
isCodonMutation (x, y)
    | (codon2aa . codon $ x) == (codon2aa . codon $ y) = False
    | otherwise                                        = True

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

-- Takes a list of mutations and returns the mutations (or stable) that are
-- actual mutations (the mutation is not just the same character, which is
-- stable) with no gaps, Codon version.
filterCodonMutStab :: (Mutation Septamer -> Bool) ->
                      [Mutation Septamer]         ->
                      [Mutation Septamer]
filterCodonMutStab isWhat = filter filterRules
  where
    filterRules x    = isWhat x            &&
                       not (inTuple '-' x) &&
                       not (inTuple '.' x) &&
                       not (inTuple '~' x) &&
                       not (inTuple 'N' x)
    inTuple c (x, y)
        | c `elem` codon x || c `elem` codon y = True
        | otherwise                            = False

-- Filters the cloneMap to remove clones with 30 mutations or greater
filterCloneMap :: CloneMap AminoAcid -> CloneMap AminoAcid
filterCloneMap = M.mapWithKey filterMutated
  where
    filterMutated k xs  = filter (not . isHighlyMutated (snd k)) xs
    isHighlyMutated k x = ((genericLength x / 3) :: Double)
                       <= ((genericLength . realMutations k $ x) :: Double)
    realMutations k x   = filterMutStab isMutation
                        . map snd
                        . countMutations k
                        $ x

-- Filters the cloneMap to remove clones with 30 mutations or greater, this
-- is the nucleotide version. There is too much adhoc going on here,
-- I don't like the state of this file.
filterCodonCloneMap :: CloneMap Septamer -> CloneMap Septamer
filterCodonCloneMap = M.mapWithKey filterMutated
  where
    filterMutated k xs  = filter (not . isHighlyMutated (snd k)) xs
    isHighlyMutated k x = ((genericLength x :: Double) / 3)
                       <= ((genericLength . realMutations k $ x) :: Double)
    realMutations k x   = filterCodonMutStab isCodonMutation
                        . map snd
                        . countMutations k
                        $ x

-- Join together mutation lists of unique mutations per clone
joinMutations :: (Eq a) => [[(Position, Mutation a)]] -> MutationMap a
joinMutations = M.map nub . groupedMutations
  where
    groupedMutations = M.fromListWith (++) . map (\(x, y) -> (x, [y])) . concat

-- Generate a CloneMutMap which will then be printed to save files.
generateCloneMutMap :: (Eq a) => CloneMap a -> CloneMutMap a
generateCloneMutMap = M.mapWithKey gatherMutations
  where
    gatherMutations k xs = joinMutations . map (countMutations (snd k)) $ xs

-- Convert a list of fasta sequences to a map of septamers at each
-- position.
generateFastaSepMap :: [FastaSequence] -> SeptamerMap
generateFastaSepMap = M.fromListWith (++) . concatMap (zip [1..] . fastaToSep)
  where
    fastaToSep = map (:[]) . septize . Split.chunksOf 3 . fastaSeq

-- Take out gaps from the MutationMap Codon
filterCodonMutationMap :: MutationMap Septamer -> MutationMap Septamer
filterCodonMutationMap = M.filter (not . null) . M.map (filter removeGaps)
  where
    removeGaps (x, y)
        | (elem '-' . codon $ x) || (elem '-' . codon $ y) = False
        | otherwise                                        = True

-- Take out gaps from a map of septamers at each position
filterFastaSepMap :: SeptamerMap -> SeptamerMap
filterFastaSepMap = M.filter (not . null) . M.map (filter removeGaps)
  where
    removeGaps x
        | (elem '-' . codon $ x) = False
        | otherwise              = True

-- Take out gaps from the MutationMap AminoAcid
filterAminoAcidMutationMap :: MutationMap AminoAcid -> MutationMap AminoAcid
filterAminoAcidMutationMap = M.filter (not . null) . M.map (filter removeGaps)
  where
    removeGaps (x, y)
        | x == '-' || y == '-' = False
        | otherwise            = True

-- Generate a ChangedAAMap which contains all of the amino acids a certain
-- amino acid at a certain diversity goes to. Also supports what a position
-- can go to as well using a flag in the input. Also keeps track of a lack
-- of mutations.
generateChangedAAMap :: DivPos
                     -> [Int]
                     -> ( DiversityMap
                       -> Position
                       -> [Mutation Septamer]
                       -> [Mutation Septamer] )
                     -> DiversityMap
                     -> MutationMap Septamer
                     -> ChangedAAMap
generateChangedAAMap isPos
                     viablePos
                     important
                     germDivMap = aaDivPosMap            .
                                  M.filter (not . null)  .
                                  M.fromListWith (++)    .
                                  makeDiversityMap isPos .
                                  filterNonviablePos     .
                                  realMutMap
  where
    aaDivPosMap                = M.map (\xs -> map numMut $ xs)
    realMutMap                 = M.map (filterCodonMutStab (\_ -> True))
    numMut (x, y)              = ChangedAA { germlineAA      = c2aaSept x
                                           , cloneAA         = c2aaSept y
                                           , germlineCodon   = codon x
                                           , cloneCodon      = codon y
                                           , numMutations    = hamming (codon x)
                                                                       (codon y)
                                           , mutPositions    = findMutPos
                                                               (codon x)
                                                               (codon y)
                                           , sortFormatAA    = formAA x y
                                           , sortFormatCodon = formCodon x y
                                           , germlineBefore  = before x
                                           , cloneBefore     = before y
                                           , germlineAfter   = after x
                                           , cloneAfter      = after y
                                           }
    formAA x y                 = ( c2aaSept x
                                 , c2aaSept y
                                 , hamming (codon x) (codon y)
                                 , before x
                                 , before y
                                 , after x
                                 , after y )
    formCodon x y              = ( codon x
                                 , codon y
                                 , hamming (codon x) (codon y)
                                 , before x
                                 , before y
                                 , after x
                                 , after y )
    filterNonviablePos         = M.filterWithKey (\k _ -> elem k viablePos)
    makeDiversityMap Diversity = map keysToDiversity .  M.toAscList
    makeDiversityMap Position  = map keysToPos . M.toAscList
    keysToDiversity (x, xs)    = (getDiversity x, important germDivMap x xs)
    keysToPos (x, xs)          = (x, important germDivMap x xs)
    getDiversity x             = extractMaybe . M.lookup x $ germDivMap
    findMutPos x y             = concatMap (\(p, _) -> show p)
                               . filter (\(p, (a, b)) -> a /= b)
                               . zip [1..]
                               . zip x
                               $ y
    c2aaSept                   = codon2aa . codon
    extractMaybe (Just x)      = x
    extractMaybe Nothing       = error "Clone position not found"

-- | Generate a ChangedAAMap which contains all of the amino acids a certain
-- amino acid at a certain diversity goes to. Also supports what a position
-- can go to as well using a flag in the input. This version is for the no
-- mutations pathway.
generateSepAAMap :: DivPos
                 -> [Int]
                 -> ( DiversityMap
                   -> Position
                   -> [Septamer]
                   -> [Septamer] )
                 -> DiversityMap
                 -> SeptamerMap
                 -> SepAAMap
generateSepAAMap isPos
                 viablePos
                 important
                 germDivMap = aaDivPosMap
                            . M.filter (not . null)
                            . M.fromListWith (++)
                            . makeDiversityMap isPos
                            . filterNonviablePos
  where
    aaDivPosMap                = M.map (map numSep)
    numSep x                   = SepAA { seqAA     = c2aaSept x
                                       , seqCodon  = codon x
                                       , seqBefore = before x
                                       , seqAfter  = after x
                                       , sortFormatSepAA    = formAA x
                                       , sortFormatSepCodon = formCodon x
                                       }
    formAA x                   = ( c2aaSept x
                                 , before x
                                 , after x
                                 )
    formCodon x                = ( codon x
                                 , before x
                                 , after x
                                 )
    filterNonviablePos         = M.filterWithKey (\k _ -> elem k viablePos)
    makeDiversityMap Diversity = map keysToDiversity .  M.toAscList
    makeDiversityMap Position  = map keysToPos . M.toAscList
    keysToDiversity (x, xs)    = (getDiversity x, important germDivMap x xs)
    keysToPos (x, xs)          = (x, important germDivMap x xs)
    getDiversity x             = extractMaybe . M.lookup x $ germDivMap
    c2aaSept                   = codon2aa . codon
    extractMaybe (Just x)      = x
    extractMaybe Nothing       = error "Clone position not found"

-- Classify a list of amino acids as a certain kind of hydrophobicity
classifyAA :: Char -> String
classifyAA aa
    | aa `elem` "IVLFCMW" = "Hydrophobic"
    | aa `elem` "AGTSYPH" = "Neutral"
    | aa `elem` "NDQEKR"  = "Hydrophilic"
    | aa == '*'           = "Stop"
    | otherwise           = error ("Amino acid not found: " ++ [aa])

-- Classify a position based on the numbers of types of hydrophobicity
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

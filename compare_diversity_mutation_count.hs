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
import Options.Applicative

import qualified Data.List.Split as Split

-- Command line arguments
data Options = Options { inputOrder                :: Double
                       , inputFasta                :: String
                       , inputDiversity            :: String
                       , unitFlag                  :: GeneticUnit
                       , outputMutCounts           :: String
                       , outputStabCounts          :: String
                       , outputMutDiversityCounts  :: String
                       , outputStabDiversityCounts :: String
                       , outputMutAAUse            :: String
                       , outputStabAAUse           :: String
                       , outputRarefaction         :: String
                       , outputChangedAAMap        :: String
                       }

-- Algebraic
data GeneticUnit = AminoAcid | Codon

-- Basic
type AminoAcid         = Char
type Codon             = String
type ID                = Int
type Sequence a        = [a]
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

-- Command line options
options :: Parser Options
options = Options
      <$> option
          ( long "inputOrder"
         <> short 'o'
         <> metavar "ORDER"
         <> value 1
         <> help "The order of true diversity" )
      <*> strOption
          ( long "inputFasta"
         <> short 'i'
         <> metavar "FILE"
         <> value ""
         <> help "The fasta file containing the germlines and clones" )
      <*> strOption
          ( long "inputDiversity"
         <> short 'd'
         <> metavar "FILE"
         <> value ""
         <> help "The csv file containing the diversities at each position\
                 \ (must be generated into a specific format" )
      <*> flag AminoAcid Codon
          ( long "nucleotides"
         <> short 'u'
         <> help "Whether these sequences are of nucleotides (Codon) or\
                 \ amino acids (AminoAcid)" )
      <*> strOption
          ( long "outputMutCounts"
         <> short 'm'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the changed amino acid counts" )
      <*> strOption
          ( long "outputStabCounts"
         <> short 's'
         <> metavar "File"
         <> value ""
         <> help "The output file for the maintained amino acid counts" )
      <*> strOption
          ( long "outputMutDiversityCounts"
         <> short 'M'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the hanged amino acid diversities" )
      <*> strOption
          ( long "outputStabDiversityCounts"
         <> short 'S'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the maintained amino acid diversities")
      <*> strOption
          ( long "outputMutAAUse"
         <> short 'j'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the specific changed amino acids used" )
      <*> strOption
          ( long "outputStabAAUse"
         <> short 'k'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the specific maintained amino acids used" )
      <*> strOption
          ( long "outputRarefaction"
         <> short 'r'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the rarefaction curves" )
      <*> strOption
          ( long "outputChangedAAMap"
         <> short 'c'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the map of changed amino acids" )

-- Converts a codon to an amino acid
-- Remember, if there is an "N" in that DNA sequence, then it is invalid
codon2aa :: Codon -> AminoAcid
codon2aa x
    | codon `elem` ["GCT", "GCC", "GCA", "GCG"]               = 'A'
    | codon `elem` ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"] = 'R'
    | codon `elem` ["AAT", "AAC"]                             = 'N'
    | codon `elem` ["GAT", "GAC"]                             = 'D'
    | codon `elem` ["TGT", "TGC"]                             = 'C'
    | codon `elem` ["CAA", "CAG"]                             = 'Q'
    | codon `elem` ["GAA", "GAG"]                             = 'E'
    | codon `elem` ["GGT", "GGC", "GGA", "GGG"]               = 'G'
    | codon `elem` ["CAT", "CAC"]                             = 'H'
    | codon `elem` ["ATT", "ATC", "ATA"]                      = 'I'
    | codon `elem` ["ATG"]                                    = 'M'
    | codon `elem` ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"] = 'L'
    | codon `elem` ["AAA", "AAG"]                             = 'K'
    | codon `elem` ["TTT", "TTC"]                             = 'F'
    | codon `elem` ["CCT", "CCC", "CCA", "CCG"]               = 'P'
    | codon `elem` ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"] = 'S'
    | codon `elem` ["ACT", "ACC", "ACA", "ACG"]               = 'T'
    | codon `elem` ["TGG"]                                    = 'W'
    | codon `elem` ["TAT", "TAC"]                             = 'Y'
    | codon `elem` ["GTT", "GTC", "GTA", "GTG"]               = 'V'
    | codon `elem` ["TAA", "TGA", "TAG"]                      = '*'
    | codon `elem` ["---", "..."]                             = '-'
    | codon == "~~~"                                          = '~'
    | 'N' `elem` codon                                        = '~'
    | '-' `elem` codon                                        = '~'
    | otherwise                                               = error errorMsg
  where
    codon    = map toUpper x
    errorMsg = "Unidentified codon: " ++ codon

-- Translates a string of nucleotides
translate :: String -> Sequence AminoAcid
translate = map codon2aa . filter ((== 3) . length) . Split.chunksOf 3

-- Takes two strings, returns Hamming distance
hamming :: String -> String -> Int
hamming xs ys = length $ filter not $ zipWith (==) xs ys

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

-- Generate  DiversityMap which contains the germline diversity at each
-- position.
generateDiversityMap :: String -> DiversityMap
generateDiversityMap = M.fromList . map diverseTuple . csvFilter . csvParse
  where
    diverseTuple x = (read (splitComma x !! 2) :: Int
                     , round' (read (splitComma x !! 3) :: Double))
    csvFilter      = filter (\x -> isHeavy x && isOrder x && isWindow x)
    isHeavy        = (==) "IGH" . flip (!!) 0 . splitComma
    isOrder        = (==) "1" . flip (!!) 1 . splitComma
    isWindow       = (==) "1" . flip (!!) 5 . splitComma
    splitComma     = Split.splitOn ","
    csvParse       = drop 1 . lines
    round' x
        | round x == 0 = 1
        | otherwise    = round x

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

geneticUnitBranch :: GeneticUnit -> Options -> IO ()
geneticUnitBranch AminoAcid opts = do
    contents <- readFile . inputFasta $ opts
    let order = inputOrder opts

    let unfilteredCloneMap  = generateCloneMap contents
    let cloneMap            = filterCloneMap unfilteredCloneMap
    let cloneMutMap         = generateCloneMutMap cloneMap
    let combinedCloneMutMap = M.unionsWith (++) .
                              map snd           .
                              M.toAscList       $
                              cloneMutMap

    writeFile (outputMutCounts opts) $ printMutStabCounts True combinedCloneMutMap
    writeFile (outputStabCounts opts) $ printMutStabCounts False combinedCloneMutMap
    writeFile (outputMutDiversityCounts opts) $ printMutStabTypeCounts True order combinedCloneMutMap
    writeFile (outputStabDiversityCounts opts) $ printMutStabTypeCounts False order combinedCloneMutMap
    writeFile (outputMutAAUse opts) $ printMutStabAAUse True combinedCloneMutMap
    writeFile (outputStabAAUse opts) $ printMutStabAAUse False combinedCloneMutMap
    writeFile (outputRarefaction opts) $ printRarefaction combinedCloneMutMap
geneticUnitBranch Codon opts = do
    contents <- readFile . inputFasta $ opts
    diversityContents <- readFile . inputDiversity $ opts

    let viablePos     = [25..30] ++ [35..59] ++ [63..72] ++ [74..106]
    let divMap        = generateDiversityMap diversityContents
    let unfilteredCloneMap  = generateCodonCloneMap contents
    let cloneMap            = filterCodonCloneMap unfilteredCloneMap
    let cloneMutMap         = generateCloneMutMap cloneMap
    let combinedCloneMutMap = M.unionsWith (++) .
                              map snd           .
                              M.toAscList       $
                              cloneMutMap

    let changedAAMap = generateChangedAAMap viablePos divMap combinedCloneMutMap

    writeFile (outputChangedAAMap opts) $ printChangedAAMap changedAAMap

compareDiversityMutationCounts :: Options -> IO ()
compareDiversityMutationCounts opts = do
    geneticUnitBranch (unitFlag opts) opts

main :: IO ()
main = execParser opts >>= compareDiversityMutationCounts
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Return various information about the relationship between\
                 \ the germline and the clones"
     <> header "Germline and Clone Comparison, Gregory W. Schwartz")

-- FastaDiversity module.
-- By G.W. Schwartz
--
-- Collection of functions for the parsing of an IMGT fasta file and the
-- collection of amino acids for the diversity calculations.

module FastaDiversity where

-- Built in
import qualified Data.Map as M
import Data.Ord
import Data.List

-- Cabal
import qualified Data.List.Split as Split

-- Local
import Types
import Diversity
import Translation

fastaParser :: String -> [FastaSequence]
fastaParser = map makeFastaSequence . Split.splitOn ">"
  where
    makeFastaSequence x = FastaSequence { fastaInfo = head . lines $ x
                                        , fastaSeq  = getSeq x
                                        }
    getSeq              = concat . drop 1 . lines

generatePositionMap :: [FastaSequence] -> PositionMap
generatePositionMap = M.fromListWith (++) . posSeqList
  where
    posSeqList    = map toList . concatMap (\x -> zip [1..] . fastaSeq $ x)
    toList (x, y) = (x, [y])

-- Generate  DiversityMap which contains the germline diversity at each
-- position.
generateDiversityMap :: String -> DiversityMap
generateDiversityMap = M.fromList . map diverseTuple . csvFilter . csvParse
  where
    diverseTuple x = ( read (splitComma x !! 2) :: Int
                     , round' (read (splitComma x !! 3) :: Double) )
    csvFilter      = filter (\x -> isHeavy x && isOrder x && isWindow x)
    isHeavy        = (==) "IGH" . flip (!!) 0 . splitComma
    isOrder        = (==) "1" . flip (!!) 1 . splitComma
    isWindow       = (==) "1" . flip (!!) 5 . splitComma
    splitComma     = Split.splitOn ","
    csvParse       = drop 1 . lines
    round' x
        | round x == 0 = 1
        | otherwise    = round x

-- Return a PositionMap where the values are a list of the most important
-- germline amino acids
importantAADiversityMap :: DiversityMap -> PositionMap -> PositionMap
importantAADiversityMap divMap = M.mapWithKey important
  where
    important k a = getGermImportantAA divMap k a

-- Returns a list of stuff ordered by how many of a type of stuff there
-- are, like [1,1,2,2,3,3,3,3,4,4,4,4,4] -> [[4,4,4,4,4], [3,3,3,3], [2,2],
-- [1,1]]
groupLengthSort :: (Ord a) => [a] -> [[a]]
groupLengthSort = reverse . sortBy (comparing length) . group . sort

-- Return the Diversity Order 1 most important GeneticUnits
getImportantAA :: (Ord a) => Sequence a -> Sequence a
getImportantAA mutList = concat . takeWhile biggerLength $ sortedMutList
  where
    biggerLength x = length x >= weight
    weight         = if shannonDiv == 0
                     then length . head $ sortedMutList
                     else length $ sortedMutList !! (shannonDiv - 1)
    shannonDiv     = round . diversity 1 $ mutList
    sortedMutList  = groupLengthSort mutList

-- Return the Diversity Order 1 most important amino acids from the
-- existing germline list only
getGermImportantAA :: (Ord a) => DiversityMap ->
                                 Position     ->
                                 Sequence a   ->
                                 Sequence a
getGermImportantAA germDivMap pos genUnitList = concat                 .
                                                takeWhile biggerLength $
                                                sortedGenUnitList
  where
    biggerLength x        = length x >= weight
    weight                = if shannonDiv == 0
                            then length . head $ sortedGenUnitList
                            else length $ sortedGenUnitList !! (shannonDiv - 1)
    shannonDiv            = getDiversity pos
    getDiversity x        = extractMaybe . M.lookup x $ germDivMap
    sortedGenUnitList     = groupLengthSort genUnitList
    extractMaybe (Just x) = x
    extractMaybe Nothing  = error "Diversity position not in germline positions"

mostImportantCodons :: DiversityMap     ->
                       Position         ->
                       [Mutation Codon] ->
                       [Mutation Codon]
mostImportantCodons germDivMap pos l = filter isImportant l
  where
    isImportant (x, y) = importantGerm x && importantClone y
    importantGerm x    = codon2aa x `elem` importantGermlines
    importantClone y   = codon2aa y `elem` getImportantAA clones
    importantGermlines = getGermImportantAA germDivMap pos $ germlines
    germlines          = map (codon2aa . fst) l
    clones             = map (codon2aa . snd) l

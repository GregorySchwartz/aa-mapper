-- Compare Diversity with Mutation Counts
-- By G.W. Schwartz

-- Takes a DW2 fasta file and generates the mutation counts of the clones
-- from each listed germline in the file per position in a dataframe type
-- format. For use with comparing the counts with the diversities generated
-- already.

-- Built in
import qualified Data.Map as M

-- Cabal
import Options.Applicative

-- Local
import Types
import CompareDiversityMutationCount
import FastaDiversity

-- Command line arguments
data Options = Options { inputOrder                    :: Double
                       , inputFasta                    :: String
                       , inputDiversity                :: String
                       , inputAAMapType                :: DivPos
                       , unitFlag                      :: GeneticUnit
                       , outputMutCounts               :: String
                       , outputStabCounts              :: String
                       , outputMutDiversityCounts      :: String
                       , outputStabDiversityCounts     :: String
                       , outputMutAAUse                :: String
                       , outputStabAAUse               :: String
                       , outputRarefaction             :: String
                       , outputAllChangedAAMap         :: String
                       , outputImportantChangedAAMap   :: String
                       , outputUnimportantChangedAAMap :: String
                       }

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
                 \ (must be generated into a specific format)" )
      <*> option
          ( long "inputAAMapType"
         <> short 't'
         <> metavar "DIVERSITY | POSITION"
         <> value Diversity
         <> help "Whether to split the amino acid map by position\
                 \ or diversity" )
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
          ( long "outputAllChangedAAMap"
         <> short 'c'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the map of all changed amino acids" )
      <*> strOption
          ( long "outputImportantChangedAAMap"
         <> short 'y'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the map of important changed\
                 \ amino acids" )
      <*> strOption
          ( long "outputUnimportantChangedAAMap"
         <> short 'z'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the map of unimportant changed\
                 \ amino acids" )

geneticUnitBranch :: GeneticUnit -> Options -> IO ()
geneticUnitBranch AminoAcid opts = do
    contents <- readFile . inputFasta $ opts
    let contentsFormatted = joinSeq contents
    let order = inputOrder opts

    let unfilteredCloneMap  = generateCloneMap contentsFormatted
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
    let divPos        = inputAAMapType opts
    let divMap        = generateDiversityMap diversityContents
    let unfilteredCloneMap  = generateCodonCloneMap contents
    let cloneMap            = filterCodonCloneMap unfilteredCloneMap
    let cloneMutMap         = generateCloneMutMap cloneMap
    let combinedCloneMutMap = M.unionsWith (++) .
                              map snd           .
                              M.toAscList       $
                              cloneMutMap

    let allImportant d p l   = l
    let important         = mostImportantCodons
    let unimportant d p l = filter (\(x, y) ->
                                    (elem x . map fst . important d p $ l) &&
                                    (notElem y . map snd . important d p $ l)) l
    let changedAAMap = generateChangedAAMap divPos
                                            viablePos
                                            allImportant
                                            divMap
                                            combinedCloneMutMap
    let importantChangedAAMap = generateChangedAAMap divPos
                                                     viablePos
                                                     important
                                                     divMap
                                                     combinedCloneMutMap
    let unimportantChangedAAMap = generateChangedAAMap divPos
                                                       viablePos
                                                       unimportant
                                                       divMap
                                                       combinedCloneMutMap

    writeFile (outputAllChangedAAMap opts) $
              printChangedAAMap divPos changedAAMap
    writeFile (outputImportantChangedAAMap opts) $
              printChangedAAMap divPos importantChangedAAMap
    writeFile (outputUnimportantChangedAAMap opts) $
              printChangedAAMap divPos unimportantChangedAAMap

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
     <> header "Germline and Clone Comparison, Gregory W. Schwartz" )

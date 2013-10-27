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
          ( long "outputAllChangedAAMap"
         <> short 'c'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the map of all changed amino acids" )
      <*> strOption
          ( long "outputImportantChangedAAMap"
         <> short 'c'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the map of important changed\
                 \ amino acids" )
      <*> strOption
          ( long "outputUnimportantChangedAAMap"
         <> short 'c'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the map of unimportant changed\
                 \ amino acids" )

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

    let allOrtant d p l   = l
    let important         = mostImportantCodons
    let unimportant d p l = filter (\(x, y) ->
                                    (notElem x . map fst . important d p $ l) &&
                                    (notElem y . map snd . important d p $ l)) l
    let changedAAMap = generateChangedAAMap
                       viablePos allOrtant divMap combinedCloneMutMap
    let importantChangedAAMap = generateChangedAAMap
                                viablePos important divMap combinedCloneMutMap
    let unimportantChangedAAMap = generateChangedAAMap
                                  viablePos unimportant divMap combinedCloneMutMap

    writeFile (outputAllChangedAAMap opts) $
              printChangedAAMap changedAAMap
    writeFile (outputImportantChangedAAMap opts) $
              printChangedAAMap importantChangedAAMap
    writeFile (outputUnimportantChangedAAMap opts) $
              printChangedAAMap unimportantChangedAAMap

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

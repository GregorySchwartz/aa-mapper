-- aa-mapper
-- By Gregory W. Schwartz

-- Takes a DW2 fasta file and generates the amino acid maps from each
-- listed germline in the file per position in a dataframe type format. For
-- use with comparing the counts with the diversities generated already.

-- Built in
import qualified Data.Map as M

-- Cabal
import Options.Applicative

-- Local
import Types
import FastaDiversity
import CompareDiversityMutationCount
import Print

-- Command line arguments
data Options = Options { inputOrder                    :: Double
                       , inputFasta                    :: String
                       , inputDiversity                :: String
                       , inputAAMapType                :: DivPos
                       , unitFlag                      :: GeneticUnit
                       , noMutations                   :: Bool
                       , outputMutCounts               :: String
                       , outputStabCounts              :: String
                       , outputMutDiversityCounts      :: String
                       , outputStabDiversityCounts     :: String
                       , outputMutAAUse                :: String
                       , outputStabAAUse               :: String
                       , outputRarefaction             :: String
                       , outputAllAAMap                :: String
                       , outputImportantAAMap          :: String
                       , outputUnimportantAAMap        :: String
                       }

-- Command line options
options :: Parser Options
options = Options
      <$> option
          ( long "input-order"
         <> short 'o'
         <> metavar "ORDER"
         <> value 1
         <> help "The order of true diversity" )
      <*> strOption
          ( long "input-fasta"
         <> short 'i'
         <> metavar "FILE"
         <> value ""
         <> help "The fasta file containing the germlines and clones" )
      <*> strOption
          ( long "input-diversity"
         <> short 'd'
         <> metavar "FILE"
         <> value ""
         <> help "The csv file containing the diversities at each position\
                 \ (must be generated into a specific format)" )
      <*> option
          ( long "input-AA-map-type"
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
      <*> switch
          ( long "no-mutations"
         <> short 'n'
         <> help "Whether to look at the codons from a fasta file, not\
                 \ from a germline to a clone sequence of mutations\
                 \ but rather (ideally) from a germline only" )
      <*> strOption
          ( long "output-mut-counts"
         <> short 'm'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the changed amino acid counts" )
      <*> strOption
          ( long "output-stab-counts"
         <> short 's'
         <> metavar "File"
         <> value ""
         <> help "The output file for the maintained amino acid counts" )
      <*> strOption
          ( long "output-mut-diversity-counts"
         <> short 'M'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the hanged amino acid diversities" )
      <*> strOption
          ( long "output-stab-diversity-counts"
         <> short 'S'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the maintained amino acid diversities")
      <*> strOption
          ( long "output-mut-AA-use"
         <> short 'j'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the specific changed amino acids used" )
      <*> strOption
          ( long "output-stab-AA-use"
         <> short 'k'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the specific maintained amino acids used" )
      <*> strOption
          ( long "output-rarefaction"
         <> short 'r'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the rarefaction curves" )
      <*> strOption
          ( long "output-all-changed-AA-map"
         <> short 'c'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the map of all changed amino acids" )
      <*> strOption
          ( long "output-important-changedAA-map"
         <> short 'y'
         <> metavar "FILE"
         <> value ""
         <> help "The output file for the map of important changed\
                 \ amino acids" )
      <*> strOption
          ( long "output-unimportant-changed-AA-map"
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
    let unfilteredCombinedCloneMutMap = M.unionsWith (++)
                                      . map snd
                                      . M.toAscList
                                      $ cloneMutMap
    let combinedCloneMutMap = filterAminoAcidMutationMap
                              unfilteredCombinedCloneMutMap

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

--    let viablePos     = [25..30] ++ [35..59] ++ [63..72] ++ [74..106]
    let viablePos     = [1..30] ++ [35..59] ++ [63..72] ++ [74..106]
    let divPos        = inputAAMapType opts
    let divMap        = generateDiversityMap diversityContents
    let unfilteredCloneMap  = generateCodonCloneMap contents
    let cloneMap            = filterCodonCloneMap unfilteredCloneMap
    let cloneMutMap         = generateCloneMutMap cloneMap
    let unfilteredCombinedCloneMutMap = M.unionsWith (++)
                                      . map snd
                                      . M.toAscList
                                      $ cloneMutMap
    let combinedCloneMutMap = filterCodonMutationMap
                              unfilteredCombinedCloneMutMap

    let allImportant _ _ l  = l
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

    writeFile (outputAllAAMap opts) $
              printChangedAAMap divPos changedAAMap
    writeFile (outputImportantAAMap opts) $
              printChangedAAMap divPos importantChangedAAMap
    writeFile (outputUnimportantAAMap opts) $
              printChangedAAMap divPos unimportantChangedAAMap

aaMapperNoMutations opts = do
    contents <- readFile . inputFasta $ opts
    diversityContents <- readFile . inputDiversity $ opts

--    let viablePos     = [25..30] ++ [35..59] ++ [63..72] ++ [74..106]
    let viablePos = [1..30] ++ [35..59] ++ [63..72] ++ [74..106]
    let divPos                = inputAAMapType opts
    let divMap                = generateDiversityMap diversityContents
    let fastaList             = fastaParser contents
    let unfilteredFastaSepMap = generateFastaSepMap fastaList
    let fastaSepMap           = filterFastaSepMap unfilteredFastaSepMap

    let allImportant _ _ l = l
    let important         = mostImportantCodonsSep
    let unimportant d p l = filter (\x -> notElem x . important d p $ l) l
    let sepAAMap = generateSepAAMap divPos
                                    viablePos
                                    allImportant
                                    divMap
                                    fastaSepMap
    let importantSepAAMap = generateSepAAMap divPos
                                             viablePos
                                             important
                                             divMap
                                             fastaSepMap
    let unimportantSepAAMap = generateSepAAMap divPos
                                               viablePos
                                               unimportant
                                               divMap
                                               fastaSepMap

    writeFile (outputAllAAMap opts) $
              printSepAAMap divPos sepAAMap
    writeFile (outputImportantAAMap opts) $
              printSepAAMap divPos importantSepAAMap
    writeFile (outputUnimportantAAMap opts) $
              printSepAAMap divPos unimportantSepAAMap

compareDiversityMutationCounts :: Options -> IO ()
compareDiversityMutationCounts opts = do
    if (noMutations opts)
        then aaMapperNoMutations opts
        else geneticUnitBranch (unitFlag opts) opts

main :: IO ()
main = execParser opts >>= compareDiversityMutationCounts
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Return various information about the relationship between\
                 \ the germline and the clones, most importantly the amino\
                 \ acid maps (nucleotide sequences only)"
     <> header "aa-mapper, Gregory W. Schwartz" )

-- Types module.
-- By G.W. Schwartz
--
-- Collects all application specific types.

module Types where

import qualified Data.Map as M

-- Algebraic
data GeneticUnit   = AminoAcid | Codon
data DivPos        = Diversity | Position deriving (Eq, Show, Read)
data FastaSequence = FastaSequence { fastaInfo :: String
                                   , fastaSeq  :: String
                                   } deriving (Eq, Ord, Show)
data ChangedAA     = ChangedAA { germlineAA      :: AminoAcid
                               , cloneAA         :: AminoAcid
                               , germlineCodon   :: Codon
                               , cloneCodon      :: Codon
                               , numMutations    :: Int
                               , mutPositions    :: String
                               , sortFormatAA    :: ChangedAASort AminoAcid
                               , sortFormatCodon :: ChangedAASort Codon
                               , germlineBefore  :: String
                               , cloneBefore     :: String
                               , germlineAfter   :: String
                               , cloneAfter      :: String
                               } deriving (Eq, Ord, Show)
data SepAA         = SepAA { seqAA              :: AminoAcid
                           , seqCodon           :: Codon
                           , seqBefore          :: String
                           , seqAfter           :: String
                           , sortFormatSepAA    :: SepAASort AminoAcid
                           , sortFormatSepCodon :: SepAASort Codon
                           } deriving (Eq, Ord, Show)
-- Septamer is for Codon only, in order to find the two nucleotides before
-- and the two after the codon for use with the mutability calculation.
data Septamer      = Septamer { before :: String
                              , after  :: String
                              , codon  :: Codon
                              } deriving (Eq, Ord, Show)

-- Basic
type AminoAcid         = Char
type Codon             = String
type Sequence a        = [a]
type ID                = Int
type Clone a           = Sequence a
type Germline a        = Sequence a
type Position          = Int
type Diversity         = Int
type Size              = Int

-- Advanced
type Mutation a      = (a, a)
type CloneMap a      = M.Map (ID, Germline a) [Clone a]
type SeptamerMap     = M.Map Position [Septamer]
type MutationMap a   = M.Map Position [Mutation a]
type CloneMutMap a   = M.Map (ID, Germline a) (MutationMap a)
type PositionMap     = M.Map Position [AminoAcid]
type DiversityMap    = M.Map Position Diversity
-- ChangedAAMap keys can be positions or diversities
type ChangedAAMap    = M.Map Int [ChangedAA]
type SepAAMap        = M.Map Int [SepAA]
type ChangedAASort a = (a, a, Int, String, String, String, String)
type SepAASort a     = (a, String, String)

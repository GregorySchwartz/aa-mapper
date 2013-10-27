-- Types module.
-- By G.W. Schwartz
--
-- Collects all application specific types.

module Types where

import qualified Data.Map as M

-- Algebraic
data GeneticUnit   = AminoAcid | Codon
data FastaSequence = FastaSequence { fastaInfo :: String
                                   , fastaSeq  :: String
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
type Mutation a    = (a, a)
type CloneMap a    = M.Map (ID, Germline a) [Clone a]
type MutationMap a = M.Map Position [Mutation a]
type CloneMutMap a = M.Map (ID, Germline a) (MutationMap a)
type PositionMap   = M.Map Position [AminoAcid]
type DiversityMap  = M.Map Position Diversity
type ChangedAAMap  = M.Map Diversity [[(AminoAcid, AminoAcid, Size)]] 

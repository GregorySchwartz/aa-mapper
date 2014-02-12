#aa-mapper

**Gregory W. Schwartz**

```
Usage: aa-mapper [-o|--inputOrder ORDER] [-i|--inputFasta FILE] [-d|--inputDiversity FILE] [-t|--inputAAMapType DIVERSITY | POSITION] [-u|--nucleotides] [-m|--outputMutCounts FILE] [-s|--outputStabCounts File] [-M|--outputMutDiversityCounts FILE] [-S|--outputStabDiversityCounts FILE] [-j|--outputMutAAUse FILE] [-k|--outputStabAAUse FILE] [-r|--outputRarefaction FILE] [-c|--outputAllChangedAAMap FILE] [-y|--outputImportantChangedAAMap FILE] [-z|--outputUnimportantChangedAAMap FILE]
  Return various information about the relationship between the germline and the clones, most importantly the amino acid maps (nucleotide sequences only)

Available options:
  -h,--help                Show this help text
  -o,--inputOrder ORDER    The order of true diversity
  -i,--inputFasta FILE     The fasta file containing the germlines and clones
  -d,--inputDiversity FILE The csv file containing the diversities at each position (must be generated into a specific format)
  -t,--inputAAMapType DIVERSITY | POSITION Whether to split the amino acid map by position or diversity
  -u,--nucleotides         Whether these sequences are of nucleotides (Codon) or amino acids (AminoAcid)
  -m,--outputMutCounts FILE The output file for the changed amino acid counts
  -s,--outputStabCounts File The output file for the maintained amino acid counts
  -M,--outputMutDiversityCounts FILE The output file for the hanged amino acid diversities
  -S,--outputStabDiversityCounts FILE The output file for the maintained amino acid diversities
  -j,--outputMutAAUse FILE The output file for the specific changed amino acids used
  -k,--outputStabAAUse FILE The output file for the specific maintained amino acids used
  -r,--outputRarefaction FILE The output file for the rarefaction curves
  -c,--outputAllChangedAAMap FILE The output file for the map of all changed amino acids
  -y,--outputImportantChangedAAMap FILE The output file for the map of important changed amino acids
  -z,--outputUnimportantChangedAAMap FILE The output file for the map of unimportant changed amino acids
```

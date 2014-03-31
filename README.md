#aa-mapper

**Gregory W. Schwartz**

This program will take a CLIP fasta file and find the amino acid usage for
mutations from the germline to the clone sequences.

To install:
```
cabal configure
cabal build
cabal install
```

```
aa-mapper, Gregory W. Schwartz

Usage: aa-mapper [-o|--input-order ORDER] [-i|--input-fasta FILE] [-d|--input-diversity FILE] [-t|--input-AA-map-type DIVERSITY | POSITION] [-u|--nucleotides] [-n|--no-mutations] [-m|--output-mut-counts FILE] [-s|--output-stab-counts File] [-M|--output-mut-diversity-counts FILE] [-S|--output-stab-diversity-counts FILE] [-j|--output-mut-AA-use FILE] [-k|--output-stab-AA-use FILE] [-r|--output-rarefaction FILE] [-c|--output-all-changed-AA-map FILE] [-y|--output-important-changedAA-map FILE] [-z|--output-unimportant-changed-AA-map FILE]
  Return various information about the relationship between the germline and the clones, most importantly the amino acid maps (nucleotide sequences only)

Available options:
  -h,--help                Show this help text
  -o,--input-order ORDER   The order of true diversity
  -i,--input-fasta FILE    The fasta file containing the germlines and clones
  -d,--input-diversity FILE The csv file containing the diversities at each position (must be generated into a specific format)
  -t,--input-AA-map-type DIVERSITY | POSITION Whether to split the amino acid map by position or diversity
  -u,--nucleotides         Whether these sequences are of nucleotides (Codon) or amino acids (AminoAcid)
  -n,--no-mutations        Whether to look at the codons from a fasta file, not from a germline to a clone sequence of mutations but rather (ideally) from a germline only
  -m,--output-mut-counts FILE The output file for the changed amino acid counts
  -s,--output-stab-counts File The output file for the maintained amino acid counts
  -M,--output-mut-diversity-counts FILE The output file for the hanged amino acid diversities
  -S,--output-stab-diversity-counts FILE The output file for the maintained amino acid diversities
  -j,--output-mut-AA-use FILE The output file for the specific changed amino acids used
  -k,--output-stab-AA-use FILE The output file for the specific maintained amino acids used
  -r,--output-rarefaction FILE The output file for the rarefaction curves
  -c,--output-all-changed-AA-map FILE The output file for the map of all changed amino acids
  -y,--output-important-changedAA-map FILE The output file for the map of important changed amino acids
  -z,--output-unimportant-changed-AA-map FILE The output file for the map of unimportant changed amino acids
```

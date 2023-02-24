# alnwrapper
Contains various programs that create sequence alignments. 

## Mode

### Pair-wise sequence alignment

#### Input

- single sequence 1 (FASTA)
- single sequence 2 (FASTA)

Sequences can be nucleotide or protein. Each FASTA file should includes only one sequence.

#### Output

- alignment (FASTA)

## How to run

```shell
python3 alnwrapper.py [aligner] [input1] [input2] [output] [args]
```

## Pair-wise sequence alignment (Batch mode)

### Input

- many sequences 1 (FASTA)
- many sequences 2 (FASTA)

This mode does all vs all comparison. Each entry in input 1 is aligned to each in input 2. 

### Output

- N x M alignments (FASTA)

### How to run

```shell
python3 alnwrapper.py [aligner] [input1] [input2] [output] [args]
```

## Available aligner

- MAFFT
- Clustal Omega
- MUSCLE (v5.1)


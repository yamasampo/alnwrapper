# alnwrapper
Contains various programs that create sequence alignments. 

## Mode

### Pair-wise sequence alignment

#### Input

- single sequence 1 (FASTA)
- single sequence 2 (FASTA)

Sequences can be nucleotide or protein. Each FASTA file should includes only one 
sequence.

#### Output

- alignment (FASTA)

### How to run

```shell
python3 alnwrapper.py [aligner] [input1] [input2] [output] [args]
```

## Pair-wise sequence alignment (Batch mode)

### Input

- FASTA file for input 1
- FASTA file for input 2

These FASTA files can contain many sequences (N and M sequences for input 1 and 
2, respectively). This mode does all vs all comparison. Each entry in input 1 is 
aligned to each in input 2. 

### Output

- folder path

The folder will contain N (# of sequences in input 1) x M (# of sequences in 
input 2) alignment files (FASTA). 

### How to run

```shell
python3 alnwrapper.py [aligner] [input1] [input2] [output] [args]
```

## Available aligner

- MAFFT
- Clustal Omega
- MUSCLE (v5.1)


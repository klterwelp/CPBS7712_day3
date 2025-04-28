# CPBS7712_day3

This repository holds the code for the CPBS7712 programming assignment. The goal of this programming assignment is to identify the sequence context of a given query sequence through targeted assembly.

## Description

The algorithm constructs a directed graph from the input reads, identifies the k-mers associated with the query sequence, and assembles the longest contig that includes the query sequence. Then the reads are aligned to the contig using a simple seed-and-extend approach.

## Installation Instructions

Before you can use this repository, you'll need to have conda installed. Follow [conda installation instructions if you need help setting up conda.](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

1. Clone github repository.

  ```sh
  git clone https://github.com/klterwelp/CPBS7712_day3.git
  ```

2. Set up conda environment

```sh
# move into repo environment folder
cd CPBS7712_day3/envs
# create conda environment
conda env create -f dbg.yaml 
# activate conda environment
conda activate dbg_env
```

## Example Usage

### Assembly and Alignment

1. Place your sequence reads and the query sequence into the `data` folder, ensuring they are both fasta files. You can only search one query sequence at a time.

```sh
# example of move commands
mv reads.fasta data/
mv query.fasta data/
```

2. Run the assembly and alignment script.

```sh
python src/assemble_align.py \
--read_path data/reads.fasta \
--query_path data/query.fasta \
--output_assembly_path results/contig.fasta \
--output_alignment_path results/alignments.tsv \
--k_assembly 31 \
--k_alignment 20 \
--interval 1 \
--min_count 1
```

This will run the assembly and alignment pipeline. The reads in `read_path` will be used to construct a de Bruijn Graph, while the query sequence in `query_path` will be used to identify the final assembly. The reads will be aligned against the longest constructed contig and the output of the assembly and alignments will be written to the specified output paths.

### Parameters

- `--read_path`: Path to the input reads file in fasta format.
- `--query_path`: Path to the input query sequence file in fasta format.
- `--output_assembly_path`: Path to the output contig file in fasta format.
- `--output_alignment_path`: Path to the output alignment file in TSV format.
- `--k_assembly`: Length of the k-mers used for assembly. Should be an odd number to avoid palindromic k-mers. Default is 31.
- `--k_alignment`: Length of the k-mers used for alignment. Default is 20.
- `--interval`: Interval used to generate the seed k-mers for alignment. Default is to calculate based on the read length.
- `--min_count`: Minimum number of times a given k-mer must appear in the reads to be considered for assembly. Default is 1.

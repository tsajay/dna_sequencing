# DNA Sequencing Toolkit

A collection of algorithms for DNA sequence alignment, approximate pattern matching, and *de novo* genome assembly — built alongside Ben Langmead's [Algorithms for DNA Sequencing](https://www.coursera.org/learn/dna-sequencing) course.

## Algorithms

### Pattern Matching

| Algorithm | Module | Description |
|-----------|--------|-------------|
| Naive | `naive_sequencer.py` | Brute-force O(nm) exact and approximate matching, with optional reverse complement search |
| Boyer-Moore | `bm_preproc.py` | Efficient exact matching using the bad character rule and good suffix rule (Z-algorithm preprocessing) |
| Pigeonhole | `pigeon_hole.py` | Approximate matching allowing up to *m* mismatches by splitting the pattern into segments and querying a k-mer index |

### Indexing

| Structure | Module | Description |
|-----------|--------|-------------|
| K-mer Index | `kmer_index.py` | Sorted list of (k-mer, offset) pairs with binary search lookup |
| Subsequence Index | `kmer_index.py` | Indexes spaced subsequences (every *ival*-th character) for sparser, more selective hits |

### Edit Distance & Approximate Matching

| Algorithm | Module | Description |
|-----------|--------|-------------|
| Edit Distance | `overlap_finder.py` | Classic DP edit distance between two strings |
| Approximate Match | `overlap_finder.py` | Finds the minimum edit distance of a pattern against all substrings of a reference |

### Genome Assembly
~~~~
| Algorithm | Module | Description |
|-----------|--------|-------------|
| Overlap Graph | `overlap_finder.py` | Builds a k-mer dictionary across all reads and computes pairwise suffix-prefix overlaps |
| Greedy SCS | `overlap_finder.py` | Greedy shortest common superstring — repeatedly merges the read pair with the largest overlap |
| De Bruijn Graph | `overlap_finder.py` | Constructs a De Bruijn graph ((k-1)-mer nodes, k-mer edges) from reads |

## Project Structure

```
dna_sequencing/
├── dna_sequencer.py            # Main CLI driver
├── sequence_file_reader.py     # FASTA/FASTQ parsing and reverse complement
├── naive_sequencer.py          # Naive exact and approximate matching
├── bm_preproc.py               # Boyer-Moore preprocessing and matching
├── kmer_index.py               # K-mer and subsequence indexing
├── pigeon_hole.py              # Pigeonhole approximate matching
├── overlap_finder.py           # Edit distance, overlap graphs, and assembly
│
├── chr1.GRCh38.excerpt.fasta   # Human chromosome 1 excerpt (GRCh38)
├── lambda_virus.fa             # Lambda phage genome
├── phix.fa                     # PhiX174 bacteriophage genome
├── ERR037900_1.first1000.fastq # First 1000 Illumina reads
├── ERR266411_1.for_asm.fastq   # Reads for assembly exercises
├── ads1_week4_reads.fq         # Course week 4 assignment reads
├── ref.fasta                   # Small test reference
├── sample.fasta                # Small test sample
└── test.fastq                  # Small test FASTQ
```

## Usage

All algorithms are accessed through `dna_sequencer.py`:

```bash
# Naive exact matching
python dna_sequencer.py -r ref.fasta -s sample.fasta

# Naive matching with up to 2 mismatches and reverse complement
python dna_sequencer.py -r ref.fasta -s sample.fasta -m 2 -c

# Boyer-Moore exact matching
python dna_sequencer.py -r lambda_virus.fa -s sample.fasta -b

# Pigeonhole approximate matching (k-mer size 8, up to 2 mismatches)
python dna_sequencer.py -r chr1.GRCh38.excerpt.fasta -s sample.fasta -p -k 8 -m 2

# Pigeonhole with subsequence index (ival=3)
python dna_sequencer.py -r chr1.GRCh38.excerpt.fasta -s sample.fasta -p -k 8 -m 2 -i 3

# Edit distance / approximate match
python dna_sequencer.py -r ref.fasta -s sample.fasta -e

# Find overlapping reads
python dna_sequencer.py -r ref.fasta -s ERR266411_1.for_asm.fastq -o -k 30

# Greedy shortest common superstring assembly
python dna_sequencer.py -r ref.fasta -s ERR266411_1.for_asm.fastq -g -k 30

# De Bruijn graph assembly
python dna_sequencer.py -r ref.fasta -s ERR266411_1.for_asm.fastq -d -k 30
```

### CLI Flags

| Flag | Description |
|------|-------------|
| `-r` / `--ref` | Reference genome file (FASTA) |
| `-s` / `--samp` | Sample/pattern file (FASTA or FASTQ for assembly) |
| `-m` / `--mismatchesAllowed` | Max mismatches for naive or pigeonhole matching |
| `-k` / `--kmerSize` | K-mer size (required for pigeonhole, overlap, and assembly) |
| `-i` / `--ival` | Subsequence interval for pigeonhole |
| `-b` / `--boyerMoore` | Use Boyer-Moore matching |
| `-p` / `--pigeonHole` | Use pigeonhole approximate matching |
| `-e` / `--editDistance` | Compute approximate edit distance |
| `-c` / `--reverseComplement` | Include reverse complement (naive only) |
| `-t` / `--textAlphabet` | Use lowercase alphabet for Boyer-Moore |
| `-o` / `--overlapReadFinder` | Find pairwise read overlaps |
| `-g` / `--greedySCS` | Greedy shortest common superstring assembly |
| `-d` / `--deBruijnSCS` | De Bruijn graph assembly |

## Requirements

Python 2.7+ or Python 3.x (standard library only — no external dependencies).

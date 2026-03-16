# DNA Sequencing Toolkit — Detailed Module Reference

A detailed walkthrough of every Python file in this project, covering purpose, functions, classes, algorithms, complexity, and how they interconnect.

---

## Table of Contents

1. [sequence_file_reader.py](#1-sequence_file_readerpy)
2. [naive_sequencer.py](#2-naive_sequencerpy)
3. [naive_sequencer_.py](#3-naive_sequencer_py)
4. [kmer_index.py](#4-kmer_indexpy)
5. [bm_preproc.py](#5-bm_preprocpy)
6. [pigeon_hole.py](#6-pigeon_holepy)
7. [overlap_finder.py](#7-overlap_finderpy)
8. [dna_sequencer.py](#8-dna_sequencerpy)
9. [Module Dependency Graph](#9-module-dependency-graph)

---

## 1. `sequence_file_reader.py`

**Purpose:** Low-level I/O utilities for reading genomic data files and computing reverse complements. Imported by most other modules.

### Functions

#### `reverseComplement(s)`

Computes the reverse complement of a DNA string. Each base is replaced by its Watson-Crick complement (A↔T, C↔G), and the string is reversed.

- **Input:** A string of characters from `{A, C, G, T, N}`
- **Output:** The reverse complement string
- **Example:** `reverseComplement("AGTC")` → `"GACT"`
- **Complexity:** O(n) where n = len(s)

#### `readGenome(filename)`

Reads a FASTA file and concatenates all sequence lines into a single string, skipping header lines (lines starting with `>`).

- **Input:** Path to a `.fasta` or `.fa` file
- **Output:** A single string containing the entire genome
- **Format expected:**

```
>header line (ignored)
ATCGATCGATCG
ATCGATCGATCG
```

#### `readFastq(filename)`

Reads a FASTQ file and returns separate lists of sequences and quality strings. FASTQ files have a 4-line-per-record structure:

```
@read_name        (line 1 — skipped)
ATCGATCG          (line 2 — sequence)
+                 (line 3 — skipped)
IIIIIIII          (line 4 — quality scores)
```

- **Returns:** `(sequences, qualities)` — two parallel lists of strings

---

## 2. `naive_sequencer.py`

**Purpose:** Modular naive pattern matching algorithms. This is the refactored version that imports `sequence_file_reader` for file I/O.

### Functions

#### `naive(p, t)`

Exact pattern matching by brute force. Slides the pattern `p` across every position in text `t` and checks character by character.

- **Complexity:** O(n \* m) where n = len(t), m = len(p)
- **Returns:** List of match positions (0-indexed offsets into `t`)

#### `naiveWithNMismatches(p, t, n)`

Same sliding-window approach but tolerates up to `n` mismatches. A position counts as a match if the number of character differences is ≤ `n`.

- **Returns:** List of match positions

#### `naiveWithReverseComplement(p, t)`

Searches for both the pattern and its reverse complement. Runs `naive()` on `p` and (if different) on `reverseComplement(p)`, then merges and sorts results.

- **Returns:** Sorted list of all match positions

#### `naiveWithReverseComplement2(p, t)`

An alternative single-pass implementation that checks both the pattern and its reverse complement within the same loop iteration rather than making two separate passes.

- **Returns:** List of match positions

---

## 3. `naive_sequencer_.py`

**Purpose:** An earlier, self-contained version of the naive sequencer. Contains its own copies of `reverseComplement`, `readGenome`, `readFastq`, argument parsing, and a `main()` entry point — all in one file.

### Key Differences from `naive_sequencer.py`

- **Self-contained:** Does not import `sequence_file_reader`; duplicates those functions locally.
- **`naiveWithNMismatches`** returns a 3-tuple `(occurrences, alignments, charComparisons)` — tracking alignment and comparison counts for performance analysis.
- **`naiveWithReverseComplement(p, t, m=0)`** accepts an optional mismatch parameter `m` and delegates to `naiveWithNMismatches` instead of `naive`.
- **Has its own CLI** via `parseArgs()` and `main()` — accepts `-r` (reference), `-s` (sample), and `-m` (mismatches).

This file is largely superseded by the modular design (`naive_sequencer.py` + `sequence_file_reader.py` + `dna_sequencer.py`).

---

## 4. `kmer_index.py`

**Purpose:** Provides two index data structures for fast substring/subsequence lookup in a text. Used by `pigeon_hole.py` to accelerate approximate matching. Authored by Ben Langmead.

### Class: `Index`

A **k-mer index** — a sorted array of `(k-mer, offset)` tuples with binary search for queries.

#### `__init__(self, t, k)`

Builds the index from text `t`:
1. Extracts every consecutive k-mer: `t[i:i+k]` for all valid `i`
2. Stores as `(k-mer_string, offset)` pairs
3. Sorts alphabetically by k-mer

- **Space:** O(n) entries where n = len(t) - k + 1
- **Build time:** O(n \* k + n log n) — extraction + sort

#### `query(self, p, idx=0)`

Looks up the `idx`-th k-mer of pattern `p` in the index using `bisect.bisect_left` (binary search), then collects all entries with a matching k-mer.

- **Input:** Pattern `p`, partition index `idx` (which k-mer-sized chunk of `p` to query)
- **Returns:** List of offsets in `t` where that k-mer occurs
- **Query time:** O(k log n + hits)

### Class: `SubseqIndex`

A **subsequence index** — instead of indexing consecutive characters, it indexes characters spaced `ival` positions apart.

#### `__init__(self, t, k, ival)`

- Extracts subsequences: `t[i : i+span : ival]` where `span = 1 + ival * (k - 1)`
- **Example:** `SubseqIndex("ATAT", 2, 2)` extracts `("AA", 0)` and `("TT", 1)`
- Fewer index hits than `Index` (more selective), but the index covers wider spans of the text

#### `query(self, p, idx=0)`

Same binary search approach but queries a spaced subsequence of `p` starting at offset `idx`.

---

## 5. `bm_preproc.py`

**Purpose:** Full implementation of the **Boyer-Moore** string matching algorithm, including all preprocessing tables and the matching loop. Authored by Ben Langmead. Includes unit tests.

### Constants

- `LOWERCASE_ALPHABET = 'abcdefghijklmnopqrstuvwxyz '` — for matching general English text
- `DNA_SEQUENCING_ALPHABET = 'ACGT'` — for matching DNA sequences

### Preprocessing Functions

These implement the theory from Gusfield's *"Algorithms on Strings, Trees, and Sequences"*.

#### `z_array(s)` — Z-Algorithm (Gusfield Theorem 1.4.1)

Computes the Z-array for string `s`. `Z[i]` is the length of the longest substring starting at position `i` that matches a prefix of `s`.

- **Complexity:** O(n) amortized
- **Used by:** `n_array` (and transitively by all good suffix preprocessing)

#### `n_array(s)` — N-Array (Gusfield Theorem 2.2.2)

`N[j]` = length of the longest suffix of `s[0..j]` that is also a suffix of `s`. Computed by running the Z-algorithm on the reversed string.

#### `big_l_prime_array(p, n)` — L' Array

`L'[i]` = largest index `j < n` such that `N[j] = |P[i:]|`. Used in the strong good suffix rule.

#### `big_l_array(p, lp)` — L Array

`L[i]` = largest index `j < n` such that `N[j] >= |P[i:]|`. The weak version of the good suffix rule.

#### `small_l_prime_array(n)` — l' Array (Gusfield Theorem 2.2.4)

`l'[i]` = length of the largest suffix of `P[i:]` that is also a prefix of `P`. Handles the case when no good suffix shift via L/L' is available.

#### `good_suffix_table(p)`

Convenience function that computes and returns all three tables: `(L', L, l')`.

#### `good_suffix_mismatch(i, big_l_prime, small_l_prime)`

Given a mismatch at offset `i`, returns the shift amount dictated by the good suffix rule.

#### `good_suffix_match(small_l_prime)`

Returns the shift amount when the full pattern matches (to continue searching for more occurrences).

#### `dense_bad_char_tab(p, amap)`

Builds the **bad character table** — a 2D array indexed by `[position][character]` giving the rightmost occurrence of each character to the left of that position.

### Class: `BoyerMoore`

Encapsulates all preprocessing for a given pattern.

#### `__init__(self, p, alphabet='ACGT')`

Builds the bad character table and good suffix tables for pattern `p` over the given alphabet.

#### `bad_character_rule(self, i, c)`

Returns how many positions to skip when character `c` in the text causes a mismatch at pattern position `i`.

#### `good_suffix_rule(self, i)`

Returns how many positions to skip based on the good suffix rule when a mismatch occurs at position `i`.

#### `match_skip(self)`

Returns how many positions to skip after a full match is found.

### Matching Function

#### `boyer_moore(p, p_bm, t)`

Runs Boyer-Moore matching with pattern `p` against text `t` using the preprocessed `BoyerMoore` object `p_bm`.

- Compares characters **right-to-left** within each alignment
- On mismatch, shifts by `max(1, bad_character_skip, good_suffix_skip)`
- On full match, shifts by the match skip amount
- **Returns:** `(occurrences, alignments, comparisons)` — match positions plus performance counters
- **Best case:** O(n/m) — sublinear, skips large portions of text
- **Worst case:** O(n \* m) — but rare in practice

### Unit Tests: `TestBoyerMoorePreproc`

Validates Z-array, N-array, L', L, l', and good suffix shift computations against known examples like `"abb"`, `"abababab"`, `"abracadabra"`, and `"GGTAGGT"`.

---

## 6. `pigeon_hole.py`

**Purpose:** Implements the **pigeonhole principle** for approximate pattern matching with up to `m` mismatches. Uses `kmer_index.py` for fast candidate lookup.

### Pigeonhole Principle

If a pattern of length `n` has at most `m` mismatches with some substring of the text, and we split the pattern into `m + 1` non-overlapping segments, then **at least one segment must match exactly**. This converts an approximate matching problem into exact index lookups + verification.

### Class: `PigeonHole`

#### `__init__(self, t, p, k, m, i=0)`

- `t` — text (reference genome)
- `p` — pattern to search for
- `k` — k-mer length for the index
- `m` — maximum mismatches allowed
- `i` — subsequence interval (`0` = use `Index`, `>0` = use `SubseqIndex`)

Builds a k-mer or subsequence index on `t` at construction time.

#### `naiveWithNMisMatches(self, p, t)`

Verifies whether pattern `p` matches text region `t` with at most `m` mismatches. Used to verify candidates found via the index.

- **Returns:** `(is_match, num_comparisons)`

#### `getMatches(self)`

The main method. Splits the pattern into slices, queries the index for each slice, then verifies each candidate position:

1. For each of the `numSlices` partitions of `p`, query the index
2. For each hit, compute the start position in `t` where the full pattern would begin
3. If not already verified, run `naiveWithNMisMatches` to confirm
4. Collect confirmed matches in a set (deduplication)

- **Returns:** `(sorted_matches, num_queries, num_comparisons)`

---

## 7. `overlap_finder.py`

**Purpose:** The largest module. Implements edit distance, approximate matching via DP, overlap graph construction, and multiple genome assembly strategies.

### Edit Distance Functions

#### `editDistance(x, y)`

Classic dynamic programming edit distance (Levenshtein distance). Computes the minimum number of single-character insertions, deletions, and substitutions to transform `x` into `y`.

- Builds a `(len(x)+1) x (len(y)+1)` matrix `D`
- `D[0][j] = j` and `D[i][0] = i` (both rows/columns initialized with gap penalties)
- Recurrence: `D[i][j] = min(D[i][j-1]+1, D[i-1][j]+1, D[i-1][j-1] + (0 if match else 1))`
- **Returns:** `D[-1][-1]` — the edit distance
- **Complexity:** O(n \* m) time and space

#### `minEditDistance(ref, samp)`

Slides a window of `len(samp)` across `ref` and returns the minimum `editDistance` found. A brute-force approach to approximate matching.

- **Complexity:** O(n \* m^2) where n = len(ref), m = len(samp)

#### `approximateMatch(ref, samp)`

A DP-based approximate matching algorithm. The key difference from `editDistance` is that the **first row is initialized to all zeros** (`D[0][j] = 0` instead of `D[0][j] = j`). This means the pattern can start matching at any position in the reference without penalty.

- **Returns:** `min(D[-1])` — the minimum edit distance of the pattern against any substring of the reference
- **Complexity:** O(n \* m) — much more efficient than `minEditDistance`

### Overlap Functions

#### `overlap(a, b, min_length=3)`

Computes the length of the longest suffix of read `a` that matches a prefix of read `b`, requiring at least `min_length` characters.

- Uses `str.find()` to locate occurrences of `b[:min_length]` within `a`
- For each hit, checks if `b.startswith(a[start:])`
- **Returns:** Overlap length, or `0` if none found

#### `overlap_opt(a, b, readLength, min_length=3)`

Optimized version of `overlap` that accepts a precomputed `readLength` to avoid recomputing `len(a)`.

#### `overlapMinCheck(a, b, min_length=3)`

A quick boolean check — returns `1` if any suffix-prefix overlap of at least `min_length` exists, `0` otherwise. Only checks the first occurrence (no loop).

### Overlap Graph Construction

#### `addKmersForRead(globalKmers, sample, k)`

Extracts all k-mers from a read and adds them to a global dictionary mapping each k-mer to the set of reads that contain it.

#### `findOverlapsForReads(reads, k)`

Builds a complete overlap graph from a list of reads:

1. **Build global k-mer dictionary:** For every read, extract all k-mers and record which reads contain each k-mer
2. **Build candidate sets:** For each read, collect all other reads that share at least one k-mer (potential overlaps)
3. **Compute overlaps:** For each read and each of its candidates, compute the actual suffix-prefix overlap using `overlap_opt`

- **Returns:** List of `(read_a, read_b, overlap_length)` triples
- Logs timing for each phase

#### `findSourceNodes(overlaps)`

Extracts the set of all reads that appear as a source (left side) in an overlap pair.

### Genome Assembly Algorithms

#### `pickMaximalOverlap(reads, k)`

Finds the pair of reads with the largest suffix-prefix overlap by calling `findOverlapsForReads` and scanning for the maximum.

- **Returns:** `(all_overlaps, best_read_a, best_read_b, best_overlap_length)`

#### `greedy_scs(reads, k)`

**Greedy Shortest Common Superstring (SCS)** assembly:

1. Find the pair with the largest overlap
2. Merge them: `reada + readb[olen:]`
3. Replace the two reads with the merged read
4. Repeat until no overlaps remain
5. Concatenate all remaining reads

This is a greedy heuristic — it does not guarantee the true shortest common superstring (which is NP-hard), but produces a reasonable approximation.

- **Returns:** The assembled genome string

#### `greedy_scs_experminetal(reads, k)`

An experimental variant that avoids recomputing the full overlap graph after each merge. Instead, it updates the overlap list in-place by:
- Removing edges involving the merged reads
- Replacing references to the old reads with the new merged read

- **Returns:** `(scs, num_merged)`

#### `scs(ss)`

**Brute-force SCS** via all permutations. Tries every permutation of the input strings, computes the superstring for each by greedily overlapping adjacent strings, and returns all shortest ones.

- **Complexity:** O(n! \* n \* m) — only feasible for very small inputs (< ~10 strings)
- **Returns:** List of all shortest common superstrings

#### `de_bruijn_graph(reads, k)`

Constructs a **De Bruijn graph** from reads:
- **Nodes:** All distinct (k-1)-mers
- **Edges:** Each k-mer becomes a directed edge from its (k-1)-length prefix to its (k-1)-length suffix

- **Returns:** `(nodes_set, edges_list)`

#### `de_bruin_graph_scs(reads, k)`

Partially implemented De Bruijn graph assembly. Currently builds the graph and writes the edge list to `debruijn.out` but does not yet traverse the graph to find an Eulerian path.

- **Returns:** Empty string (assembly not completed)

---

## 8. `dna_sequencer.py`

**Purpose:** The main CLI entry point that ties all modules together. Parses command-line arguments, loads data, dispatches to the selected algorithm, and prints results.

### Imports

```
naive_sequencer  → naive
sequence_file_reader → reader
bm_preproc → bm
pigeon_hole → ph
overlap_finder → of
```

### `parseArgs()`

Defines and validates all CLI flags:

| Flag | Type | Description |
|------|------|-------------|
| `-r` / `--ref` | string | Reference genome file (default: `ref.fasta`) |
| `-s` / `--samp` | string | Sample/pattern file (default: `sample.fasta`) |
| `-m` / `--mismatchesAllowed` | int | Max mismatches (default: 0) |
| `-k` / `--kmerSize` | int | K-mer size for index-based algorithms |
| `-i` / `--ival` | int | Subsequence interval for pigeonhole |
| `-b` / `--boyerMoore` | flag | Use Boyer-Moore |
| `-p` / `--pigeonHole` | flag | Use pigeonhole |
| `-e` / `--editDistance` | flag | Compute edit distance |
| `-c` / `--reverseComplement` | flag | Include reverse complement (naive only) |
| `-t` / `--textAlphabet` | flag | Use lowercase alphabet (Boyer-Moore only) |
| `-o` / `--overlapReadFinder` | flag | Find read overlaps |
| `-g` / `--greedySCS` | flag | Greedy SCS assembly |
| `-d` / `--deBruijnSCS` | flag | De Bruijn assembly |

**Validation rules:**
- Reference and sample files must exist
- Boyer-Moore cannot be combined with mismatches or reverse complement
- Text alphabet flag only applies to Boyer-Moore
- Pigeonhole, overlap, greedy SCS, and De Bruijn all require `kmerSize > 0`

### `main()`

Dispatches to the selected algorithm based on flags:

- **No algorithm flags** → Naive matching (with optional reverse complement and mismatches)
- **`-b`** → Boyer-Moore exact matching
- **`-p`** → Pigeonhole approximate matching
- **`-e`** → Approximate match via edit distance DP
- **`-o`** → Overlap read finder (reads sample as FASTQ)
- **`-g`** → Greedy SCS assembly (reads sample as FASTQ, writes output to `greedy_assembly_genome.fa`)
- **`-d`** → De Bruijn graph assembly (reads sample as FASTQ)

Prints a summary of all metrics: match positions, alignment count, query count, character comparisons, duration, edit distance, overlap count, source count, and SCS length.

---

## 9. Module Dependency Graph

```
dna_sequencer.py
├── sequence_file_reader.py    (file I/O)
├── naive_sequencer.py         (naive matching)
│   └── sequence_file_reader.py
├── bm_preproc.py              (Boyer-Moore)
├── pigeon_hole.py             (pigeonhole)
│   └── kmer_index.py          (k-mer / subsequence indexing)
└── overlap_finder.py          (edit distance, overlaps, assembly)
```

`dna_sequencer.py` is the orchestrator. Each algorithm module is self-contained (apart from shared I/O via `sequence_file_reader`) and can also be imported independently in other scripts or notebooks.

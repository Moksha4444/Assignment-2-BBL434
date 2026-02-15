# Assignment-2-BBL434
# Affine Gapped sequence alignment tool

A computational tool for pairwise global alignment of DNA/RNA sequences using the Needleman-Wunsch algorithm with affine gap penalties.

**Author:** Tamme Mokshagna (2023BB11055)  
**Course:** BBL434-Assignment 2

---
## 📋 Overview

This Python-based tool performs global sequence alignment by:

Reading sequences from standard FASTA format files
Computing optimal alignment using dynamic programming with affine gap penalties
Displaying results with detailed statistics including alignment score, identity percentage, and gap analysis

The tool uses classical bioinformatics algorithms including the Needleman-Wunsch algorithm with a three-matrix approach for biologically accurate gap modeling.

---
## 🧬 Biological Background

### What is Sequence Alignment?

Sequence alignment is the arrangement of DNA, RNA, or protein sequences to identify regions of similarity. These similarities may indicate:

- **Evolutionary relationships** between organisms
- **Functional regions** in genes
- **Structural similarities** in proteins
- **Point mutations and indels** (insertions/deletions)

### Global vs Local Alignment

**Global Alignment (Needleman-Wunsch):** Aligns sequences end-to-end
- Best for: Sequences of similar length and high similarity
- Use cases: Comparing orthologous genes, closely related species

**Local Alignment (Smith-Waterman):** Finds best matching regions
- Best for: Sequences of different lengths or distantly related
- Use cases: Database searches, finding conserved domains

**This tool implements Global Alignment.**

### How This Tool Works

The tool uses a sophisticated three-matrix dynamic programming approach:
- **M matrix:** Tracks match/mismatch alignments
- **X matrix:** Tracks gaps in sequence 2 (horizontal gaps)
- **Y matrix:** Tracks gaps in sequence 1 (vertical gaps)

This approach models **affine gap penalties**, where opening a gap is more costly than extending it—matching biological reality where insertions/deletions typically occur in clusters.

---

## 🛠️ Features

### Core Capabilities

- **Needleman-Wunsch Algorithm:** Industry-standard global alignment
- **Affine Gap Penalties:** Biologically accurate gap scoring (gap_open + gap_extend)
- **FASTA Format Support:** Reads standard bioinformatics file format
- **Customizable Scoring:** User-defined match, mismatch, and gap parameters
- **Detailed Statistics:** Identity percentage, gap count, alignment score
- **Visual Alignment:** Easy-to-read output with match indicators
- **Error Handling:** Robust file validation and error reporting

### Scoring System

The tool uses a flexible scoring matrix:

| Parameter | Default | Description |
|-----------|---------|-------------|
| Match | +5 | Reward for identical nucleotides |
| Mismatch | -4 | Penalty for non-matching nucleotides |
| Gap Open | -12 | Penalty for starting a new gap |
| Gap Extend | -4 | Penalty for each position in a gap |

**Gap Cost Formula:**
```
Total gap cost = gap_open + (gap_length × gap_extend)

Example:
- Gap of length 1: -12 + (1 × -4) = -16
- Gap of length 5: -12 + (5 × -4) = -32
```

---

## 🚀 Usage

### Basic Command

```bash
python3 sequence_align_final.py
```

---

## 🔬 Algorithm Details

### 1. Needleman-Wunsch with Affine Gaps

The algorithm uses three dynamic programming matrices:

**Matrix Definitions:**
- **M[i,j]:** Best score where seq1[i] aligns to seq2[j] (match/mismatch)
- **X[i,j]:** Best score ending with a gap in seq2 (insertion in seq1)
- **Y[i,j]:** Best score ending with a gap in seq1 (insertion in seq2)

**Recurrence Relations:**

```
M[i,j] = max(
    M[i-1,j-1] + s(seq1[i], seq2[j]),
    X[i-1,j-1] + s(seq1[i], seq2[j]),
    Y[i-1,j-1] + s(seq1[i], seq2[j])
)

X[i,j] = max(
    M[i-1,j] + gap_open + gap_extend,
    X[i-1,j] + gap_extend
)

Y[i,j] = max(
    M[i,j-1] + gap_open + gap_extend,
    Y[i,j-1] + gap_extend
)
```

Where `s(a,b)` is the match score if a=b, else mismatch penalty.

**Initialization:**
```
M[0,0] = 0
X[i,0] = gap_open + i × gap_extend  (for all i > 0)
Y[0,j] = gap_open + j × gap_extend  (for all j > 0)
```

### 2. Traceback

Starting from the highest scoring position at [n,m], trace back through matrices:

1. If in **M state**: Consume both sequences, move diagonally
2. If in **X state**: Consume seq1, add gap to seq2, move up
3. If in **Y state**: Consume seq2, add gap to seq1, move left
4. Continue until reaching [0,0]

The traceback produces the aligned sequences with gaps represented as `-`.

### 3. Scoring

**Alignment Score Formula:**
```
Score = Σ(matches × match_score) + 
        Σ(mismatches × mismatch_penalty) + 
        Σ(gap_costs)
```

**Statistics Calculated:**
- **Identities:** Number of matching positions
- **Gaps:** Number of gap characters in alignment
- **Identity Percentage:** (Identities / Alignment Length) × 100

---

## 📈 Output Format

### Alignment Display

```
======================================================================
Alignment Score: 226.00
Length: 116 bp
Identities: 78/116 (67.2%)
Gaps: 38/116 (32.8%)
======================================================================

Seq1: ATGCGATACGCTTGAACCTAGCTAAGCTTAAGGCTTAGCGTATCGATCGATCGTAGCTAG
      ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Seq2: ATGCGATACGCTTGAACCTAGCTAAGCTTAAGGCTTAGCGTATCGATCGATCGTAGCTAG

Seq1: CTAGCTAGCTGACTGACT--------------------------------------
      ||||||||||||||||||                                      
Seq2: CTAGCTAGCTGACTGACTAAAATTTTTCCCCCGGGGGTATATATATGCGCGCGCGC
```

**Match Indicators:**
- `|` = Match (identical nucleotides)
- `.` = Mismatch (different nucleotides aligned)
- ` ` = Gap (space for insertion/deletion)

**Sequences are displayed in 60-character chunks for readability.**

---

## 🏗️ Code Architecture

### Module Structure

```
sequence_align_final.py
│
├── read_fasta(filepath)
│   └── Reads FASTA files and returns uppercase sequence
│
├── needleman_wunsch(seq1, seq2, match, mismatch, gap_open, gap_extend)
│   ├── Initialize three matrices (M, X, Y)
│   ├── Initialize traceback matrices
│   ├── Fill matrices using dynamic programming
│   └── Return matrices, score, and end position
│
├── traceback(seq1, seq2, M, X, Y, trace_M, trace_X, trace_Y, end_info)
│   ├── Start from end position
│   ├── Follow traceback pointers
│   └── Reconstruct aligned sequences
│
├── print_alignment(seq1, seq2, score)
│   ├── Calculate statistics
│   ├── Generate match line
│   └── Display formatted alignment
│
└── main()
    ├── Get user input
    ├── Load sequences
    ├── Get scoring parameters
    ├── Perform alignment
    └── Display results
```

### Complexity Analysis

**Time Complexity:** O(n × m)
- Where n = length of sequence 1
- Where m = length of sequence 2
- Must fill three matrices of size (n+1) × (m+1)

**Space Complexity:** O(n × m)
- Three score matrices: M, X, Y
- Three traceback matrices

**Example:**
- Aligning 1000 bp × 1000 bp sequences:
  - ~6 million matrix cells to compute
  - ~24 MB memory for matrices
  - Runtime: <1 second on modern hardware

---

## ⚙️ Configuration

### Adjustable Parameters

#### Scoring Parameters

**Match Score** (default: 5)
- Higher values emphasize identity
- Typical range: 1-10

**Mismatch Penalty** (default: -4)
- More negative = stronger penalty for mismatches
- Typical range: -1 to -10

**Gap Opening Penalty** (default: -12)
- Cost to start a new gap
- Should be more negative than gap extend
- Typical range: -10 to -20

**Gap Extension Penalty** (default: -4)
- Cost per position in a gap
- Should be less negative than gap open
- Typical range: -1 to -5

#### When to Adjust Parameters

**High Similarity Expected:**
```
Match: +10, Mismatch: -6, Gap Open: -15, Gap Extend: -2
```

**Distant Sequences:**
```
Match: +3, Mismatch: -2, Gap Open: -8, Gap Extend: -1
```

**Minimize Gaps:**
```
Match: +5, Mismatch: -4, Gap Open: -20, Gap Extend: -5
```

---

## 🧪 Testing

### Test Files Provided

**test_seq1.fasta** (78 bp)
```
>Test_Sequence_1
ATGCGATACGCTTGAACCTAGCTAAGCTTAAGGCTTAGCGTATCGATCGATCGTAGCTAG
CTAGCTAGCTGACTGACT
```

**test_seq2.fasta** (116 bp)
```
>Test_Sequence_2
ATGCGATACGCTTGAACCTAGCTAAGCTTAAGGCTTAGCGTATCGATCGATCGTAGCTAG
CTAGCTAGCTGACTGACTAAAATTTTTCCCCCGGGGGTATATATATGCGCGCGCGC
```

### Expected Output

With default parameters:
```
Alignment Score: 226.00
Length: 116 bp
Identities: 78/116 (67.2%)
Gaps: 38/116 (32.8%)
```

### Validation Against EBI Tool

Compare results with the European Bioinformatics Institute's alignment tool:

1. Visit: https://www.ebi.ac.uk/jdispatcher/psa/emboss_needle
2. Upload your FASTA files
3. Set parameters: Match=5, Mismatch=-4, Gap Open=12, Gap Extend=4
   - **Note:** EBI uses positive values for penalties
4. Compare alignment score and identity percentage

**Your results should match EBI's output.**

---

## ⚠️ Limitations

### Algorithmic Limitations

1. **Global Alignment Only**
   - Aligns sequences end-to-end
   - Not suitable for finding local similarities
   - Use Smith-Waterman for local alignment

2. **No Multiple Alignment**
   - Only aligns two sequences at a time
   - For multiple sequences, use tools like MUSCLE or ClustalW

3. **Simple Scoring**
   - Uses uniform match/mismatch scoring
   - Real-world applications may need substitution matrices (e.g., BLOSUM)

### Computational Limitations

1. **Memory Usage**
   - O(n×m) space requirement
   - Large sequences (>10,000 bp) require significant RAM
   - Example: 10,000 × 10,000 = 600 MB for matrices

2. **Processing Time**
   - Quadratic time complexity
   - Very long sequences may take minutes to align
   - Consider specialized tools (BLAST) for large-scale comparisons

3. **No Parallelization**
   - Single-threaded execution
   - Could be optimized with parallel processing

### Biological Limitations

1. **DNA/RNA Only**
   - Not designed for protein sequences
   - Proteins require amino acid substitution matrices

2. **No Structural Information**
   - Ignores secondary structure
   - Ignores functional domains

3. **Linear Sequences**
   - Does not handle circular DNA (plasmids)
   - No consideration for reverse complement

---

## 📝 File Structure

```
sequence-alignment-tool/
├── sequence_align_final.py    # Main program (180 lines)
├── test_seq1.fasta            # Test sequence 1 (78 bp)
├── test_seq2.fasta            # Test sequence 2 (116 bp)
├── README.md                  # This file
└── Output.txt                 # Example output (optional)
```

---




##  Acknowledgments

- Course instructor for assignment specifications
- Needleman & Wunsch for the foundational algorithm
- Gotoh for the affine gap penalty model
- EBI for providing validation tools
- Biopython community for bioinformatics resources

---


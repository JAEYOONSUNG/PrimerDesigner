# PrimerDesigner <img src="man/figures/logo.png" align="right" height="139" />

[![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-CrisprVerse-green)](https://bioconductor.org/packages/crisprVerse/)

End-to-end CRISPR primer design pipeline for non-model organisms.

Just provide a GenBank accession and a vector `.dna` file — BSgenome building, Bowtie indexing, gRNA library generation, off-target analysis, Golden Gate / Gibson cloning primer design, deletion arm primers, and SnapGene-compatible GenBank output are all handled **in a single function call**.

---

## Features

- **Any organism** — Works with any NCBI genome accession. BSgenome + Bowtie index are **auto-built** on first run (no separate setup)
- **Single function** — `design_grna_and_deletion()` handles everything: genome download → BSgenome → Bowtie → gRNA → primers → Excel → GenBank
- **Genome-wide gRNA library** — Composite scoring (GC%, Tm, position, off-target), methylation site filtering with IUPAC support
- **Off-target analysis** — CrisprVerse `crisprBowtie::runCrisprBowtie` integration, verified for custom BSgenome
- **Vector-based oligo annealing** — Reads actual vector sequence, scans ±20bp for Type IIS enzyme cut sites on both strands, computes real overhangs + fill sequences (not hardcoded)
- **9 Type IIS enzymes** — BbsI, BpiI, BsaI, Eco31I, BsmBI, Esp3I, SapI, BspQI, PaqCI with NEB-verified cut patterns
- **Cut pattern polarity** — Handles both 5' overhang (a < b, e.g. BbsI 2/6) and 3' overhang (a > b) with automatic oligo structure adjustment
- **Gibson Assembly** — Tm-optimized overlap primers for inverse PCR amplification
- **Deletion arm primers** — 4-primer Gibson design for upstream/downstream homology arms, circular genome support
- **Combined construct** — All-in-one vector: gRNA spacer + deletion arms in a single construct
- **GenBank output** — Annotated `.gbk` files for SnapGene (spacer, arms, primer binding sites with Tm)
- **Excel output** — Multi-sheet `.xlsx` with unified primer tables and conditional formatting
- **SnapGene `.dna` support** — Direct reading of SnapGene binary vector files via `reticulate`
- **Resume mode** — Interrupt and restart long genome-wide runs without losing progress

---

## Installation

```r
# 1. Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "crisprVerse", "crisprBase", "crisprDesign", "crisprBowtie",
  "BSgenome", "Biostrings", "GenomicFeatures", "GenomeInfoDb",
  "GenomicRanges", "Rsamtools", "rtracklayer", "Rbowtie"
))

# 2. Install from GitHub
devtools::install_github("JAEYOONSUNG/PrimerDesigner")
```

**Python dependency** (for SnapGene `.dna` files):

```bash
pip install snapgene_reader
```

---

## Quick Start

Just provide a GenBank accession and a vector `.dna` file. Genome download, BSgenome package build, and Bowtie index creation are all handled automatically.

```r
library(PrimerDesigner)

# Single call does everything — includes auto BSgenome/Bowtie build
result <- design_grna_and_deletion(
  genbank_accession    = "GCF_030376765.1",
  locus_tags           = c("QT235_RS00005", "QT235_RS00010"),
  nuclease             = "GeoCas9",
  grna_vector_file     = "~/vectors/pG1Kt-GeoCas9EF.dna",
  grna_start           = 7823,
  grna_end             = 7850,
  deletion_start = 3977,
  deletion_end   = 3987,
  upstream_bp = 500,
  downstream_bp = 500,
  grna_cloning_method = "gibson",
  methylation_patterns = c("GCCAT","TACNNNNNNCTC","RTAYNNNNNCTC","GAGNNNNNNNTGG"),
  output_file          = "~/crispr_project/all_primers.xlsx",
  output_dir           = "~/crispr_project/constructs",
  dir                  = "~/crispr_project"
)

# Access results
result$grna_primers      # gRNA cloning primers (2 per gRNA)
result$deletion_primers  # deletion arm primers (4 per locus_tag)
result$unified_primers   # all primers in unified long-format table
```

**What happens automatically on first run:**
1. Download `.fna` / `.gbff` from NCBI (if not present)
2. Build + load BSgenome package
3. Build Bowtie index
4. Generate gRNA library + off-target analysis
5. Select best gRNA(s) + design cloning primers
6. Design deletion arm primers
7. Export Excel + GenBank output

On subsequent runs, pre-built BSgenome/Bowtie are reused automatically.

---

## Pipeline Overview

```
                    design_grna_and_deletion()
                              │
         ┌────────────────────┼────────────────────┐
         │      All steps run automatically         │
         │                    │                     │
         ▼                    ▼                     ▼
  ┌──────────────┐  ┌─────────────────┐  ┌────────────────┐
  │ Genome Setup │  │ gRNA Library    │  │ Primer Design  │
  │ (auto-build) │  │                 │  │                │
  ├──────────────┤  ├─────────────────┤  ├────────────────┤
  │ .fna download│  │ findSpacers     │  │ gRNA cloning   │
  │ .gbff download│ │ Off-target n0/n1│  │  ├─ Golden Gate│
  │ BSgenome build│ │ Composite score │  │  └─ Gibson     │
  │ Bowtie index │  │ Methylation     │  │ Deletion arms  │
  └──────────────┘  │ filtering       │  │ (4-primer)     │
                    └─────────────────┘  │ Combined GenBank│
                                         │ Excel output   │
                                         └────────────────┘
```

**Internal call chain:**

```
design_grna_and_deletion()
  └─ design_grna_construct()
       └─ run_gRNA_list_generator()
            ├─ download_genbank_fna()       # auto-download if missing
            ├─ download_genbank_gbff()      # auto-download if missing
            ├─ run_bowtie_build()           # auto-build if missing
            ├─ build_bsgenome_from_accession()  # auto-build if missing
            └─ findSpacers + off-target + scoring
  └─ design_deletion_primers()              # 4-primer Gibson deletion arms
  └─ write_combined_construct_genbank()     # SnapGene-compatible GenBank output
```

---

## Golden Gate Oligo Annealing — Vector-Based Overhang Detection

Instead of using hardcoded overhangs, PrimerDesigner scans the actual vector sequence to locate Type IIS enzyme recognition sites and designs oligo annealing primers with real overhangs and fill sequences.

### How it works

1. Load vector file (`.dna` / `.gb`) → scan ±20bp flanking regions around the stuffer
2. Search **both strands** → use only sites that cut toward the stuffer (warns on wrong-orientation sites)
3. Automatically determine overhang type from cut pattern `(a/b)`:
   - `a < b` (e.g. BbsI 2/6) → **5' overhang**, length = b - a
   - `a > b` (hypothetical 6/2) → **3' overhang**, length = a - b, oligo structure auto-adjusted
   - `a == b` → blunt end → incompatible with Golden Gate (warning)
4. Automatically include **fill sequences** between enzyme cut site and stuffer boundary

### Supported Type IIS Enzymes

| Enzyme | Recognition | Cut Pattern | Overhang | Source |
|--------|-------------|-------------|----------|--------|
| BbsI   | `GAAGAC`    | (2/6)       | 4nt 5'   | NEB R0539 |
| BpiI   | `GAAGAC`    | (2/6)       | 4nt 5'   | Thermo ER1011 |
| BsaI   | `GGTCTC`    | (1/5)       | 4nt 5'   | NEB R0535 |
| Eco31I | `GGTCTC`    | (1/5)       | 4nt 5'   | Thermo ER0291 |
| BsmBI  | `CGTCTC`    | (1/5)       | 4nt 5'   | NEB R0739 |
| Esp3I  | `CGTCTC`    | (1/5)       | 4nt 5'   | NEB R0734 |
| SapI   | `GCTCTTC`   | (1/4)       | **3nt** 5' | NEB R0569 |
| BspQI  | `GCTCTTC`   | (1/4)       | **3nt** 5' | NEB R0712 |
| PaqCI  | `CACCTGC`   | (4/8)       | 4nt 5'   | NEB R0745 |

### Oligo design formula

```
F oligo = [LEFT_OH] + [fill_L] + [G] + spacer + [fill_R]
R oligo = [RIGHT_OH] + RC(fill_R) + RC(spacer) + [C] + RC(fill_L)
```

Where `LEFT_OH`, `RIGHT_OH`, `fill_L`, `fill_R` are all derived from the actual vector sequence around the enzyme cut sites. For 3' overhang enzymes, the RC logic is automatically adjusted.

### Usage modes

```r
# Mode 1: Vector-based (recommended) — real overhangs from vector sequence
design_cloning_primers(gRNA_df, vector_file = "vector.dna",
                       start = 8811, end = 8840,
                       cloning_method = "golden_gate", enzyme = "BbsI")

# Mode 2: Default overhangs (no vector needed, e.g. BbsI = CACC/AAAC)
design_cloning_primers(gRNA_df, cloning_method = "golden_gate", enzyme = "BbsI")

# Mode 3: Custom overhangs
design_cloning_primers(gRNA_df, cloning_method = "golden_gate",
                       custom_overhangs = list(F_5prime = "CACC", R_5prime = "AAAC"))
```

---

## Supported Nucleases

| Nuclease | PAM | PAM Side | Spacer Length |
|----------|-----|----------|---------------|
| GeoCas9  | NNNNCAAA | 3' | 21 bp |
| SpCas9   | NGG | 3' | 20 bp |
| FnCas12a | TTTV | 5' | 23 bp |

Custom nucleases can be added via `crisprBase::CrisprNuclease()`.

---

## Key Functions

### Main Pipeline

| Function | Description |
|----------|-------------|
| `design_grna_and_deletion()` | **Full pipeline**: accession → BSgenome → gRNA → cloning primers → deletion arms → Excel → GenBank. All steps automated |

### Step-by-Step (call individually if needed)

| Function | Description |
|----------|-------------|
| `design_grna_construct()` | gRNA selection + cloning primers only (no deletion) |
| `run_gRNA_list_generator()` | gRNA library only (includes auto BSgenome/Bowtie build) |
| `design_cloning_primers()` | Golden Gate (OA) / Gibson (IA) primers only |
| `design_deletion_primers()` | 4-primer deletion arm design (single gene) |
| `batch_deletion_primers()` | Multi-gene deletion primer batch |

### gRNA Library Utilities

| Function | Description |
|----------|-------------|
| `generate_gRNA_for_sequence()` | Generate gRNAs from a DNA sequence string |
| `generate_gRNA_for_locus_tags()` | Generate gRNAs for specific locus tags |
| `build_gRNA_library()` | Build gRNA library (scoring + GenBank merge) |
| `calculate_composite_score()` | Weighted scoring: GC%, Tm, position, off-target |
| `generate_gRNA_names()` | Standardized gRNA naming (e.g. `Geo_dnaA_g1`) |
| `filter_methylation_sites()` | IUPAC methylation motif filtering |
| `nuclease_with_parameter()` | Return nuclease object + PAM/spacer parameters |

### Genome Setup (called automatically — no manual call needed)

| Function | Description |
|----------|-------------|
| `build_bsgenome_from_accession()` | Build BSgenome (auto-called internally, bypasses UCSC validation) |
| `run_bowtie_build()` | Build Bowtie index (auto-called internally) |
| `download_genbank_fna()` | Download .fna (auto-called internally) |
| `download_genbank_gbff()` | Download .gbff (auto-called internally) |
| `forge_BSgenome()` | Low-level BSgenome package builder |

### File I/O

| Function | Description |
|----------|-------------|
| `read_vector_file()` | Read `.dna` / `.gb` / `.fasta` vector files |
| `read_genome_genbank()` | Parse genome `.gbff` → sequence + feature table |
| `write_grna_vector_genbank()` | GenBank output: spacer inserted in vector |
| `write_deletion_genbank()` | GenBank output: deletion construct |
| `format_primer_name()` | Customizable primer naming pattern |

---

## Cloning Methods

### Golden Gate (Oligo Annealing)

Anneals two oligos and ligates directly into the sticky ends created by Type IIS restriction enzyme digestion. When a vector file is provided, overhangs and fill sequences are automatically computed from the actual cut pattern.

```r
result <- design_grna_and_deletion(
  genbank_accession   = "GCF_030376765.1",
  locus_tags          = "QT235_RS00005",
  grna_vector_file    = "~/vectors/pG1Kt-GeoCas9EF.dna",
  grna_start          = 8811,
  grna_end            = 8840,
  grna_cloning_method = "golden_gate",
  grna_enzyme         = "BbsI",
  output_file         = "primers.xlsx",
  output_dir          = "constructs/",
  dir                 = "~/crispr_project"
)
```

### Gibson Assembly (Inverse PCR)

Inserts the spacer via inverse PCR amplification of the vector using Tm-optimized overlap primers.

```r
result <- design_grna_and_deletion(
  genbank_accession   = "GCF_030376765.1",
  locus_tags          = "QT235_RS00005",
  grna_vector_file    = "~/vectors/pG1Kt-GeoCas9EF.dna",
  grna_start          = 8811,
  grna_end            = 8840,
  grna_cloning_method = "gibson",
  tm_target           = 60,
  output_file         = "primers.xlsx",
  output_dir          = "constructs/",
  dir                 = "~/crispr_project"
)
```

### Deletion (4-Primer Gibson)

Clones upstream/downstream homology arms into the vector via Gibson Assembly. Supports circular genomes. Automatically included in `design_grna_and_deletion()`.

---

## Output Files

```
output_dir/
├── all_primers.xlsx                          # Unified primer table (multi-sheet)
│   ├── Sheet: all_primers                    #   Long-format: gRNA + deletion primers
│   ├── Sheet: gRNA_library                   #   Full gRNA library with scores
│   └── Sheet: deletion_detail                #   Wide-format deletion arms
├── QT235_RS00005_g1_combined_construct.gbk   # GenBank: gRNA + arms in one vector
├── QT235_RS00010_g1_combined_construct.gbk
└── .grna_primers_cache.rds                   # Resume cache (auto-generated)
```

Open `.gbk` files in SnapGene to visually inspect spacer, homology arms, and primer binding sites.

---

## Requirements

| Category | Packages |
|----------|----------|
| **R** | >= 4.0 |
| **Bioconductor** | crisprVerse, crisprBase, crisprDesign, crisprBowtie, BSgenome, Biostrings, GenomicFeatures, GenomeInfoDb, GenomicRanges, Rsamtools, rtracklayer, Rbowtie |
| **CRAN** | openxlsx, TmCalculator, glue, dplyr, reticulate, data.table, stringr |
| **System** | Bowtie (installed via `Rbowtie`) |
| **Python** | `snapgene_reader` (optional, for `.dna` file support) |

---

## Citation

If you use PrimerDesigner in your research, please cite:

```
Sung, J.-Y. (2025). PrimerDesigner: End-to-end CRISPR primer design pipeline
for non-model organisms. R package version 1.0.0.
https://github.com/your-username/PrimerDesigner
```

---

## License

MIT © [Jae-Yoon Sung](https://github.com/your-username)

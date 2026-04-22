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
- **Shared primer design across genomes** — `shared_primer_design()` picks one UF/UR/DF/DR oligo set that works across multiple closely related target genomes, falling back to subgroup-splitting per strain only when no shared candidate exists

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

Just provide a GenBank accession (or a local `.gbff` file) and a vector `.dna` file. Genome download, BSgenome package build, and Bowtie index creation are all handled automatically.

```r
library(PrimerDesigner)

result <- design_grna_and_deletion(
  genbank_accession    = "GCF_030376765.1",        # NCBI accession or local .gbff/.gb/.gbk file path
  locus_tags           = c("QT235_RS00005", "QT235_RS00010"),  # single, multiple, or "all"
  nuclease             = "GeoCas9",
  grna_vector_file     = "~/vectors/pG1Kt-GeoCas9EF.dna",
  grna_start           = 8811,
  grna_end             = 8840,
  grna_cloning_method  = "golden_gate",
  grna_enzyme          = "BbsI",
  methylation_patterns = c("GCCAT", "CCANNNNNTTG"),
  output_file          = "all_primers.xlsx",
  output_dir           = "constructs"
)

# Access results
result$grna_primers      # gRNA cloning primers (2 per gRNA)
result$deletion_primers  # deletion arm primers (4 per locus_tag)
result$unified_primers   # all primers in unified long-format table
```

> **`genbank_accession`** — accepts an NCBI accession (e.g. `"GCF_030376765.1"`, auto-downloads) or a local GenBank file path (e.g. `"~/genomes/my_genome.gbff"`). Local `.gbff` / `.gb` / `.gbk` files are auto-detected; FASTA is extracted internally and multi-contig genomes are fully supported.
>
> **`locus_tags`** — `"QT235_RS00005"` for a single gene, `c("QT235_RS00005", "QT235_RS00010")` for multiple genes, or `"all"` to run genome-wide gRNA + deletion primer design for every CDS.

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

## Multi-Genome Shared Primer Design

For knocking out the same gene / ortholog across multiple closely related strains, `shared_primer_design()` picks primers that bind every target genome simultaneously and falls back to per-strain subgroups only where cross-genome sharing is biologically impossible. All outputs land in a single 3-sheet Excel workbook plus an optional single-page summary PDF.

```r
# 1) Discover the target CDS across every genome in a folder.
targets <- find_target_across_genomes(
  genbank_dirs = "path/to/gbk",
  query        = "TIGR02679",        # locus_tag, gene name, or product keyword
  query_type   = "product"
)

# 2) End-to-end design: shared UF/UR/DF/DR + gRNA cloning + check primers +
#    GenBank constructs + Excel workbook.
res <- shared_primer_design(
  target_table   = targets[, c("genome_id", "genbank_file", "locus_tag")],
  nuclease       = "GeoCas9",
  overlap_policy = "strict",         # inner primers pinned at 0 bp
  upstream_bp    = 500, downstream_bp = 500,
  min_arm_bp     = 150L, max_arm_bp   = 3000L,
  grna_vector_file  = "pJET.gb",
  grna_start        = 7823, grna_end  = 7850,
  combined_vector_file          = "pJET.gb",
  combined_grna_start           = 7823, combined_grna_end      = 7850,
  combined_deletion_start       = 3977, combined_deletion_end  = 4986,
  combined_construct_output_dir = "out/constructs",
  output_file                   = "out/design.xlsx"
)

# 3) One-page visual summary (arm alignment heatmaps, shared-primer bands,
#    flanking gene context, construct overview).
visualize_shared_design(
  result            = res,
  genbank_dir       = "path/to/gbk",
  construct_gbk_dir = "out/constructs",
  output_file       = "out/summary.pdf"
)
```

**Key behaviours:**

- **Strict overlap policy** — Under `overlap_policy = "strict"` the inner primers (UR, DF) are hard-pinned to the deletion boundary (0 bp tolerance). If no shared candidate exists at the boundary, the selector subgroup-splits into per-strain clusters rather than silently drifting off the stop codon.
- **Shared core + subgroup split** — Outer primers (UF, DR) and check primers progressively share across as many genomes as possible; strain-specific clusters are emitted only when biologically forced (e.g., transposon-adjacent targets).
- **Cross-role RC collision guard** — The selector rejects any primer whose reverse complement equals its partner on the same arm (prevents the Gibson pair from collapsing the effective homology arm).
- **Colony-PCR check primers** — cF / cR pairs are placed just outside the effective homology arm (`check_outer_pad` + `check_search_window`, default 50 bp + 800 bp) with a **progressive window-widening loop**: the window doubles up to `check_search_window_max` (default 6000 bp) until a single shared pair covers every genome, falling back to per-strain cR / cF clusters only when the surrounding context genuinely diverges. Tm is enforced at 55 ± 3 °C with a pair ΔTm ≤ 2 °C, strict genome-uniqueness, and a 6 bp 3'-complementarity heterodimer guard.
- **Per-genome WT / deletion PCR band sizes** — For every (genome, cF, cR) triple the expected wild-type amplicon and post-deletion amplicon length are reported alongside `delta_bp` (equal to the deletion span), so colony PCR gel interpretation is one table lookup.
- **3-sheet Excel output**:
  1. `primer_order` — single combined list of every oligo to order (arm primers with Gibson overhangs, gRNA cloning primers with Type IIS overhangs, cF / cR check primers). A `group` column separates arm / gRNA / check, rows are shaded by `used_by` cluster, and thin borders automatically outline each role block.
  2. `check_primers` — per-genome diagnostic PCR band sizes (WT, deletion, Δ) plus pair ΔTm and dimer flags.
  3. `final_construct_groups` — canonical construct list with gRNA protospacer and the full UF/UR/DF/DR primer cluster assignment for each build.
- **Single-page summary PDF** — `visualize_shared_design()` renders a 7-panel overview: construct maps, insert-site zoom, sharing matrix, oligo order list, upstream / downstream arm alignment heatmaps with cluster-banded primer lanes, and the flanking gene context (which CDSes occupy the homology arms on each genome).

The legacy lower-level `design_shared_grna_and_deletion()` is kept exported for back-compat; `shared_primer_design()` is the recommended entry point.

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
| `design_grna_and_deletion()` | **Full pipeline (single genome)**: accession → BSgenome → gRNA → cloning primers → deletion arms → Excel → GenBank. All steps automated |
| `shared_primer_design()` | **Full pipeline (multi-genome)**: shared UF/UR/DF/DR across strains + gRNA cloning + colony-PCR check primers + combined construct GenBank + 3-sheet Excel. See [Multi-Genome Shared Primer Design](#multi-genome-shared-primer-design) |
| `find_target_across_genomes()` | Resolve the same ortholog (by locus_tag / gene / product keyword) across a folder of genome files, with interactive ambiguity handling |
| `visualize_shared_design()` | Single-page 7-panel PDF summary for a `shared_primer_design()` result |

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
  output_dir          = "constructs/"
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
  output_dir          = "constructs/"
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

## Troubleshooting

**UCSC seqlengths warning on first run**

When building a BSgenome for non-model organisms, you may see warnings like:

```
Error in evaluating the argument 'x' in selecting a method for function 'seqlengths':
  UCSC library operation failed
BSgenome getSeq validation: FAILED -- will rebuild
```

This is **expected and harmless**. Non-UCSC genomes (i.e., most non-model organisms) do not have UCSC-style chromosome naming, so the initial seqlengths validation fails. PrimerDesigner automatically detects this and rebuilds the BSgenome package, which then works correctly. You can safely ignore this warning.

---

## Citation

If you use PrimerDesigner in your research, please cite:

```
Sung, J.-Y. (2025). PrimerDesigner: End-to-end CRISPR primer design pipeline
for non-model organisms. R package version 1.0.0.
https://github.com/JAEYOONSUNG/PrimerDesigner
```

---

## License

MIT © [Jae-Yoon Sung](https://github.com/JAEYOONSUNG)

# PrimerDesigner <img src="man/figures/logo.png" align="right" height="139" />

[![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-CrisprVerse-green)](https://bioconductor.org/packages/crisprVerse/)

End-to-end CRISPR primer design pipeline for non-model organisms.

GenBank accession과 벡터 `.dna` 파일만 넣으면 BSgenome 구축, Bowtie 인덱스 생성, gRNA 라이브러리, off-target 분석, Golden Gate / Gibson 클로닝 프라이머 설계, deletion arm 프라이머, SnapGene용 GenBank 출력까지 **한 함수 호출로** 전부 처리합니다.

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
devtools::install_github("your-username/PrimerDesigner")
```

**Python dependency** (for SnapGene `.dna` files):

```bash
pip install snapgene_reader
```

---

## Quick Start

GenBank accession + 벡터 `.dna` 파일만 주면 됩니다. 게놈 다운로드, BSgenome 패키지 빌드, Bowtie 인덱스 생성은 모두 자동으로 처리됩니다.

```r
library(PrimerDesigner)

# 이것만 호출하면 끝 — BSgenome/Bowtie 자동 빌드 포함
result <- design_grna_and_deletion(
  genbank_accession    = "GCF_030376765.1",
  locus_tags           = c("QT235_RS00005", "QT235_RS00010"),
  nuclease             = "GeoCas9",
  grna_vector_file     = "~/vectors/pG1Kt-GeoCas9EF.dna",
  grna_start           = 8811,
  grna_end             = 8840,
  grna_cloning_method  = "golden_gate",
  grna_enzyme          = "BbsI",
  methylation_patterns = c("GCCAT", "CCANNNNNTTG"),
  output_file          = "~/crispr_project/all_primers.xlsx",
  output_dir           = "~/crispr_project/constructs",
  dir                  = "~/crispr_project"
)

# 결과 확인
result$grna_primers      # gRNA cloning primers (2 per gRNA)
result$deletion_primers  # deletion arm primers (4 per locus_tag)
result$unified_primers   # 모든 프라이머 통합 long-format table
```

**첫 실행 시 자동으로 수행되는 작업:**
1. NCBI에서 `.fna` / `.gbff` 다운로드 (없으면)
2. BSgenome 패키지 빌드 + 로드
3. Bowtie 인덱스 빌드
4. gRNA 라이브러리 생성 + off-target 분석
5. 최적 gRNA 선별 + 클로닝 프라이머 설계
6. Deletion arm 프라이머 설계
7. Excel + GenBank 출력

두 번째 실행부터는 이미 빌드된 BSgenome/Bowtie를 재사용합니다.

---

## Pipeline Overview

```
                    design_grna_and_deletion()
                              │
         ┌────────────────────┼────────────────────┐
         │         자동 실행 (내부 호출)              │
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

**내부 호출 체인:**

```
design_grna_and_deletion()
  └─ design_grna_construct()
       └─ run_gRNA_list_generator()
            ├─ download_genbank_fna()       # .fna 없으면 자동 다운로드
            ├─ download_genbank_gbff()      # .gbff 없으면 자동 다운로드
            ├─ run_bowtie_build()           # Bowtie 인덱스 없으면 자동 빌드
            ├─ build_bsgenome_from_accession()  # BSgenome 없으면 자동 빌드
            └─ findSpacers + off-target + scoring
  └─ design_deletion_primers()              # 4-primer Gibson deletion arms
  └─ write_combined_construct_genbank()     # SnapGene용 GenBank 출력
```

---

## Golden Gate Oligo Annealing — Vector-Based Overhang Detection

PrimerDesigner는 하드코딩된 overhang 대신, 실제 벡터 서열에서 제한효소 인식서열을 찾아 올리고 어닐링 프라이머를 설계합니다.

### How it works

1. 벡터 파일 (`.dna` / `.gb`) 로드 → stuffer 양쪽 ±20bp 영역 스캔
2. **양쪽 strand** 모두 검색 → stuffer 방향으로 자르는 올바른 site만 사용 (잘못된 방향은 warning)
3. 컷패턴 `(a/b)` 에 따라 overhang 타입 자동 판별:
   - `a < b` (예: BbsI 2/6) → **5' overhang**, 길이 = b - a
   - `a > b` (가상 6/2) → **3' overhang**, 길이 = a - b, 올리고 구조 자동 조정
   - `a == b` → blunt end → Golden Gate 사용 불가 (warning)
4. Enzyme cut 과 stuffer 사이 **fill 서열** 자동 포함

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

### Main Pipeline (이것만 쓰면 됩니다)

| Function | Description |
|----------|-------------|
| `design_grna_and_deletion()` | **Full pipeline**: accession → BSgenome → gRNA → cloning primers → deletion arms → Excel → GenBank. 모든 단계 자동 실행 |

### Step-by-Step (필요 시 개별 호출)

| Function | Description |
|----------|-------------|
| `design_grna_construct()` | gRNA 선별 + cloning primer만 (deletion 제외) |
| `run_gRNA_list_generator()` | gRNA 라이브러리만 (BSgenome/Bowtie 자동 빌드 포함) |
| `design_cloning_primers()` | Golden Gate (OA) / Gibson (IA) primer만 |
| `design_deletion_primers()` | 4-primer deletion arm (single gene) |
| `batch_deletion_primers()` | Multi-gene deletion primer batch |

### gRNA Library Utilities

| Function | Description |
|----------|-------------|
| `generate_gRNA_for_sequence()` | DNA 서열 문자열에서 gRNA 생성 |
| `generate_gRNA_for_locus_tags()` | Locus tag 기반 gRNA 생성 |
| `build_gRNA_library()` | gRNA 라이브러리 구축 (scoring + merge) |
| `calculate_composite_score()` | 가중 점수 계산: GC%, Tm, position, off-target |
| `generate_gRNA_names()` | 표준화된 gRNA 이름 생성 (e.g. `Geo_dnaA_g1`) |
| `filter_methylation_sites()` | IUPAC methylation motif 필터링 |
| `nuclease_with_parameter()` | Nuclease 객체 + PAM/spacer 파라미터 반환 |

### Genome Setup (자동 호출됨 — 수동 호출 불필요)

| Function | Description |
|----------|-------------|
| `build_bsgenome_from_accession()` | BSgenome 빌드 (내부에서 자동 호출, UCSC validation 우회) |
| `run_bowtie_build()` | Bowtie 인덱스 빌드 (내부에서 자동 호출) |
| `download_genbank_fna()` | .fna 다운로드 (내부에서 자동 호출) |
| `download_genbank_gbff()` | .gbff 다운로드 (내부에서 자동 호출) |
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

두 올리고를 어닐링하여 Type IIS 제한효소가 만든 sticky end에 직접 ligation. 벡터 파일 제공 시 실제 컷 패턴 기반으로 overhang + fill 서열 자동 계산.

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

Tm 최적화된 overlap 프라이머로 벡터를 PCR 증폭하여 spacer 삽입.

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

유전자 상하류 homology arm을 Gibson Assembly로 벡터에 클로닝. Circular genome 지원. `design_grna_and_deletion()`에서 자동으로 포함됩니다.

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

GenBank 파일은 SnapGene에서 열어 spacer, homology arm, primer binding site를 시각적으로 확인할 수 있습니다.

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

# PrimerDesginer package
## Usage 
-----
### Requirements
- Python Packages:
  - `pandas`
  - `biopython`
  - `snapgene_reader`
- SnapGene: Requires SnapGene CLI for converting `.dna` files to GenBank format.
- Operating System: macOS (SnapGene path set to `/Applications/SnapGene.app/Contents/MacOS/SnapGene`).

### Installation
Install the required packages:
```bash
pip install pandas biopython snapgene_reader
```

# 1.PrimerDesigner for Gibson Assembly
----------------------------------

This Python script is designed to create Gibson Assembly primers and modify GenBank files by integrating target gene information extracted from an `xlsx` file generated by DNMB (DNA Manipulation and Molecular Biology). It leverages SnapGene software to convert `.dna` files into GenBank format and supports both single genes and multiple genes, either concatenated (e.g., `QT234_RS00005-00010`) or processed separately (e.g., `QT234_RS00005,QT234_RS00010`).

## Key Features
------------

1. Extract Target Gene Data from `xlsx`:
   - Retrieves nucleotide sequences (`rearranged_nt_seq`), direction (`direction`), product descriptions (`product`), translated amino acid sequences (`translation`), and protein IDs (`protein_id`) from an `xlsx` file produced by DNMB.
   - Supports single genes (e.g., `QT234_RS00005`), concatenated genes with `-` (e.g., `QT234_RS00005-00010`), or multiple separate genes with `,` (e.g., `QT234_RS00005,QT234_RS00010`).

2. Gene Sequence Concatenation:
   - Concatenates sequences of multiple genes when provided as a combined `locus_tag` with `-` (e.g., `QT234_RS00005-00010`), adjusting for directionality.
   - Combines `/translation` fields into a single amino acid sequence with quotes for concatenated genes.

3. Gibson Assembly Primer Design:
   - Designs Gibson Assembly primers based on the specified vector and target sequences.
   - Calculates melting temperatures (Tm) to ensure primer suitability.

4. GenBank File Modification:
   - Converts SnapGene `.dna` files to GenBank format.
   - Inserts target genes into the vector sequence and adds primer information as `primer` features.
   - Preserves existing GenBank features without quotes, applying quotes only to the target gene’s `/translation`.

5. Log Generation:
   - Saves designed primers and Tm values in an `xlsx` log file, with a single file for all processed `locus_tag`s.

Command-Line Execution Examples

1. Single Gene
For designing primers and generating a GenBank file for a single gene:
```bash
python PrimerDesigner_for_Gibson.py --genbank_table path/to/genbank_table.xlsx --locus_tag QT234_RS00005 --vector_file path/to/vector.dna --start 1000 --end 2000 --output_dir cloning_results
```
- Output: `QT234_RS00005_vector.gbk` and `log_QT234_RS00005.xlsx`.

2. Multiple Separate Genes (Comma-Separated)
For processing multiple genes independently in a loop:
```bash
python script.py --genbank_table path/to/genbank_table.xlsx --locus_tag QT234_RS00005,QT234_RS00010 --vector_file path/to/vector.dna --start 1000 --end 2000 --output_dir cloning_results
```
- Output: `QT234_RS00005_vector.gbk`, `QT234_RS00010_vector.gbk`, and `log_QT234_RS00005_QT234_RS00010.xlsx`.

3. Concatenated Operon (Hyphen-Connected)
For concatenating genes into a single operon-like target:
```bash
python script.py --genbank_table path/to/genbank_table.xlsx --locus_tag QT234_RS00005-00010 --vector_file path/to/vector.dna --start 1000 --end 2000 --output_dir cloning_results
```

### Arguments
---------
--genbank_table: Path to the `xlsx` or `tsv` file from DNMB (required).
--locus_tag: Target gene identifier (single: `QT234_RS00005`, concatenated: `QT234_RS00005-00010`) (required).
--vector_file: Path to the SnapGene `.dna` vector file (required).
--start: Insertion start position in the vector (1-based, required).
--end: Insertion end position in the vector (1-based, required).
--output_dir: Output directory (default: `cloning_results`).
--tm: Target Tm for primers (default: 60°C).

### Input Requirements
------------------
#### `xlsx` File Format
The `xlsx` file from DNMB must include the following columns:

Column Name        | Description                       | Example
-------------------|-----------------------------------|-----------------------
`locus_tag`        | Gene identifier                  | `QT234_RS00005`
`rearranged_nt_seq`| Nucleotide sequence              | `GTGGAA...`
`direction`        | Gene direction (`+` or `-`)      | `+`
`product`          | Gene product description         | `chromosomal replication initiator protein DnaA`
`translation`      | Translated amino acid sequence   | `MENIHDLWDRVL...`
`protein_id`       | Protein ID                       | `WP_289651751.1`

- For concatenated genes (e.g., `QT234_RS00005-00010`), each `locus_tag` must exist as a separate row in the `xlsx`.

#### Vector File
- A SnapGene `.dna` file is required.
- Insertion positions (`start`, `end`) must be valid within the vector sequence.


# 2.PrimerDesigner for Deletion vector
----------------------------------
## Overview
This Python script is designed to automate the design of Gibson Assembly primers for constructing deletion vectors in molecular biology. It supports both single-gene deletions and multi-gene region deletions by specifying a range of locus tags. The tool reads genomic data from a GenBank file and a corresponding table (Excel or TSV), designs primers with optimized melting temperatures (Tm) and overlap regions, and generates an updated GenBank file with the deletion construct.

## Features
- **Single Gene Deletion**: Delete a single gene by specifying its `locus_tag`.
- **Multi-Gene Region Deletion**: Delete a range of genes by specifying the start and end `locus_tag` in the format `START-END`.
- **Customizable Deletion Arms**: Specify the number of base pairs (bp) upstream and downstream of the target region.
- **Primer Design**: 
  - Target sequence Tm: 55–65°C (optimal 60°C).
  - Target sequence length: 18–50 bp.
  - Total primer length (including overlap): ≤70 bp.
  - Overlap region: 20–24 bp (even numbers only) with Tm ≥55°C.
- **Output**: 
  - Updated GenBank file with the deletion construct.
  - Excel log file with primer sequences and Tm values.

```bash
python PrimerDesigner_for_Deletion.py --genbank_file <GENBANK_FILE> --genbank_table <TABLE_FILE> --locus_tag <LOCUS_TAG> --upstream_bp <UPSTREAM_BP> --downstream_bp <DOWNSTREAM_BP> --vector_file <VECTOR_FILE> --start <START_POS> --end <END_POS> [--output_dir <OUTPUT_DIR>] [--tm <TM>]
```
### Arguments
---------
--genbank_file: Path to the GenBank file containing the target genome sequence (required).
--genbank_table: Path to the table file (Excel .xlsx or TSV) with locus tag positions and metadata (required).
--locus_tag: Target locus tag to delete. Use START-END for a range (e.g., QT235_RS00080-QT235_RS00090) (required).
--upstream_bp: Number of base pairs upstream of the target (required).
--downstream_bp: Number of base pairs downstream of the target (required).
--vector_file: Path to the SnapGene .dna vector file (required).
--start: Start position in the vector for insertion (1-based, required).
--end: End position in the vector for insertion (1-based, required).
--output_dir: Directory to save output files (default: deletion_results).
--tm: Target Tm for primers (default: 60°C, minimum 55°C).

### examples
Single Gene Deletion
```bash
python script.py --genbank_file genome.gb --genbank_table loci.xlsx --locus_tag QT235_RS00080 --upstream_bp 500 --downstream_bp 500 --vector_file vector.dna --start 3977 --end 4976 --output_dir deletion_output
```
Multi-Gene Region Deletion
```bash
python script.py --genbank_file genome.gb --genbank_table loci.xlsx --locus_tag QT235_RS00080-QT235_RS00090 --upstream_bp 500 --downstream_bp 500 --vector_file vector.dna --start 3977 --end 4976 --output_dir deletion_output
```

Output
Updated GenBank File: <output_dir>/<locus_tag>_deletion_vector.gbk
Log File: <output_dir>/log_<locus_tag>_deletion.xlsx
The log file contains primer sequences, target Tm, and full Tm values for each designed primer.

Key Functions
get_deletion_arms: Extracts upstream and downstream sequences based on the specified locus_tag(s) and base pair lengths.
design_deletion_primers: Designs four primers (Upstream_Forward, Upstream_Reverse, Downstream_Forward, Downstream_Reverse) with optimized Tm and overlap.
replace_sequence: Inserts the deletion arms into the vector sequence at the specified position.
modify_genbank: Updates the GenBank file with the new sequence and annotations.

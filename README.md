# NKG2A Ectodomain Comparison: Step 1 Starter

## Goal

Build a traceable sequence panel for NKG2A and NKG2C in:

- **Homo sapiens** (human)
- **Macaca mulatta** (rhesus macaque)
- **Macaca nemestrina** (pig-tailed macaque)

The later design objective is an NKG2A-selective VHH/nanobody that avoids NKG2C. Step 1 therefore prioritizes accession/locus correctness, retention of all credible NKG2C paralogs, explicit ectodomain boundaries, and source provenance before any RFantibody run.

## Data separation rule

Do not silently merge records from different sources.

- `02_raw_data/sequences/seed_sequences.fasta` contains the four sequences from JC's RTF plus two explicitly labelled bootstrap comparison records.
- `02_raw_data/sequences/reference_snapshots/reviewed_uniprot_2026-08-12.full_length.fasta` contains a separately preserved reviewed UniProt snapshot for human NKG2A/NKG2C and rhesus NKG2A/NKG2C.
- `03_metadata/source_comparison_pairs.tsv` declares which records may be compared. The comparison script reports substitutions and gaps but never replaces a record.

The rhesus JC records and the reviewed UniProt records are not identical. Their accession history, allele/transcript choice, and exact locus/paralog mapping must be resolved before the final panel is frozen.

## Quick start 

From the project root:

```bash
make step1-offline
```

This run:

1. checks the local environment;
2. rebuilds the bundled preliminary full-length panel;
3. validates amino-acid content, length, membrane-protein plausibility, cysteines, and N-X-S/T motifs;
4. extracts only ectodomains with explicit metadata boundaries;
5. computes preliminary pairwise identity tables;
6. rebuilds the reviewed-reference ectodomain panel;
7. compares JC and reviewed sequence sources without merging them;
8. runs unit tests; and
9. writes a SHA-256 project inventory.

Key outputs:

- `09_results/qc/sequence_qc.tsv`
- `09_results/tables/ectodomain_manifest.tsv`
- `09_results/tables/ectodomain_pairwise_identity_matrix.tsv`
- `09_results/tables/reviewed_reference_pairwise_identity_matrix.tsv`
- `09_results/tables/source_sequence_comparison.tsv`
- `09_results/tables/source_sequence_differences.tsv`
- `00_admin/PROJECT_FILE_MANIFEST.tsv`

The dependency-free pairwise alignment is for bootstrap QC only. The final analysis must use a manually reviewed multiple-sequence alignment.

## Optional curation environment

The bootstrap uses only Python's standard library. MAFFT and other convenience tools can be installed with:

```bash
conda env create -f environment-curation.yml
conda activate nkg2a-curation
```

Then run:

```bash
make alignment
```

The MAFFT wrapper stops if MAFFT is absent; it does not silently substitute a different aligner.

## Online sequence and candidate collection

Use an internet-connected workstation. NCBI requests require a contact email:

```bash
export NCBI_EMAIL='your.name@institution.edu'
make step1-online
```

An optional `NCBI_API_KEY` environment variable is supported.

The online workflow:

- refetches accessions already listed in `03_metadata/sequence_manifest.tsv`;
- records whether fetched sequences match bundled fallbacks;
- searches NCBI Protein using `03_metadata/ncbi_search_queries.tsv`;
- writes raw candidate-search JSON and summary TSV files;
- downloads the listed reference structures; and
- reruns QC and tests.

Search hits are **unreviewed candidates**, not accepted NKG2A/NKG2C assignments. Review gene/locus identity, synteny, transcript completeness, splice form, and topology before editing the master manifest.

## Current unresolved records

- Pig-tailed macaque NKG2A
- Every credible pig-tailed macaque NKG2C-like paralog/splice form
- Exact mapping of the complete rhesus NKG2C-1/NKG2C-2/NKG2C-3 set
- Species-matched rhesus and pig-tailed CD94 records
- RefSeq/UniProt reconciliation where source sequences differ

These gaps are intentionally visible in `03_metadata/sequence_manifest.tsv` and `00_admin/PROJECT_STATUS.md`.

## RFantibody later, on the GPU host

RFantibody is not required for Step 1. The guarded installer is:

```bash
bash 04_scripts/90_install_rfantibody_later.sh --help
```

Example on the intended Linux/NVIDIA workstation:

```bash
bash 04_scripts/90_install_rfantibody_later.sh \
  --ref 8fe3114 \
  --install-uv \
  --with-weights
```

The script verifies Linux, Git, `nvidia-smi`, and `uv`; records the exact Git commit; downloads model weights only when explicitly requested; runs `uv sync`; and checks the RFdiffusion command-line interface. It never uses `sudo`.

## Directory layout

```text
00_admin/                  status, decisions, provenance, checksums
01_references/uploaded/    original RTF, PowerPoint, and paper
01_references/external/    captured RFantibody setup files
02_raw_data/sequences/     source FASTA, fetched FASTA, ectodomains
02_raw_data/structures/    downloaded experimental structures
02_raw_data/genomes/       future KLRC locus/synteny extracts
03_metadata/               accessions, boundaries, search and review tables
04_scripts/                reproducible curation and setup scripts
05_alignments/             final MSA files
06_structures/             prepared CD94/NKG2 target/off-target models
07_analysis/               later sequence/surface analyses
08_rfantibody/             later configs, inputs, checkout, outputs
09_results/                QC tables, reports, figures
logs/                      execution logs
```

## Step 1 scientific completion criteria

Step 1 is complete only after:

1. canonical human NKG2A/NKG2C are reconciled across RefSeq and UniProt;
2. rhesus NKG2A source differences are resolved;
3. all expressed rhesus NKG2C paralogs and splice forms have verified accessions and locus assignments;
4. pig-tailed NKG2A and all plausible NKG2C-like records have accession plus synteny support;
5. species-matched CD94 records are curated;
6. every accepted ectodomain boundary is documented;
7. a final full-length and ectodomain FASTA are frozen with checksums; and
8. the MSA, identity matrix, and variable-position table are manually reviewed.

Monalizumab/Z199 residue annotations and RFantibody hotspot selection come after the target sequence panel and numbering scheme are accepted.

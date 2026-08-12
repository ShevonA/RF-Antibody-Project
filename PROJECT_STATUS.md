# Project status — Step 1 sequence curation

**Status date:** 2026-08-12

## Completed in this starter

- [x] Stable project directory tree created
- [x] Uploaded RTF, PowerPoint, and Nature paper preserved with SHA-256 checksums
- [x] Key RFantibody setup files captured for reference
- [x] JC's four protein sequences preserved without silent replacement
- [x] Separate reviewed UniProt snapshot preserved for human and rhesus NKG2A/NKG2C
- [x] Historical 246-aa rhesus NKG2-C2 record retained as a candidate, not treated as the complete negative panel
- [x] Offline sequence workflow implemented and executed
- [x] Online NCBI/UniProt retrieval and source-discrepancy logging implemented
- [x] NCBI candidate searches prepared for pig-tailed macaque, all three rhesus NKG2C loci, and CD94
- [x] Sequence QC, ectodomain extraction, pairwise identity, PDB download, MAFFT wrapper, source comparison, tests, checksums, and packaging implemented
- [x] RFantibody installation deferred behind explicit Linux/GPU/weight checks

## Current offline run

- Six bootstrap records were written to the preliminary working panel.
- Six explicit ectodomains were extracted from that panel.
- Four reviewed-reference records and four reviewed-reference ectodomains were preserved separately.
- Sequence QC: 6 PASS, 0 WARN, 0 FAIL.
- Eight unit tests passed.
- The reviewed human NKG2C sequence is identical to the JC human NKG2C seed in the captured records.
- The captured reviewed rhesus NKG2A and NKG2C records differ from JC's rhesus records and remain separate.
- Pairwise identities are bootstrap global-alignment QC values, not final MSA statistics.

## Required before Step 1 is scientifically complete

- [ ] Fetch and compare RefSeq NP_002250.2 with UniProt P26715
- [ ] Map NP_001098647.2 and Q9MZK6 to exact rhesus NKG2C locus/paralog assignments
- [ ] Curate rhesus NKG2C-1, NKG2C-2, and NKG2C-3 separately, including splice forms
- [ ] Identify pig-tailed macaque NKG2A with accession and locus/synteny support
- [ ] Identify every credible pig-tailed macaque NKG2C-like paralog and splice form
- [ ] Curate species-matched CD94 records
- [ ] Confirm ectodomain boundaries for every accepted record
- [ ] Generate and manually review the final multiple-sequence alignment
- [ ] Freeze the final panel with accessions, version dates, and checksums

## Runtime limitation

The build host has no visible NVIDIA GPU and its shell cannot resolve NCBI, UniProt, GitHub, or RCSB hosts. Live database searches, full Git cloning, structure downloads, and RFantibody model-weight downloads are therefore staged in scripts rather than claimed as complete. The reviewed FASTA snapshot is explicitly labelled by source and capture date and is not a substitute for final accession-history/locus curation.

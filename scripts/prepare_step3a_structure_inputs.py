#!/usr/bin/env python3
"""
STEP 3A.1
Prepare exact species-specific NKG2A sequences for structural prediction.

Purpose
-------
Prepare reproducible full-length NKG2A sequence inputs that will later be used
to generate species-specific structural models for RFantibody target
preparation.

This script DOES NOT:
    - predict structures
    - create RFantibody targets
    - alter Step 2V candidate rankings
    - introduce new structural evidence
    - use cross-species conservation/reactivity for ranking

Primary sequence authority
--------------------------
data/curated/primary_sequence_panel.tsv

Required positive targets
-------------------------
Rhesus macaque:
    accession: NP_001028001.3
    expected hotspot identities:
        113 A
        148 A

Pig-tailed macaque:
    accession: XP_070928357.1
    expected hotspot identities:
        113 A
        121 K
        148 T

Outputs
-------
structures/models/step3a/inputs/
    rhesus_NKG2A_NP_001028001_3.fasta
    pigtail_NKG2A_XP_070928357_1.fasta

results/tables/step3a/
    step3a_structure_input_manifest.tsv
    step3a_hotspot_sequence_validation.tsv

The script fails rather than writing final outputs if sequence provenance or
hotspot identity validation fails.
"""

from __future__ import annotations

import csv
import hashlib
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple


# ============================================================================
# PROJECT PATHS
# ============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

PRIMARY_PANEL = (
    PROJECT_ROOT
    / "data"
    / "curated"
    / "primary_sequence_panel.tsv"
)

RAW_DIR = PROJECT_ROOT / "data" / "raw"

OUTPUT_FASTA_DIR = (
    PROJECT_ROOT
    / "structures"
    / "models"
    / "step3a"
    / "inputs"
)

OUTPUT_TABLE_DIR = (
    PROJECT_ROOT
    / "results"
    / "tables"
    / "step3a"
)

MANIFEST_OUT = OUTPUT_TABLE_DIR / "step3a_structure_input_manifest.tsv"

HOTSPOT_OUT = OUTPUT_TABLE_DIR / "step3a_hotspot_sequence_validation.tsv"


# ============================================================================
# LOCKED STEP 3A TARGET DEFINITIONS
# ============================================================================

TARGETS = {
    "rhesus": {
        "species_label": "rhesus_macaque",
        "receptor": "NKG2A",
        "record_id": "rhesus_NKG2A",
        "accession": "NP_001028001.3",
        "raw_fasta": "rhesus_NKG2A_refseq_NP_001028001_3.fasta",
        "output_fasta": "rhesus_NKG2A_NP_001028001_3.fasta",
        "hotspots": {
            113: "A",
            148: "A",
        },
        "design_jobs": [
            "RH113",
            "RH148",
        ],
    },

    "pigtail": {
        "species_label": "pig_tailed_macaque",
        "receptor": "NKG2A",
        "record_id": "pigtail_NKG2A",
        "accession": "XP_070928357.1",
        "raw_fasta": "pigtail_NKG2A_X1_XP_070928357_1.fasta",
        "output_fasta": "pigtail_NKG2A_XP_070928357_1.fasta",
        "hotspots": {
            113: "A",
            121: "K",
            148: "T",
        },
        "design_jobs": [
            "PT113_121",
            "PT148",
        ],
    },
}


# ============================================================================
# BASIC UTILITIES
# ============================================================================

def fail(message: str) -> None:
    """Stop execution with a clearly labeled error."""
    print()
    print("=" * 78)
    print("STEP 3A VALIDATION FAILED")
    print("=" * 78)
    print(message)
    print()
    sys.exit(1)


def normalize_accession(value: str) -> str:
    """
    Normalize accession text for exact comparison.

    We intentionally preserve accession version numbers.
    """
    return (value or "").strip()


def sequence_sha256(sequence: str) -> str:
    """Return SHA-256 digest of an amino-acid sequence."""
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def validate_protein_sequence(sequence: str, label: str) -> None:
    """
    Validate that the sequence contains ordinary protein letters only.

    X is allowed because predicted/reference proteins can theoretically contain
    unknown residues. Stop characters and alignment gaps are not allowed in
    the structure-prediction FASTA.
    """
    allowed = set("ACDEFGHIKLMNPQRSTVWYXBZUOJ")

    bad = sorted(set(sequence.upper()) - allowed)

    if bad:
        fail(
            f"{label}: invalid characters found in protein sequence: "
            + ", ".join(repr(x) for x in bad)
        )

    if not sequence:
        fail(f"{label}: sequence is empty.")


# ============================================================================
# TSV READING
# ============================================================================

def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        fail(f"Required TSV does not exist:\n{path}")

    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            if reader.fieldnames is None:
                fail(f"TSV has no header:\n{path}")

            return [dict(row) for row in reader]

    except OSError as exc:
        fail(
            "Could not read required TSV.\n"
            f"Path: {path}\n"
            f"Error: {exc}\n\n"
            "If this is a Box Drive project, confirm Box Drive is running "
            "and that the file is available locally."
        )


def find_primary_panel_record(
    rows: List[Dict[str, str]],
    expected_record_id: str,
    expected_accession: str,
) -> Dict[str, str]:

    matches = [
        row
        for row in rows
        if (row.get("record_id") or "").strip() == expected_record_id
    ]

    if len(matches) != 1:
        fail(
            f"Expected exactly one primary-panel row for "
            f"{expected_record_id!r}; found {len(matches)}."
        )

    row = matches[0]

    accession = normalize_accession(row.get("accession", ""))

    if accession != expected_accession:
        fail(
            f"Primary-panel accession mismatch for {expected_record_id}.\n"
            f"Expected: {expected_accession}\n"
            f"Observed: {accession}"
        )

    status = (row.get("status") or "").strip().lower()

    if status != "primary":
        fail(
            f"{expected_record_id} is not marked primary in "
            f"primary_sequence_panel.tsv.\n"
            f"Observed status: {row.get('status', '')}"
        )

    receptor = (row.get("receptor") or "").strip().upper()

    if receptor != "NKG2A":
        fail(
            f"{expected_record_id} does not have receptor=NKG2A.\n"
            f"Observed receptor: {row.get('receptor', '')}"
        )

    return row


# ============================================================================
# FASTA READING
# ============================================================================

def read_single_fasta(path: Path) -> Tuple[str, str]:
    """
    Read exactly one FASTA record.

    Returns
    -------
    header
    sequence
    """

    if not path.exists():
        fail(f"Required FASTA does not exist:\n{path}")

    try:
        text = path.read_text(encoding="utf-8-sig")
    except OSError as exc:
        fail(
            "Could not read required FASTA.\n"
            f"Path: {path}\n"
            f"Error: {exc}\n\n"
            "If this is a Box Drive project, confirm Box Drive is running "
            "and that the file is available locally."
        )

    lines = [line.strip() for line in text.splitlines() if line.strip()]

    if not lines:
        fail(f"FASTA is empty:\n{path}")

    headers = [line for line in lines if line.startswith(">")]

    if len(headers) != 1:
        fail(
            f"Expected exactly one FASTA record in:\n{path}\n"
            f"Observed FASTA headers: {len(headers)}"
        )

    if not lines[0].startswith(">"):
        fail(f"First non-empty FASTA line is not a header:\n{path}")

    header = lines[0][1:].strip()

    sequence_lines = []

    for line in lines[1:]:
        if line.startswith(">"):
            fail(f"Multiple FASTA records detected in:\n{path}")

        sequence_lines.append(re.sub(r"\s+", "", line))

    sequence = "".join(sequence_lines).upper()

    validate_protein_sequence(sequence, str(path))

    return header, sequence


def accession_from_fasta_header(header: str) -> str:
    """
    Extract the first whitespace-delimited token from an NCBI-style FASTA
    header.
    """
    return header.split()[0].strip()


# ============================================================================
# HOTSPOT VALIDATION
# ============================================================================

def validate_hotspots(
    species_key: str,
    sequence: str,
    hotspots: Dict[int, str],
) -> List[Dict[str, str]]:

    rows: List[Dict[str, str]] = []

    for position in sorted(hotspots):
        expected = hotspots[position]

        if position < 1 or position > len(sequence):
            observed = ""
            passed = False
            reason = "position_outside_sequence"
        else:
            observed = sequence[position - 1]
            passed = observed == expected

            if passed:
                reason = "expected_residue_confirmed"
            else:
                reason = "residue_identity_mismatch"

        rows.append(
            {
                "species": species_key,
                "position_1based_full_length": str(position),
                "expected_nkg2a_residue": expected,
                "observed_nkg2a_residue": observed,
                "validation_pass": "yes" if passed else "no",
                "validation_reason": reason,
            }
        )

    return rows


# ============================================================================
# OUTPUT WRITING
# ============================================================================

def write_fasta(path: Path, header: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{header}\n")

        width = 70
        for start in range(0, len(sequence), width):
            handle.write(sequence[start:start + width] + "\n")


def write_tsv(
    path: Path,
    rows: List[Dict[str, str]],
    fieldnames: List[str],
) -> None:

    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(row)


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print("STEP 3A.1 - PREPARE NKG2A STRUCTURE-PREDICTION INPUTS")
    print("=" * 78)
    print()
    print("Primary sequence authority:")
    print(f"  {PRIMARY_PANEL}")
    print()
    print("This step:")
    print("  - locks exact macaque NKG2A sequences")
    print("  - verifies primary accessions")
    print("  - verifies intended hotspot residue identities")
    print("  - creates clean full-length structure-prediction FASTAs")
    print("  - records SHA-256 sequence hashes")
    print()
    print("This step does NOT:")
    print("  - predict structures")
    print("  - create RFantibody targets")
    print("  - change Step 2V rankings")
    print("  - add new structural evidence")
    print()

    panel_rows = read_tsv(PRIMARY_PANEL)

    manifest_rows: List[Dict[str, str]] = []
    hotspot_rows: List[Dict[str, str]] = []

    prepared_sequences: Dict[str, Tuple[str, str]] = {}

    # ----------------------------------------------------------------------
    # Validate each locked target
    # ----------------------------------------------------------------------

    for species_key, cfg in TARGETS.items():

        print("-" * 78)
        print(species_key.upper())
        print("-" * 78)

        panel_row = find_primary_panel_record(
            panel_rows,
            expected_record_id=cfg["record_id"],
            expected_accession=cfg["accession"],
        )

        fasta_path = RAW_DIR / cfg["raw_fasta"]

        header, sequence = read_single_fasta(fasta_path)

        fasta_accession = accession_from_fasta_header(header)

        if fasta_accession != cfg["accession"]:
            fail(
                f"{species_key}: FASTA accession does not match locked "
                f"primary accession.\n"
                f"Expected: {cfg['accession']}\n"
                f"Observed FASTA header accession: {fasta_accession}\n"
                f"FASTA: {fasta_path}"
            )

        hotspot_validation = validate_hotspots(
            species_key,
            sequence,
            cfg["hotspots"],
        )

        failed_hotspots = [
            row
            for row in hotspot_validation
            if row["validation_pass"] != "yes"
        ]

        if failed_hotspots:
            details = "\n".join(
                (
                    f"  position "
                    f"{row['position_1based_full_length']}: "
                    f"expected {row['expected_nkg2a_residue']}, "
                    f"observed {row['observed_nkg2a_residue'] or '<missing>'}"
                )
                for row in failed_hotspots
            )

            fail(
                f"{species_key}: hotspot sequence validation failed.\n"
                f"{details}\n\n"
                "No final Step 3A structure-prediction inputs were written."
            )

        hotspot_rows.extend(hotspot_validation)

        digest = sequence_sha256(sequence)

        output_fasta = OUTPUT_FASTA_DIR / cfg["output_fasta"]

        output_header = (
            f"{cfg['accession']} "
            f"{cfg['species_label']} "
            f"NKG2A "
            f"Step3A_exact_primary_sequence"
        )

        prepared_sequences[species_key] = (
            output_header,
            sequence,
        )

        hotspot_text = ",".join(
            f"{position}{residue}"
            for position, residue in sorted(cfg["hotspots"].items())
        )

        manifest_rows.append(
            {
                "species": species_key,
                "species_label": cfg["species_label"],
                "receptor": cfg["receptor"],
                "primary_record_id": cfg["record_id"],
                "primary_accession": cfg["accession"],
                "primary_panel_status": (
                    panel_row.get("status") or ""
                ).strip(),
                "primary_panel_role": (
                    panel_row.get("role") or ""
                ).strip(),
                "source_fasta": str(
                    fasta_path.relative_to(PROJECT_ROOT)
                ).replace("\\", "/"),
                "source_fasta_header": header,
                "sequence_length_aa": str(len(sequence)),
                "sequence_sha256": digest,
                "intended_hotspots_full_length": hotspot_text,
                "hotspot_validation": "pass",
                "design_jobs": ";".join(cfg["design_jobs"]),
                "prediction_input_fasta": str(
                    output_fasta.relative_to(PROJECT_ROOT)
                ).replace("\\", "/"),
                "sequence_scope": "full_length",
                "numbering_scheme": "full_length_1based",
                "step3a_status": "validated_for_structure_prediction",
            }
        )

        print(f"Primary accession: {cfg['accession']}")
        print(f"Source FASTA:      {fasta_path}")
        print(f"Length:            {len(sequence)} aa")
        print(f"SHA-256:           {digest}")
        print(f"Hotspots:          {hotspot_text}")
        print("Hotspot validation: PASS")

        for row in hotspot_validation:
            print(
                "  "
                f"{row['position_1based_full_length']}"
                f"{row['observed_nkg2a_residue']}"
                "  PASS"
            )

        print()

    # ----------------------------------------------------------------------
    # Only write final outputs after ALL targets pass validation
    # ----------------------------------------------------------------------

    OUTPUT_FASTA_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_TABLE_DIR.mkdir(parents=True, exist_ok=True)

    for species_key, cfg in TARGETS.items():

        header, sequence = prepared_sequences[species_key]

        output_fasta = OUTPUT_FASTA_DIR / cfg["output_fasta"]

        write_fasta(
            output_fasta,
            header,
            sequence,
        )

    write_tsv(
        MANIFEST_OUT,
        manifest_rows,
        [
            "species",
            "species_label",
            "receptor",
            "primary_record_id",
            "primary_accession",
            "primary_panel_status",
            "primary_panel_role",
            "source_fasta",
            "source_fasta_header",
            "sequence_length_aa",
            "sequence_sha256",
            "intended_hotspots_full_length",
            "hotspot_validation",
            "design_jobs",
            "prediction_input_fasta",
            "sequence_scope",
            "numbering_scheme",
            "step3a_status",
        ],
    )

    write_tsv(
        HOTSPOT_OUT,
        hotspot_rows,
        [
            "species",
            "position_1based_full_length",
            "expected_nkg2a_residue",
            "observed_nkg2a_residue",
            "validation_pass",
            "validation_reason",
        ],
    )

    # ----------------------------------------------------------------------
    # Final summary
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("STEP 3A.1 VALIDATION SUMMARY")
    print("=" * 78)
    print()

    for row in manifest_rows:
        print(
            f"{row['species']:<9} "
            f"{row['primary_accession']:<18} "
            f"length={row['sequence_length_aa']:<4} "
            f"hotspots={row['intended_hotspots_full_length']}"
        )

    print()
    print("All locked primary sequences passed validation.")
    print()

    print("=" * 78)
    print("STRUCTURE-PREDICTION FASTA INPUTS")
    print("=" * 78)
    print()

    for cfg in TARGETS.values():
        print(OUTPUT_FASTA_DIR / cfg["output_fasta"])

    print()
    print("=" * 78)
    print("QC OUTPUTS")
    print("=" * 78)
    print()
    print(MANIFEST_OUT)
    print(HOTSPOT_OUT)

    print()
    print("=" * 78)
    print("IMPORTANT")
    print("=" * 78)
    print()
    print(
        "These FASTAs contain the exact validated full-length primary "
        "NKG2A sequences."
    )
    print(
        "Residue numbering is full-length 1-based numbering so Step 2 "
        "candidate positions are preserved."
    )
    print(
        "No structure has yet been predicted and no RFantibody target "
        "PDB has yet been created."
    )
    print(
        "The next stage is species-specific structural prediction followed "
        "by sequence/numbering QC before RFantibody target preparation."
    )


if __name__ == "__main__":
    main()
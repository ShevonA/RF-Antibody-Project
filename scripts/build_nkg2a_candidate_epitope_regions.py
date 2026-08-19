#!/usr/bin/env python3

"""
STEP 2S - BUILD WITHIN-SPECIES NKG2A-SPECIFIC CANDIDATE EPITOPE REGIONS

Purpose
-------
Integrate the Step 2R within-species NKG2A-vs-NKG2C specificity evidence
into species-specific candidate epitope regions.

The analysis is deliberately restricted to Axis 1:

    human:
        human NKG2A vs human NKG2C

    rhesus macaque:
        rhesus NKG2A vs BOTH rhesus NKG2C isoforms

    pigtail macaque:
        pigtail NKG2A vs pigtail NKG2C

Cross-species NKG2A conservation and antibody cross-reactivity are NOT used
as ranking criteria in this step.

Important structural interpretation rules
-----------------------------------------
1. Experimental human NKG2A coordinates are used where available.

2. Human residues in the experimentally unresolved N-terminal region are
   retained as a SEQUENCE-DEFINED candidate region.

3. AlphaFold coordinates are NOT used to infer a fixed 3D epitope geometry
   for the unresolved N-terminal region. Step 2Q showed low confidence and
   high relative-position PAE for that region.

4. Rhesus and pigtail structural annotations are homologous-position
   projections from human NKG2A / 3CDG, not direct macaque structures.

5. Structurally resolved discriminatory residues are grouped into candidate
   regions using minimum heavy-atom distances from the human NKG2A structure.

6. A connected structural region means that residues are linked through
   pairwise discriminatory-residue distances <= STRUCTURAL_PATCH_CUTOFF_A.
   The residues do not all need to be directly within the cutoff of every
   other residue.

Inputs
------
results/tables/structure/
    human_NKG2A_specificity_candidates.tsv
    rhesus_NKG2A_specificity_candidates.tsv
    pigtail_NKG2A_specificity_candidates.tsv
    nkg2a_candidate_spatial_distances.tsv
    nkg2_ectodomain_n_glycosylation_sites.tsv
    alphafold_nkg2a_model_validation.tsv

Outputs
-------
results/tables/structure/
    human_NKG2A_candidate_epitope_regions.tsv
    rhesus_NKG2A_candidate_epitope_regions.tsv
    pigtail_NKG2A_candidate_epitope_regions.tsv
    nkg2a_within_species_candidate_epitope_regions.tsv
"""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple


# =============================================================================
# CONFIGURATION
# =============================================================================

ROOT = Path(__file__).resolve().parents[1]

STRUCTURE_DIR = ROOT / "results" / "tables" / "structure"

HUMAN_INPUT = STRUCTURE_DIR / "human_NKG2A_specificity_candidates.tsv"
RHESUS_INPUT = STRUCTURE_DIR / "rhesus_NKG2A_specificity_candidates.tsv"
PIGTAIL_INPUT = STRUCTURE_DIR / "pigtail_NKG2A_specificity_candidates.tsv"

DISTANCE_INPUT = STRUCTURE_DIR / "nkg2a_candidate_spatial_distances.tsv"
GLYCOSYLATION_INPUT = STRUCTURE_DIR / "nkg2_ectodomain_n_glycosylation_sites.tsv"
ALPHAFOLD_VALIDATION_INPUT = (
    STRUCTURE_DIR / "alphafold_nkg2a_model_validation.tsv"
)

HUMAN_OUTPUT = STRUCTURE_DIR / "human_NKG2A_candidate_epitope_regions.tsv"
RHESUS_OUTPUT = STRUCTURE_DIR / "rhesus_NKG2A_candidate_epitope_regions.tsv"
PIGTAIL_OUTPUT = STRUCTURE_DIR / "pigtail_NKG2A_candidate_epitope_regions.tsv"

COMBINED_OUTPUT = (
    STRUCTURE_DIR / "nkg2a_within_species_candidate_epitope_regions.tsv"
)


# Structural candidate residues are connected if a discriminatory residue
# pair has minimum heavy-atom distance <= this cutoff.
STRUCTURAL_PATCH_CUTOFF_A = 10.0

# Unresolved discriminatory residues separated by no more than this many
# full-length sequence positions are combined into one sequence-defined region.
#
# Example:
# 95,96,97,98,99,100,101,106
#
# A gap allowance of 5 keeps these as a single N-terminal candidate region,
# while still preventing unrelated distant sequence positions from being
# merged indiscriminately.
SEQUENCE_REGION_MAX_GAP = 5

# Canonical N-linked glycosylation sequons within this sequence distance of
# a candidate region are reported as nearby glycosylation context.
GLYCOSYLATION_CONTEXT_WINDOW = 5


# =============================================================================
# GENERAL HELPERS
# =============================================================================

def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Required input file not found:\n{path}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def write_tsv(
    path: Path,
    rows: Sequence[Dict[str, Any]],
    fieldnames: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in rows:
            out = {}
            for field in fieldnames:
                value = row.get(field, "")
                if value is None:
                    value = ""
                out[field] = value
            writer.writerow(out)


def clean(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def first_nonempty(row: Dict[str, str], names: Iterable[str]) -> str:
    for name in names:
        value = clean(row.get(name, ""))
        if value != "":
            return value
    return ""


def parse_int(value: Any) -> Optional[int]:
    text = clean(value)
    if text == "":
        return None

    try:
        return int(float(text))
    except ValueError:
        return None


def parse_float(value: Any) -> Optional[float]:
    text = clean(value)
    if text == "":
        return None

    try:
        value_float = float(text)
    except ValueError:
        return None

    if not math.isfinite(value_float):
        return None

    return value_float


def yes(value: Any) -> bool:
    return clean(value).lower() in {
        "yes",
        "y",
        "true",
        "1",
        "t",
    }


def mean(values: Iterable[Optional[float]]) -> Optional[float]:
    usable = [x for x in values if x is not None]

    if not usable:
        return None

    return sum(usable) / len(usable)


def fmt_float(
    value: Optional[float],
    digits: int = 4,
) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def print_banner(title: str) -> None:
    print("=" * 78)
    print(title)
    print("=" * 78)


def print_subbanner(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def natural_region_sort_key(region_id: str) -> Tuple[str, int]:
    match = re.match(r"([A-Za-z]+)(\d+)$", region_id)

    if not match:
        return region_id, 999999

    return match.group(1), int(match.group(2))


# =============================================================================
# COLUMN NORMALIZATION
# =============================================================================

SPECIES_CONFIG = {
    "human": {
        "input": HUMAN_INPUT,
        "output": HUMAN_OUTPUT,
        "prefix": "H",
        "nkg2a_residue_columns": [
            "human_NKG2A_residue",
            "full_length_residue",
        ],
        "nkg2a_aa_columns": [
            "human_NKG2A_aa",
            "nkg2a_aa",
        ],
        "nkg2c_aa_columns": [
            "human_NKG2C_aa",
        ],
        "tier_columns": [
            "human_evidence_tier",
            "evidence_tier",
        ],
    },
    "rhesus": {
        "input": RHESUS_INPUT,
        "output": RHESUS_OUTPUT,
        "prefix": "R",
        "nkg2a_residue_columns": [
            "rhesus_NKG2A_residue",
            "human_NKG2A_residue",
            "full_length_residue",
        ],
        "nkg2a_aa_columns": [
            "rhesus_NKG2A_aa",
            "nkg2a_aa",
        ],
        "nkg2c1_aa_columns": [
            "rhesus_NKG2C1_aa",
        ],
        "nkg2c2_aa_columns": [
            "rhesus_NKG2C2_aa",
        ],
        "tier_columns": [
            "rhesus_evidence_tier",
            "evidence_tier",
        ],
    },
    "pigtail": {
        "input": PIGTAIL_INPUT,
        "output": PIGTAIL_OUTPUT,
        "prefix": "P",
        "nkg2a_residue_columns": [
            "pigtail_NKG2A_residue",
            "human_NKG2A_residue",
            "full_length_residue",
        ],
        "nkg2a_aa_columns": [
            "pigtail_NKG2A_aa",
            "nkg2a_aa",
        ],
        "nkg2c_aa_columns": [
            "pigtail_NKG2C_aa",
        ],
        "tier_columns": [
            "pigtail_evidence_tier",
            "evidence_tier",
        ],
    },
}


def normalize_candidate_row(
    species: str,
    row: Dict[str, str],
) -> Optional[Dict[str, Any]]:
    config = SPECIES_CONFIG[species]

    residue = parse_int(
        first_nonempty(
            row,
            config["nkg2a_residue_columns"],
        )
    )

    if residue is None:
        return None

    nkg2a_aa = first_nonempty(
        row,
        config["nkg2a_aa_columns"],
    )

    tier = first_nonempty(
        row,
        config["tier_columns"],
    )

    if species == "rhesus":
        nkg2c1 = first_nonempty(
            row,
            config["nkg2c1_aa_columns"],
        )
        nkg2c2 = first_nonempty(
            row,
            config["nkg2c2_aa_columns"],
        )

        nkg2c_display = f"{nkg2c1}/{nkg2c2}"

        discriminatory = (
            nkg2a_aa != ""
            and nkg2c1 != ""
            and nkg2c2 != ""
            and nkg2a_aa != nkg2c1
            and nkg2a_aa != nkg2c2
        )

    else:
        nkg2c = first_nonempty(
            row,
            config["nkg2c_aa_columns"],
        )

        nkg2c_display = nkg2c

        discriminatory = (
            nkg2a_aa != ""
            and nkg2c != ""
            and nkg2a_aa != nkg2c
        )

    experimental_structure_class = first_nonempty(
        row,
        [
            "experimental_structure_class",
            "structure_class",
        ],
    )

    complex_rsa = parse_float(
        first_nonempty(
            row,
            [
                "complex_rsa",
                "human_3CDG_complex_rsa",
            ],
        )
    )

    interface = yes(
        first_nonempty(
            row,
            [
                "any_interface_contact",
                "interface_contact",
            ],
        )
    )

    compact_footprint_ids = first_nonempty(
        row,
        [
            "compact_footprint_ids",
            "footprint_ids",
        ],
    )

    alphafold_plddt = parse_float(
        first_nonempty(
            row,
            [
                "alphafold_plddt",
                "plddt",
            ],
        )
    )

    alphafold_geometry_evidence = first_nonempty(
        row,
        [
            "alphafold_geometry_evidence",
        ],
    )

    tier_lower = tier.lower()
    structure_lower = experimental_structure_class.lower()

    unresolved = (
        "unresolved" in tier_lower
        or "no_experiment" in structure_lower
        or "unresolved" in structure_lower
    )

    resolved = not unresolved

    exposed_noninterface = (
        resolved
        and complex_rsa is not None
        and complex_rsa >= 0.25
        and not interface
    )

    exposed_interface = (
        resolved
        and complex_rsa is not None
        and complex_rsa >= 0.25
        and interface
    )

    partially_exposed_noninterface = (
        resolved
        and complex_rsa is not None
        and 0.10 <= complex_rsa < 0.25
        and not interface
    )

    partially_exposed_interface = (
        resolved
        and complex_rsa is not None
        and 0.10 <= complex_rsa < 0.25
        and interface
    )

    buried = (
        resolved
        and complex_rsa is not None
        and complex_rsa < 0.10
    )

    return {
        "species": species,
        "residue": residue,
        "nkg2a_aa": nkg2a_aa,
        "nkg2c_display": nkg2c_display,
        "tier": tier,
        "experimental_structure_class": experimental_structure_class,
        "complex_rsa": complex_rsa,
        "interface": interface,
        "compact_footprint_ids": compact_footprint_ids,
        "alphafold_plddt": alphafold_plddt,
        "alphafold_geometry_evidence": alphafold_geometry_evidence,
        "unresolved": unresolved,
        "resolved": resolved,
        "exposed_noninterface": exposed_noninterface,
        "exposed_interface": exposed_interface,
        "partially_exposed_noninterface": partially_exposed_noninterface,
        "partially_exposed_interface": partially_exposed_interface,
        "buried": buried,
        "discriminatory": discriminatory,
        "raw": row,
    }


# =============================================================================
# DISTANCE GRAPH
# =============================================================================

def load_distance_map(
    path: Path,
) -> Dict[Tuple[int, int], float]:
    rows = read_tsv(path)

    distances: Dict[Tuple[int, int], float] = {}

    for row in rows:
        residue_1 = parse_int(
            first_nonempty(
                row,
                [
                    "residue_1",
                    "full_length_residue_1",
                ],
            )
        )

        residue_2 = parse_int(
            first_nonempty(
                row,
                [
                    "residue_2",
                    "full_length_residue_2",
                ],
            )
        )

        distance = parse_float(
            first_nonempty(
                row,
                [
                    "minimum_heavy_atom_distance_A",
                    "minimum_distance_A",
                ],
            )
        )

        if (
            residue_1 is None
            or residue_2 is None
            or distance is None
        ):
            continue

        key = tuple(sorted((residue_1, residue_2)))
        distances[key] = distance

    return distances


def connected_components(
    residues: Sequence[int],
    distance_map: Dict[Tuple[int, int], float],
    cutoff: float,
) -> List[List[int]]:
    residue_set = set(residues)

    adjacency: Dict[int, Set[int]] = {
        residue: set()
        for residue in residue_set
    }

    for residue_1 in residue_set:
        for residue_2 in residue_set:
            if residue_2 <= residue_1:
                continue

            key = tuple(sorted((residue_1, residue_2)))
            distance = distance_map.get(key)

            if distance is not None and distance <= cutoff:
                adjacency[residue_1].add(residue_2)
                adjacency[residue_2].add(residue_1)

    components: List[List[int]] = []
    visited: Set[int] = set()

    for start in sorted(residue_set):
        if start in visited:
            continue

        stack = [start]
        component: List[int] = []

        while stack:
            current = stack.pop()

            if current in visited:
                continue

            visited.add(current)
            component.append(current)

            for neighbor in sorted(
                adjacency[current],
                reverse=True,
            ):
                if neighbor not in visited:
                    stack.append(neighbor)

        components.append(sorted(component))

    return components


# =============================================================================
# SEQUENCE-DEFINED REGIONS
# =============================================================================

def sequence_components(
    residues: Sequence[int],
    max_gap: int,
) -> List[List[int]]:
    if not residues:
        return []

    ordered = sorted(set(residues))

    components: List[List[int]] = [[ordered[0]]]

    for residue in ordered[1:]:
        previous = components[-1][-1]

        if residue - previous <= max_gap:
            components[-1].append(residue)
        else:
            components.append([residue])

    return components


# =============================================================================
# GLYCOSYLATION CONTEXT
# =============================================================================

GLYCO_RECORD_IDS = {
    "human": "human_NKG2A",
    "rhesus": "rhesus_NKG2A",
    "pigtail": "pigtail_NKG2A",
}


def load_glycosylation_sites(
    path: Path,
) -> Dict[str, List[Tuple[int, str]]]:
    rows = read_tsv(path)

    by_species: Dict[str, List[Tuple[int, str]]] = defaultdict(list)

    record_to_species = {
        record_id: species
        for species, record_id in GLYCO_RECORD_IDS.items()
    }

    for row in rows:
        record_id = clean(row.get("record_id", ""))

        if record_id not in record_to_species:
            continue

        residue = parse_int(
            first_nonempty(
                row,
                [
                    "sequon_full_length_residue",
                    "full_length_residue",
                ],
            )
        )

        motif = first_nonempty(
            row,
            [
                "motif",
                "sequon",
            ],
        )

        if residue is None:
            continue

        species = record_to_species[record_id]
        by_species[species].append((residue, motif))

    for species in by_species:
        by_species[species].sort()

    return by_species


def glycosylation_context(
    species: str,
    residues: Sequence[int],
    glyco_sites: Dict[str, List[Tuple[int, str]]],
) -> Tuple[str, str]:
    if not residues:
        return "none", ""

    region_min = min(residues)
    region_max = max(residues)

    overlapping: List[str] = []
    nearby: List[str] = []

    for sequon_residue, motif in glyco_sites.get(species, []):
        label = f"{sequon_residue}{motif}"

        if region_min <= sequon_residue <= region_max:
            overlapping.append(label)

        elif (
            region_min - GLYCOSYLATION_CONTEXT_WINDOW
            <= sequon_residue
            <= region_max + GLYCOSYLATION_CONTEXT_WINDOW
        ):
            nearby.append(label)

    if overlapping:
        context = "overlapping_canonical_N_glycosylation_sequon"
        labels = ",".join(overlapping)

        if nearby:
            labels += ";nearby=" + ",".join(nearby)

        return context, labels

    if nearby:
        return (
            "nearby_canonical_N_glycosylation_sequon",
            ",".join(nearby),
        )

    return "none_detected_near_region", ""


# =============================================================================
# ALPHAFOLD INTERPRETATION
# =============================================================================

def load_alphafold_decision(path: Path) -> str:
    rows = read_tsv(path)

    if not rows:
        return ""

    return first_nonempty(
        rows[0],
        [
            "model_use_decision",
            "alphafold_model_use_decision",
        ],
    )


# =============================================================================
# REGION SCORING
# =============================================================================

def region_priority(
    region_type: str,
    members: Sequence[Dict[str, Any]],
) -> str:
    exposed_noninterface_count = sum(
        1
        for member in members
        if member["exposed_noninterface"]
    )

    exposed_interface_count = sum(
        1
        for member in members
        if member["exposed_interface"]
    )

    partially_exposed_noninterface_count = sum(
        1
        for member in members
        if member["partially_exposed_noninterface"]
    )

    partially_exposed_interface_count = sum(
        1
        for member in members
        if member["partially_exposed_interface"]
    )

    buried_count = sum(
        1
        for member in members
        if member["buried"]
    )

    unresolved_count = sum(
        1
        for member in members
        if member["unresolved"]
    )

    residue_count = len(members)

    if region_type == "sequence_defined_unresolved":
        if residue_count >= 3:
            return "priority_sequence_defined_unresolved_region"
        return "secondary_sequence_defined_unresolved_region"

    if exposed_noninterface_count >= 2:
        return "priority_resolved_exposed_noninterface_patch"

    if exposed_noninterface_count == 1 and residue_count == 1:
        return "priority_resolved_exposed_noninterface_site"

    if exposed_noninterface_count == 1:
        return "priority_resolved_mixed_surface_patch"

    if exposed_interface_count >= 1:
        return "secondary_resolved_interface_patch"

    if partially_exposed_noninterface_count >= 1:
        return "secondary_partially_exposed_noninterface_patch"

    if partially_exposed_interface_count >= 1:
        return "secondary_partially_exposed_interface_patch"

    if buried_count == residue_count and residue_count > 0:
        return "low_priority_buried_region"

    if unresolved_count == residue_count and residue_count > 0:
        return "sequence_defined_unresolved_region"

    return "lower_structural_priority"


def region_interpretation(
    species: str,
    region_type: str,
    members: Sequence[Dict[str, Any]],
) -> str:
    labels = ",".join(
        f"{member['residue']}{member['nkg2a_aa']}"
        for member in members
    )

    exposed_noninterface_count = sum(
        1
        for member in members
        if member["exposed_noninterface"]
    )

    interface_count = sum(
        1
        for member in members
        if member["interface"]
    )

    buried_count = sum(
        1
        for member in members
        if member["buried"]
    )

    if region_type == "sequence_defined_unresolved":
        if species == "human":
            return (
                f"{labels} is a contiguous/near-contiguous human NKG2A-vs-NKG2C "
                "sequence-discriminatory region lacking experimental NKG2A "
                "coordinates. AlphaFold may provide local residue confidence, "
                "but Step 2Q does not support treating its orientation as a "
                "fixed antibody epitope geometry."
            )

        return (
            f"{labels} is a contiguous/near-contiguous {species} "
            "NKG2A-vs-NKG2C sequence-discriminatory region corresponding to "
            "the experimentally unresolved N-terminal portion of human NKG2A. "
            "It is retained as sequence evidence rather than a resolved "
            "three-dimensional epitope patch."
        )

    if exposed_noninterface_count >= 2 and interface_count == 0:
        return (
            f"{labels} forms a resolved discriminatory surface neighborhood "
            "containing multiple exposed non-interface positions and is a "
            "strong structural candidate for within-species NKG2A-specific "
            "recognition."
        )

    if exposed_noninterface_count >= 1 and interface_count == 0:
        return (
            f"{labels} contains an exposed non-interface discriminatory "
            "position and is a favorable within-species NKG2A-specific "
            "surface candidate."
        )

    if exposed_noninterface_count >= 1 and interface_count > 0:
        return (
            f"{labels} contains accessible discriminatory residue(s), but "
            "the local region also contains receptor/ligand-interface "
            "involvement. It remains a candidate, with mechanism-dependent "
            "interpretation."
        )

    if interface_count > 0:
        return (
            f"{labels} is a discriminatory structural region with interface "
            "involvement. It is retained because interface-directed antibodies "
            "may still be useful, but it is not a clean non-interface "
            "specificity surface."
        )

    if buried_count == len(members):
        return (
            f"{labels} is sequence-discriminatory but buried in the human "
            "3CDG structural context and is therefore a low-priority "
            "standalone antibody-accessible region."
        )

    return (
        f"{labels} is retained as a within-species NKG2A-vs-NKG2C "
        "discriminatory region, but current structural evidence does not "
        "place it among the strongest exposed non-interface candidates."
    )


# =============================================================================
# REGION CONSTRUCTION
# =============================================================================

def build_region_row(
    species: str,
    region_id: str,
    region_type: str,
    members: Sequence[Dict[str, Any]],
    glyco_sites: Dict[str, List[Tuple[int, str]]],
    alphafold_decision: str,
) -> Dict[str, Any]:
    members = sorted(
        members,
        key=lambda member: member["residue"],
    )

    residues = [member["residue"] for member in members]

    residue_labels = ",".join(
        f"{member['residue']}{member['nkg2a_aa']}"
        for member in members
    )

    nkg2a_sequence = "|".join(
        f"{member['residue']}{member['nkg2a_aa']}"
        for member in members
    )

    nkg2c_sequence = "|".join(
        f"{member['residue']}{member['nkg2c_display']}"
        for member in members
    )

    comparison_labels = "|".join(
        f"{member['residue']}{member['nkg2a_aa']}>{member['nkg2c_display']}"
        for member in members
    )

    resolved_count = sum(
        1
        for member in members
        if member["resolved"]
    )

    unresolved_count = sum(
        1
        for member in members
        if member["unresolved"]
    )

    exposed_noninterface_count = sum(
        1
        for member in members
        if member["exposed_noninterface"]
    )

    exposed_interface_count = sum(
        1
        for member in members
        if member["exposed_interface"]
    )

    partially_exposed_noninterface_count = sum(
        1
        for member in members
        if member["partially_exposed_noninterface"]
    )

    partially_exposed_interface_count = sum(
        1
        for member in members
        if member["partially_exposed_interface"]
    )

    interface_count = sum(
        1
        for member in members
        if member["interface"]
    )

    buried_count = sum(
        1
        for member in members
        if member["buried"]
    )

    mean_rsa = mean(
        member["complex_rsa"]
        for member in members
    )

    mean_af_plddt = mean(
        member["alphafold_plddt"]
        for member in members
    )

    glyco_context, glyco_labels = glycosylation_context(
        species,
        residues,
        glyco_sites,
    )

    footprint_ids: Set[str] = set()

    for member in members:
        text = member["compact_footprint_ids"]

        if not text:
            continue

        for token in re.split(r"[,;| ]+", text):
            token = token.strip()
            if token:
                footprint_ids.add(token)

    if region_type == "sequence_defined_unresolved":
        structural_evidence_basis = (
            "sequence-defined; no experimental NKG2A coordinates for this "
            "region"
        )

        if species == "human" and alphafold_decision:
            structural_evidence_basis += (
                f"; AlphaFold Step 2Q={alphafold_decision}"
            )

    elif species == "human":
        structural_evidence_basis = (
            "direct homologous human NKG2A experimental structural evidence "
            "from 3CDG"
        )

    else:
        structural_evidence_basis = (
            "homologous-position structural annotation projected from human "
            "NKG2A 3CDG; not a direct macaque structure"
        )

    return {
        "species": species,
        "region_id": region_id,
        "region_type": region_type,
        "region_start": min(residues),
        "region_end": max(residues),
        "residue_count": len(members),
        "residue_labels": residue_labels,
        "within_species_comparisons": comparison_labels,
        "NKG2A_sequence": nkg2a_sequence,
        "NKG2C_sequence": nkg2c_sequence,
        "discriminatory_residue_count": len(members),
        "resolved_residue_count": resolved_count,
        "unresolved_residue_count": unresolved_count,
        "exposed_noninterface_count": exposed_noninterface_count,
        "exposed_interface_count": exposed_interface_count,
        "partially_exposed_noninterface_count": (
            partially_exposed_noninterface_count
        ),
        "partially_exposed_interface_count": (
            partially_exposed_interface_count
        ),
        "interface_residue_count": interface_count,
        "buried_residue_count": buried_count,
        "mean_complex_rsa": fmt_float(mean_rsa, 4),
        "mean_alphafold_plddt": fmt_float(mean_af_plddt, 3),
        "glycosylation_context": glyco_context,
        "glycosylation_sequons": glyco_labels,
        "compact_footprint_ids": ",".join(sorted(footprint_ids)),
        "structural_evidence_basis": structural_evidence_basis,
        "region_priority": region_priority(
            region_type,
            members,
        ),
        "interpretation": region_interpretation(
            species,
            region_type,
            members,
        ),
    }


def build_species_regions(
    species: str,
    candidates: Sequence[Dict[str, Any]],
    distance_map: Dict[Tuple[int, int], float],
    glyco_sites: Dict[str, List[Tuple[int, str]]],
    alphafold_decision: str,
) -> List[Dict[str, Any]]:
    config = SPECIES_CONFIG[species]

    # Safety check: keep only actual within-species discriminators.
    candidates = [
        candidate
        for candidate in candidates
        if candidate["discriminatory"]
    ]

    unresolved = [
        candidate
        for candidate in candidates
        if candidate["unresolved"]
    ]

    resolved = [
        candidate
        for candidate in candidates
        if candidate["resolved"]
    ]

    unresolved_by_residue = {
        candidate["residue"]: candidate
        for candidate in unresolved
    }

    resolved_by_residue = {
        candidate["residue"]: candidate
        for candidate in resolved
    }

    unresolved_components = sequence_components(
        list(unresolved_by_residue),
        SEQUENCE_REGION_MAX_GAP,
    )

    resolved_components = connected_components(
        list(resolved_by_residue),
        distance_map,
        STRUCTURAL_PATCH_CUTOFF_A,
    )

    provisional: List[
        Tuple[str, List[Dict[str, Any]]]
    ] = []

    for component in unresolved_components:
        members = [
            unresolved_by_residue[residue]
            for residue in component
        ]

        provisional.append(
            (
                "sequence_defined_unresolved",
                members,
            )
        )

    for component in resolved_components:
        members = [
            resolved_by_residue[residue]
            for residue in component
        ]

        if len(component) == 1:
            region_type = "resolved_structural_site"
        else:
            region_type = "resolved_structural_patch"

        provisional.append(
            (
                region_type,
                members,
            )
        )

    # Sort regions by a biologically useful hierarchy:
    #   1. strong resolved exposed regions
    #   2. sequence-defined unresolved regions
    #   3. other resolved regions
    #
    # Final region IDs are assigned after this sort.
    def provisional_sort_key(
        item: Tuple[str, List[Dict[str, Any]]],
    ) -> Tuple[int, int, int]:
        region_type, members = item

        priority = region_priority(
            region_type,
            members,
        )

        priority_order = {
            "priority_resolved_exposed_noninterface_patch": 0,
            "priority_resolved_exposed_noninterface_site": 1,
            "priority_resolved_mixed_surface_patch": 2,
            "priority_sequence_defined_unresolved_region": 3,
            "secondary_sequence_defined_unresolved_region": 4,
            "secondary_resolved_interface_patch": 5,
            "secondary_partially_exposed_noninterface_patch": 6,
            "secondary_partially_exposed_interface_patch": 7,
            "low_priority_buried_region": 8,
            "sequence_defined_unresolved_region": 9,
            "lower_structural_priority": 10,
        }

        first_residue = min(
            member["residue"]
            for member in members
        )

        return (
            priority_order.get(priority, 99),
            -len(members),
            first_residue,
        )

    provisional.sort(key=provisional_sort_key)

    regions: List[Dict[str, Any]] = []

    for index, (region_type, members) in enumerate(
        provisional,
        start=1,
    ):
        region_id = f"{config['prefix']}{index}"

        region = build_region_row(
            species=species,
            region_id=region_id,
            region_type=region_type,
            members=members,
            glyco_sites=glyco_sites,
            alphafold_decision=alphafold_decision,
        )

        regions.append(region)

    return regions


# =============================================================================
# REPORTING
# =============================================================================

OUTPUT_FIELDS = [
    "species",
    "region_id",
    "region_type",
    "region_start",
    "region_end",
    "residue_count",
    "residue_labels",
    "within_species_comparisons",
    "NKG2A_sequence",
    "NKG2C_sequence",
    "discriminatory_residue_count",
    "resolved_residue_count",
    "unresolved_residue_count",
    "exposed_noninterface_count",
    "exposed_interface_count",
    "partially_exposed_noninterface_count",
    "partially_exposed_interface_count",
    "interface_residue_count",
    "buried_residue_count",
    "mean_complex_rsa",
    "mean_alphafold_plddt",
    "glycosylation_context",
    "glycosylation_sequons",
    "compact_footprint_ids",
    "structural_evidence_basis",
    "region_priority",
    "interpretation",
]


def print_species_report(
    species: str,
    regions: Sequence[Dict[str, Any]],
) -> None:
    print_subbanner(
        f"{species.upper()} WITHIN-SPECIES CANDIDATE EPITOPE REGIONS"
    )

    if not regions:
        print("No candidate regions.")
        return

    for region in regions:
        print()
        print(
            f"{region['region_id']}  "
            f"{region['residue_labels']}"
        )

        print(
            f"  Type:       {region['region_type']}"
        )

        print(
            f"  Comparison: {region['within_species_comparisons']}"
        )

        print(
            "  Evidence:   "
            f"resolved={region['resolved_residue_count']}  "
            f"unresolved={region['unresolved_residue_count']}  "
            f"exposed_noninterface="
            f"{region['exposed_noninterface_count']}  "
            f"interface={region['interface_residue_count']}  "
            f"buried={region['buried_residue_count']}"
        )

        rsa = region["mean_complex_rsa"]
        if rsa:
            print(f"  Mean RSA:   {rsa}")

        af = region["mean_alphafold_plddt"]
        if af and region["unresolved_residue_count"]:
            print(f"  Mean AF:    {af}")

        print(
            f"  Glyco:      {region['glycosylation_context']}"
            + (
                f" ({region['glycosylation_sequons']})"
                if region["glycosylation_sequons"]
                else ""
            )
        )

        print(
            f"  Priority:   {region['region_priority']}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print_banner(
        "STEP 2S - BUILD WITHIN-SPECIES NKG2A-SPECIFIC CANDIDATE EPITOPE REGIONS"
    )

    print()
    print("Analysis axis:")
    print("  human   NKG2A vs human NKG2C")
    print("  rhesus  NKG2A vs BOTH rhesus NKG2C isoforms")
    print("  pigtail NKG2A vs pigtail NKG2C")
    print()
    print(
        "Cross-species conservation/reactivity is NOT used as a ranking criterion."
    )
    print(
        f"Resolved structural patch cutoff: {STRUCTURAL_PATCH_CUTOFF_A:.1f} A"
    )
    print(
        f"Unresolved sequence-region maximum gap: {SEQUENCE_REGION_MAX_GAP} residues"
    )

    # -------------------------------------------------------------------------
    # Load supporting evidence.
    # -------------------------------------------------------------------------

    distance_map = load_distance_map(DISTANCE_INPUT)

    print()
    print(
        f"Structural candidate-pair distances loaded: {len(distance_map)}"
    )

    glyco_sites = load_glycosylation_sites(
        GLYCOSYLATION_INPUT
    )

    print(
        "Canonical NKG2A glycosylation sequons loaded:"
    )
    for species in ("human", "rhesus", "pigtail"):
        print(
            f"  {species:<8} {len(glyco_sites.get(species, []))}"
        )

    alphafold_decision = load_alphafold_decision(
        ALPHAFOLD_VALIDATION_INPUT
    )

    print()
    print("AlphaFold Step 2Q interpretation:")
    print(
        f"  {alphafold_decision if alphafold_decision else 'not available'}"
    )

    # -------------------------------------------------------------------------
    # Load and normalize species-specific Step 2R candidate tables.
    # -------------------------------------------------------------------------

    species_candidates: Dict[
        str,
        List[Dict[str, Any]],
    ] = {}

    for species, config in SPECIES_CONFIG.items():
        raw_rows = read_tsv(config["input"])

        normalized: List[Dict[str, Any]] = []

        for row in raw_rows:
            candidate = normalize_candidate_row(
                species,
                row,
            )

            if candidate is None:
                continue

            normalized.append(candidate)

        normalized.sort(
            key=lambda candidate: candidate["residue"]
        )

        species_candidates[species] = normalized

    print()
    print("Step 2R candidate rows loaded:")
    for species in ("human", "rhesus", "pigtail"):
        candidates = species_candidates[species]

        discriminatory = sum(
            1
            for candidate in candidates
            if candidate["discriminatory"]
        )

        print(
            f"  {species:<8} "
            f"{len(candidates):>2} rows; "
            f"{discriminatory:>2} within-species discriminators"
        )

    # -------------------------------------------------------------------------
    # Build regions.
    # -------------------------------------------------------------------------

    all_regions: List[Dict[str, Any]] = []

    species_regions: Dict[
        str,
        List[Dict[str, Any]],
    ] = {}

    for species in ("human", "rhesus", "pigtail"):
        regions = build_species_regions(
            species=species,
            candidates=species_candidates[species],
            distance_map=distance_map,
            glyco_sites=glyco_sites,
            alphafold_decision=alphafold_decision,
        )

        species_regions[species] = regions
        all_regions.extend(regions)

    # -------------------------------------------------------------------------
    # Write outputs.
    # -------------------------------------------------------------------------

    for species, config in SPECIES_CONFIG.items():
        write_tsv(
            config["output"],
            species_regions[species],
            OUTPUT_FIELDS,
        )

    write_tsv(
        COMBINED_OUTPUT,
        all_regions,
        OUTPUT_FIELDS,
    )

    # -------------------------------------------------------------------------
    # Console report.
    # -------------------------------------------------------------------------

    for species in ("human", "rhesus", "pigtail"):
        print_species_report(
            species,
            species_regions[species],
        )

    # -------------------------------------------------------------------------
    # Summary.
    # -------------------------------------------------------------------------

    print_subbanner("STEP 2S SUMMARY")

    for species in ("human", "rhesus", "pigtail"):
        regions = species_regions[species]

        resolved_regions = sum(
            1
            for region in regions
            if region["resolved_residue_count"] > 0
        )

        unresolved_regions = sum(
            1
            for region in regions
            if region["region_type"]
            == "sequence_defined_unresolved"
        )

        priority_regions = sum(
            1
            for region in regions
            if clean(region["region_priority"]).startswith(
                "priority_"
            )
        )

        print(
            f"{species:<8} "
            f"{len(regions):>2} total region(s); "
            f"{resolved_regions:>2} structurally resolved; "
            f"{unresolved_regions:>2} sequence-defined unresolved; "
            f"{priority_regions:>2} priority"
        )

    print_subbanner("OUTPUTS")

    print()
    print(HUMAN_OUTPUT)
    print(RHESUS_OUTPUT)
    print(PIGTAIL_OUTPUT)
    print(COMBINED_OUTPUT)

    print()
    print(
        "NOTE: Step 2S ranks candidate epitope regions independently for "
        "within-species NKG2A-vs-NKG2C specificity."
    )
    print(
        "No cross-species NKG2A conservation or antibody cross-reactivity "
        "criterion is used."
    )
    print(
        "Human experimentally unresolved N-terminal residues are treated as "
        "sequence-defined regions rather than AlphaFold-defined 3D epitopes."
    )
    print(
        "Rhesus and pigtail structural accessibility/interface annotations "
        "are homologous-position projections from human NKG2A 3CDG."
    )
    print(
        "Candidate regions are screening hypotheses and are not experimentally "
        "validated antibody epitopes."
    )


if __name__ == "__main__":
    main()
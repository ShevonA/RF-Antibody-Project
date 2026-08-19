#!/usr/bin/env python3

"""
STEP 2T - BUILD LOCAL NKG2A-SPECIFIC EPITOPE CORES

Purpose
-------
Refine Step 2S within-species candidate epitope regions into compact local
epitope cores.

Axis 1 only:

    human:
        human NKG2A vs human NKG2C

    rhesus macaque:
        rhesus NKG2A vs BOTH rhesus NKG2C isoforms

    pigtail macaque:
        pigtail NKG2A vs pigtail NKG2C

Cross-species NKG2A conservation and antibody cross-reactivity are NOT used.

Important difference from Step 2S
---------------------------------
Step 2S used connected structural components. A chain such as:

    residue A within 10 A of B
    residue B within 10 A of C

can place A, B, and C in the same connected component even if A and C are
far apart.

Step 2T avoids that problem.

For every resolved discriminatory residue, this step constructs a LOCAL CORE
containing only discriminatory residues that lie within a fixed heavy-atom
distance of that same center residue.

Thus every residue in a local core must be directly close to its defining
center.

Unresolved N-terminal regions are preserved as sequence-defined cores and
are NOT clustered using AlphaFold coordinates.

Inputs
------
results/tables/structure/
    nkg2a_within_species_candidate_epitope_regions.tsv
    human_NKG2A_specificity_candidates.tsv
    rhesus_NKG2A_specificity_candidates.tsv
    pigtail_NKG2A_specificity_candidates.tsv
    nkg2a_candidate_spatial_distances.tsv
    nkg2_ectodomain_n_glycosylation_sites.tsv
    alphafold_nkg2a_model_validation.tsv

Outputs
-------
results/tables/structure/
    human_NKG2A_local_epitope_cores.tsv
    rhesus_NKG2A_local_epitope_cores.tsv
    pigtail_NKG2A_local_epitope_cores.tsv
    nkg2a_within_species_local_epitope_cores.tsv
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

STRUCTURE_DIR = (
    ROOT
    / "results"
    / "tables"
    / "structure"
)

REGION_INPUT = (
    STRUCTURE_DIR
    / "nkg2a_within_species_candidate_epitope_regions.tsv"
)

HUMAN_CANDIDATE_INPUT = (
    STRUCTURE_DIR
    / "human_NKG2A_specificity_candidates.tsv"
)

RHESUS_CANDIDATE_INPUT = (
    STRUCTURE_DIR
    / "rhesus_NKG2A_specificity_candidates.tsv"
)

PIGTAIL_CANDIDATE_INPUT = (
    STRUCTURE_DIR
    / "pigtail_NKG2A_specificity_candidates.tsv"
)

DISTANCE_INPUT = (
    STRUCTURE_DIR
    / "nkg2a_candidate_spatial_distances.tsv"
)

GLYCOSYLATION_INPUT = (
    STRUCTURE_DIR
    / "nkg2_ectodomain_n_glycosylation_sites.tsv"
)

ALPHAFOLD_VALIDATION_INPUT = (
    STRUCTURE_DIR
    / "alphafold_nkg2a_model_validation.tsv"
)

HUMAN_OUTPUT = (
    STRUCTURE_DIR
    / "human_NKG2A_local_epitope_cores.tsv"
)

RHESUS_OUTPUT = (
    STRUCTURE_DIR
    / "rhesus_NKG2A_local_epitope_cores.tsv"
)

PIGTAIL_OUTPUT = (
    STRUCTURE_DIR
    / "pigtail_NKG2A_local_epitope_cores.tsv"
)

COMBINED_OUTPUT = (
    STRUCTURE_DIR
    / "nkg2a_within_species_local_epitope_cores.tsv"
)


# Every resolved member in a core must lie within this distance of the
# defining center residue.
LOCAL_CORE_CUTOFF_A = 8.0

# Nearby glycosylation sequons are reported within this number of sequence
# residues from the minimum/maximum residue number represented by a core.
GLYCOSYLATION_CONTEXT_WINDOW = 5


# =============================================================================
# GENERAL HELPERS
# =============================================================================

def clean(value: Any) -> str:
    if value is None:
        return ""

    return str(value).strip()


def yes(value: Any) -> bool:
    return clean(value).lower() in {
        "yes",
        "y",
        "true",
        "1",
        "t",
    }


def parse_int(value: Any) -> Optional[int]:
    text = clean(value)

    if not text:
        return None

    try:
        return int(float(text))

    except ValueError:
        return None


def parse_float(value: Any) -> Optional[float]:
    text = clean(value)

    if not text:
        return None

    try:
        result = float(text)

    except ValueError:
        return None

    if not math.isfinite(result):
        return None

    return result


def mean(
    values: Iterable[Optional[float]],
) -> Optional[float]:

    usable = [
        value
        for value in values
        if value is not None
    ]

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


def first_nonempty(
    row: Dict[str, str],
    names: Sequence[str],
) -> str:

    for name in names:

        value = clean(
            row.get(
                name,
                "",
            )
        )

        if value:
            return value

    return ""


def read_tsv(
    path: Path,
) -> List[Dict[str, str]]:

    if not path.exists():

        raise FileNotFoundError(
            f"Required input file not found:\n{path}"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        return list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )


def write_tsv(
    path: Path,
    rows: Sequence[Dict[str, Any]],
    fieldnames: Sequence[str],
) -> None:

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )

        writer.writeheader()

        for row in rows:

            output = {}

            for field in fieldnames:

                value = row.get(
                    field,
                    "",
                )

                if value is None:
                    value = ""

                output[field] = value

            writer.writerow(
                output
            )


def print_banner(
    title: str,
) -> None:

    print("=" * 78)
    print(title)
    print("=" * 78)


def print_section(
    title: str,
) -> None:

    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


# =============================================================================
# SPECIES CONFIGURATION
# =============================================================================

SPECIES_CONFIG = {

    "human": {
        "prefix": "HC",
        "candidate_file":
            HUMAN_CANDIDATE_INPUT,

        "residue_columns": [
            "human_NKG2A_residue",
        ],

        "nkg2a_aa_columns": [
            "human_NKG2A_aa",
        ],

        "nkg2c_columns": [
            "human_NKG2C_aa",
        ],

        "tier_columns": [
            "human_evidence_tier",
        ],
    },

    "rhesus": {
        "prefix": "RC",
        "candidate_file":
            RHESUS_CANDIDATE_INPUT,

        "residue_columns": [
            "rhesus_NKG2A_residue",
        ],

        "nkg2a_aa_columns": [
            "rhesus_NKG2A_aa",
        ],

        "nkg2c1_columns": [
            "rhesus_NKG2C1_aa",
        ],

        "nkg2c2_columns": [
            "rhesus_NKG2C2_aa",
        ],

        "tier_columns": [
            "rhesus_evidence_tier",
        ],
    },

    "pigtail": {
        "prefix": "PC",
        "candidate_file":
            PIGTAIL_CANDIDATE_INPUT,

        "residue_columns": [
            "pigtail_NKG2A_residue",
        ],

        "nkg2a_aa_columns": [
            "pigtail_NKG2A_aa",
        ],

        "nkg2c_columns": [
            "pigtail_NKG2C_aa",
        ],

        "tier_columns": [
            "pigtail_evidence_tier",
        ],
    },
}


# =============================================================================
# RESIDUE LABEL PARSING
# =============================================================================

def parse_residue_labels(
    value: Any,
) -> List[int]:

    residues = []

    for label in clean(value).split(","):

        label = label.strip()

        if not label:
            continue

        digits = ""

        for character in label:

            if character.isdigit():
                digits += character

            else:
                break

        if digits:

            residues.append(
                int(digits)
            )

    return sorted(
        set(residues)
    )


# =============================================================================
# DISTANCE MAP
# =============================================================================

def load_distance_map(
    path: Path,
) -> Dict[Tuple[int, int], float]:

    rows = read_tsv(
        path
    )

    distances = {}

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

        key = tuple(
            sorted(
                (
                    residue_1,
                    residue_2,
                )
            )
        )

        distances[
            key
        ] = distance

    return distances


def residue_distance(
    residue_1: int,
    residue_2: int,
    distance_map: Dict[Tuple[int, int], float],
) -> Optional[float]:

    if residue_1 == residue_2:
        return 0.0

    key = tuple(
        sorted(
            (
                residue_1,
                residue_2,
            )
        )
    )

    return distance_map.get(
        key
    )


# =============================================================================
# CANDIDATE NORMALIZATION
# =============================================================================

def normalize_candidate(
    species: str,
    row: Dict[str, str],
) -> Optional[Dict[str, Any]]:

    config = SPECIES_CONFIG[
        species
    ]

    residue = parse_int(
        first_nonempty(
            row,
            config[
                "residue_columns"
            ],
        )
    )

    if residue is None:
        return None

    nkg2a_aa = first_nonempty(
        row,
        config[
            "nkg2a_aa_columns"
        ],
    )

    if species == "rhesus":

        nkg2c1 = first_nonempty(
            row,
            config[
                "nkg2c1_columns"
            ],
        )

        nkg2c2 = first_nonempty(
            row,
            config[
                "nkg2c2_columns"
            ],
        )

        nkg2c_display = (
            f"{nkg2c1}/{nkg2c2}"
        )

        discriminatory = (
            bool(nkg2a_aa)
            and bool(nkg2c1)
            and bool(nkg2c2)
            and nkg2a_aa != nkg2c1
            and nkg2a_aa != nkg2c2
        )

    else:

        nkg2c = first_nonempty(
            row,
            config[
                "nkg2c_columns"
            ],
        )

        nkg2c_display = (
            nkg2c
        )

        discriminatory = (
            bool(nkg2a_aa)
            and bool(nkg2c)
            and nkg2a_aa != nkg2c
        )

    tier = first_nonempty(
        row,
        config[
            "tier_columns"
        ],
    )

    structure_class = first_nonempty(
        row,
        [
            "experimental_structure_class",
        ],
    )

    rsa = parse_float(
        first_nonempty(
            row,
            [
                "complex_rsa",
            ],
        )
    )

    interface = yes(
        first_nonempty(
            row,
            [
                "any_interface_contact",
            ],
        )
    )

    contact_cd94 = yes(
        row.get(
            "contact_CD94",
            "",
        )
    )

    contact_hla_e = yes(
        row.get(
            "contact_HLA_E",
            "",
        )
    )

    contact_peptide = yes(
        row.get(
            "contact_peptide",
            "",
        )
    )

    contact_b2m = yes(
        row.get(
            "contact_B2M",
            "",
        )
    )

    unresolved = (
        "unresolved"
        in tier.lower()
        or
        "no_experiment"
        in structure_class.lower()
    )

    resolved = not unresolved

    exposed = (
        resolved
        and rsa is not None
        and rsa >= 0.25
    )

    partially_exposed = (
        resolved
        and rsa is not None
        and 0.10 <= rsa < 0.25
    )

    buried = (
        resolved
        and rsa is not None
        and rsa < 0.10
    )

    exposed_noninterface = (
        exposed
        and not interface
    )

    exposed_interface = (
        exposed
        and interface
    )

    partially_exposed_noninterface = (
        partially_exposed
        and not interface
    )

    partially_exposed_interface = (
        partially_exposed
        and interface
    )

    plddt = parse_float(
        first_nonempty(
            row,
            [
                "alphafold_plddt",
            ],
        )
    )

    return {
        "species":
            species,

        "residue":
            residue,

        "nkg2a_aa":
            nkg2a_aa,

        "nkg2c_display":
            nkg2c_display,

        "discriminatory":
            discriminatory,

        "tier":
            tier,

        "structure_class":
            structure_class,

        "rsa":
            rsa,

        "interface":
            interface,

        "contact_cd94":
            contact_cd94,

        "contact_hla_e":
            contact_hla_e,

        "contact_peptide":
            contact_peptide,

        "contact_b2m":
            contact_b2m,

        "resolved":
            resolved,

        "unresolved":
            unresolved,

        "exposed":
            exposed,

        "partially_exposed":
            partially_exposed,

        "buried":
            buried,

        "exposed_noninterface":
            exposed_noninterface,

        "exposed_interface":
            exposed_interface,

        "partially_exposed_noninterface":
            partially_exposed_noninterface,

        "partially_exposed_interface":
            partially_exposed_interface,

        "plddt":
            plddt,

        "raw":
            row,
    }


def load_species_candidates(
    species: str,
) -> Dict[int, Dict[str, Any]]:

    path = SPECIES_CONFIG[
        species
    ][
        "candidate_file"
    ]

    rows = read_tsv(
        path
    )

    lookup = {}

    for row in rows:

        candidate = normalize_candidate(
            species,
            row,
        )

        if candidate is None:
            continue

        if not candidate[
            "discriminatory"
        ]:
            continue

        lookup[
            candidate[
                "residue"
            ]
        ] = candidate

    return lookup


# =============================================================================
# GLYCOSYLATION
# =============================================================================

GLYCO_RECORD_IDS = {
    "human":
        "human_NKG2A",

    "rhesus":
        "rhesus_NKG2A",

    "pigtail":
        "pigtail_NKG2A",
}


def load_glycosylation_sites(
    path: Path,
) -> Dict[str, List[Tuple[int, str]]]:

    rows = read_tsv(
        path
    )

    reverse = {
        record_id: species
        for species, record_id
        in GLYCO_RECORD_IDS.items()
    }

    result = defaultdict(
        list
    )

    for row in rows:

        record_id = clean(
            row.get(
                "record_id"
            )
        )

        species = reverse.get(
            record_id
        )

        if species is None:
            continue

        residue = parse_int(
            row.get(
                "sequon_full_length_residue"
            )
        )

        motif = clean(
            row.get(
                "motif"
            )
        )

        if residue is None:
            continue

        result[
            species
        ].append(
            (
                residue,
                motif,
            )
        )

    for species in result:

        result[
            species
        ].sort()

    return result


def glycosylation_context(
    species: str,
    residues: Sequence[int],
    sites: Dict[str, List[Tuple[int, str]]],
) -> Tuple[str, str]:

    if not residues:

        return (
            "none",
            "",
        )

    minimum = min(
        residues
    )

    maximum = max(
        residues
    )

    overlapping = []

    nearby = []

    for position, motif in sites.get(
        species,
        [],
    ):

        label = (
            f"{position}{motif}"
        )

        if (
            minimum
            <= position
            <= maximum
        ):

            overlapping.append(
                label
            )

        elif (
            minimum
            - GLYCOSYLATION_CONTEXT_WINDOW
            <= position
            <= maximum
            + GLYCOSYLATION_CONTEXT_WINDOW
        ):

            nearby.append(
                label
            )

    if overlapping:

        description = (
            ",".join(
                overlapping
            )
        )

        if nearby:

            description += (
                ";nearby="
                + ",".join(
                    nearby
                )
            )

        return (
            "overlapping_canonical_N_glycosylation_sequon",
            description,
        )

    if nearby:

        return (
            "nearby_canonical_N_glycosylation_sequon",
            ",".join(
                nearby
            ),
        )

    return (
        "none_detected_near_core",
        "",
    )


# =============================================================================
# ALPHAFOLD
# =============================================================================

def load_alphafold_decision(
    path: Path,
) -> str:

    rows = read_tsv(
        path
    )

    if not rows:
        return ""

    return first_nonempty(
        rows[0],
        [
            "model_use_decision",
        ],
    )


# =============================================================================
# STEP 2S REGION LOOKUP
# =============================================================================

def load_regions(
    path: Path,
) -> Dict[str, List[Dict[str, str]]]:

    rows = read_tsv(
        path
    )

    result = defaultdict(
        list
    )

    for row in rows:

        species = clean(
            row.get(
                "species"
            )
        ).lower()

        if species not in (
            "human",
            "rhesus",
            "pigtail",
        ):
            continue

        result[
            species
        ].append(
            row
        )

    return result


# =============================================================================
# LOCAL CORE CREATION
# =============================================================================

def local_members_for_center(
    center: int,
    allowed_residues: Set[int],
    distance_map: Dict[Tuple[int, int], float],
) -> List[int]:

    members = [
        center
    ]

    for residue in sorted(
        allowed_residues
    ):

        if residue == center:
            continue

        distance = residue_distance(
            center,
            residue,
            distance_map,
        )

        if (
            distance is not None
            and distance
            <= LOCAL_CORE_CUTOFF_A
        ):

            members.append(
                residue
            )

    return sorted(
        set(members)
    )


def max_center_distance(
    center: int,
    members: Sequence[int],
    distance_map: Dict[Tuple[int, int], float],
) -> Optional[float]:

    distances = []

    for residue in members:

        distance = residue_distance(
            center,
            residue,
            distance_map,
        )

        if distance is not None:

            distances.append(
                distance
            )

    if not distances:
        return None

    return max(
        distances
    )


def pairwise_core_diameter(
    members: Sequence[int],
    distance_map: Dict[Tuple[int, int], float],
) -> Optional[float]:
    """
    Maximum known pairwise minimum-heavy-atom distance among members.

    This is descriptive only.

    Missing pair distances are ignored because the original distance table
    contains candidate pairs rather than a guaranteed complete all-residue
    matrix.
    """

    distances = []

    for index, residue_1 in enumerate(
        members
    ):

        for residue_2 in members[
            index + 1:
        ]:

            distance = residue_distance(
                residue_1,
                residue_2,
                distance_map,
            )

            if distance is not None:

                distances.append(
                    distance
                )

    if not distances:

        return 0.0

    return max(
        distances
    )


# =============================================================================
# CORE CLASSIFICATION
# =============================================================================

def core_priority(
    core_type: str,
    members: Sequence[Dict[str, Any]],
) -> str:

    if core_type == "sequence_defined_unresolved":

        if len(members) >= 3:

            return (
                "priority_sequence_defined_unresolved_core"
            )

        return (
            "secondary_sequence_defined_unresolved_core"
        )

    exposed_noninterface = sum(
        member[
            "exposed_noninterface"
        ]
        for member in members
    )

    exposed_interface = sum(
        member[
            "exposed_interface"
        ]
        for member in members
    )

    partial_noninterface = sum(
        member[
            "partially_exposed_noninterface"
        ]
        for member in members
    )

    partial_interface = sum(
        member[
            "partially_exposed_interface"
        ]
        for member in members
    )

    buried = sum(
        member[
            "buried"
        ]
        for member in members
    )

    interface = sum(
        member[
            "interface"
        ]
        for member in members
    )

    if (
        exposed_noninterface >= 2
        and interface == 0
    ):

        return (
            "priority_compact_exposed_noninterface_core"
        )

    if (
        exposed_noninterface == 1
        and interface == 0
        and len(members) == 1
    ):

        return (
            "priority_exposed_noninterface_single_site"
        )

    if (
        exposed_noninterface >= 1
        and interface == 0
    ):

        return (
            "priority_mixed_accessible_noninterface_core"
        )

    if (
        exposed_noninterface >= 1
        and interface >= 1
    ):

        return (
            "secondary_mixed_interface_core"
        )

    if exposed_interface >= 1:

        return (
            "secondary_exposed_interface_core"
        )

    if partial_noninterface >= 1:

        return (
            "secondary_partially_exposed_noninterface_core"
        )

    if partial_interface >= 1:

        return (
            "secondary_partially_exposed_interface_core"
        )

    if (
        buried == len(members)
        and len(members) > 0
    ):

        return (
            "low_priority_buried_core"
        )

    return (
        "lower_structural_priority"
    )


def core_interpretation(
    species: str,
    core_type: str,
    labels: str,
    members: Sequence[Dict[str, Any]],
) -> str:

    if (
        core_type
        == "sequence_defined_unresolved"
    ):

        if species == "human":

            return (
                f"{labels} is retained as a human NKG2A-vs-NKG2C "
                "sequence-defined specificity core. No experimental "
                "NKG2A coordinates exist for this region, and Step 2Q "
                "does not support treating the AlphaFold conformation "
                "as fixed epitope geometry."
            )

        return (
            f"{labels} is retained as a {species} NKG2A-vs-NKG2C "
            "sequence-defined specificity core corresponding to the "
            "experimentally unresolved N-terminal region in human NKG2A. "
            "No direct macaque structural geometry is inferred."
        )

    exposed_noninterface = sum(
        member[
            "exposed_noninterface"
        ]
        for member in members
    )

    interface = sum(
        member[
            "interface"
        ]
        for member in members
    )

    buried = sum(
        member[
            "buried"
        ]
        for member in members
    )

    if (
        exposed_noninterface >= 2
        and interface == 0
    ):

        return (
            f"{labels} is a compact local discriminatory surface "
            "containing multiple exposed non-interface positions and is "
            "a strong within-species NKG2A-specific antibody core candidate."
        )

    if (
        exposed_noninterface >= 1
        and interface == 0
    ):

        return (
            f"{labels} contains an exposed non-interface discriminatory "
            "position in a compact local neighborhood and remains a favorable "
            "within-species specificity core."
        )

    if (
        exposed_noninterface >= 1
        and interface > 0
    ):

        return (
            f"{labels} combines accessible discriminatory residue(s) with "
            "receptor/ligand-interface involvement. It is a plausible local "
            "core but has greater mechanism-dependent liability."
        )

    if interface > 0:

        return (
            f"{labels} is a compact discriminatory region dominated by "
            "interface-associated residues. It is retained as a secondary "
            "candidate rather than a clean non-interface specificity core."
        )

    if (
        buried
        == len(members)
        and members
    ):

        return (
            f"{labels} is discriminatory but buried in the human 3CDG "
            "structural context and is a low-priority standalone "
            "antibody-accessible core."
        )

    return (
        f"{labels} is a compact within-species NKG2A-vs-NKG2C "
        "discriminatory neighborhood with intermediate structural evidence."
    )


# =============================================================================
# CORE OUTPUT ROW
# =============================================================================

def build_core_row(
    species: str,
    core_id: str,
    source_region_id: str,
    core_type: str,
    center_residue: Optional[int],
    member_residues: Sequence[int],
    candidate_lookup: Dict[int, Dict[str, Any]],
    distance_map: Dict[Tuple[int, int], float],
    glyco_sites: Dict[str, List[Tuple[int, str]]],
    alphafold_decision: str,
) -> Dict[str, Any]:

    members = [
        candidate_lookup[
            residue
        ]
        for residue
        in member_residues
        if residue
        in candidate_lookup
    ]

    members.sort(
        key=lambda member:
        member[
            "residue"
        ]
    )

    residues = [
        member[
            "residue"
        ]
        for member
        in members
    ]

    residue_labels = ",".join(
        f"{member['residue']}"
        f"{member['nkg2a_aa']}"
        for member
        in members
    )

    comparisons = "|".join(
        f"{member['residue']}"
        f"{member['nkg2a_aa']}"
        f">{member['nkg2c_display']}"
        for member
        in members
    )

    resolved_count = sum(
        member[
            "resolved"
        ]
        for member
        in members
    )

    unresolved_count = sum(
        member[
            "unresolved"
        ]
        for member
        in members
    )

    exposed_noninterface_count = sum(
        member[
            "exposed_noninterface"
        ]
        for member
        in members
    )

    exposed_interface_count = sum(
        member[
            "exposed_interface"
        ]
        for member
        in members
    )

    partial_noninterface_count = sum(
        member[
            "partially_exposed_noninterface"
        ]
        for member
        in members
    )

    partial_interface_count = sum(
        member[
            "partially_exposed_interface"
        ]
        for member
        in members
    )

    interface_count = sum(
        member[
            "interface"
        ]
        for member
        in members
    )

    buried_count = sum(
        member[
            "buried"
        ]
        for member
        in members
    )

    cd94_count = sum(
        member[
            "contact_cd94"
        ]
        for member
        in members
    )

    hla_e_count = sum(
        member[
            "contact_hla_e"
        ]
        for member
        in members
    )

    peptide_count = sum(
        member[
            "contact_peptide"
        ]
        for member
        in members
    )

    b2m_count = sum(
        member[
            "contact_b2m"
        ]
        for member
        in members
    )

    mean_rsa = mean(
        member[
            "rsa"
        ]
        for member
        in members
    )

    mean_plddt = mean(
        member[
            "plddt"
        ]
        for member
        in members
    )

    glyco_context, glyco_labels = (
        glycosylation_context(
            species,
            residues,
            glyco_sites,
        )
    )

    if (
        center_residue
        is not None
    ):

        center_candidate = (
            candidate_lookup[
                center_residue
            ]
        )

        center_label = (
            f"{center_residue}"
            f"{center_candidate['nkg2a_aa']}"
        )

        maximum_center_distance = (
            max_center_distance(
                center_residue,
                residues,
                distance_map,
            )
        )

        diameter = (
            pairwise_core_diameter(
                residues,
                distance_map,
            )
        )

    else:

        center_label = ""

        maximum_center_distance = None

        diameter = None

    if (
        core_type
        == "sequence_defined_unresolved"
    ):

        geometry_basis = (
            "sequence_defined_only;"
            "no_experimental_fixed_geometry"
        )

        if (
            species
            == "human"
            and alphafold_decision
        ):

            geometry_basis += (
                ";AlphaFold="
                + alphafold_decision
            )

    elif species == "human":

        geometry_basis = (
            "human_NKG2A_3CDG_experimental_geometry"
        )

    else:

        geometry_basis = (
            "human_3CDG_homologous_position_geometry;"
            "not_direct_macaque_structure"
        )

    priority = core_priority(
        core_type,
        members,
    )

    interpretation = (
        core_interpretation(
            species,
            core_type,
            residue_labels,
            members,
        )
    )

    return {
        "species":
            species,

        "core_id":
            core_id,

        "source_region_id":
            source_region_id,

        "core_type":
            core_type,

        "center_residue":
            (
                center_residue
                if center_residue
                is not None
                else ""
            ),

        "center_label":
            center_label,

        "local_distance_cutoff_A":
            (
                f"{LOCAL_CORE_CUTOFF_A:.1f}"
                if center_residue
                is not None
                else ""
            ),

        "maximum_center_to_member_distance_A":
            fmt_float(
                maximum_center_distance,
                3,
            ),

        "known_pairwise_core_diameter_A":
            fmt_float(
                diameter,
                3,
            ),

        "residue_count":
            len(
                residues
            ),

        "residue_labels":
            residue_labels,

        "within_species_comparisons":
            comparisons,

        "resolved_residue_count":
            resolved_count,

        "unresolved_residue_count":
            unresolved_count,

        "exposed_noninterface_count":
            exposed_noninterface_count,

        "exposed_interface_count":
            exposed_interface_count,

        "partially_exposed_noninterface_count":
            partial_noninterface_count,

        "partially_exposed_interface_count":
            partial_interface_count,

        "interface_residue_count":
            interface_count,

        "CD94_contact_residue_count":
            cd94_count,

        "HLA_E_contact_residue_count":
            hla_e_count,

        "peptide_contact_residue_count":
            peptide_count,

        "B2M_contact_residue_count":
            b2m_count,

        "buried_residue_count":
            buried_count,

        "mean_complex_rsa":
            fmt_float(
                mean_rsa,
                4,
            ),

        "mean_alphafold_plddt":
            fmt_float(
                mean_plddt,
                3,
            ),

        "glycosylation_context":
            glyco_context,

        "glycosylation_sequons":
            glyco_labels,

        "geometry_evidence_basis":
            geometry_basis,

        "core_priority":
            priority,

        "interpretation":
            interpretation,
    }


# =============================================================================
# BUILD SPECIES CORES
# =============================================================================

def build_species_cores(
    species: str,
    regions: Sequence[Dict[str, str]],
    candidate_lookup: Dict[int, Dict[str, Any]],
    distance_map: Dict[Tuple[int, int], float],
    glyco_sites: Dict[str, List[Tuple[int, str]]],
    alphafold_decision: str,
) -> List[Dict[str, Any]]:

    provisional = []

    # -------------------------------------------------------------------------
    # Process every Step 2S region independently.
    # -------------------------------------------------------------------------

    for region in regions:

        region_id = clean(
            region.get(
                "region_id"
            )
        )

        region_type = clean(
            region.get(
                "region_type"
            )
        )

        region_residues = [
            residue
            for residue
            in parse_residue_labels(
                region.get(
                    "residue_labels"
                )
            )
            if residue
            in candidate_lookup
        ]

        if not region_residues:
            continue

        # ---------------------------------------------------------------------
        # Preserve unresolved regions exactly as sequence-defined cores.
        # ---------------------------------------------------------------------

        if (
            region_type
            == "sequence_defined_unresolved"
        ):

            provisional.append(
                {
                    "source_region_id":
                        region_id,

                    "core_type":
                        "sequence_defined_unresolved",

                    "center_residue":
                        None,

                    "member_residues":
                        tuple(
                            sorted(
                                region_residues
                            )
                        ),
                }
            )

            continue

        # ---------------------------------------------------------------------
        # Resolved region:
        # construct a direct local neighborhood around every center residue.
        # ---------------------------------------------------------------------

        allowed = set(
            region_residues
        )

        seen_member_sets = {}

        for center in sorted(
            region_residues
        ):

            center_candidate = (
                candidate_lookup.get(
                    center
                )
            )

            if (
                center_candidate
                is None
                or not center_candidate[
                    "resolved"
                ]
            ):

                continue

            members = (
                local_members_for_center(
                    center,
                    allowed,
                    distance_map,
                )
            )

            key = tuple(
                members
            )

            # If multiple centers yield the exact same local member set,
            # keep the structurally more attractive center.
            if key in seen_member_sets:

                old_center = (
                    seen_member_sets[
                        key
                    ]
                )

                old_candidate = (
                    candidate_lookup[
                        old_center
                    ]
                )

                def center_rank(
                    candidate:
                    Dict[str, Any],
                ) -> Tuple[int, int, float]:

                    return (
                        0
                        if candidate[
                            "exposed_noninterface"
                        ]
                        else 1,

                        0
                        if not candidate[
                            "interface"
                        ]
                        else 1,

                        -(
                            candidate[
                                "rsa"
                            ]
                            if candidate[
                                "rsa"
                            ]
                            is not None
                            else -1
                        ),
                    )

                if (
                    center_rank(
                        center_candidate
                    )
                    <
                    center_rank(
                        old_candidate
                    )
                ):

                    seen_member_sets[
                        key
                    ] = center

                continue

            seen_member_sets[
                key
            ] = center

        for (
            member_tuple,
            center,
        ) in seen_member_sets.items():

            provisional.append(
                {
                    "source_region_id":
                        region_id,

                    "core_type":
                        (
                            "resolved_local_site"
                            if len(
                                member_tuple
                            ) == 1
                            else
                            "resolved_local_core"
                        ),

                    "center_residue":
                        center,

                    "member_residues":
                        member_tuple,
                }
            )

    # -------------------------------------------------------------------------
    # Convert provisional objects to temporary core rows so we can rank them.
    # -------------------------------------------------------------------------

    temporary = []

    for item in provisional:

        temporary_row = (
            build_core_row(
                species=species,

                core_id="",

                source_region_id=item[
                    "source_region_id"
                ],

                core_type=item[
                    "core_type"
                ],

                center_residue=item[
                    "center_residue"
                ],

                member_residues=item[
                    "member_residues"
                ],

                candidate_lookup=
                    candidate_lookup,

                distance_map=
                    distance_map,

                glyco_sites=
                    glyco_sites,

                alphafold_decision=
                    alphafold_decision,
            )
        )

        temporary.append(
            temporary_row
        )

    # -------------------------------------------------------------------------
    # Sort by evidence priority.
    # -------------------------------------------------------------------------

    priority_order = {

        "priority_compact_exposed_noninterface_core":
            0,

        "priority_exposed_noninterface_single_site":
            1,

        "priority_mixed_accessible_noninterface_core":
            2,

        "priority_sequence_defined_unresolved_core":
            3,

        "secondary_sequence_defined_unresolved_core":
            4,

        "secondary_mixed_interface_core":
            5,

        "secondary_exposed_interface_core":
            6,

        "secondary_partially_exposed_noninterface_core":
            7,

        "secondary_partially_exposed_interface_core":
            8,

        "low_priority_buried_core":
            9,

        "lower_structural_priority":
            10,
    }

    temporary.sort(
        key=lambda row: (
            priority_order.get(
                row[
                    "core_priority"
                ],
                99,
            ),

            -int(
                row[
                    "exposed_noninterface_count"
                ]
            ),

            int(
                row[
                    "interface_residue_count"
                ]
            ),

            -int(
                row[
                    "residue_count"
                ]
            ),

            (
                int(
                    row[
                        "center_residue"
                    ]
                )
                if clean(
                    row[
                        "center_residue"
                    ]
                )
                else
                9999
            ),
        )
    )

    # -------------------------------------------------------------------------
    # Assign final species core IDs.
    # -------------------------------------------------------------------------

    prefix = (
        SPECIES_CONFIG[
            species
        ][
            "prefix"
        ]
    )

    for index, row in enumerate(
        temporary,
        start=1,
    ):

        row[
            "core_id"
        ] = (
            f"{prefix}{index}"
        )

    return temporary


# =============================================================================
# OUTPUT FIELDS
# =============================================================================

OUTPUT_FIELDS = [

    "species",

    "core_id",

    "source_region_id",

    "core_type",

    "center_residue",

    "center_label",

    "local_distance_cutoff_A",

    "maximum_center_to_member_distance_A",

    "known_pairwise_core_diameter_A",

    "residue_count",

    "residue_labels",

    "within_species_comparisons",

    "resolved_residue_count",

    "unresolved_residue_count",

    "exposed_noninterface_count",

    "exposed_interface_count",

    "partially_exposed_noninterface_count",

    "partially_exposed_interface_count",

    "interface_residue_count",

    "CD94_contact_residue_count",

    "HLA_E_contact_residue_count",

    "peptide_contact_residue_count",

    "B2M_contact_residue_count",

    "buried_residue_count",

    "mean_complex_rsa",

    "mean_alphafold_plddt",

    "glycosylation_context",

    "glycosylation_sequons",

    "geometry_evidence_basis",

    "core_priority",

    "interpretation",
]


# =============================================================================
# REPORT
# =============================================================================

def print_species_report(
    species: str,
    cores: Sequence[
        Dict[str, Any]
    ],
) -> None:

    print_section(
        f"{species.upper()} LOCAL NKG2A-SPECIFIC EPITOPE CORES"
    )

    if not cores:

        print(
            "No cores generated."
        )

        return

    for row in cores:

        print()

        print(
            f"{row['core_id']}  "
            f"{row['residue_labels']}"
        )

        print(
            f"  Source:     "
            f"{row['source_region_id']}"
        )

        print(
            f"  Type:       "
            f"{row['core_type']}"
        )

        if row[
            "center_label"
        ]:

            print(
                f"  Center:     "
                f"{row['center_label']}"
            )

            print(
                "  Radius:     "
                f"max center-member="
                f"{row['maximum_center_to_member_distance_A']} A"
            )

        print(
            f"  Comparison: "
            f"{row['within_species_comparisons']}"
        )

        print(
            "  Evidence:   "
            f"exposed_noninterface="
            f"{row['exposed_noninterface_count']}  "
            f"interface="
            f"{row['interface_residue_count']}  "
            f"buried="
            f"{row['buried_residue_count']}"
        )

        if row[
            "mean_complex_rsa"
        ]:

            print(
                f"  Mean RSA:   "
                f"{row['mean_complex_rsa']}"
            )

        if (
            row[
                "unresolved_residue_count"
            ]
            and row[
                "mean_alphafold_plddt"
            ]
        ):

            print(
                f"  Mean AF:    "
                f"{row['mean_alphafold_plddt']}"
            )

        print(
            f"  Glyco:      "
            f"{row['glycosylation_context']}"
            + (
                " ("
                + row[
                    "glycosylation_sequons"
                ]
                + ")"
                if row[
                    "glycosylation_sequons"
                ]
                else ""
            )
        )

        print(
            f"  Priority:   "
            f"{row['core_priority']}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    print_banner(
        "STEP 2T - BUILD LOCAL WITHIN-SPECIES NKG2A-SPECIFIC EPITOPE CORES"
    )

    print()

    print(
        "Analysis axis:"
    )

    print(
        "  human   NKG2A vs human NKG2C"
    )

    print(
        "  rhesus  NKG2A vs BOTH rhesus NKG2C isoforms"
    )

    print(
        "  pigtail NKG2A vs pigtail NKG2C"
    )

    print()

    print(
        "Cross-species conservation/reactivity is NOT used."
    )

    print(
        f"Local resolved-core cutoff: "
        f"{LOCAL_CORE_CUTOFF_A:.1f} A"
    )

    print(
        "Every resolved member must lie within the cutoff "
        "of the SAME defining center residue."
    )

    # =========================================================================
    # Inputs
    # =========================================================================

    regions_by_species = (
        load_regions(
            REGION_INPUT
        )
    )

    distance_map = (
        load_distance_map(
            DISTANCE_INPUT
        )
    )

    glyco_sites = (
        load_glycosylation_sites(
            GLYCOSYLATION_INPUT
        )
    )

    alphafold_decision = (
        load_alphafold_decision(
            ALPHAFOLD_VALIDATION_INPUT
        )
    )

    print()

    print(
        f"Structural pair distances loaded: "
        f"{len(distance_map)}"
    )

    print()

    print(
        "Step 2S regions loaded:"
    )

    for species in (
        "human",
        "rhesus",
        "pigtail",
    ):

        print(
            f"  {species:<8} "
            f"{len(regions_by_species.get(species, []))}"
        )

    print()

    print(
        "AlphaFold Step 2Q interpretation:"
    )

    print(
        f"  {alphafold_decision}"
    )

    # =========================================================================
    # Candidate lookups
    # =========================================================================

    candidate_lookups = {}

    for species in (
        "human",
        "rhesus",
        "pigtail",
    ):

        candidate_lookups[
            species
        ] = (
            load_species_candidates(
                species
            )
        )

    print()

    print(
        "Within-species discriminatory residues loaded:"
    )

    for species in (
        "human",
        "rhesus",
        "pigtail",
    ):

        print(
            f"  {species:<8} "
            f"{len(candidate_lookups[species])}"
        )

    # =========================================================================
    # Build cores
    # =========================================================================

    cores_by_species = {}

    all_cores = []

    for species in (
        "human",
        "rhesus",
        "pigtail",
    ):

        cores = (
            build_species_cores(
                species=species,

                regions=
                    regions_by_species.get(
                        species,
                        [],
                    ),

                candidate_lookup=
                    candidate_lookups[
                        species
                    ],

                distance_map=
                    distance_map,

                glyco_sites=
                    glyco_sites,

                alphafold_decision=
                    alphafold_decision,
            )
        )

        cores_by_species[
            species
        ] = cores

        all_cores.extend(
            cores
        )

    # =========================================================================
    # Write outputs
    # =========================================================================

    write_tsv(
        HUMAN_OUTPUT,
        cores_by_species[
            "human"
        ],
        OUTPUT_FIELDS,
    )

    write_tsv(
        RHESUS_OUTPUT,
        cores_by_species[
            "rhesus"
        ],
        OUTPUT_FIELDS,
    )

    write_tsv(
        PIGTAIL_OUTPUT,
        cores_by_species[
            "pigtail"
        ],
        OUTPUT_FIELDS,
    )

    write_tsv(
        COMBINED_OUTPUT,
        all_cores,
        OUTPUT_FIELDS,
    )

    # =========================================================================
    # Report
    # =========================================================================

    for species in (
        "human",
        "rhesus",
        "pigtail",
    ):

        print_species_report(
            species,
            cores_by_species[
                species
            ],
        )

    # =========================================================================
    # Summary
    # =========================================================================

    print_section(
        "STEP 2T SUMMARY"
    )

    for species in (
        "human",
        "rhesus",
        "pigtail",
    ):

        cores = (
            cores_by_species[
                species
            ]
        )

        sequence_defined = sum(
            1
            for core in cores
            if core[
                "core_type"
            ]
            == "sequence_defined_unresolved"
        )

        resolved = (
            len(cores)
            - sequence_defined
        )

        priority = sum(
            1
            for core in cores
            if clean(
                core[
                    "core_priority"
                ]
            ).startswith(
                "priority_"
            )
        )

        clean_noninterface = sum(
            1
            for core in cores
            if core[
                "core_priority"
            ]
            in {
                "priority_compact_exposed_noninterface_core",
                "priority_exposed_noninterface_single_site",
                "priority_mixed_accessible_noninterface_core",
            }
        )

        print(
            f"{species:<8} "
            f"{len(cores):>2} core(s); "
            f"{resolved:>2} resolved; "
            f"{sequence_defined:>2} unresolved sequence-defined; "
            f"{priority:>2} priority; "
            f"{clean_noninterface:>2} clean non-interface priority"
        )

    print_section(
        "OUTPUTS"
    )

    print()

    print(
        HUMAN_OUTPUT
    )

    print(
        RHESUS_OUTPUT
    )

    print(
        PIGTAIL_OUTPUT
    )

    print(
        COMBINED_OUTPUT
    )

    print()

    print(
        "NOTE: Step 2T uses center-defined local neighborhoods rather than "
        "connected-component chaining."
    )

    print(
        "Every resolved residue in a core must lie directly within "
        f"{LOCAL_CORE_CUTOFF_A:.1f} A of the same center residue."
    )

    print(
        "Unresolved N-terminal regions remain sequence-defined and are not "
        "assigned AlphaFold-derived fixed geometry."
    )

    print(
        "Human, rhesus, and pigtail remain independent within-species "
        "NKG2A-vs-NKG2C analyses."
    )

    print(
        "Rhesus and pigtail structural geometry remains homologous-position "
        "evidence from human 3CDG, not direct macaque structural evidence."
    )


if __name__ == "__main__":
    main()
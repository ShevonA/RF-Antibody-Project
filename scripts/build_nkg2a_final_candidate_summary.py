#!/usr/bin/env python3
"""
STEP 2V - FINAL WITHIN-SPECIES NKG2A-SPECIFIC CANDIDATE SUMMARY

Purpose
-------
Convert the Step 2U ranked NKG2A-specific epitope hypotheses into a concise,
species-specific final candidate selection and evidence summary.

Analysis remains strictly Axis 1:

    human:
        human NKG2A vs human NKG2C

    rhesus macaque:
        rhesus NKG2A vs BOTH rhesus NKG2C isoforms

    pigtail macaque:
        pigtail NKG2A vs pigtail NKG2C

Cross-species conservation and antibody cross-reactivity are NOT used as
selection criteria in this step.

This step does not generate new structural geometry. It summarizes the
evidence already established in Steps 2S-2U.

Evidence classes
----------------
Resolved structural hypotheses:
    tier_1_resolved_exposed_noninterface
    tier_2_resolved_mixed_accessibility
    tier_3_resolved_interface_associated
    tier_4_resolved_low_accessibility

Sequence-defined unresolved hypotheses:
    tier_U_sequence_defined_unresolved

Final selection categories
--------------------------
primary_resolved_candidate
    Tier 1 resolved exposed/non-interface hypothesis.

secondary_resolved_candidate
    Tier 2 resolved hypothesis with useful but less favorable accessibility.

primary_sequence_defined_candidate
    Strong sequence-discriminatory unresolved hypothesis retained because
    absence of experimental coordinates represents uncertainty rather than
    evidence of burial or interface involvement.

interface_associated_candidate
    Tier 3 resolved hypothesis with receptor/ligand-interface involvement.

low_accessibility_candidate
    Tier 4 resolved hypothesis with burial or poor accessibility.

Important interpretation rule
-----------------------------
A sequence-defined unresolved candidate may receive a high experimental
selection priority, but it MUST remain explicitly labeled unresolved.
No RSA, interface state, diameter, or fixed epitope geometry is inferred
for that candidate.

Inputs
------
results/tables/structure/
    nkg2a_within_species_ranked_epitope_hypotheses.tsv

Outputs
-------
results/tables/structure/
    human_NKG2A_final_candidate_summary.tsv
    rhesus_NKG2A_final_candidate_summary.tsv
    pigtail_NKG2A_final_candidate_summary.tsv
    nkg2a_within_species_final_candidate_summary.tsv
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


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

INPUT_FILE = (
    STRUCTURE_DIR
    / "nkg2a_within_species_ranked_epitope_hypotheses.tsv"
)

OUTPUT_FILES = {
    "human": (
        STRUCTURE_DIR
        / "human_NKG2A_final_candidate_summary.tsv"
    ),
    "rhesus": (
        STRUCTURE_DIR
        / "rhesus_NKG2A_final_candidate_summary.tsv"
    ),
    "pigtail": (
        STRUCTURE_DIR
        / "pigtail_NKG2A_final_candidate_summary.tsv"
    ),
}

COMBINED_OUTPUT = (
    STRUCTURE_DIR
    / "nkg2a_within_species_final_candidate_summary.tsv"
)

SPECIES = (
    "human",
    "rhesus",
    "pigtail",
)

SPECIES_ORDER = {
    "human": 0,
    "rhesus": 1,
    "pigtail": 2,
}


# =============================================================================
# BASIC HELPERS
# =============================================================================

def clean(value) -> str:
    if value is None:
        return ""

    return str(value).strip()


def lower(value) -> str:
    return clean(value).lower()


def safe_int(
    value,
    default: int = 0,
) -> int:
    text = clean(value)

    if not text:
        return default

    try:
        return int(float(text))

    except (TypeError, ValueError):
        return default


def safe_float(
    value,
) -> Optional[float]:
    text = clean(value)

    if not text:
        return None

    try:
        value_float = float(text)

    except (TypeError, ValueError):
        return None

    if not math.isfinite(value_float):
        return None

    return value_float


def fmt_float(
    value: Optional[float],
    digits: int = 4,
) -> str:
    if value is None:
        return ""

    return f"{value:.{digits}f}"


def read_tsv(
    path: Path,
) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(
            f"Required Step 2U input file not found:\n{path}"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise RuntimeError(
                f"No TSV header found in:\n{path}"
            )

        return [
            {
                clean(key): clean(value)
                for key, value in row.items()
            }
            for row in reader
        ]


def write_tsv(
    path: Path,
    rows: Sequence[Dict[str, object]],
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


def banner(
    text: str,
) -> None:
    print("=" * 78)
    print(text)
    print("=" * 78)


def section(
    text: str,
) -> None:
    print()
    print("=" * 78)
    print(text)
    print("=" * 78)
    print()


def unique_join(
    values: Iterable[str],
    separator: str = ";",
) -> str:
    output = []
    seen = set()

    for value in values:

        value = clean(value)

        if not value:
            continue

        if value in seen:
            continue

        seen.add(value)
        output.append(value)

    return separator.join(
        output
    )


# =============================================================================
# VALIDATION
# =============================================================================

ALLOWED_EVIDENCE_TIERS = {
    "tier_1_resolved_exposed_noninterface",
    "tier_2_resolved_mixed_accessibility",
    "tier_3_resolved_interface_associated",
    "tier_4_resolved_low_accessibility",
    "tier_U_sequence_defined_unresolved",
}

ALLOWED_HYPOTHESIS_TYPES = {
    "resolved_structural_hypothesis",
    "sequence_defined_unresolved_hypothesis",
}


def validate_input(
    rows: Sequence[Dict[str, str]],
) -> None:

    if not rows:
        raise RuntimeError(
            "Step 2U ranked hypothesis table contains no rows."
        )

    required_columns = {
        "species",
        "hypothesis_rank",
        "hypothesis_id",
        "hypothesis_type",
        "parent_region_id",
        "representative_residue_labels",
        "evidence_tier",
        "discriminator_count",
    }

    missing_columns = sorted(
        required_columns
        - set(rows[0])
    )

    if missing_columns:
        raise RuntimeError(
            "Step 2U input is missing required column(s):\n  "
            + "\n  ".join(
                missing_columns
            )
        )

    errors = []

    seen_ids = set()

    for row in rows:

        species = lower(
            row.get("species")
        )

        hypothesis_id = clean(
            row.get("hypothesis_id")
        )

        hypothesis_type = clean(
            row.get("hypothesis_type")
        )

        tier = clean(
            row.get("evidence_tier")
        )

        if species not in SPECIES:
            errors.append(
                f"{hypothesis_id or '?'}: "
                f"unsupported species {species!r}"
            )

        if not hypothesis_id:
            errors.append(
                f"{species or '?'}: missing hypothesis_id"
            )

        if hypothesis_id in seen_ids:
            errors.append(
                f"duplicate hypothesis_id: {hypothesis_id}"
            )

        seen_ids.add(
            hypothesis_id
        )

        if (
            hypothesis_type
            not in ALLOWED_HYPOTHESIS_TYPES
        ):
            errors.append(
                f"{species} {hypothesis_id}: "
                f"unsupported hypothesis_type "
                f"{hypothesis_type!r}"
            )

        if tier not in ALLOWED_EVIDENCE_TIERS:
            errors.append(
                f"{species} {hypothesis_id}: "
                f"unsupported evidence_tier {tier!r}"
            )

        if (
            hypothesis_type
            == "sequence_defined_unresolved_hypothesis"
        ):
            if (
                tier
                != "tier_U_sequence_defined_unresolved"
            ):
                errors.append(
                    f"{species} {hypothesis_id}: "
                    "unresolved hypothesis has structural tier"
                )

            # These must remain blank for unresolved hypotheses.
            for field in (
                "mean_complex_rsa",
                "diameter_A",
                "interface_residue_count",
                "buried_residue_count",
            ):
                if clean(
                    row.get(field)
                ):
                    errors.append(
                        f"{species} {hypothesis_id}: "
                        f"unresolved hypothesis unexpectedly "
                        f"contains {field}"
                    )

        if (
            hypothesis_type
            == "resolved_structural_hypothesis"
        ):
            if (
                tier
                == "tier_U_sequence_defined_unresolved"
            ):
                errors.append(
                    f"{species} {hypothesis_id}: "
                    "resolved structural hypothesis has "
                    "unresolved evidence tier"
                )

    missing_species = [
        species
        for species in SPECIES
        if not any(
            lower(row.get("species"))
            == species
            for row in rows
        )
    ]

    if missing_species:
        errors.append(
            "Missing species from Step 2U table: "
            + ", ".join(
                missing_species
            )
        )

    if errors:
        raise RuntimeError(
            "Step 2V input validation failed:\n\n  "
            + "\n  ".join(
                errors
            )
        )


# =============================================================================
# FINAL SELECTION CLASSIFICATION
# =============================================================================

def selection_category(
    row: Dict[str, str],
) -> str:

    tier = clean(
        row.get("evidence_tier")
    )

    if (
        tier
        == "tier_1_resolved_exposed_noninterface"
    ):
        return "primary_resolved_candidate"

    if (
        tier
        == "tier_2_resolved_mixed_accessibility"
    ):
        return "secondary_resolved_candidate"

    if (
        tier
        == "tier_U_sequence_defined_unresolved"
    ):
        return "primary_sequence_defined_candidate"

    if (
        tier
        == "tier_3_resolved_interface_associated"
    ):
        return "interface_associated_candidate"

    if (
        tier
        == "tier_4_resolved_low_accessibility"
    ):
        return "low_accessibility_candidate"

    raise ValueError(
        f"Unsupported evidence tier: {tier}"
    )


SELECTION_ORDER = {
    "primary_resolved_candidate": 1,
    "primary_sequence_defined_candidate": 2,
    "secondary_resolved_candidate": 3,
    "interface_associated_candidate": 4,
    "low_accessibility_candidate": 5,
}


# =============================================================================
# STRUCTURAL CONFIDENCE
# =============================================================================

def structural_confidence_class(
    row: Dict[str, str],
) -> str:

    hypothesis_type = clean(
        row.get("hypothesis_type")
    )

    species = lower(
        row.get("species")
    )

    if (
        hypothesis_type
        == "sequence_defined_unresolved_hypothesis"
    ):
        return "unresolved_no_fixed_geometry"

    if species == "human":
        return "direct_human_experimental_structure"

    if species in {
        "rhesus",
        "pigtail",
    }:
        return "human_structure_homology_projection"

    return "unknown"


# =============================================================================
# TARGET ATTRACTIVENESS
# =============================================================================

def target_attractiveness_class(
    row: Dict[str, str],
) -> str:

    category = selection_category(
        row
    )

    if category == "primary_resolved_candidate":
        return "high"

    if category == "primary_sequence_defined_candidate":
        return "high_but_structurally_unresolved"

    if category == "secondary_resolved_candidate":
        return "moderate"

    if category == "interface_associated_candidate":
        return "mechanism_dependent"

    if category == "low_accessibility_candidate":
        return "low"

    return "unclassified"


# =============================================================================
# EXPERIMENTAL PRIORITY
# =============================================================================

def experimental_priority_class(
    row: Dict[str, str],
) -> str:

    category = selection_category(
        row
    )

    if category in {
        "primary_resolved_candidate",
        "primary_sequence_defined_candidate",
    }:
        return "priority_1"

    if category == "secondary_resolved_candidate":
        return "priority_2"

    if category == "interface_associated_candidate":
        return "priority_3"

    if category == "low_accessibility_candidate":
        return "priority_4"

    return "unclassified"


# =============================================================================
# FLAGS / CAVEATS
# =============================================================================

def glycosylation_flag(
    row: Dict[str, str],
) -> str:

    context = lower(
        row.get("glycosylation_context")
    )

    if not context:
        return "unknown"

    if "overlapping" in context:
        return "overlapping"

    if "nearby" in context:
        return "nearby"

    if (
        "none_detected" in context
        or context == "none"
    ):
        return "none_detected"

    return "reported"


def build_strength_summary(
    row: Dict[str, str],
) -> str:

    strengths = []

    discriminator_n = safe_int(
        row.get("discriminator_count")
    )

    if discriminator_n > 1:
        strengths.append(
            f"{discriminator_n}_within_species_discriminators"
        )

    elif discriminator_n == 1:
        strengths.append(
            "within_species_discriminator"
        )

    exposed = safe_int(
        row.get("exposed_noninterface_count")
    )

    if exposed > 1:
        strengths.append(
            f"{exposed}_exposed_noninterface_residues"
        )

    elif exposed == 1:
        strengths.append(
            "exposed_noninterface_residue"
        )

    if (
        clean(row.get("hypothesis_type"))
        == "sequence_defined_unresolved_hypothesis"
    ):
        strengths.append(
            "strong_sequence_discrimination_retained_despite_missing_structure"
        )

    rsa = safe_float(
        row.get("mean_complex_rsa")
    )

    if rsa is not None:

        if rsa >= 0.50:
            strengths.append(
                "high_mean_accessibility"
            )

        elif rsa >= 0.25:
            strengths.append(
                "exposed_mean_accessibility"
            )

    diameter = safe_float(
        row.get("diameter_A")
    )

    if diameter is not None:

        if diameter <= 8.0:
            strengths.append(
                "compact_structural_neighborhood"
            )

    interface = safe_int(
        row.get("interface_residue_count")
    )

    if (
        clean(row.get("hypothesis_type"))
        == "resolved_structural_hypothesis"
        and interface == 0
    ):
        strengths.append(
            "no_detected_interface_residue"
        )

    glyco = glycosylation_flag(
        row
    )

    if glyco == "none_detected":
        strengths.append(
            "no_detected_nearby_canonical_N_glycosylation_sequon"
        )

    # Retain Step 2U strengths as provenance.
    prior = clean(
        row.get("hypothesis_strengths")
    )

    if prior:
        strengths.append(
            prior
        )

    return unique_join(
        strengths
    )


def build_limitation_summary(
    row: Dict[str, str],
) -> str:

    limitations = []

    hypothesis_type = clean(
        row.get("hypothesis_type")
    )

    species = lower(
        row.get("species")
    )

    if (
        hypothesis_type
        == "sequence_defined_unresolved_hypothesis"
    ):
        limitations.extend(
            [
                "no_experimental_structural_coordinates",
                "surface_accessibility_not_established",
                "fixed_3D_geometry_not_assigned",
            ]
        )

        if species == "human":
            limitations.append(
                "AlphaFold_94_112_geometry_not_used_for_fixed_epitope_interpretation"
            )

        else:
            limitations.append(
                "no_direct_macaque_structural_geometry"
            )

    interface = safe_int(
        row.get("interface_residue_count")
    )

    if interface > 0:
        limitations.append(
            f"{interface}_interface_associated_residue"
            + (
                "s"
                if interface != 1
                else ""
            )
        )

    buried = safe_int(
        row.get("buried_residue_count")
    )

    if buried > 0:
        limitations.append(
            f"{buried}_buried_residue"
            + (
                "s"
                if buried != 1
                else ""
            )
        )

    rsa = safe_float(
        row.get("mean_complex_rsa")
    )

    if rsa is not None:

        if rsa < 0.10:
            limitations.append(
                "very_low_mean_accessibility"
            )

        elif rsa < 0.25:
            limitations.append(
                "limited_mean_accessibility"
            )

    glyco = glycosylation_flag(
        row
    )

    if glyco == "overlapping":
        limitations.append(
            "overlapping_canonical_N_glycosylation_sequon"
        )

    elif glyco == "nearby":
        limitations.append(
            "nearby_canonical_N_glycosylation_sequon"
        )

    if (
        hypothesis_type
        == "resolved_structural_hypothesis"
        and species in {
            "rhesus",
            "pigtail",
        }
    ):
        limitations.append(
            "structural_evidence_projected_from_human_3CDG"
        )

        limitations.append(
            "not_direct_macaque_structure"
        )

    prior = clean(
        row.get("hypothesis_limitations")
    )

    if prior:
        limitations.append(
            prior
        )

    return unique_join(
        limitations
    )


# =============================================================================
# FINAL RECOMMENDATION TEXT
# =============================================================================

def recommendation_text(
    row: Dict[str, str],
) -> str:

    category = selection_category(
        row
    )

    labels = clean(
        row.get(
            "representative_residue_labels"
        )
    )

    if category == "primary_resolved_candidate":
        return (
            f"{labels} is retained as a primary resolved within-species "
            "NKG2A-specific candidate because it combines NKG2A-vs-NKG2C "
            "sequence discrimination with favorable exposed non-interface "
            "structural evidence."
        )

    if category == "primary_sequence_defined_candidate":
        return (
            f"{labels} is retained as a primary sequence-defined "
            "NKG2A-specific candidate because of strong within-species "
            "NKG2A-vs-NKG2C sequence discrimination. Structural accessibility "
            "remains unresolved and no fixed 3D epitope geometry is assigned."
        )

    if category == "secondary_resolved_candidate":
        return (
            f"{labels} is retained as a secondary resolved NKG2A-specific "
            "candidate. Sequence discrimination is supported, but accessibility "
            "is less favorable than the primary resolved candidates."
        )

    if category == "interface_associated_candidate":
        return (
            f"{labels} is retained as an interface-associated NKG2A-specific "
            "candidate. Its usefulness may depend on the intended antibody "
            "mechanism because the discriminatory surface participates in "
            "receptor/ligand contacts."
        )

    if category == "low_accessibility_candidate":
        return (
            f"{labels} remains sequence-discriminatory but is a low-priority "
            "antibody-accessible candidate because the available structure "
            "indicates poor accessibility or burial."
        )

    return (
        f"{labels} remains an unclassified within-species NKG2A-specific "
        "candidate."
    )


# =============================================================================
# FINAL ROW CONSTRUCTION
# =============================================================================

def build_final_row(
    row: Dict[str, str],
) -> Dict[str, object]:

    species = lower(
        row.get("species")
    )

    category = selection_category(
        row
    )

    structural_confidence = (
        structural_confidence_class(
            row
        )
    )

    attractiveness = (
        target_attractiveness_class(
            row
        )
    )

    experimental_priority = (
        experimental_priority_class(
            row
        )
    )

    output = {
        "species":
            species,

        "step2u_hypothesis_rank":
            clean(
                row.get(
                    "hypothesis_rank"
                )
            ),

        "hypothesis_id":
            clean(
                row.get(
                    "hypothesis_id"
                )
            ),

        "parent_region_id":
            clean(
                row.get(
                    "parent_region_id"
                )
            ),

        "representative_core":
            clean(
                row.get(
                    "representative_core"
                )
            ),

        "candidate_residues":
            clean(
                row.get(
                    "representative_residue_labels"
                )
            ),

        "within_species_comparisons":
            clean(
                row.get(
                    "within_species_comparisons"
                )
            ),

        "hypothesis_type":
            clean(
                row.get(
                    "hypothesis_type"
                )
            ),

        "evidence_tier":
            clean(
                row.get(
                    "evidence_tier"
                )
            ),

        "final_selection_category":
            category,

        "experimental_priority":
            experimental_priority,

        "target_attractiveness":
            attractiveness,

        "structural_confidence":
            structural_confidence,

        "discriminator_count":
            clean(
                row.get(
                    "discriminator_count"
                )
            ),

        "exposed_noninterface_count":
            clean(
                row.get(
                    "exposed_noninterface_count"
                )
            ),

        "interface_residue_count":
            clean(
                row.get(
                    "interface_residue_count"
                )
            ),

        "buried_residue_count":
            clean(
                row.get(
                    "buried_residue_count"
                )
            ),

        "mean_complex_rsa":
            clean(
                row.get(
                    "mean_complex_rsa"
                )
            ),

        "diameter_A":
            clean(
                row.get(
                    "diameter_A"
                )
            ),

        "glycosylation_context":
            clean(
                row.get(
                    "glycosylation_context"
                )
            ),

        "glycosylation_flag":
            glycosylation_flag(
                row
            ),

        "strength_summary":
            build_strength_summary(
                row
            ),

        "limitation_summary":
            build_limitation_summary(
                row
            ),

        "final_interpretation":
            recommendation_text(
                row
            ),

        "step2u_interpretation":
            clean(
                row.get(
                    "interpretation"
                )
            ),
    }

    return output


# =============================================================================
# FINAL SORTING
# =============================================================================

def final_sort_key(
    row: Dict[str, object],
) -> Tuple:

    species = clean(
        row.get("species")
    )

    category = clean(
        row.get(
            "final_selection_category"
        )
    )

    discriminator_count = safe_int(
        row.get(
            "discriminator_count"
        )
    )

    exposed = safe_int(
        row.get(
            "exposed_noninterface_count"
        )
    )

    rsa = safe_float(
        row.get(
            "mean_complex_rsa"
        )
    )

    if rsa is None:
        rsa_sort = -1.0

    else:
        rsa_sort = rsa

    step2u_rank = safe_int(
        row.get(
            "step2u_hypothesis_rank"
        ),
        default=999,
    )

    return (
        SPECIES_ORDER.get(
            species,
            99,
        ),
        SELECTION_ORDER.get(
            category,
            99,
        ),
        -exposed,
        -discriminator_count,
        -rsa_sort,
        step2u_rank,
        clean(
            row.get(
                "hypothesis_id"
            )
        ),
    )


# =============================================================================
# FINAL RANK WITHIN SPECIES
# =============================================================================

def assign_final_species_rank(
    rows: Sequence[Dict[str, object]],
) -> List[Dict[str, object]]:

    by_species = {
        species: []
        for species in SPECIES
    }

    for row in rows:

        species = clean(
            row.get("species")
        )

        by_species[
            species
        ].append(
            dict(row)
        )

    output = []

    for species in SPECIES:

        ordered = sorted(
            by_species[
                species
            ],
            key=final_sort_key,
        )

        for index, row in enumerate(
            ordered,
            start=1,
        ):

            row[
                "final_species_rank"
            ] = index

            output.append(
                row
            )

    return output


# =============================================================================
# OUTPUT FIELDS
# =============================================================================

OUTPUT_FIELDS = [
    "species",
    "final_species_rank",
    "step2u_hypothesis_rank",
    "hypothesis_id",
    "parent_region_id",
    "representative_core",
    "candidate_residues",
    "within_species_comparisons",
    "hypothesis_type",
    "evidence_tier",
    "final_selection_category",
    "experimental_priority",
    "target_attractiveness",
    "structural_confidence",
    "discriminator_count",
    "exposed_noninterface_count",
    "interface_residue_count",
    "buried_residue_count",
    "mean_complex_rsa",
    "diameter_A",
    "glycosylation_context",
    "glycosylation_flag",
    "strength_summary",
    "limitation_summary",
    "final_interpretation",
    "step2u_interpretation",
]


# =============================================================================
# CONSOLE REPORT
# =============================================================================

def print_species_summary(
    species: str,
    rows: Sequence[Dict[str, object]],
) -> None:

    section(
        f"{species.upper()} FINAL NKG2A-SPECIFIC CANDIDATES"
    )

    species_rows = [
        row
        for row in rows
        if clean(
            row.get(
                "species"
            )
        ) == species
    ]

    species_rows.sort(
        key=lambda row:
        safe_int(
            row.get(
                "final_species_rank"
            )
        )
    )

    for row in species_rows:

        print(
            f"{row['hypothesis_id']:<7} "
            f"final_rank={row['final_species_rank']:<2} "
            f"{row['candidate_residues']}"
        )

        print(
            f"  Selection:   "
            f"{row['final_selection_category']}"
        )

        print(
            f"  Experimental:"
            f" {row['experimental_priority']}"
        )

        print(
            f"  Attractiveness: "
            f"{row['target_attractiveness']}"
        )

        print(
            f"  Evidence:    "
            f"{row['evidence_tier']}"
        )

        print(
            f"  Structure:   "
            f"{row['structural_confidence']}"
        )

        print(
            f"  Discriminators: "
            f"{row['discriminator_count']}"
        )

        if clean(
            row.get(
                "mean_complex_rsa"
            )
        ):

            print(
                f"  Mean RSA:    "
                f"{row['mean_complex_rsa']}"
            )

        if clean(
            row.get(
                "diameter_A"
            )
        ):

            print(
                f"  Diameter:    "
                f"{row['diameter_A']} A"
            )

        print(
            f"  Glyco:       "
            f"{row['glycosylation_flag']}"
        )

        print(
            f"  Interpretation:"
            f" {row['final_interpretation']}"
        )

        print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    banner(
        "STEP 2V - FINAL WITHIN-SPECIES NKG2A-SPECIFIC CANDIDATE SUMMARY"
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
        "Step 2V adds no new structural geometry."
    )

    print(
        "Selection priority and structural confidence remain separate concepts."
    )

    # -------------------------------------------------------------------------
    # Load and validate Step 2U.
    # -------------------------------------------------------------------------

    rows = read_tsv(
        INPUT_FILE
    )

    validate_input(
        rows
    )

    print()

    print(
        f"Step 2U hypotheses loaded: "
        f"{len(rows)}"
    )

    for species in SPECIES:

        count = sum(
            1
            for row in rows
            if lower(
                row.get(
                    "species"
                )
            ) == species
        )

        print(
            f"  {species:<8} {count}"
        )

    # -------------------------------------------------------------------------
    # Construct final evidence summaries.
    # -------------------------------------------------------------------------

    final_rows = [
        build_final_row(
            row
        )
        for row in rows
    ]

    final_rows = (
        assign_final_species_rank(
            final_rows
        )
    )

    # -------------------------------------------------------------------------
    # Write species outputs.
    # -------------------------------------------------------------------------

    for species in SPECIES:

        species_rows = [
            row
            for row in final_rows
            if row[
                "species"
            ] == species
        ]

        species_rows.sort(
            key=lambda row:
            safe_int(
                row[
                    "final_species_rank"
                ]
            )
        )

        write_tsv(
            OUTPUT_FILES[
                species
            ],
            species_rows,
            OUTPUT_FIELDS,
        )

    # -------------------------------------------------------------------------
    # Combined output.
    # -------------------------------------------------------------------------

    combined = sorted(
        final_rows,
        key=lambda row: (
            SPECIES_ORDER[
                row[
                    "species"
                ]
            ],
            safe_int(
                row[
                    "final_species_rank"
                ]
            ),
        ),
    )

    write_tsv(
        COMBINED_OUTPUT,
        combined,
        OUTPUT_FIELDS,
    )

    # -------------------------------------------------------------------------
    # Console report.
    # -------------------------------------------------------------------------

    for species in SPECIES:

        print_species_summary(
            species,
            combined,
        )

    # -------------------------------------------------------------------------
    # Final Axis-1 summary.
    # -------------------------------------------------------------------------

    section(
        "STEP 2V AXIS-1 SUMMARY"
    )

    for species in SPECIES:

        species_rows = [
            row
            for row in combined
            if row[
                "species"
            ] == species
        ]

        primary_resolved = sum(
            1
            for row in species_rows
            if row[
                "final_selection_category"
            ]
            == "primary_resolved_candidate"
        )

        primary_sequence = sum(
            1
            for row in species_rows
            if row[
                "final_selection_category"
            ]
            == "primary_sequence_defined_candidate"
        )

        secondary = sum(
            1
            for row in species_rows
            if row[
                "final_selection_category"
            ]
            == "secondary_resolved_candidate"
        )

        interface = sum(
            1
            for row in species_rows
            if row[
                "final_selection_category"
            ]
            == "interface_associated_candidate"
        )

        low = sum(
            1
            for row in species_rows
            if row[
                "final_selection_category"
            ]
            == "low_accessibility_candidate"
        )

        print(
            f"{species:<8} "
            f"primary_resolved={primary_resolved}  "
            f"primary_sequence={primary_sequence}  "
            f"secondary={secondary}  "
            f"interface={interface}  "
            f"low_accessibility={low}"
        )

    section(
        "OUTPUTS"
    )

    print(
        OUTPUT_FILES[
            "human"
        ]
    )

    print(
        OUTPUT_FILES[
            "rhesus"
        ]
    )

    print(
        OUTPUT_FILES[
            "pigtail"
        ]
    )

    print(
        COMBINED_OUTPUT
    )

    print()

    print(
        "NOTE: Step 2V is the final evidence-summary step for Axis 1."
    )

    print(
        "Human, rhesus, and pigtail remain independent within-species "
        "NKG2A-vs-NKG2C analyses."
    )

    print(
        "A primary_sequence_defined_candidate is not structurally validated; "
        "its priority reflects sequence discrimination and experimental "
        "interest while structural uncertainty remains explicit."
    )

    print(
        "Resolved rhesus and pigtail structural evidence remains a homologous-"
        "position projection from human NKG2A 3CDG rather than direct macaque "
        "structural measurement."
    )

    print(
        "Cross-species antibody reactivity has not yet been introduced."
    )

    print(
        "All final candidates remain computational screening hypotheses and "
        "are not experimentally validated antibody epitopes."
    )


if __name__ == "__main__":
    main()
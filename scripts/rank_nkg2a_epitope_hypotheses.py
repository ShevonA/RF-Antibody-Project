#!/usr/bin/env python3
"""
STEP 2U
Rank within-species NKG2A-specific epitope hypotheses.

PURPOSE
-------
Integrate Step 2S candidate epitope regions with Step 2T local epitope
cores and produce a final within-species ranking of NKG2A-specific
antibody epitope hypotheses.

ANALYSIS AXIS
-------------
Human:
    human NKG2A vs human NKG2C

Rhesus macaque:
    rhesus NKG2A vs BOTH rhesus NKG2C isoforms

Pigtail macaque:
    pigtail NKG2A vs pigtail NKG2C

IMPORTANT
---------
Cross-species conservation/reactivity is NOT used as a ranking criterion.

Step 2T contains TWO conceptually different evidence classes:

1. Structural cores
       resolved_local_core
       resolved_local_site

2. Sequence-defined unresolved cores
       sequence_defined_unresolved

Sequence-defined unresolved cores are legitimate Step 2T outputs.
They MUST NOT be subjected to structural validation, structural
accessibility ranking, interface ranking, structural diameter ranking,
or fixed-geometry interpretation.

They are instead carried forward as:
       sequence_defined_unresolved_hypothesis

This distinction prevents experimentally unresolved N-terminal regions
from being accidentally promoted to resolved structural hypotheses.

INPUTS
------
results/tables/structure/
    nkg2a_within_species_candidate_epitope_regions.tsv
    nkg2a_within_species_local_epitope_cores.tsv

OUTPUTS
-------
results/tables/structure/
    human_NKG2A_ranked_epitope_hypotheses.tsv
    rhesus_NKG2A_ranked_epitope_hypotheses.tsv
    pigtail_NKG2A_ranked_epitope_hypotheses.tsv
    nkg2a_within_species_ranked_epitope_hypotheses.tsv

DESIGN PRINCIPLES
-----------------
* Within-species NKG2A-vs-NKG2C discrimination is mandatory.
* Structural hypotheses and unresolved sequence hypotheses are kept
  conceptually separate.
* Structural ranking uses categorical evidence tiers rather than an
  arbitrary additive numerical score.
* Redundant structural cores from the same Step 2S region are
  consolidated.
* Sequence-defined unresolved regions are represented once and are
  not duplicated as structural hypotheses.
* Human experimental structure is the direct structural evidence.
* Rhesus and pigtail structural annotations are homologous-position
  projections from human NKG2A 3CDG and are labeled accordingly.
* AlphaFold geometry for residues 94-112 is NOT used as fixed epitope
  geometry.
"""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# =============================================================================
# CONFIGURATION
# =============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

STRUCTURE_TABLE_DIR = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "structure"
)

REGION_FILE = (
    STRUCTURE_TABLE_DIR
    / "nkg2a_within_species_candidate_epitope_regions.tsv"
)

CORE_FILE = (
    STRUCTURE_TABLE_DIR
    / "nkg2a_within_species_local_epitope_cores.tsv"
)

OUTPUT_COMBINED = (
    STRUCTURE_TABLE_DIR
    / "nkg2a_within_species_ranked_epitope_hypotheses.tsv"
)

OUTPUT_BY_SPECIES = {
    "human": (
        STRUCTURE_TABLE_DIR
        / "human_NKG2A_ranked_epitope_hypotheses.tsv"
    ),
    "rhesus": (
        STRUCTURE_TABLE_DIR
        / "rhesus_NKG2A_ranked_epitope_hypotheses.tsv"
    ),
    "pigtail": (
        STRUCTURE_TABLE_DIR
        / "pigtail_NKG2A_ranked_epitope_hypotheses.tsv"
    ),
}

SPECIES_ORDER = {
    "human": 0,
    "rhesus": 1,
    "pigtail": 2,
}

HYPOTHESIS_PREFIX = {
    "human": "HHYP",
    "rhesus": "RHYP",
    "pigtail": "PHYP",
}

STRUCTURAL_CORE_TYPES = {
    "resolved_local_core",
    "resolved_local_site",
}

UNRESOLVED_CORE_TYPES = {
    "sequence_defined_unresolved",
}

ALLOWED_CORE_TYPES = STRUCTURAL_CORE_TYPES | UNRESOLVED_CORE_TYPES

STRUCTURAL_REGION_TYPES = {
    "resolved_structural_patch",
    "resolved_structural_site",
}

UNRESOLVED_REGION_TYPES = {
    "sequence_defined_unresolved",
}

ALLOWED_REGION_TYPES = STRUCTURAL_REGION_TYPES | UNRESOLVED_REGION_TYPES


# =============================================================================
# OUTPUT COLUMN ORDER
# =============================================================================

OUTPUT_FIELDS = [
    "species",
    "hypothesis_rank",
    "hypothesis_id",
    "hypothesis_type",
    "parent_region_id",
    "representative_core",
    "representative_residue_labels",
    "within_species_comparisons",
    "evidence_tier",
    "discriminator_count",
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
    "diameter_A",
    "glycosylation_context",
    "glycosylation_sequons",
    "geometry_evidence_basis",
    "alternative_overlapping_cores",
    "redundant_core_count",
    "hypothesis_strengths",
    "hypothesis_limitations",
    "interpretation",
]


# =============================================================================
# GENERAL HELPERS
# =============================================================================

def banner(title: str) -> None:
    print("=" * 78)
    print(title)
    print("=" * 78)


def subsection(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)
    print()


def clean(value) -> str:
    if value is None:
        return ""
    return str(value).strip()


def lower(value) -> str:
    return clean(value).lower()


def safe_int(value, default: int = 0) -> int:
    text = clean(value)
    if not text:
        return default

    try:
        return int(float(text))
    except (TypeError, ValueError):
        return default


def safe_float(value) -> Optional[float]:
    text = clean(value)

    if not text:
        return None

    try:
        value_float = float(text)
    except (TypeError, ValueError):
        return None

    if math.isnan(value_float):
        return None

    return value_float


def fmt_float(
    value: Optional[float],
    decimals: int = 4,
) -> str:
    if value is None:
        return ""

    return f"{value:.{decimals}f}"


def split_labels(value: str) -> List[str]:
    text = clean(value)

    if not text:
        return []

    return [
        item.strip()
        for item in text.split(",")
        if item.strip()
    ]


def split_semicolon(value: str) -> List[str]:
    text = clean(value)

    if not text:
        return []

    return [
        item.strip()
        for item in text.split(";")
        if item.strip()
    ]


def unique_preserve_order(values: Iterable[str]) -> List[str]:
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

    return output


def join_semicolon(values: Iterable[str]) -> str:
    return ";".join(unique_preserve_order(values))


def require_file(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"{description} not found:\n  {path}"
        )


def read_tsv(path: Path) -> List[Dict[str, str]]:
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
                f"No TSV header found in:\n  {path}"
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
    rows: Sequence[Dict[str, str]],
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
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    field: clean(row.get(field, ""))
                    for field in OUTPUT_FIELDS
                }
            )


def get_first(
    row: Dict[str, str],
    *names: str,
    default: str = "",
) -> str:
    for name in names:
        if name in row:
            value = clean(row.get(name))

            if value:
                return value

    return default


# =============================================================================
# VALIDATION
# =============================================================================

def validate_region_table(
    regions: Sequence[Dict[str, str]],
) -> None:
    if not regions:
        raise RuntimeError(
            "Step 2S region table contains no rows."
        )

    required = {
        "species",
        "region_id",
        "region_type",
        "residue_labels",
    }

    missing_messages = []

    for index, row in enumerate(regions, start=2):
        missing = [
            field
            for field in required
            if field not in row
        ]

        if missing:
            missing_messages.append(
                f"row {index}: missing column(s): "
                + ", ".join(sorted(missing))
            )

    if missing_messages:
        raise RuntimeError(
            "Step 2S table validation failed:\n\n  "
            + "\n  ".join(missing_messages)
        )

    errors = []

    seen = set()

    for row in regions:
        species = lower(row.get("species"))
        region_id = clean(row.get("region_id"))
        region_type = lower(row.get("region_type"))

        if species not in SPECIES_ORDER:
            errors.append(
                f"{region_id or '?'}: unsupported species "
                f"{species!r}"
            )

        if not region_id:
            errors.append(
                f"{species or '?'}: blank region_id"
            )

        key = (species, region_id)

        if key in seen:
            errors.append(
                f"{species} {region_id}: duplicate Step 2S region"
            )

        seen.add(key)

        if region_type not in ALLOWED_REGION_TYPES:
            errors.append(
                f"{species} {region_id}: unsupported region_type "
                f"{region_type!r}"
            )

    if errors:
        raise RuntimeError(
            "Step 2S region validation failed:\n\n  "
            + "\n  ".join(errors)
        )


def validate_core_table(
    cores: Sequence[Dict[str, str]],
) -> None:
    if not cores:
        raise RuntimeError(
            "Step 2T core table contains no rows."
        )

    required = {
        "species",
        "core_id",
        "source_region_id",
        "core_type",
        "residue_labels",
    }

    errors = []

    seen = set()

    for row in cores:
        species = lower(row.get("species"))
        core_id = clean(row.get("core_id"))
        region_id = clean(row.get("source_region_id"))
        core_type = lower(row.get("core_type"))

        if species not in SPECIES_ORDER:
            errors.append(
                f"{core_id or '?'}: unsupported species "
                f"{species!r}"
            )

        if not core_id:
            errors.append(
                f"{species or '?'}: blank core_id"
            )

        if not region_id:
            errors.append(
                f"{species} {core_id}: blank source_region_id"
            )

        key = (species, core_id)

        if key in seen:
            errors.append(
                f"{species} {core_id}: duplicate Step 2T core"
            )

        seen.add(key)

        if core_type not in ALLOWED_CORE_TYPES:
            errors.append(
                f"{species} {core_id}: unsupported core_type "
                f"{core_type!r}"
            )

    if errors:
        raise RuntimeError(
            "Step 2T core-table validation failed:\n\n  "
            + "\n  ".join(errors)
        )


def build_region_lookup(
    regions: Sequence[Dict[str, str]],
) -> Dict[Tuple[str, str], Dict[str, str]]:
    return {
        (
            lower(row["species"]),
            clean(row["region_id"]),
        ): row
        for row in regions
    }


def validate_core_region_relationships(
    cores: Sequence[Dict[str, str]],
    regions: Sequence[Dict[str, str]],
) -> None:
    """
    This is the critical Step 2U safeguard.

    Structural Step 2T cores must originate from structural Step 2S regions.

    Sequence-defined unresolved Step 2T cores must originate from
    sequence-defined unresolved Step 2S regions.

    Crucially, sequence_defined_unresolved is NOT considered an invalid
    Step 2T row. It is a valid evidence class and is routed separately.
    """

    region_lookup = build_region_lookup(regions)

    errors = []

    for core in cores:
        species = lower(core["species"])
        core_id = clean(core["core_id"])
        source_region_id = clean(core["source_region_id"])
        core_type = lower(core["core_type"])

        key = (
            species,
            source_region_id,
        )

        if key not in region_lookup:
            errors.append(
                f"{species} {core_id}: source region "
                f"{source_region_id!r} does not exist in Step 2S"
            )
            continue

        region = region_lookup[key]
        region_type = lower(region["region_type"])

        if core_type in STRUCTURAL_CORE_TYPES:
            if region_type not in STRUCTURAL_REGION_TYPES:
                errors.append(
                    f"{species} {core_id}: structural core type "
                    f"{core_type!r} points to Step 2S region "
                    f"{source_region_id} with non-structural type "
                    f"{region_type!r}"
                )

        elif core_type in UNRESOLVED_CORE_TYPES:
            if region_type not in UNRESOLVED_REGION_TYPES:
                errors.append(
                    f"{species} {core_id}: unresolved core type "
                    f"{core_type!r} points to Step 2S region "
                    f"{source_region_id} with structural type "
                    f"{region_type!r}"
                )

    if errors:
        raise RuntimeError(
            "Step 2T-to-Step 2S evidence-type validation failed:\n\n  "
            + "\n  ".join(errors)
        )


def validate_structural_cores(
    structural_cores: Sequence[Dict[str, str]],
    regions: Sequence[Dict[str, str]],
) -> None:
    """
    Structural validation is intentionally applied ONLY to resolved
    structural cores/sites.

    sequence_defined_unresolved rows never enter this function.
    """

    region_lookup = build_region_lookup(regions)

    errors = []

    for core in structural_cores:
        species = lower(core["species"])
        core_id = clean(core["core_id"])
        source_region_id = clean(core["source_region_id"])
        core_type = lower(core["core_type"])

        if core_type not in STRUCTURAL_CORE_TYPES:
            errors.append(
                f"{species} {core_id}: non-structural core entered "
                "structural validation"
            )
            continue

        key = (
            species,
            source_region_id,
        )

        region = region_lookup.get(key)

        if region is None:
            errors.append(
                f"{species} {core_id}: parent region "
                f"{source_region_id} missing"
            )
            continue

        region_type = lower(region["region_type"])

        if region_type not in STRUCTURAL_REGION_TYPES:
            errors.append(
                f"{species} {core_id}: parent region "
                f"{source_region_id} is {region_type!r}; "
                "structural cores may only originate from resolved "
                "structural Step 2S regions"
            )

        resolved = safe_int(
            core.get("resolved_residue_count"),
            default=0,
        )

        unresolved = safe_int(
            core.get("unresolved_residue_count"),
            default=0,
        )

        if resolved <= 0:
            errors.append(
                f"{species} {core_id}: structural core has "
                f"resolved_residue_count={resolved}"
            )

        if unresolved > 0:
            errors.append(
                f"{species} {core_id}: structural core unexpectedly "
                f"contains {unresolved} unresolved residue(s)"
            )

        diameter = safe_float(
            get_first(
                core,
                "known_pairwise_core_diameter_A",
                "diameter_A",
            )
        )

        if core_type == "resolved_local_core":
            if diameter is None:
                errors.append(
                    f"{species} {core_id}: resolved_local_core "
                    "has no structural diameter"
                )

    if errors:
        raise RuntimeError(
            "Step 2T structural-core validation failed:\n\n  "
            + "\n  ".join(errors)
        )


def validate_unresolved_cores(
    unresolved_cores: Sequence[Dict[str, str]],
    regions: Sequence[Dict[str, str]],
) -> None:
    region_lookup = build_region_lookup(regions)

    errors = []

    for core in unresolved_cores:
        species = lower(core["species"])
        core_id = clean(core["core_id"])
        source_region_id = clean(core["source_region_id"])
        core_type = lower(core["core_type"])

        if core_type not in UNRESOLVED_CORE_TYPES:
            errors.append(
                f"{species} {core_id}: structural core entered "
                "unresolved-core validation"
            )
            continue

        region = region_lookup.get(
            (
                species,
                source_region_id,
            )
        )

        if region is None:
            errors.append(
                f"{species} {core_id}: parent region "
                f"{source_region_id} missing"
            )
            continue

        region_type = lower(region["region_type"])

        if region_type not in UNRESOLVED_REGION_TYPES:
            errors.append(
                f"{species} {core_id}: unresolved core points to "
                f"Step 2S region type {region_type!r}"
            )

        resolved = safe_int(
            core.get("resolved_residue_count"),
            default=0,
        )

        if resolved != 0:
            errors.append(
                f"{species} {core_id}: sequence-defined unresolved "
                f"core unexpectedly reports {resolved} resolved "
                "residue(s)"
            )

    if errors:
        raise RuntimeError(
            "Step 2T unresolved-core validation failed:\n\n  "
            + "\n  ".join(errors)
        )


# =============================================================================
# EVIDENCE EXTRACTION
# =============================================================================

def discriminator_count(row: Dict[str, str]) -> int:
    comparisons = clean(
        row.get("within_species_comparisons")
    )

    if comparisons:
        return len(
            [
                item
                for item in comparisons.split("|")
                if item.strip()
            ]
        )

    return safe_int(
        row.get("residue_count"),
        default=len(
            split_labels(
                row.get("residue_labels", "")
            )
        ),
    )


def structural_counts(
    row: Dict[str, str],
) -> Dict[str, int]:
    return {
        "resolved": safe_int(
            row.get("resolved_residue_count")
        ),
        "unresolved": safe_int(
            row.get("unresolved_residue_count")
        ),
        "exposed_noninterface": safe_int(
            row.get("exposed_noninterface_count")
        ),
        "exposed_interface": safe_int(
            row.get("exposed_interface_count")
        ),
        "partially_exposed_noninterface": safe_int(
            row.get(
                "partially_exposed_noninterface_count"
            )
        ),
        "partially_exposed_interface": safe_int(
            row.get(
                "partially_exposed_interface_count"
            )
        ),
        "interface": safe_int(
            row.get("interface_residue_count")
        ),
        "CD94": safe_int(
            row.get("CD94_contact_residue_count")
        ),
        "HLA_E": safe_int(
            row.get("HLA_E_contact_residue_count")
        ),
        "peptide": safe_int(
            row.get("peptide_contact_residue_count")
        ),
        "B2M": safe_int(
            row.get("B2M_contact_residue_count")
        ),
        "buried": safe_int(
            row.get("buried_residue_count")
        ),
    }


def glyco_category(row: Dict[str, str]) -> str:
    value = lower(
        row.get("glycosylation_context")
    )

    if "overlapping" in value:
        return "overlapping"

    if "nearby" in value:
        return "nearby"

    if not value or value in {
        "none",
        "none_detected",
        "none_detected_near_core",
        "none_detected_near_region",
    }:
        return "none"

    return "other"


# =============================================================================
# STRUCTURAL EVIDENCE TIERS
# =============================================================================

STRUCTURAL_TIER_ORDER = {
    "tier_1_resolved_exposed_noninterface": 1,
    "tier_2_resolved_mixed_accessibility": 2,
    "tier_3_resolved_interface_associated": 3,
    "tier_4_resolved_low_accessibility": 4,
}


def assign_structural_tier(
    core: Dict[str, str],
) -> str:
    """
    Categorical structural ranking.

    Tier 1:
        At least one exposed non-interface discriminator and no
        detected interface residue.

    Tier 2:
        Accessible/partially accessible non-interface evidence but
        not a completely clean exposed site.

    Tier 3:
        Interface-associated structural candidate.

    Tier 4:
        Low-accessibility / buried candidate.
    """

    counts = structural_counts(core)

    exposed_noninterface = counts[
        "exposed_noninterface"
    ]

    partial_noninterface = counts[
        "partially_exposed_noninterface"
    ]

    interface = counts["interface"]
    buried = counts["buried"]

    resolved = max(
        counts["resolved"],
        safe_int(core.get("residue_count")),
    )

    if (
        exposed_noninterface > 0
        and interface == 0
    ):
        return "tier_1_resolved_exposed_noninterface"

    if (
        interface == 0
        and (
            partial_noninterface > 0
            or (
                exposed_noninterface > 0
                and buried > 0
            )
        )
    ):
        return "tier_2_resolved_mixed_accessibility"

    if interface > 0:
        return "tier_3_resolved_interface_associated"

    if (
        buried >= resolved
        or exposed_noninterface == 0
    ):
        return "tier_4_resolved_low_accessibility"

    return "tier_4_resolved_low_accessibility"


# =============================================================================
# STRUCTURAL CORE RANKING / CONSOLIDATION
# =============================================================================

def structural_core_sort_key(
    core: Dict[str, str],
) -> Tuple:
    tier = assign_structural_tier(core)
    tier_rank = STRUCTURAL_TIER_ORDER[tier]

    counts = structural_counts(core)

    rsa = safe_float(
        core.get("mean_complex_rsa")
    )

    if rsa is None:
        rsa = -1.0

    diameter = safe_float(
        get_first(
            core,
            "known_pairwise_core_diameter_A",
            "diameter_A",
        )
    )

    if diameter is None:
        diameter = 9999.0

    discriminators = discriminator_count(core)

    glyco = glyco_category(core)

    glyco_rank = {
        "none": 0,
        "nearby": 1,
        "other": 2,
        "overlapping": 3,
    }.get(
        glyco,
        2,
    )

    # Lower tuple = better.
    return (
        tier_rank,
        -counts["exposed_noninterface"],
        -counts["partially_exposed_noninterface"],
        counts["interface"],
        counts["buried"],
        glyco_rank,
        -discriminators,
        -rsa,
        diameter,
        clean(core.get("core_id")),
    )


def choose_representative_core(
    cores: Sequence[Dict[str, str]],
) -> Dict[str, str]:
    if not cores:
        raise ValueError(
            "Cannot choose representative from empty core list."
        )

    return sorted(
        cores,
        key=structural_core_sort_key,
    )[0]


def group_structural_cores_by_region(
    structural_cores: Sequence[Dict[str, str]],
) -> Dict[Tuple[str, str], List[Dict[str, str]]]:
    groups = defaultdict(list)

    for core in structural_cores:
        key = (
            lower(core["species"]),
            clean(core["source_region_id"]),
        )

        groups[key].append(core)

    return dict(groups)


# =============================================================================
# HYPOTHESIS CONSTRUCTION
# =============================================================================

def structural_strengths(
    core: Dict[str, str],
) -> str:
    counts = structural_counts(core)
    strengths = []

    n_discriminators = discriminator_count(core)

    if n_discriminators == 1:
        strengths.append(
            "within_species_NKG2A_vs_NKG2C_discriminator"
        )
    else:
        strengths.append(
            f"{n_discriminators}_within_species_NKG2A_vs_NKG2C_discriminators"
        )

    if counts["exposed_noninterface"] > 0:
        strengths.append(
            f"{counts['exposed_noninterface']}_exposed_noninterface_discriminatory_residue"
            + (
                "s"
                if counts["exposed_noninterface"] != 1
                else ""
            )
        )

    if (
        counts["interface"] == 0
        and counts["CD94"] == 0
        and counts["HLA_E"] == 0
        and counts["peptide"] == 0
        and counts["B2M"] == 0
    ):
        strengths.append(
            "no_detected_CD94_HLAE_B2M_peptide_interface_contact"
        )

    if (
        lower(core.get("core_type"))
        == "resolved_local_core"
    ):
        strengths.append(
            "compact_structural_neighborhood"
        )

    return join_semicolon(strengths)


def structural_limitations(
    core: Dict[str, str],
) -> str:
    counts = structural_counts(core)
    limitations = []

    species = lower(core.get("species"))

    if counts["buried"] > 0:
        limitations.append(
            f"{counts['buried']}_buried_discriminatory_residue"
            + (
                "s"
                if counts["buried"] != 1
                else ""
            )
        )

    if counts["interface"] > 0:
        limitations.append(
            f"{counts['interface']}_interface_associated_discriminatory_residue"
            + (
                "s"
                if counts["interface"] != 1
                else ""
            )
        )

    glyco = glyco_category(core)

    if glyco == "overlapping":
        limitations.append(
            "overlaps_canonical_N_glycosylation_sequon"
        )

    elif glyco == "nearby":
        limitations.append(
            "nearby_canonical_N_glycosylation_sequon"
        )

    rsa = safe_float(
        core.get("mean_complex_rsa")
    )

    if rsa is not None and rsa < 0.10:
        limitations.append(
            "low_mean_structural_accessibility"
        )

    elif rsa is not None and rsa < 0.25:
        limitations.append(
            "limited_mean_structural_accessibility"
        )

    if species in {"rhesus", "pigtail"}:
        limitations.append(
            "structural_accessibility_projected_from_human_3CDG_homologous_position"
        )
        limitations.append(
            "not_direct_macaque_structure"
        )

    return join_semicolon(limitations)


def unresolved_strengths(
    row: Dict[str, str],
) -> str:
    n_discriminators = discriminator_count(row)

    strengths = [
        "within_species_NKG2A_vs_NKG2C_sequence_discrimination",
        "retained_despite_missing_experimental_coordinates",
    ]

    if n_discriminators > 1:
        strengths.append(
            f"{n_discriminators}_sequence_discriminators_in_local_region"
        )

    return join_semicolon(strengths)


def unresolved_limitations(
    row: Dict[str, str],
) -> str:
    limitations = [
        "no_experimental_structural_coordinates",
        "surface_accessibility_not_established",
        "fixed_3D_geometry_not_assigned",
    ]

    glyco = glyco_category(row)

    if glyco == "overlapping":
        limitations.append(
            "overlaps_canonical_N_glycosylation_sequon"
        )

    elif glyco == "nearby":
        limitations.append(
            "nearby_canonical_N_glycosylation_sequon"
        )

    species = lower(row.get("species"))

    if species == "human":
        limitations.append(
            "AlphaFold_94_112_geometry_low_confidence_and_high_relative_PAE"
        )

    else:
        limitations.append(
            "no_direct_macaque_structural_geometry_inferred"
        )

    return join_semicolon(limitations)


def structural_interpretation(
    species: str,
    core: Dict[str, str],
    tier: str,
) -> str:
    labels = clean(core.get("residue_labels"))

    if tier == "tier_1_resolved_exposed_noninterface":
        text = (
            f"{labels} contains within-species NKG2A-vs-NKG2C "
            "discriminatory residue(s) with favorable exposed "
            "non-interface structural evidence."
        )

    elif tier == "tier_2_resolved_mixed_accessibility":
        text = (
            f"{labels} is a resolved discriminatory structural "
            "hypothesis with intermediate accessibility evidence."
        )

    elif tier == "tier_3_resolved_interface_associated":
        text = (
            f"{labels} is structurally resolved and discriminatory "
            "but contains receptor/ligand-interface involvement. "
            "Its antibody value may therefore depend on the intended "
            "mechanism."
        )

    else:
        text = (
            f"{labels} is structurally resolved and discriminatory "
            "but has low accessibility or burial evidence and is "
            "therefore a lower-priority antibody-accessible hypothesis."
        )

    if species in {"rhesus", "pigtail"}:
        text += (
            " Structural accessibility/interface annotations are "
            "homologous-position projections from human NKG2A 3CDG "
            "rather than direct macaque structural measurements."
        )

    return text


def unresolved_interpretation(
    species: str,
    row: Dict[str, str],
) -> str:
    labels = clean(row.get("residue_labels"))

    if species == "human":
        return (
            f"{labels} is retained as a human NKG2A-vs-NKG2C "
            "sequence-defined specificity hypothesis. Experimental "
            "human NKG2A structures do not resolve this region, and "
            "the AlphaFold 94-112 model does not support assigning "
            "a fixed antibody epitope geometry."
        )

    return (
        f"{labels} is retained as a {species} NKG2A-vs-NKG2C "
        "sequence-defined specificity hypothesis corresponding to "
        "an experimentally unresolved region in human NKG2A. "
        "No direct macaque fixed structural geometry is inferred."
    )


def build_structural_hypothesis(
    species: str,
    region_id: str,
    representative: Dict[str, str],
    all_region_cores: Sequence[Dict[str, str]],
) -> Dict[str, str]:
    tier = assign_structural_tier(
        representative
    )

    counts = structural_counts(
        representative
    )

    alternatives = [
        clean(core.get("core_id"))
        for core in sorted(
            all_region_cores,
            key=structural_core_sort_key,
        )
        if clean(core.get("core_id"))
        != clean(representative.get("core_id"))
    ]

    diameter = safe_float(
        get_first(
            representative,
            "known_pairwise_core_diameter_A",
            "diameter_A",
        )
    )

    rsa = safe_float(
        representative.get("mean_complex_rsa")
    )

    plddt = safe_float(
        representative.get("mean_alphafold_plddt")
    )

    return {
        "species": species,
        "hypothesis_rank": "",
        "hypothesis_id": "",
        "hypothesis_type": "resolved_structural_hypothesis",
        "parent_region_id": region_id,
        "representative_core": clean(
            representative.get("core_id")
        ),
        "representative_residue_labels": clean(
            representative.get("residue_labels")
        ),
        "within_species_comparisons": clean(
            representative.get(
                "within_species_comparisons"
            )
        ),
        "evidence_tier": tier,
        "discriminator_count": str(
            discriminator_count(representative)
        ),
        "resolved_residue_count": str(
            counts["resolved"]
        ),
        "unresolved_residue_count": str(
            counts["unresolved"]
        ),
        "exposed_noninterface_count": str(
            counts["exposed_noninterface"]
        ),
        "exposed_interface_count": str(
            counts["exposed_interface"]
        ),
        "partially_exposed_noninterface_count": str(
            counts["partially_exposed_noninterface"]
        ),
        "partially_exposed_interface_count": str(
            counts["partially_exposed_interface"]
        ),
        "interface_residue_count": str(
            counts["interface"]
        ),
        "CD94_contact_residue_count": str(
            counts["CD94"]
        ),
        "HLA_E_contact_residue_count": str(
            counts["HLA_E"]
        ),
        "peptide_contact_residue_count": str(
            counts["peptide"]
        ),
        "B2M_contact_residue_count": str(
            counts["B2M"]
        ),
        "buried_residue_count": str(
            counts["buried"]
        ),
        "mean_complex_rsa": fmt_float(
            rsa,
            decimals=4,
        ),
        "mean_alphafold_plddt": fmt_float(
            plddt,
            decimals=3,
        ),
        "diameter_A": fmt_float(
            diameter,
            decimals=3,
        ),
        "glycosylation_context": clean(
            representative.get(
                "glycosylation_context"
            )
        ),
        "glycosylation_sequons": clean(
            representative.get(
                "glycosylation_sequons"
            )
        ),
        "geometry_evidence_basis": clean(
            representative.get(
                "geometry_evidence_basis"
            )
        ),
        "alternative_overlapping_cores": ",".join(
            alternatives
        ),
        "redundant_core_count": str(
            len(all_region_cores)
        ),
        "hypothesis_strengths": structural_strengths(
            representative
        ),
        "hypothesis_limitations": structural_limitations(
            representative
        ),
        "interpretation": structural_interpretation(
            species,
            representative,
            tier,
        ),
    }


def choose_unresolved_evidence_row(
    species: str,
    region_id: str,
    region: Dict[str, str],
    unresolved_cores_by_region: Dict[
        Tuple[str, str],
        List[Dict[str, str]],
    ],
) -> Tuple[Dict[str, str], str]:
    """
    Prefer the Step 2T unresolved core because it carries the refined
    Step 2T evidence fields.

    Fall back to Step 2S region information if no unresolved Step 2T
    core exists. This keeps Step 2U robust without pretending that
    the region is structural.
    """

    matching = unresolved_cores_by_region.get(
        (
            species,
            region_id,
        ),
        [],
    )

    if matching:
        matching = sorted(
            matching,
            key=lambda row: (
                -discriminator_count(row),
                clean(row.get("core_id")),
            ),
        )

        return matching[0], clean(
            matching[0].get("core_id")
        )

    return region, ""


def build_unresolved_hypothesis(
    species: str,
    region_id: str,
    region: Dict[str, str],
    evidence_row: Dict[str, str],
    unresolved_core_id: str,
) -> Dict[str, str]:
    labels = get_first(
        evidence_row,
        "residue_labels",
        default=clean(region.get("residue_labels")),
    )

    comparisons = get_first(
        evidence_row,
        "within_species_comparisons",
        default=clean(
            region.get("within_species_comparisons")
        ),
    )

    glyco_context = get_first(
        evidence_row,
        "glycosylation_context",
        default=clean(
            region.get("glycosylation_context")
        ),
    )

    glyco_sequons = get_first(
        evidence_row,
        "glycosylation_sequons",
        default=clean(
            region.get("glycosylation_sequons")
        ),
    )

    plddt = safe_float(
        get_first(
            evidence_row,
            "mean_alphafold_plddt",
            default=clean(
                region.get("mean_alphafold_plddt")
            ),
        )
    )

    geometry = get_first(
        evidence_row,
        "geometry_evidence_basis",
        default=(
            "sequence_defined_only;"
            "no_experimental_fixed_geometry"
        ),
    )

    temp_row = dict(region)
    temp_row.update(evidence_row)
    temp_row["species"] = species
    temp_row["residue_labels"] = labels
    temp_row["within_species_comparisons"] = comparisons
    temp_row["glycosylation_context"] = glyco_context

    n_discriminators = discriminator_count(
        temp_row
    )

    unresolved_count = safe_int(
        evidence_row.get("unresolved_residue_count"),
        default=len(split_labels(labels)),
    )

    return {
        "species": species,
        "hypothesis_rank": "",
        "hypothesis_id": "",
        "hypothesis_type": "sequence_defined_unresolved_hypothesis",
        "parent_region_id": region_id,

        # Important:
        # The Step 2T unresolved core may be named here for traceability,
        # but it is NOT treated as a structural representative.
        "representative_core": unresolved_core_id,

        "representative_residue_labels": labels,
        "within_species_comparisons": comparisons,
        "evidence_tier": "tier_U_sequence_defined_unresolved",
        "discriminator_count": str(
            n_discriminators
        ),
        "resolved_residue_count": "0",
        "unresolved_residue_count": str(
            unresolved_count
        ),
        "exposed_noninterface_count": "",
        "exposed_interface_count": "",
        "partially_exposed_noninterface_count": "",
        "partially_exposed_interface_count": "",
        "interface_residue_count": "",
        "CD94_contact_residue_count": "",
        "HLA_E_contact_residue_count": "",
        "peptide_contact_residue_count": "",
        "B2M_contact_residue_count": "",
        "buried_residue_count": "",
        "mean_complex_rsa": "",
        "mean_alphafold_plddt": fmt_float(
            plddt,
            decimals=3,
        ),
        "diameter_A": "",
        "glycosylation_context": glyco_context,
        "glycosylation_sequons": glyco_sequons,
        "geometry_evidence_basis": geometry,
        "alternative_overlapping_cores": "",
        "redundant_core_count": "1" if unresolved_core_id else "0",
        "hypothesis_strengths": unresolved_strengths(
            temp_row
        ),
        "hypothesis_limitations": unresolved_limitations(
            temp_row
        ),
        "interpretation": unresolved_interpretation(
            species,
            temp_row,
        ),
    }


# =============================================================================
# FINAL HYPOTHESIS RANKING
# =============================================================================

def hypothesis_sort_key(
    row: Dict[str, str],
) -> Tuple:
    """
    Structural hypotheses and unresolved sequence hypotheses are not
    made artificially commensurate through an additive score.

    Experimental-priority ordering:
      1. Tier-1 resolved exposed non-interface
      2. Tier-2 resolved mixed-accessibility
      3. Sequence-defined unresolved
      4. Tier-3 resolved interface-associated
      5. Tier-4 resolved low-accessibility/buried

    IMPORTANT
    ---------
    Tier U is not being assigned stronger structural evidence than
    Tier 3. It has no resolved structural evidence.

    Instead, this ordering reflects experimental antibody-discovery
    priority. An unresolved but strongly discriminatory sequence
    region is retained ahead of a resolved interface-associated
    candidate because interface involvement is known negative/caveat
    evidence for a clean antibody-accessible specificity surface,
    whereas lack of coordinates represents uncertainty rather than
    demonstrated interface involvement or burial.

    Tier U therefore remains explicitly labeled:
        tier_U_sequence_defined_unresolved

    and is never assigned RSA, interface, diameter, or fixed 3D
    geometry.
    """

    tier = clean(
        row.get("evidence_tier")
    )

    tier_order = {
        "tier_1_resolved_exposed_noninterface": 1,
        "tier_2_resolved_mixed_accessibility": 2,
        "tier_U_sequence_defined_unresolved": 3,
        "tier_3_resolved_interface_associated": 4,
        "tier_4_resolved_low_accessibility": 5,
    }

    rank = tier_order.get(
        tier,
        99,
    )

    discriminators = safe_int(
        row.get("discriminator_count")
    )

    exposed = safe_int(
        row.get("exposed_noninterface_count")
    )

    interface = safe_int(
        row.get("interface_residue_count")
    )

    buried = safe_int(
        row.get("buried_residue_count")
    )

    rsa = safe_float(
        row.get("mean_complex_rsa")
    )

    if rsa is None:
        rsa = -1.0

    diameter = safe_float(
        row.get("diameter_A")
    )

    if diameter is None:
        diameter = 9999.0

    return (
        rank,
        -exposed,
        interface,
        buried,
        -discriminators,
        -rsa,
        diameter,
        clean(row.get("parent_region_id")),
    )


def assign_hypothesis_ids_and_ranks(
    species: str,
    hypotheses: Sequence[Dict[str, str]],
) -> List[Dict[str, str]]:
    ordered = sorted(
        hypotheses,
        key=hypothesis_sort_key,
    )

    prefix = HYPOTHESIS_PREFIX[species]

    output = []

    for index, hypothesis in enumerate(
        ordered,
        start=1,
    ):
        row = dict(hypothesis)

        row["hypothesis_rank"] = str(index)
        row["hypothesis_id"] = f"{prefix}{index}"

        output.append(row)

    return output


# =============================================================================
# BUILD SPECIES HYPOTHESES
# =============================================================================

def build_species_hypotheses(
    species: str,
    regions: Sequence[Dict[str, str]],
    structural_cores: Sequence[Dict[str, str]],
    unresolved_cores: Sequence[Dict[str, str]],
) -> List[Dict[str, str]]:
    species_regions = [
        row
        for row in regions
        if lower(row.get("species")) == species
    ]

    species_structural = [
        row
        for row in structural_cores
        if lower(row.get("species")) == species
    ]

    species_unresolved = [
        row
        for row in unresolved_cores
        if lower(row.get("species")) == species
    ]

    structural_by_region = group_structural_cores_by_region(
        species_structural
    )

    unresolved_by_region = defaultdict(list)

    for core in species_unresolved:
        unresolved_by_region[
            (
                species,
                clean(core.get("source_region_id")),
            )
        ].append(core)

    hypotheses = []

    for region in species_regions:
        region_id = clean(
            region.get("region_id")
        )

        region_type = lower(
            region.get("region_type")
        )

        if region_type in STRUCTURAL_REGION_TYPES:
            region_cores = structural_by_region.get(
                (
                    species,
                    region_id,
                ),
                [],
            )

            if not region_cores:
                raise RuntimeError(
                    f"{species} Step 2S structural region "
                    f"{region_id} has no corresponding resolved "
                    "Step 2T structural core/site."
                )

            representative = choose_representative_core(
                region_cores
            )

            hypotheses.append(
                build_structural_hypothesis(
                    species=species,
                    region_id=region_id,
                    representative=representative,
                    all_region_cores=region_cores,
                )
            )

        elif region_type in UNRESOLVED_REGION_TYPES:
            evidence_row, unresolved_core_id = (
                choose_unresolved_evidence_row(
                    species=species,
                    region_id=region_id,
                    region=region,
                    unresolved_cores_by_region=unresolved_by_region,
                )
            )

            hypotheses.append(
                build_unresolved_hypothesis(
                    species=species,
                    region_id=region_id,
                    region=region,
                    evidence_row=evidence_row,
                    unresolved_core_id=unresolved_core_id,
                )
            )

        else:
            raise RuntimeError(
                f"{species} {region_id}: unsupported Step 2S "
                f"region_type {region_type!r}"
            )

    return assign_hypothesis_ids_and_ranks(
        species,
        hypotheses,
    )


# =============================================================================
# REPORTING
# =============================================================================

def print_input_summary(
    regions: Sequence[Dict[str, str]],
    cores: Sequence[Dict[str, str]],
    structural_cores: Sequence[Dict[str, str]],
    unresolved_cores: Sequence[Dict[str, str]],
) -> None:
    print(f"Step 2S regions loaded: {len(regions)}")

    for species in SPECIES_ORDER:
        count = sum(
            1
            for row in regions
            if lower(row.get("species")) == species
        )

        print(f"  {species:<8} {count}")

    print()

    print("Step 2T combined compact-core file:")
    print(f"  {CORE_FILE}")
    print()

    print(f"Step 2T core rows loaded: {len(cores)}")

    for species in SPECIES_ORDER:
        total = sum(
            1
            for row in cores
            if lower(row.get("species")) == species
        )

        structural = sum(
            1
            for row in structural_cores
            if lower(row.get("species")) == species
        )

        unresolved = sum(
            1
            for row in unresolved_cores
            if lower(row.get("species")) == species
        )

        print(
            f"  {species:<8} "
            f"{total:2d} total; "
            f"{structural:2d} structural; "
            f"{unresolved:2d} sequence-defined unresolved"
        )


def print_species_hypotheses(
    species: str,
    rows: Sequence[Dict[str, str]],
) -> None:
    subsection(
        f"{species.upper()} WITHIN-SPECIES NKG2A EPITOPE HYPOTHESES"
    )

    for row in rows:
        hypothesis_id = clean(
            row.get("hypothesis_id")
        )

        rank = clean(
            row.get("hypothesis_rank")
        )

        labels = clean(
            row.get("representative_residue_labels")
        )

        tier = clean(
            row.get("evidence_tier")
        )

        hypothesis_type = clean(
            row.get("hypothesis_type")
        )

        parent = clean(
            row.get("parent_region_id")
        )

        core = clean(
            row.get("representative_core")
        )

        discriminator_n = clean(
            row.get("discriminator_count")
        )

        exposed = clean(
            row.get("exposed_noninterface_count")
        )

        interface = clean(
            row.get("interface_residue_count")
        )

        buried = clean(
            row.get("buried_residue_count")
        )

        rsa = clean(
            row.get("mean_complex_rsa")
        )

        diameter = clean(
            row.get("diameter_A")
        )

        glyco = clean(
            row.get("glycosylation_context")
        )

        alternatives = clean(
            row.get("alternative_overlapping_cores")
        )

        print(
            f"{hypothesis_id:<7} "
            f"rank={rank:<2} "
            f"{labels:<50} "
            f"{tier}"
        )

        print(
            f"  Type:       {hypothesis_type}"
        )

        print(
            f"  Parent:     {parent}"
        )

        if core:
            if hypothesis_type == "sequence_defined_unresolved_hypothesis":
                print(
                    f"  Step 2T:    {core} "
                    "(sequence-defined unresolved; not structural)"
                )
            else:
                print(
                    f"  Core:       {core}"
                )

        print(
            "  Evidence:   "
            f"discriminators={discriminator_n}  "
            f"exposed_noninterface={exposed}  "
            f"interface={interface}  "
            f"buried={buried}"
        )

        if rsa:
            print(
                f"  Mean RSA:   {rsa}"
            )

        if diameter:
            print(
                f"  Diameter:   {diameter} A"
            )

        if glyco:
            print(
                f"  Glyco:      {glyco}"
            )

        if alternatives:
            print(
                f"  Alternatives: {alternatives}"
            )

        strengths = clean(
            row.get("hypothesis_strengths")
        )

        limitations = clean(
            row.get("hypothesis_limitations")
        )

        if strengths:
            print(
                f"  Strengths:  {strengths}"
            )

        if limitations:
            print(
                f"  Limits:     {limitations}"
            )

        print()


def print_summary(
    all_species_rows: Dict[str, List[Dict[str, str]]],
) -> None:
    subsection("STEP 2U SUMMARY")

    for species in SPECIES_ORDER:
        rows = all_species_rows.get(
            species,
            [],
        )

        structural = sum(
            1
            for row in rows
            if clean(row.get("hypothesis_type"))
            == "resolved_structural_hypothesis"
        )

        unresolved = sum(
            1
            for row in rows
            if clean(row.get("hypothesis_type"))
            == "sequence_defined_unresolved_hypothesis"
        )

        tier1 = sum(
            1
            for row in rows
            if clean(row.get("evidence_tier"))
            == "tier_1_resolved_exposed_noninterface"
        )

        print(
            f"{species:<8} "
            f"{len(rows):2d} hypothesis/hypotheses; "
            f"{structural:2d} resolved structural; "
            f"{unresolved:2d} sequence-defined unresolved; "
            f"{tier1:2d} tier-1 structural"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    banner(
        "STEP 2U - RANK WITHIN-SPECIES NKG2A-SPECIFIC EPITOPE HYPOTHESES"
    )

    print()
    print("Analysis axis:")
    print("  human   NKG2A vs human NKG2C")
    print("  rhesus  NKG2A vs BOTH rhesus NKG2C isoforms")
    print("  pigtail NKG2A vs pigtail NKG2C")
    print()
    print(
        "Cross-species conservation/reactivity is NOT used as a "
        "ranking criterion."
    )
    print(
        "Ranking uses categorical evidence tiers rather than an "
        "arbitrary additive score."
    )
    print(
        "Resolved structural cores are consolidated by Step 2S "
        "parent region before final ranking."
    )
    print(
        "Step 2T sequence_defined_unresolved cores are retained as "
        "unresolved hypotheses and are NOT structurally validated."
    )
    print()

    require_file(
        REGION_FILE,
        "Step 2S combined candidate-region table",
    )

    require_file(
        CORE_FILE,
        "Step 2T combined local-core table",
    )

    regions = read_tsv(
        REGION_FILE
    )

    cores = read_tsv(
        CORE_FILE
    )

    validate_region_table(
        regions
    )

    validate_core_table(
        cores
    )

    validate_core_region_relationships(
        cores,
        regions,
    )

    structural_cores = [
        row
        for row in cores
        if lower(row.get("core_type"))
        in STRUCTURAL_CORE_TYPES
    ]

    unresolved_cores = [
        row
        for row in cores
        if lower(row.get("core_type"))
        in UNRESOLVED_CORE_TYPES
    ]

    validate_structural_cores(
        structural_cores,
        regions,
    )

    validate_unresolved_cores(
        unresolved_cores,
        regions,
    )

    print_input_summary(
        regions,
        cores,
        structural_cores,
        unresolved_cores,
    )

    all_species_rows = {}

    for species in [
        "human",
        "rhesus",
        "pigtail",
    ]:
        rows = build_species_hypotheses(
            species=species,
            regions=regions,
            structural_cores=structural_cores,
            unresolved_cores=unresolved_cores,
        )

        all_species_rows[species] = rows

        print_species_hypotheses(
            species,
            rows,
        )

    combined = []

    for species in [
        "human",
        "rhesus",
        "pigtail",
    ]:
        combined.extend(
            all_species_rows[species]
        )

    # Final sanity checks.
    errors = []

    for row in combined:
        hypothesis_type = clean(
            row.get("hypothesis_type")
        )

        tier = clean(
            row.get("evidence_tier")
        )

        parent = clean(
            row.get("parent_region_id")
        )

        species = clean(
            row.get("species")
        )

        if (
            hypothesis_type
            == "sequence_defined_unresolved_hypothesis"
        ):
            if tier != "tier_U_sequence_defined_unresolved":
                errors.append(
                    f"{species} {parent}: unresolved hypothesis "
                    f"has structural tier {tier!r}"
                )

            if clean(row.get("mean_complex_rsa")):
                errors.append(
                    f"{species} {parent}: unresolved hypothesis "
                    "unexpectedly has mean_complex_rsa"
                )

            if clean(row.get("diameter_A")):
                errors.append(
                    f"{species} {parent}: unresolved hypothesis "
                    "unexpectedly has structural diameter"
                )

        elif (
            hypothesis_type
            == "resolved_structural_hypothesis"
        ):
            if tier == "tier_U_sequence_defined_unresolved":
                errors.append(
                    f"{species} {parent}: structural hypothesis "
                    "has unresolved evidence tier"
                )

            if not clean(
                row.get("representative_core")
            ):
                errors.append(
                    f"{species} {parent}: structural hypothesis "
                    "has no representative structural core"
                )

        else:
            errors.append(
                f"{species} {parent}: unsupported final "
                f"hypothesis_type {hypothesis_type!r}"
            )

    if errors:
        raise RuntimeError(
            "Final Step 2U hypothesis validation failed:\n\n  "
            + "\n  ".join(errors)
        )

    # Write species-specific tables.
    for species, rows in all_species_rows.items():
        write_tsv(
            OUTPUT_BY_SPECIES[species],
            rows,
        )

    # Write combined table.
    write_tsv(
        OUTPUT_COMBINED,
        combined,
    )

    print_summary(
        all_species_rows
    )

    subsection("OUTPUTS")

    for species in [
        "human",
        "rhesus",
        "pigtail",
    ]:
        print(
            OUTPUT_BY_SPECIES[species]
        )

    print(
        OUTPUT_COMBINED
    )

    print()
    print(
        "NOTE: Step 2U ranks NKG2A-specific antibody epitope "
        "hypotheses independently"
    )
    print(
        "for human, rhesus macaque, and pigtail macaque using "
        "within-species"
    )
    print(
        "NKG2A-vs-NKG2C separation only."
    )
    print(
        "Cross-species conservation and antibody cross-reactivity "
        "are not ranking criteria."
    )
    print()
    print(
        "Resolved Step 2T local cores/sites are structurally "
        "validated and ranked using"
    )
    print(
        "categorical evidence tiers based on accessibility, "
        "interface context, burial,"
    )
    print(
        "glycosylation context, discrimination, and compactness."
    )
    print()
    print(
        "Step 2T sequence_defined_unresolved cores are explicitly "
        "recognized as a separate"
    )
    print(
        "evidence class. They are retained as sequence-defined "
        "unresolved hypotheses and"
    )
    print(
        "are NOT assigned RSA, interface, diameter, or fixed 3D "
        "epitope geometry."
    )
    print()
    print(
        "Human experimentally unresolved N-terminal hypotheses "
        "remain sequence-defined because"
    )
    print(
        "Step 2Q showed low-confidence/high-PAE AlphaFold geometry "
        "for residues 94-112."
    )
    print(
        "Rhesus and pigtail structural annotations remain "
        "homologous-position projections from"
    )
    print(
        "human NKG2A 3CDG rather than direct macaque structural "
        "measurements."
    )
    print()
    print(
        "These hypotheses are screening priorities and are not "
        "experimentally validated epitopes."
    )


if __name__ == "__main__":
    try:
        main()

    except Exception as exc:
        print(
            f"\nERROR: {exc}",
            file=sys.stderr,
        )
        raise
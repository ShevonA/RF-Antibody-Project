#!/usr/bin/env python3
"""
STEP 2W - BUILD NKG2A EXPERIMENTAL VALIDATION MATRIX
====================================================

Purpose
-------
Translate the final Step 2V within-species NKG2A-specific candidate shortlist
into an experimental-validation evidence matrix and a mutation/sequence-mapping
plan.

Step 2W DOES NOT rerank candidates.

Step 2V remains the computational candidate-selection endpoint.

Analysis axis
-------------
Axis 1 only:

    human:
        human NKG2A vs human NKG2C

    rhesus:
        rhesus NKG2A vs BOTH rhesus NKG2C isoforms

    pigtail:
        pigtail NKG2A vs pigtail NKG2C

Cross-species conservation/reactivity is NOT used to rank candidates here.

Key Step 2W rules
-----------------
1. Step 2V final_species_rank is preserved exactly.

2. Resolved structural candidates:
       Favor reciprocal NKG2A <-> NKG2C residue-swap experiments.

3. Multi-residue resolved candidates:
       Generate individual substitutions plus a combined substitution.

4. Sequence-defined unresolved candidates:
       Do NOT assign fixed structural geometry.
       Generate exploratory sequence-mapping designs only.

5. Alignment/gap differences:
       A comparison such as:

           95S>-
           95S>-/-          (rhesus two-isoform notation)

       is NOT represented as a point mutation.

       It is classified as:

           sequence_length_or_alignment_state_difference

       and retained as an alignment/indel-state experimental hypothesis.

6. Rhesus identical isoform states:
       A comparison such as:

           98T>F/F

       means:

           NKG2C isoform 1 = F
           NKG2C isoform 2 = F

       For conventional mutation naming this is collapsed to:

           T98F
           F98T

       while the original isoform state "F/F" is retained separately.

7. Rhesus discordant isoform states:
       If the two NKG2C isoforms differ at a candidate position, no single
       mutation label is silently invented. The row is flagged as requiring
       isoform-specific design.

8. Glycosylation:
       Candidates overlapping or near canonical N-X-S/T sequons are flagged
       for glycosylation-context validation.

9. Expression/folding:
       All mutagenesis experiments require matched receptor surface-expression
       or orthogonal folding/expression controls.

10. Binding and function:
       NKG2A-vs-NKG2C antibody specificity and functional blockade are treated
       as separate experimental questions.

Input
-----
results/tables/structure/
    nkg2a_within_species_final_candidate_summary.tsv

Outputs
-------
results/tables/validation/
    human_NKG2A_experimental_validation_plan.tsv
    rhesus_NKG2A_experimental_validation_plan.tsv
    pigtail_NKG2A_experimental_validation_plan.tsv
    nkg2a_within_species_experimental_validation_matrix.tsv
    nkg2a_within_species_mutagenesis_plan.tsv
"""

from __future__ import annotations

import csv
import math
import re
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

VALIDATION_DIR = (
    ROOT
    / "results"
    / "tables"
    / "validation"
)

INPUT_FILE = (
    STRUCTURE_DIR
    / "nkg2a_within_species_final_candidate_summary.tsv"
)

OUTPUT_FILES = {
    "human": (
        VALIDATION_DIR
        / "human_NKG2A_experimental_validation_plan.tsv"
    ),
    "rhesus": (
        VALIDATION_DIR
        / "rhesus_NKG2A_experimental_validation_plan.tsv"
    ),
    "pigtail": (
        VALIDATION_DIR
        / "pigtail_NKG2A_experimental_validation_plan.tsv"
    ),
}

COMBINED_OUTPUT = (
    VALIDATION_DIR
    / "nkg2a_within_species_experimental_validation_matrix.tsv"
)

MUTAGENESIS_OUTPUT = (
    VALIDATION_DIR
    / "nkg2a_within_species_mutagenesis_plan.tsv"
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


def yes_no(
    condition: bool,
) -> str:
    return "yes" if condition else "no"


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

    return separator.join(output)


def read_tsv(
    path: Path,
) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(
            f"Required Step 2V input file not found:\n{path}"
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

        rows = []

        for row in reader:

            cleaned_row = {}

            for key, value in row.items():

                if key is None:
                    continue

                cleaned_row[
                    clean(key)
                ] = clean(value)

            rows.append(
                cleaned_row
            )

        return rows


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


# =============================================================================
# INPUT VALIDATION
# =============================================================================

ALLOWED_SELECTION_CATEGORIES = {
    "primary_resolved_candidate",
    "primary_sequence_defined_candidate",
    "secondary_resolved_candidate",
    "interface_associated_candidate",
    "low_accessibility_candidate",
}

ALLOWED_STRUCTURAL_CONFIDENCE = {
    "direct_human_experimental_structure",
    "human_structure_homology_projection",
    "unresolved_no_fixed_geometry",
}


def validate_input(
    rows: Sequence[Dict[str, str]],
) -> None:

    if not rows:
        raise RuntimeError(
            "Step 2V final candidate table contains no rows."
        )

    required_columns = {
        "species",
        "final_species_rank",
        "hypothesis_id",
        "candidate_residues",
        "within_species_comparisons",
        "hypothesis_type",
        "evidence_tier",
        "final_selection_category",
        "experimental_priority",
        "target_attractiveness",
        "structural_confidence",
        "discriminator_count",
        "glycosylation_flag",
    }

    missing = sorted(
        required_columns
        - set(rows[0].keys())
    )

    if missing:
        raise RuntimeError(
            "Step 2V input is missing required column(s):\n  "
            + "\n  ".join(missing)
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

        category = clean(
            row.get(
                "final_selection_category"
            )
        )

        confidence = clean(
            row.get(
                "structural_confidence"
            )
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

        if category not in ALLOWED_SELECTION_CATEGORIES:
            errors.append(
                f"{species} {hypothesis_id}: "
                f"unsupported final_selection_category "
                f"{category!r}"
            )

        if confidence not in ALLOWED_STRUCTURAL_CONFIDENCE:
            errors.append(
                f"{species} {hypothesis_id}: "
                f"unsupported structural_confidence "
                f"{confidence!r}"
            )

        if (
            confidence
            == "unresolved_no_fixed_geometry"
        ):

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
                        f"unresolved candidate unexpectedly contains "
                        f"{field}"
                    )

    missing_species = [
        species
        for species in SPECIES
        if not any(
            lower(
                row.get("species")
            ) == species
            for row in rows
        )
    ]

    if missing_species:
        errors.append(
            "Missing species from Step 2V table: "
            + ", ".join(missing_species)
        )

    if errors:
        raise RuntimeError(
            "Step 2W input validation failed:\n\n  "
            + "\n  ".join(errors)
        )


# =============================================================================
# RESIDUE LABEL PARSING
# =============================================================================

RESIDUE_LABEL_PATTERN = re.compile(
    r"^(\d+)([A-Za-z])$"
)


def parse_residue_labels(
    text: str,
) -> List[Tuple[int, str]]:

    text = clean(text)

    if not text:
        return []

    output = []

    for token in text.split(","):

        token = token.strip()

        match = RESIDUE_LABEL_PATTERN.match(
            token
        )

        if not match:
            raise ValueError(
                f"Could not parse candidate residue label: "
                f"{token!r}"
            )

        position = int(
            match.group(1)
        )

        residue = (
            match.group(2)
            .upper()
        )

        output.append(
            (
                position,
                residue,
            )
        )

    return output


# =============================================================================
# WITHIN-SPECIES COMPARISON PARSING
# =============================================================================

def split_comparison_targets(
    text: str,
) -> List[str]:

    text = clean(text)

    if not text:
        return []

    return [
        target.strip().upper()
        for target in text.split("/")
    ]


def parse_comparison_token(
    token: str,
) -> Dict[str, object]:

    token = clean(token)

    if ">" not in token:
        raise ValueError(
            f"Comparison token lacks '>': "
            f"{token!r}"
        )

    left, right = token.split(
        ">",
        1,
    )

    left = left.strip()
    right = right.strip()

    match = RESIDUE_LABEL_PATTERN.match(
        left
    )

    if not match:
        raise ValueError(
            f"Could not parse NKG2A side of comparison: "
            f"{token!r}"
        )

    position = int(
        match.group(1)
    )

    nkg2a_residue = (
        match.group(2)
        .upper()
    )

    targets = split_comparison_targets(
        right
    )

    if not targets:
        targets = [""]

    return {
        "position":
            position,

        "nkg2a_residue":
            nkg2a_residue,

        "nkg2c_targets":
            targets,

        "raw_target_state":
            right,

        "raw":
            token,
    }


def parse_comparisons(
    text: str,
) -> List[Dict[str, object]]:

    text = clean(text)

    if not text:
        return []

    output = []

    for token in text.split("|"):

        token = token.strip()

        if not token:
            continue

        output.append(
            parse_comparison_token(
                token
            )
        )

    return output


# =============================================================================
# COMPARISON STATE NORMALIZATION
# =============================================================================

def normalized_target_states(
    targets: Sequence[str],
) -> List[str]:
    """
    Preserve the biological NKG2C target states but normalize whitespace/case.

    Examples
    --------
    ["F", "F"] -> ["F", "F"]
    ["-", "-"] -> ["-", "-"]
    ["V"]      -> ["V"]
    """

    return [
        clean(target).upper()
        for target in targets
    ]


def unique_target_states(
    targets: Sequence[str],
) -> List[str]:

    output = []

    for target in normalized_target_states(
        targets
    ):

        if target not in output:
            output.append(
                target
            )

    return output


def isoform_state_string(
    targets: Sequence[str],
) -> str:
    """
    Retain the original per-isoform states.

    Examples
    --------
    ["F", "F"] -> "F/F"
    ["-", "-"] -> "-/-"
    ["V"]      -> "V"
    """

    return "/".join(
        normalized_target_states(
            targets
        )
    )


def collapsed_target_state(
    targets: Sequence[str],
) -> str:
    """
    Collapse identical NKG2C states for conventional mutation naming.

    Examples
    --------
    ["F", "F"] -> "F"
    ["-", "-"] -> "-"
    ["V"]      -> "V"

    Discordant isoforms:
    ["F", "L"] -> ""

    Empty means no single target state should be used for a conventional
    point-mutation label.
    """

    unique_states = unique_target_states(
        targets
    )

    if len(unique_states) != 1:
        return ""

    return unique_states[0]


def all_gap_states(
    targets: Sequence[str],
) -> bool:

    normalized = normalized_target_states(
        targets
    )

    return (
        bool(normalized)
        and all(
            target == "-"
            for target in normalized
        )
    )


def contains_gap_state(
    targets: Sequence[str],
) -> bool:

    return any(
        target == "-"
        for target in normalized_target_states(
            targets
        )
    )


def comparison_event_type(
    targets: Sequence[str],
) -> str:
    """
    Classify the biological comparison event.

    amino_acid_substitution
        All NKG2C states are amino-acid residues and identical.

    isoform_specific_amino_acid_states
        NKG2C isoforms differ from each other.

    sequence_length_or_alignment_state_difference
        NKG2C comparison contains gap state(s).

    missing_or_unresolved_target_state
        No interpretable target state.
    """

    normalized = normalized_target_states(
        targets
    )

    if not normalized:
        return (
            "missing_or_unresolved_target_state"
        )

    if any(
        not target
        for target in normalized
    ):
        return (
            "missing_or_unresolved_target_state"
        )

    if contains_gap_state(
        normalized
    ):
        return (
            "sequence_length_or_alignment_state_difference"
        )

    unique_states = unique_target_states(
        normalized
    )

    if len(unique_states) == 1:
        return (
            "amino_acid_substitution"
        )

    return (
        "isoform_specific_amino_acid_states"
    )


def reciprocal_design_status(
    targets: Sequence[str],
) -> str:

    event_type = comparison_event_type(
        targets
    )

    if (
        event_type
        == "amino_acid_substitution"
    ):
        return (
            "conventional_reciprocal_point_mutation_defined"
        )

    if (
        event_type
        == "sequence_length_or_alignment_state_difference"
    ):
        return (
            "requires_indel_or_region_level_sequence_design"
        )

    if (
        event_type
        == "isoform_specific_amino_acid_states"
    ):
        return (
            "requires_isoform_specific_mutation_design"
        )

    return (
        "manual_review_required"
    )


# =============================================================================
# MUTATION LABEL CONSTRUCTION
# =============================================================================

def build_point_mutation_labels(
    position: int,
    nkg2a_residue: str,
    targets: Sequence[str],
) -> Tuple[str, str]:

    event_type = comparison_event_type(
        targets
    )

    if (
        event_type
        != "amino_acid_substitution"
    ):
        return "", ""

    target = collapsed_target_state(
        targets
    )

    if (
        not target
        or target == "-"
    ):
        return "", ""

    forward = (
        f"{nkg2a_residue}"
        f"{position}"
        f"{target}"
    )

    reciprocal = (
        f"{target}"
        f"{position}"
        f"{nkg2a_residue}"
    )

    return (
        forward,
        reciprocal,
    )


def build_alignment_state_design_label(
    position: int,
    nkg2a_residue: str,
    targets: Sequence[str],
) -> str:

    state_string = isoform_state_string(
        targets
    )

    if all_gap_states(
        targets
    ):
        return (
            f"NKG2A_{nkg2a_residue}{position}_"
            f"versus_NKG2C_gap_state_{state_string}"
        )

    return (
        f"NKG2A_{nkg2a_residue}{position}_"
        f"versus_NKG2C_alignment_states_{state_string}"
    )


def build_isoform_specific_design_label(
    position: int,
    nkg2a_residue: str,
    targets: Sequence[str],
) -> str:

    state_string = isoform_state_string(
        targets
    )

    return (
        f"NKG2A_{nkg2a_residue}{position}_"
        f"versus_NKG2C_isoform_states_{state_string}"
    )


# =============================================================================
# CANDIDATE TYPE HELPERS
# =============================================================================

def is_unresolved(
    row: Dict[str, str],
) -> bool:

    return (
        clean(
            row.get(
                "structural_confidence"
            )
        )
        == "unresolved_no_fixed_geometry"
    )


def is_resolved(
    row: Dict[str, str],
) -> bool:

    return not is_unresolved(
        row
    )


def is_interface_candidate(
    row: Dict[str, str],
) -> bool:

    return (
        clean(
            row.get(
                "final_selection_category"
            )
        )
        == "interface_associated_candidate"
    )


def is_low_accessibility_candidate(
    row: Dict[str, str],
) -> bool:

    return (
        clean(
            row.get(
                "final_selection_category"
            )
        )
        == "low_accessibility_candidate"
    )


def has_glyco_overlap(
    row: Dict[str, str],
) -> bool:

    return (
        lower(
            row.get(
                "glycosylation_flag"
            )
        )
        == "overlapping"
    )


def has_nearby_glyco(
    row: Dict[str, str],
) -> bool:

    return (
        lower(
            row.get(
                "glycosylation_flag"
            )
        )
        == "nearby"
    )


def has_glyco_context(
    row: Dict[str, str],
) -> bool:

    return (
        has_glyco_overlap(row)
        or has_nearby_glyco(row)
    )


# =============================================================================
# VALIDATION STRATEGY
# =============================================================================

def primary_validation_test(
    row: Dict[str, str],
) -> str:

    residues = parse_residue_labels(
        row.get(
            "candidate_residues",
            "",
        )
    )

    if is_unresolved(row):

        if has_glyco_context(row):
            return (
                "region_level_sequence_discrimination_test_plus_"
                "linear_peptide_screen_plus_"
                "glycosylation_context_test"
            )

        return (
            "region_level_sequence_discrimination_test_plus_"
            "linear_peptide_screen"
        )

    if len(residues) > 1:

        if has_glyco_context(row):
            return (
                "individual_and_combined_reciprocal_"
                "NKG2A_NKG2C_residue_swaps_plus_"
                "glycosylation_context_control"
            )

        return (
            "individual_and_combined_reciprocal_"
            "NKG2A_NKG2C_residue_swaps"
        )

    if len(residues) == 1:

        if has_glyco_context(row):
            return (
                "single_residue_reciprocal_"
                "NKG2A_NKG2C_swap_plus_"
                "glycosylation_context_control"
            )

        return (
            "single_residue_reciprocal_"
            "NKG2A_NKG2C_swap"
        )

    return (
        "manual_review_required"
    )


def validation_strategy(
    row: Dict[str, str],
) -> str:

    strategies = []

    residues = parse_residue_labels(
        row.get(
            "candidate_residues",
            "",
        )
    )

    if is_unresolved(row):

        strategies.extend(
            [
                "sequence_defined_region_validation",
                "linear_epitope_screening",
                "targeted_sequence_substitution_if_practical",
                "alignment_or_indel_state_mapping_if_applicable",
            ]
        )

    else:

        strategies.append(
            "reciprocal_NKG2A_NKG2C_mutagenesis"
        )

        if len(residues) > 1:

            strategies.extend(
                [
                    "individual_residue_swaps",
                    "combined_residue_swap",
                ]
            )

        else:

            strategies.append(
                "single_residue_swap"
            )

    if has_glyco_context(row):

        strategies.append(
            "glycosylation_context_validation"
        )

    strategies.extend(
        [
            "surface_expression_control",
            "antibody_binding_assay",
        ]
    )

    if is_interface_candidate(row):

        strategies.append(
            "functional_receptor_ligand_assay"
        )

    return unique_join(
        strategies
    )


# =============================================================================
# EXPERIMENTAL FLAGS
# =============================================================================

def reciprocal_swap_recommended(
    row: Dict[str, str],
) -> bool:

    return is_resolved(
        row
    )


def single_residue_swap_recommended(
    row: Dict[str, str],
) -> bool:

    residues = parse_residue_labels(
        row.get(
            "candidate_residues",
            "",
        )
    )

    return (
        is_resolved(row)
        and len(residues) >= 1
    )


def combined_swap_recommended(
    row: Dict[str, str],
) -> bool:

    residues = parse_residue_labels(
        row.get(
            "candidate_residues",
            "",
        )
    )

    return (
        is_resolved(row)
        and len(residues) > 1
    )


def linear_peptide_test_recommended(
    row: Dict[str, str],
) -> bool:

    return is_unresolved(
        row
    )


def glycosylation_control_required(
    row: Dict[str, str],
) -> bool:

    return has_glyco_context(
        row
    )


def functional_validation_priority(
    row: Dict[str, str],
) -> str:

    if is_interface_candidate(row):
        return "high"

    if is_low_accessibility_candidate(row):
        return "low"

    if is_unresolved(row):
        return (
            "after_binding_specificity_is_established"
        )

    return (
        "secondary_to_binding_specificity"
    )


# =============================================================================
# INTERPRETATION LOGIC
# =============================================================================

def interpretation_if_binding_lost(
    row: Dict[str, str],
) -> str:

    if is_unresolved(row):

        return (
            "Loss of binding after targeted sequence or region-state "
            "modification would support involvement of the sequence-defined "
            "region, provided matched surface-expression and folding controls "
            "remain acceptable."
        )

    return (
        "Loss or substantial reduction of antibody binding after converting "
        "NKG2A residue(s) toward the within-species NKG2C state would support "
        "involvement of the predicted discriminatory residue(s), provided "
        "surface expression and receptor integrity are preserved."
    )


def interpretation_if_binding_gained(
    row: Dict[str, str],
) -> str:

    if is_unresolved(row):

        return (
            "Gain of binding after introducing NKG2A-like sequence or region "
            "state into the corresponding NKG2C region would strengthen "
            "evidence that the sequence-defined region contributes to "
            "NKG2A-specific recognition."
        )

    return (
        "Gain or increased binding after introducing the NKG2A residue(s) "
        "into NKG2C would provide reciprocal evidence that the candidate "
        "residue(s) contribute to NKG2A-specific recognition."
    )


def interpretation_if_binding_retained(
    row: Dict[str, str],
) -> str:

    if is_unresolved(row):

        return (
            "Retention of binding after sequence or region-state modification "
            "would argue that the tested feature is not individually "
            "sufficient to explain specificity, although neighboring residues, "
            "glycosylation, or conformational context could still contribute."
        )

    return (
        "Retention of binding after the predicted NKG2A-to-NKG2C substitution "
        "would weaken the hypothesis that the tested residue(s) are necessary "
        "for antibody recognition, assuming expression and folding controls "
        "are satisfactory."
    )


# =============================================================================
# GLYCOSYLATION PLAN
# =============================================================================

def glycosylation_validation_plan(
    row: Dict[str, str],
) -> str:

    if has_glyco_overlap(row):

        return unique_join(
            [
                "candidate_overlaps_canonical_NXS_T_context",
                "compare_binding_under_matched_glycosylation_conditions",
                "interpret_peptide_results_cautiously",
                "consider_sequon_preserving_and_sequon_disrupting_controls",
            ]
        )

    if has_nearby_glyco(row):

        return unique_join(
            [
                "canonical_NXS_T_sequon_is_near_candidate",
                "compare_binding_with_glycosylation_context_preserved",
                "evaluate_whether_local_glycan_context_changes_candidate_accessibility",
            ]
        )

    return (
        "no_specific_glycosylation_control_triggered_"
        "by_current_sequence_screen"
    )


# =============================================================================
# EXPRESSION / FOLDING CONTROLS
# =============================================================================

def expression_control_plan(
    row: Dict[str, str],
) -> str:

    if is_unresolved(row):

        return (
            "verify_matched_surface_expression_of_WT_and_sequence_or_region_"
            "variant_receptors_with_an_independent_non_candidate_region_"
            "reagent_or_orthogonal_expression_measure"
        )

    return (
        "verify_matched_surface_expression_and_receptor_integrity_for_WT_"
        "and_each_mutant_using_an_independent_non_epitope_reagent_or_"
        "orthogonal_expression_measure"
    )


# =============================================================================
# BINDING PLAN
# =============================================================================

def binding_validation_plan(
    row: Dict[str, str],
) -> str:

    if is_unresolved(row):

        return (
            "compare_antibody_binding_to_WT_NKG2A_and_corresponding_"
            "sequence_or_region_state_variants;use_peptide_or_region_level_"
            "screen_as_supportive_evidence_not_as_a_substitute_for_native_"
            "receptor_binding"
        )

    return (
        "compare_antibody_binding_across_WT_NKG2A_WT_NKG2C_"
        "NKG2A_to_NKG2C_mutant_and_reciprocal_NKG2C_to_NKG2A_mutant"
    )


# =============================================================================
# FUNCTIONAL PLAN
# =============================================================================

def functional_validation_plan(
    row: Dict[str, str],
) -> str:

    if is_interface_candidate(row):

        return (
            "after_specific_binding_is_confirmed_evaluate_whether_antibody_"
            "binding_alters_CD94_NKG2A_ligand_interaction_or_downstream_"
            "receptor_function"
        )

    if is_low_accessibility_candidate(row):

        return (
            "functional_testing_not_prioritized_until_specific_native_"
            "surface_binding_is_demonstrated"
        )

    if is_unresolved(row):

        return (
            "functional_testing_after_native_receptor_binding_specificity_"
            "and_surface_accessibility_are_established"
        )

    return (
        "functional_testing_is_secondary_to_demonstrating_"
        "NKG2A_vs_NKG2C_binding_specificity"
    )


# =============================================================================
# VALIDATION NOTES
# =============================================================================

def validation_notes(
    row: Dict[str, str],
) -> str:

    notes = []

    species = lower(
        row.get(
            "species"
        )
    )

    residues = parse_residue_labels(
        row.get(
            "candidate_residues",
            "",
        )
    )

    if is_unresolved(row):

        notes.append(
            "do_not_assign_fixed_3D_epitope_geometry"
        )

        notes.append(
            "negative_linear_peptide_binding_does_not_exclude_a_native_"
            "conformational_epitope"
        )

        notes.append(
            "gap_states_are_alignment_or_sequence_length_features_not_"
            "conventional_point_mutations"
        )

    if (
        is_resolved(row)
        and species in {
            "rhesus",
            "pigtail",
        }
    ):

        notes.append(
            "structural_accessibility_is_projected_from_human_3CDG_"
            "homologous_position_not_direct_macaque_structure"
        )

    if (
        len(residues) > 1
        and is_resolved(row)
    ):

        notes.append(
            "test_single_mutants_before_or_alongside_combined_mutant_to_"
            "separate_dominant_from_composite_effects"
        )

    if has_glyco_overlap(row):

        notes.append(
            "candidate_overlaps_canonical_N_glycosylation_sequence_context"
        )

    if has_nearby_glyco(row):

        notes.append(
            "candidate_is_near_canonical_N_glycosylation_sequence_context"
        )

    if is_interface_candidate(row):

        notes.append(
            "binding_specificity_and_functional_blockade_are_separate_"
            "experimental_questions"
        )

    if is_low_accessibility_candidate(row):

        notes.append(
            "low_structural_accessibility_makes_native_surface_binding_"
            "less_likely_under_current_model"
        )

    notes.append(
        "computational_candidate_not_experimentally_validated_epitope"
    )

    return unique_join(
        notes
    )


# =============================================================================
# BUILD CANDIDATE-LEVEL VALIDATION ROW
# =============================================================================

def build_validation_row(
    row: Dict[str, str],
) -> Dict[str, object]:

    output = dict(
        row
    )

    residues = parse_residue_labels(
        row.get(
            "candidate_residues",
            "",
        )
    )

    comparisons = parse_comparisons(
        row.get(
            "within_species_comparisons",
            "",
        )
    )

    event_types = unique_join(
        comparison_event_type(
            comparison[
                "nkg2c_targets"
            ]
        )
        for comparison in comparisons
    )

    alignment_difference_count = sum(
        1
        for comparison in comparisons
        if comparison_event_type(
            comparison[
                "nkg2c_targets"
            ]
        )
        == "sequence_length_or_alignment_state_difference"
    )

    amino_acid_substitution_count = sum(
        1
        for comparison in comparisons
        if comparison_event_type(
            comparison[
                "nkg2c_targets"
            ]
        )
        == "amino_acid_substitution"
    )

    isoform_specific_state_count = sum(
        1
        for comparison in comparisons
        if comparison_event_type(
            comparison[
                "nkg2c_targets"
            ]
        )
        == "isoform_specific_amino_acid_states"
    )

    output.update(
        {
            "candidate_residue_count":
                len(residues),

            "parsed_comparison_count":
                len(comparisons),

            "comparison_event_types":
                event_types,

            "amino_acid_substitution_count":
                amino_acid_substitution_count,

            "alignment_state_difference_count":
                alignment_difference_count,

            "isoform_specific_state_count":
                isoform_specific_state_count,

            "validation_strategy":
                validation_strategy(row),

            "primary_validation_test":
                primary_validation_test(row),

            "reciprocal_swap_recommended":
                yes_no(
                    reciprocal_swap_recommended(
                        row
                    )
                ),

            "single_residue_swap_recommended":
                yes_no(
                    single_residue_swap_recommended(
                        row
                    )
                ),

            "combined_swap_recommended":
                yes_no(
                    combined_swap_recommended(
                        row
                    )
                ),

            "linear_peptide_test_recommended":
                yes_no(
                    linear_peptide_test_recommended(
                        row
                    )
                ),

            "glycosylation_control_required":
                yes_no(
                    glycosylation_control_required(
                        row
                    )
                ),

            "surface_expression_control_required":
                "yes",

            "binding_validation_required":
                "yes",

            "functional_validation_priority":
                functional_validation_priority(
                    row
                ),

            "glycosylation_validation_plan":
                glycosylation_validation_plan(
                    row
                ),

            "surface_expression_control_plan":
                expression_control_plan(
                    row
                ),

            "binding_validation_plan":
                binding_validation_plan(
                    row
                ),

            "functional_validation_plan":
                functional_validation_plan(
                    row
                ),

            "interpretation_if_binding_lost":
                interpretation_if_binding_lost(
                    row
                ),

            "interpretation_if_binding_gained":
                interpretation_if_binding_gained(
                    row
                ),

            "interpretation_if_binding_retained":
                interpretation_if_binding_retained(
                    row
                ),

            "validation_notes":
                validation_notes(
                    row
                ),
        }
    )

    return output


# =============================================================================
# COMPARISON LOOKUP
# =============================================================================

def comparison_lookup(
    row: Dict[str, str],
) -> Dict[int, Dict[str, object]]:

    comparisons = parse_comparisons(
        row.get(
            "within_species_comparisons",
            "",
        )
    )

    return {
        int(
            item["position"]
        ): item
        for item in comparisons
    }


# =============================================================================
# MUTAGENESIS / SEQUENCE-MAPPING PLAN
# =============================================================================

def base_mutagenesis_row(
    row: Dict[str, str],
) -> Dict[str, object]:

    return {
        "species":
            lower(
                row.get(
                    "species"
                )
            ),

        "final_species_rank":
            clean(
                row.get(
                    "final_species_rank"
                )
            ),

        "hypothesis_id":
            clean(
                row.get(
                    "hypothesis_id"
                )
            ),

        "candidate_residues":
            clean(
                row.get(
                    "candidate_residues"
                )
            ),

        "structural_status":
            clean(
                row.get(
                    "structural_confidence"
                )
            ),

        "surface_expression_control_required":
            "yes",
    }


def build_unresolved_mapping_rows(
    row: Dict[str, str],
) -> List[Dict[str, object]]:

    residues = parse_residue_labels(
        row.get(
            "candidate_residues",
            "",
        )
    )

    lookup = comparison_lookup(
        row
    )

    output = []

    for position, nkg2a_residue in residues:

        base = base_mutagenesis_row(
            row
        )

        comparison = lookup.get(
            position
        )

        if comparison is None:

            base.update(
                {
                    "mutation_class":
                        "exploratory_sequence_region_mapping",

                    "comparison_event_type":
                        "missing_or_unresolved_target_state",

                    "position":
                        position,

                    "nkg2a_residue":
                        nkg2a_residue,

                    "nkg2c_isoform_states":
                        "",

                    "collapsed_nkg2c_target_state":
                        "",

                    "forward_NKG2A_to_NKG2C_mutation":
                        "",

                    "reciprocal_NKG2C_to_NKG2A_mutation":
                        "",

                    "alignment_or_region_design":
                        "",

                    "reciprocal_design_status":
                        "manual_review_required",

                    "recommended_use":
                        "manual_review_comparison_missing",

                    "notes":
                        (
                            "No parsed within-species comparison was found "
                            "for this candidate residue."
                        ),
                }
            )

            output.append(
                base
            )

            continue

        targets = comparison[
            "nkg2c_targets"
        ]

        event_type = comparison_event_type(
            targets
        )

        state_string = isoform_state_string(
            targets
        )

        collapsed_state = collapsed_target_state(
            targets
        )

        reciprocal_status = reciprocal_design_status(
            targets
        )

        forward = ""
        reciprocal = ""
        alignment_design = ""

        if (
            event_type
            == "amino_acid_substitution"
        ):

            (
                forward,
                reciprocal,
            ) = build_point_mutation_labels(
                position,
                nkg2a_residue,
                targets,
            )

            mutation_class = (
                "exploratory_sequence_region_substitution"
            )

            recommended_use = (
                "exploratory_region_mapping_point_substitution"
            )

            notes = (
                "Candidate remains structurally unresolved; conventional "
                "point-mutation notation is used only for sequence mapping "
                "and does not imply a fixed 3D epitope."
            )

        elif (
            event_type
            == "sequence_length_or_alignment_state_difference"
        ):

            mutation_class = (
                "sequence_length_or_alignment_state_difference"
            )

            alignment_design = (
                build_alignment_state_design_label(
                    position,
                    nkg2a_residue,
                    targets,
                )
            )

            recommended_use = (
                "exploratory_region_mapping_indel_or_alignment_state_test"
            )

            notes = (
                "This is an alignment/sequence-length difference, not a "
                "conventional amino-acid point substitution. Experimental "
                "design should preserve the relevant local sequence context."
            )

        elif (
            event_type
            == "isoform_specific_amino_acid_states"
        ):

            mutation_class = (
                "isoform_specific_sequence_region_substitution"
            )

            alignment_design = (
                build_isoform_specific_design_label(
                    position,
                    nkg2a_residue,
                    targets,
                )
            )

            recommended_use = (
                "exploratory_isoform_specific_region_mapping"
            )

            notes = (
                "NKG2C isoforms have different target states at this "
                "position; design and interpret isoform-specific variants "
                "separately."
            )

        else:

            mutation_class = (
                "exploratory_sequence_region_mapping"
            )

            recommended_use = (
                "manual_review_required"
            )

            notes = (
                "Target state could not be converted into a defined "
                "experimental sequence design."
            )

        base.update(
            {
                "mutation_class":
                    mutation_class,

                "comparison_event_type":
                    event_type,

                "position":
                    position,

                "nkg2a_residue":
                    nkg2a_residue,

                "nkg2c_isoform_states":
                    state_string,

                "collapsed_nkg2c_target_state":
                    collapsed_state,

                "forward_NKG2A_to_NKG2C_mutation":
                    forward,

                "reciprocal_NKG2C_to_NKG2A_mutation":
                    reciprocal,

                "alignment_or_region_design":
                    alignment_design,

                "reciprocal_design_status":
                    reciprocal_status,

                "recommended_use":
                    recommended_use,

                "notes":
                    notes,
            }
        )

        output.append(
            base
        )

    return output


def build_resolved_mutation_rows(
    row: Dict[str, str],
) -> List[Dict[str, object]]:

    residues = parse_residue_labels(
        row.get(
            "candidate_residues",
            "",
        )
    )

    lookup = comparison_lookup(
        row
    )

    output = []

    combined_forward = []
    combined_reciprocal = []

    combined_all_defined = True

    for position, nkg2a_residue in residues:

        base = base_mutagenesis_row(
            row
        )

        comparison = lookup.get(
            position
        )

        if comparison is None:

            combined_all_defined = False

            base.update(
                {
                    "mutation_class":
                        "individual_resolved_candidate",

                    "comparison_event_type":
                        "missing_or_unresolved_target_state",

                    "position":
                        position,

                    "nkg2a_residue":
                        nkg2a_residue,

                    "nkg2c_isoform_states":
                        "",

                    "collapsed_nkg2c_target_state":
                        "",

                    "forward_NKG2A_to_NKG2C_mutation":
                        "",

                    "reciprocal_NKG2C_to_NKG2A_mutation":
                        "",

                    "alignment_or_region_design":
                        "",

                    "reciprocal_design_status":
                        "manual_review_required",

                    "recommended_use":
                        "manual_review_comparison_missing",

                    "notes":
                        (
                            "No parsed within-species comparison found."
                        ),
                }
            )

            output.append(
                base
            )

            continue

        targets = comparison[
            "nkg2c_targets"
        ]

        event_type = comparison_event_type(
            targets
        )

        state_string = isoform_state_string(
            targets
        )

        collapsed_state = collapsed_target_state(
            targets
        )

        reciprocal_status = reciprocal_design_status(
            targets
        )

        forward = ""
        reciprocal = ""
        alignment_design = ""

        if (
            event_type
            == "amino_acid_substitution"
        ):

            (
                forward,
                reciprocal,
            ) = build_point_mutation_labels(
                position,
                nkg2a_residue,
                targets,
            )

            recommended_use = (
                "primary_residue_swap_test"
            )

            notes = (
                "Interpret binding changes only when matched "
                "surface-expression/folding controls are acceptable."
            )

            if forward:
                combined_forward.append(
                    forward
                )

            else:
                combined_all_defined = False

            if reciprocal:
                combined_reciprocal.append(
                    reciprocal
                )

            else:
                combined_all_defined = False

        elif (
            event_type
            == "sequence_length_or_alignment_state_difference"
        ):

            combined_all_defined = False

            alignment_design = (
                build_alignment_state_design_label(
                    position,
                    nkg2a_residue,
                    targets,
                )
            )

            recommended_use = (
                "resolved_position_but_requires_indel_or_region_state_design"
            )

            notes = (
                "The comparison contains a gap/alignment state and therefore "
                "cannot be represented as a conventional point mutation."
            )

        elif (
            event_type
            == "isoform_specific_amino_acid_states"
        ):

            combined_all_defined = False

            alignment_design = (
                build_isoform_specific_design_label(
                    position,
                    nkg2a_residue,
                    targets,
                )
            )

            recommended_use = (
                "resolved_position_requires_isoform_specific_mutation_design"
            )

            notes = (
                "NKG2C isoforms have different target residues; evaluate "
                "isoform-specific reciprocal mutations separately."
            )

        else:

            combined_all_defined = False

            recommended_use = (
                "manual_review_required"
            )

            notes = (
                "Target state could not be converted into a defined "
                "experimental mutation."
            )

        base.update(
            {
                "mutation_class":
                    "individual_resolved_candidate",

                "comparison_event_type":
                    event_type,

                "position":
                    position,

                "nkg2a_residue":
                    nkg2a_residue,

                "nkg2c_isoform_states":
                    state_string,

                "collapsed_nkg2c_target_state":
                    collapsed_state,

                "forward_NKG2A_to_NKG2C_mutation":
                    forward,

                "reciprocal_NKG2C_to_NKG2A_mutation":
                    reciprocal,

                "alignment_or_region_design":
                    alignment_design,

                "reciprocal_design_status":
                    reciprocal_status,

                "recommended_use":
                    recommended_use,

                "notes":
                    notes,
            }
        )

        output.append(
            base
        )

    # -------------------------------------------------------------------------
    # Combined mutant only when every member is a conventional substitution.
    # -------------------------------------------------------------------------

    if len(residues) > 1:

        base = base_mutagenesis_row(
            row
        )

        if (
            combined_all_defined
            and len(combined_forward)
            == len(residues)
            and len(combined_reciprocal)
            == len(residues)
        ):

            base.update(
                {
                    "mutation_class":
                        "combined_resolved_candidate",

                    "comparison_event_type":
                        "combined_amino_acid_substitution",

                    "position":
                        "",

                    "nkg2a_residue":
                        "",

                    "nkg2c_isoform_states":
                        "",

                    "collapsed_nkg2c_target_state":
                        "",

                    "forward_NKG2A_to_NKG2C_mutation":
                        "+".join(
                            combined_forward
                        ),

                    "reciprocal_NKG2C_to_NKG2A_mutation":
                        "+".join(
                            combined_reciprocal
                        ),

                    "alignment_or_region_design":
                        "",

                    "reciprocal_design_status":
                        "conventional_combined_reciprocal_mutation_defined",

                    "recommended_use":
                        "combined_candidate_surface_test",

                    "notes":
                        (
                            "Compare the combined mutant with individual "
                            "mutants to distinguish dominant from composite "
                            "residue effects."
                        ),
                }
            )

        else:

            base.update(
                {
                    "mutation_class":
                        "combined_resolved_candidate_not_directly_defined",

                    "comparison_event_type":
                        "mixed_or_non_point_mutation_candidate",

                    "position":
                        "",

                    "nkg2a_residue":
                        "",

                    "nkg2c_isoform_states":
                        "",

                    "collapsed_nkg2c_target_state":
                        "",

                    "forward_NKG2A_to_NKG2C_mutation":
                        "",

                    "reciprocal_NKG2C_to_NKG2A_mutation":
                        "",

                    "alignment_or_region_design":
                        "",

                    "reciprocal_design_status":
                        "requires_custom_combined_design",

                    "recommended_use":
                        "manual_combined_design_after_individual_review",

                    "notes":
                        (
                            "At least one member cannot be represented as a "
                            "single conventional point substitution; do not "
                            "construct a combined mutation label automatically."
                        ),
                }
            )

        output.append(
            base
        )

    return output


def build_mutagenesis_rows_for_candidate(
    row: Dict[str, str],
) -> List[Dict[str, object]]:

    if is_unresolved(
        row
    ):
        return build_unresolved_mapping_rows(
            row
        )

    return build_resolved_mutation_rows(
        row
    )


def build_mutagenesis_plan(
    rows: Sequence[Dict[str, str]],
) -> List[Dict[str, object]]:

    output = []

    for row in rows:

        output.extend(
            build_mutagenesis_rows_for_candidate(
                row
            )
        )

    return output


# =============================================================================
# OUTPUT FIELDS
# =============================================================================

STEP2V_FIELDS = [
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

STEP2W_FIELDS = [
    "candidate_residue_count",
    "parsed_comparison_count",
    "comparison_event_types",
    "amino_acid_substitution_count",
    "alignment_state_difference_count",
    "isoform_specific_state_count",
    "validation_strategy",
    "primary_validation_test",
    "reciprocal_swap_recommended",
    "single_residue_swap_recommended",
    "combined_swap_recommended",
    "linear_peptide_test_recommended",
    "glycosylation_control_required",
    "surface_expression_control_required",
    "binding_validation_required",
    "functional_validation_priority",
    "glycosylation_validation_plan",
    "surface_expression_control_plan",
    "binding_validation_plan",
    "functional_validation_plan",
    "interpretation_if_binding_lost",
    "interpretation_if_binding_gained",
    "interpretation_if_binding_retained",
    "validation_notes",
]

OUTPUT_FIELDS = (
    STEP2V_FIELDS
    + STEP2W_FIELDS
)

MUTAGENESIS_FIELDS = [
    "species",
    "final_species_rank",
    "hypothesis_id",
    "candidate_residues",
    "mutation_class",
    "comparison_event_type",
    "position",
    "nkg2a_residue",
    "nkg2c_isoform_states",
    "collapsed_nkg2c_target_state",
    "forward_NKG2A_to_NKG2C_mutation",
    "reciprocal_NKG2C_to_NKG2A_mutation",
    "alignment_or_region_design",
    "reciprocal_design_status",
    "structural_status",
    "recommended_use",
    "surface_expression_control_required",
    "notes",
]


# =============================================================================
# SORTING
# =============================================================================

def candidate_sort_key(
    row: Dict[str, object],
) -> Tuple:

    species = lower(
        row.get(
            "species"
        )
    )

    rank = safe_int(
        row.get(
            "final_species_rank"
        ),
        default=999,
    )

    return (
        SPECIES_ORDER.get(
            species,
            99,
        ),
        rank,
        clean(
            row.get(
                "hypothesis_id"
            )
        ),
    )


def mutagenesis_sort_key(
    row: Dict[str, object],
) -> Tuple:

    species = lower(
        row.get(
            "species"
        )
    )

    rank = safe_int(
        row.get(
            "final_species_rank"
        ),
        default=999,
    )

    position = safe_int(
        row.get(
            "position"
        ),
        default=9999,
    )

    mutation_class = clean(
        row.get(
            "mutation_class"
        )
    )

    class_order = {
        "individual_resolved_candidate": 0,
        "combined_resolved_candidate": 1,
        "combined_resolved_candidate_not_directly_defined": 2,
        "exploratory_sequence_region_substitution": 3,
        "sequence_length_or_alignment_state_difference": 4,
        "isoform_specific_sequence_region_substitution": 5,
        "exploratory_sequence_region_mapping": 6,
    }

    return (
        SPECIES_ORDER.get(
            species,
            99,
        ),
        rank,
        class_order.get(
            mutation_class,
            99,
        ),
        position,
    )


# =============================================================================
# SANITY CHECKS
# =============================================================================

def sanity_checks(
    validation_rows: Sequence[Dict[str, object]],
    mutagenesis_rows: Sequence[Dict[str, object]],
) -> None:

    errors = []

    # -------------------------------------------------------------------------
    # Step 2V ranks must remain intact.
    # -------------------------------------------------------------------------

    for species in SPECIES:

        species_rows = [
            row
            for row in validation_rows
            if lower(
                row.get(
                    "species"
                )
            )
            == species
        ]

        ranks = sorted(
            safe_int(
                row.get(
                    "final_species_rank"
                )
            )
            for row in species_rows
        )

        expected = list(
            range(
                1,
                len(species_rows) + 1,
            )
        )

        if ranks != expected:

            errors.append(
                f"{species}: final species ranks are not contiguous: "
                f"{ranks}"
            )

    # -------------------------------------------------------------------------
    # Unresolved candidate-level hypotheses cannot be called structurally
    # resolved reciprocal-swap candidates.
    # -------------------------------------------------------------------------

    for row in validation_rows:

        if (
            clean(
                row.get(
                    "structural_confidence"
                )
            )
            == "unresolved_no_fixed_geometry"
            and clean(
                row.get(
                    "reciprocal_swap_recommended"
                )
            )
            == "yes"
        ):

            errors.append(
                f"{row.get('hypothesis_id')}: unresolved candidate "
                f"incorrectly marked for resolved reciprocal swap."
            )

    # -------------------------------------------------------------------------
    # Gap/alignment states must NEVER be emitted as conventional point mutation
    # strings.
    # -------------------------------------------------------------------------

    for row in mutagenesis_rows:

        event_type = clean(
            row.get(
                "comparison_event_type"
            )
        )

        forward = clean(
            row.get(
                "forward_NKG2A_to_NKG2C_mutation"
            )
        )

        reciprocal = clean(
            row.get(
                "reciprocal_NKG2C_to_NKG2A_mutation"
            )
        )

        if (
            event_type
            == "sequence_length_or_alignment_state_difference"
        ):

            if forward or reciprocal:

                errors.append(
                    f"{row.get('hypothesis_id')} position "
                    f"{row.get('position')}: gap/alignment state was "
                    f"incorrectly emitted as a point mutation."
                )

    # -------------------------------------------------------------------------
    # Identical rhesus isoform states must collapse for mutation naming.
    #
    # Example:
    #     F/F -> F
    #     T98F, not T98F/F
    # -------------------------------------------------------------------------

    for row in mutagenesis_rows:

        species = lower(
            row.get(
                "species"
            )
        )

        states = clean(
            row.get(
                "nkg2c_isoform_states"
            )
        )

        collapsed = clean(
            row.get(
                "collapsed_nkg2c_target_state"
            )
        )

        event_type = clean(
            row.get(
                "comparison_event_type"
            )
        )

        if (
            species == "rhesus"
            and "/"
            in states
            and event_type
            == "amino_acid_substitution"
        ):

            state_parts = [
                part.strip()
                for part in states.split("/")
            ]

            if (
                state_parts
                and len(
                    set(state_parts)
                )
                == 1
            ):

                expected = state_parts[0]

                if collapsed != expected:

                    errors.append(
                        f"{row.get('hypothesis_id')} position "
                        f"{row.get('position')}: identical rhesus "
                        f"isoform states {states} did not collapse "
                        f"to {expected}."
                    )

                forward = clean(
                    row.get(
                        "forward_NKG2A_to_NKG2C_mutation"
                    )
                )

                reciprocal = clean(
                    row.get(
                        "reciprocal_NKG2C_to_NKG2A_mutation"
                    )
                )

                if "/" in forward or "/" in reciprocal:

                    errors.append(
                        f"{row.get('hypothesis_id')} position "
                        f"{row.get('position')}: slash-containing "
                        f"mutation name remains after rhesus "
                        f"isoform-state collapse."
                    )

    # -------------------------------------------------------------------------
    # Multi-residue resolved candidates must have a combined design row.
    # -------------------------------------------------------------------------

    for row in validation_rows:

        if (
            clean(
                row.get(
                    "structural_confidence"
                )
            )
            == "unresolved_no_fixed_geometry"
        ):
            continue

        residues = parse_residue_labels(
            clean(
                row.get(
                    "candidate_residues"
                )
            )
        )

        if len(residues) <= 1:
            continue

        hypothesis_id = clean(
            row.get(
                "hypothesis_id"
            )
        )

        has_combined = any(
            clean(
                mutation.get(
                    "hypothesis_id"
                )
            )
            == hypothesis_id
            and clean(
                mutation.get(
                    "mutation_class"
                )
            )
            in {
                "combined_resolved_candidate",
                "combined_resolved_candidate_not_directly_defined",
            }
            for mutation in mutagenesis_rows
        )

        if not has_combined:

            errors.append(
                f"{hypothesis_id}: multi-residue resolved candidate "
                f"lacks a combined-design row."
            )

    if errors:

        raise RuntimeError(
            "Step 2W sanity checks failed:\n\n  "
            + "\n  ".join(errors)
        )


# =============================================================================
# CONSOLE REPORTING
# =============================================================================

def print_candidate(
    row: Dict[str, object],
) -> None:

    print(
        f"{row['hypothesis_id']:<7} "
        f"rank={row['final_species_rank']:<2} "
        f"{row['candidate_residues']}"
    )

    print(
        f"  Selection:      "
        f"{row['final_selection_category']}"
    )

    print(
        f"  Structure:      "
        f"{row['structural_confidence']}"
    )

    print(
        f"  Primary test:   "
        f"{row['primary_validation_test']}"
    )

    print(
        f"  Event types:    "
        f"{row['comparison_event_types']}"
    )

    print(
        f"  AA changes:     "
        f"{row['amino_acid_substitution_count']}"
    )

    print(
        f"  Alignment/gap:  "
        f"{row['alignment_state_difference_count']}"
    )

    print(
        f"  Isoform-specific: "
        f"{row['isoform_specific_state_count']}"
    )

    print(
        f"  Reciprocal swap: "
        f"{row['reciprocal_swap_recommended']}"
    )

    print(
        f"  Combined swap:   "
        f"{row['combined_swap_recommended']}"
    )

    print(
        f"  Linear peptide:  "
        f"{row['linear_peptide_test_recommended']}"
    )

    print(
        f"  Glyco control:   "
        f"{row['glycosylation_control_required']}"
    )

    print()


def print_species_summary(
    species: str,
    rows: Sequence[Dict[str, object]],
) -> None:

    section(
        f"{species.upper()} STEP 2W VALIDATION PLAN"
    )

    species_rows = [
        row
        for row in rows
        if lower(
            row.get(
                "species"
            )
        )
        == species
    ]

    species_rows.sort(
        key=candidate_sort_key
    )

    for row in species_rows:

        print_candidate(
            row
        )


def print_mutation_row(
    row: Dict[str, object],
) -> None:

    species = clean(
        row.get(
            "species"
        )
    )

    hypothesis = clean(
        row.get(
            "hypothesis_id"
        )
    )

    event_type = clean(
        row.get(
            "comparison_event_type"
        )
    )

    position = clean(
        row.get(
            "position"
        )
    )

    states = clean(
        row.get(
            "nkg2c_isoform_states"
        )
    )

    collapsed = clean(
        row.get(
            "collapsed_nkg2c_target_state"
        )
    )

    forward = clean(
        row.get(
            "forward_NKG2A_to_NKG2C_mutation"
        )
    )

    reciprocal = clean(
        row.get(
            "reciprocal_NKG2C_to_NKG2A_mutation"
        )
    )

    alignment_design = clean(
        row.get(
            "alignment_or_region_design"
        )
    )

    print(
        f"{species:<8} "
        f"{hypothesis:<7} "
        f"position={position or '-':<4} "
        f"{event_type}"
    )

    if states:

        print(
            f"  NKG2C state(s): "
            f"{states}"
        )

    if collapsed:

        print(
            f"  Collapsed state: "
            f"{collapsed}"
        )

    if forward:

        print(
            f"  Forward:         "
            f"{forward}"
        )

    if reciprocal:

        print(
            f"  Reciprocal:      "
            f"{reciprocal}"
        )

    if alignment_design:

        print(
            f"  Region design:   "
            f"{alignment_design}"
        )

    print(
        f"  Design status:   "
        f"{row['reciprocal_design_status']}"
    )

    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    banner(
        "STEP 2W - NKG2A WITHIN-SPECIES EXPERIMENTAL VALIDATION MATRIX"
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
        "Step 2V candidate rankings are preserved."
    )

    print(
        "Cross-species conservation/reactivity is NOT used."
    )

    print(
        "Step 2W adds experimental-validation logic, not new structural geometry."
    )

    print(
        "Gap/alignment states are kept separate from conventional point mutations."
    )

    print(
        "Identical rhesus NKG2C isoform states are collapsed only for mutation naming."
    )

    # -------------------------------------------------------------------------
    # Load Step 2V.
    # -------------------------------------------------------------------------

    rows = read_tsv(
        INPUT_FILE
    )

    validate_input(
        rows
    )

    print()

    print(
        f"Step 2V candidates loaded: "
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
            )
            == species
        )

        print(
            f"  {species:<8} "
            f"{count}"
        )

    # -------------------------------------------------------------------------
    # Build candidate-level validation matrix.
    # -------------------------------------------------------------------------

    validation_rows = [
        build_validation_row(
            row
        )
        for row in rows
    ]

    validation_rows.sort(
        key=candidate_sort_key
    )

    # -------------------------------------------------------------------------
    # Build mutation / sequence-mapping plan.
    # -------------------------------------------------------------------------

    mutagenesis_rows = (
        build_mutagenesis_plan(
            rows
        )
    )

    mutagenesis_rows.sort(
        key=mutagenesis_sort_key
    )

    # -------------------------------------------------------------------------
    # Sanity checks.
    # -------------------------------------------------------------------------

    sanity_checks(
        validation_rows,
        mutagenesis_rows,
    )

    # -------------------------------------------------------------------------
    # Write species-specific candidate validation plans.
    # -------------------------------------------------------------------------

    for species in SPECIES:

        species_rows = [
            row
            for row in validation_rows
            if lower(
                row.get(
                    "species"
                )
            )
            == species
        ]

        write_tsv(
            OUTPUT_FILES[
                species
            ],
            species_rows,
            OUTPUT_FIELDS,
        )

    # -------------------------------------------------------------------------
    # Write combined validation matrix.
    # -------------------------------------------------------------------------

    write_tsv(
        COMBINED_OUTPUT,
        validation_rows,
        OUTPUT_FIELDS,
    )

    # -------------------------------------------------------------------------
    # Write mutation / sequence-mapping plan.
    # -------------------------------------------------------------------------

    write_tsv(
        MUTAGENESIS_OUTPUT,
        mutagenesis_rows,
        MUTAGENESIS_FIELDS,
    )

    # -------------------------------------------------------------------------
    # Candidate-level console reports.
    # -------------------------------------------------------------------------

    for species in SPECIES:

        print_species_summary(
            species,
            validation_rows,
        )

    # -------------------------------------------------------------------------
    # Mutation-level console report.
    # -------------------------------------------------------------------------

    section(
        "MUTAGENESIS / SEQUENCE-MAPPING PLAN"
    )

    for row in mutagenesis_rows:

        print_mutation_row(
            row
        )

    # -------------------------------------------------------------------------
    # Summary statistics.
    # -------------------------------------------------------------------------

    section(
        "STEP 2W SUMMARY"
    )

    print(
        f"Candidate-level validation rows: "
        f"{len(validation_rows)}"
    )

    print(
        f"Mutation/sequence-mapping rows: "
        f"{len(mutagenesis_rows)}"
    )

    print()

    for species in SPECIES:

        species_rows = [
            row
            for row in validation_rows
            if lower(
                row.get(
                    "species"
                )
            )
            == species
        ]

        mutation_rows_species = [
            row
            for row in mutagenesis_rows
            if lower(
                row.get(
                    "species"
                )
            )
            == species
        ]

        point_mutations = sum(
            1
            for row in mutation_rows_species
            if clean(
                row.get(
                    "comparison_event_type"
                )
            )
            == "amino_acid_substitution"
        )

        alignment_states = sum(
            1
            for row in mutation_rows_species
            if clean(
                row.get(
                    "comparison_event_type"
                )
            )
            == "sequence_length_or_alignment_state_difference"
        )

        isoform_specific = sum(
            1
            for row in mutation_rows_species
            if clean(
                row.get(
                    "comparison_event_type"
                )
            )
            == "isoform_specific_amino_acid_states"
        )

        combined = sum(
            1
            for row in mutation_rows_species
            if clean(
                row.get(
                    "mutation_class"
                )
            )
            == "combined_resolved_candidate"
        )

        print(
            f"{species:<8} "
            f"candidates={len(species_rows)}  "
            f"point_mutations={point_mutations}  "
            f"alignment_states={alignment_states}  "
            f"isoform_specific={isoform_specific}  "
            f"combined_mutants={combined}"
        )

    # -------------------------------------------------------------------------
    # Key expected mutation checks.
    # -------------------------------------------------------------------------

    section(
        "KEY RESOLVED MUTATION DESIGNS"
    )

    key_hypotheses = {
        "HHYP2",
        "RHYP1",
        "RHYP2",
        "RHYP3",
        "PHYP1",
        "PHYP2",
        "PHYP3",
    }

    for row in mutagenesis_rows:

        hypothesis_id = clean(
            row.get(
                "hypothesis_id"
            )
        )

        if hypothesis_id not in key_hypotheses:
            continue

        forward = clean(
            row.get(
                "forward_NKG2A_to_NKG2C_mutation"
            )
        )

        reciprocal = clean(
            row.get(
                "reciprocal_NKG2C_to_NKG2A_mutation"
            )
        )

        if not forward:
            continue

        print(
            f"{row['species']:<8} "
            f"{hypothesis_id:<7} "
            f"{forward:<25} "
            f"<-> "
            f"{reciprocal}"
        )

    # -------------------------------------------------------------------------
    # Outputs.
    # -------------------------------------------------------------------------

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

    print(
        MUTAGENESIS_OUTPUT
    )

    print()

    print(
        "IMPORTANT INTERPRETATION RULES:"
    )

    print(
        "  1. Step 2W does NOT change Step 2V candidate rankings."
    )

    print(
        "  2. Gap/alignment states are NOT conventional amino-acid "
        "point mutations."
    )

    print(
        "  3. Rhesus F/F, L/L, E/E, Q/Q, P/P-style NKG2C states are "
        "collapsed to a single residue only for conventional mutation naming."
    )

    print(
        "  4. The original rhesus per-isoform state is retained in "
        "nkg2c_isoform_states."
    )

    print(
        "  5. Discordant rhesus isoform states, if encountered, require "
        "isoform-specific experimental designs."
    )

    print(
        "  6. Binding changes are interpretable only with matched "
        "surface-expression/folding controls."
    )

    print(
        "  7. Reciprocal gain-of-binding in NKG2C provides stronger "
        "epitope evidence than NKG2A loss-of-binding alone."
    )

    print(
        "  8. Negative peptide binding does NOT eliminate an unresolved "
        "candidate as a possible conformational/native-receptor epitope."
    )

    print(
        "  9. Glycosylation-context flags identify experimental questions; "
        "canonical sequons do not prove native glycan occupancy."
    )

    print(
        " 10. Rhesus and pigtail structural accessibility remains projected "
        "from homologous positions in human 3CDG, not directly measured "
        "macaque structures."
    )

    print(
        " 11. Binding specificity and functional blockade are separate "
        "experimental questions."
    )

    print(
        " 12. All candidates remain computational hypotheses until "
        "experimentally validated."
    )


if __name__ == "__main__":
    main()
from pathlib import Path
import csv


# =============================================================================
# STEP 2R - BUILD NKG2A WITHIN-SPECIES SPECIFICITY EVIDENCE MATRIX
# =============================================================================
#
# Goal
# ----
# Consolidate sequence and structural evidence for NKG2A-vs-NKG2C
# discrimination independently in:
#
#   human
#   rhesus macaque
#   pig-tailed macaque
#
# This step does NOT rank cross-species antibody cross-reactivity.
#
# Evidence hierarchy:
#
#   residues 113-232:
#       experimental human NKG2A structural evidence from 3CDG
#
#   residues 94-112:
#       sequence evidence retained;
#       experimental structural geometry unavailable;
#       AlphaFold confidence may be reported but must not be interpreted
#       as fixed geometry because Step 2Q found low confidence / high PAE.
#
# Rhesus structural annotations remain homologous-position evidence derived
# from human 3CDG rather than direct rhesus structural measurements.
#
# Pigtail structural annotations remain homologous-position evidence derived
# from human 3CDG rather than direct pigtail structural measurements.
# =============================================================================


ROOT = Path(__file__).resolve().parent.parent


# =============================================================================
# INPUT FILES
# =============================================================================

RANKING_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_within_species_specificity_ranking.tsv"
)

AF_CONFIDENCE_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "alphafold_nkg2a_residue_confidence.tsv"
)

AF_VALIDATION_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "alphafold_nkg2a_model_validation.tsv"
)

FOOTPRINT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_antibody_footprint_candidates.tsv"
)


# =============================================================================
# OUTPUT FILES
# =============================================================================

OUTPUT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_within_species_specificity_evidence_matrix.tsv"
)

HUMAN_OUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "human_NKG2A_specificity_candidates.tsv"
)

RHESUS_OUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "rhesus_NKG2A_specificity_candidates.tsv"
)

PIGTAIL_OUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "pigtail_NKG2A_specificity_candidates.tsv"
)


# =============================================================================
# BASIC HELPERS
# =============================================================================

def clean(value):
    """
    Normalize table values.
    """

    if value is None:
        return ""

    return str(value).strip()


def as_int(value):
    """
    Convert numeric-looking text to int.
    """

    text = clean(value)

    if not text:
        return None

    return int(float(text))


def as_float(value):
    """
    Convert numeric-looking text to float.
    """

    text = clean(value)

    if not text:
        return None

    try:
        return float(text)

    except ValueError:
        return None


def yes(value):
    """
    Interpret a yes/no field.
    """

    return clean(value).lower() == "yes"


def read_tsv(path):
    """
    Read TSV file.
    """

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
    path,
    rows,
    fieldnames,
):
    """
    Write TSV.
    """

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
        )

        writer.writeheader()
        writer.writerows(
            rows
        )


# =============================================================================
# ALPHAFOLD LOOKUPS
# =============================================================================

def build_af_confidence_lookup(rows):
    """
    full-length residue -> AlphaFold confidence row
    """

    lookup = {}

    for row in rows:

        residue = as_int(
            row.get(
                "full_length_residue"
            )
        )

        if residue is None:
            continue

        lookup[
            residue
        ] = row

    return lookup


def load_af_model_decision(rows):
    """
    Return Step 2Q model-level interpretation.
    """

    if not rows:
        return ""

    return clean(
        rows[0].get(
            "model_use_decision"
        )
    )


# =============================================================================
# FOOTPRINT LOOKUP
# =============================================================================

def parse_residue_labels(value):
    """
    Parse labels such as:
        171P,213V,214N
    into residue numbers.
    """

    output = []

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

            output.append(
                int(digits)
            )

    return output


def build_footprint_lookup(rows):
    """
    Map residue -> footprint IDs.
    """

    lookup = {}

    for row in rows:

        footprint_id = clean(
            row.get(
                "footprint_id"
            )
        )

        residues = parse_residue_labels(
            row.get(
                "residue_labels"
            )
        )

        for residue in residues:

            lookup.setdefault(
                residue,
                [],
            ).append(
                footprint_id
            )

    return lookup


# =============================================================================
# EVIDENCE CLASSIFICATION
# =============================================================================

def experimental_structure_class(row):
    """
    Summarize structural evidence from 3CDG.
    """

    status = clean(
        row.get(
            "human_3CDG_coordinate_status"
        )
    )

    if status != "resolved":

        return (
            "no_experimental_coordinate"
        )

    exposure = clean(
        row.get(
            "complex_exposure_class"
        )
    ).lower()

    interface = yes(
        row.get(
            "any_interface_contact"
        )
    )

    if (
        exposure == "exposed"
        and not interface
    ):

        return (
            "resolved_exposed_noninterface"
        )

    if (
        exposure == "exposed"
        and interface
    ):

        return (
            "resolved_exposed_interface"
        )

    if (
        exposure
        == "partially_exposed"
        and not interface
    ):

        return (
            "resolved_partially_exposed_noninterface"
        )

    if (
        exposure
        == "partially_exposed"
        and interface
    ):

        return (
            "resolved_partially_exposed_interface"
        )

    if exposure == "buried":

        return (
            "resolved_buried"
        )

    return (
        "resolved_exposure_uncertain"
    )


def alphafold_evidence_class(
    residue,
    af_row,
    af_model_decision,
):
    """
    Interpret AlphaFold evidence conservatively.

    Step 2Q established that residues 94-112 have uncertain geometry.
    Therefore pLDDT is reported as local-confidence information only.
    """

    if af_row is None:

        return (
            "not_available"
        )

    plddt = as_float(
        af_row.get(
            "plddt"
        )
    )

    confidence = clean(
        af_row.get(
            "confidence_class"
        )
    )

    if residue is None:

        return (
            "not_available"
        )

    if 94 <= residue <= 112:

        if (
            plddt is not None
            and plddt >= 70
        ):

            return (
                "local_prediction_confident_"
                "relative_geometry_uncertain"
            )

        if (
            plddt is not None
            and plddt >= 50
        ):

            return (
                "low_confidence_prediction_"
                "do_not_use_fixed_geometry"
            )

        return (
            "very_low_confidence_prediction_"
            "do_not_use_fixed_geometry"
        )

    return (
        f"{confidence}_prediction"
        if confidence
        else
        "prediction_available"
    )


# =============================================================================
# SPECIES-SPECIFIC EVIDENCE TIERS
# =============================================================================

def classify_species_candidate(
    discriminatory,
    structure_class,
    sequon_status,
):
    """
    Assign a qualitative evidence tier.

    These tiers describe candidate evidence within one species.
    They are NOT antibody validation scores.
    """

    if not discriminatory:

        return (
            "not_discriminatory"
        )

    if (
        structure_class
        == "resolved_exposed_noninterface"
    ):

        if (
            sequon_status
            in {
                "sequon_asparagine",
                "within_sequon",
            }
        ):

            return (
                "tier_1b_exposed_noninterface_"
                "glycosylation_context"
            )

        return (
            "tier_1a_exposed_noninterface"
        )

    if (
        structure_class
        == "resolved_exposed_interface"
    ):

        return (
            "tier_2_exposed_interface"
        )

    if (
        structure_class
        == "resolved_partially_exposed_noninterface"
    ):

        return (
            "tier_3_partially_exposed_noninterface"
        )

    if (
        structure_class
        == "resolved_partially_exposed_interface"
    ):

        return (
            "tier_4_partially_exposed_interface"
        )

    if (
        structure_class
        == "resolved_buried"
    ):

        return (
            "tier_5_buried"
        )

    if (
        structure_class
        == "no_experimental_coordinate"
    ):

        return (
            "unresolved_sequence_discriminator"
        )

    return (
        "structural_status_uncertain"
    )


def tier_sort_key(value):
    """
    Sorting order for reports.
    """

    order = {
        "tier_1a_exposed_noninterface": 1,
        "tier_1b_exposed_noninterface_glycosylation_context": 2,
        "tier_2_exposed_interface": 3,
        "tier_3_partially_exposed_noninterface": 4,
        "tier_4_partially_exposed_interface": 5,
        "unresolved_sequence_discriminator": 6,
        "tier_5_buried": 7,
        "structural_status_uncertain": 8,
        "not_discriminatory": 99,
    }

    return order.get(
        value,
        50,
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "STEP 2R - NKG2A WITHIN-SPECIES SPECIFICITY EVIDENCE MATRIX"
    )
    print("=" * 78)

    ranking_rows = read_tsv(
        RANKING_FILE
    )

    af_rows = read_tsv(
        AF_CONFIDENCE_FILE
    )

    af_validation_rows = read_tsv(
        AF_VALIDATION_FILE
    )

    footprint_rows = read_tsv(
        FOOTPRINT_FILE
    )

    af_lookup = (
        build_af_confidence_lookup(
            af_rows
        )
    )

    af_model_decision = (
        load_af_model_decision(
            af_validation_rows
        )
    )

    footprint_lookup = (
        build_footprint_lookup(
            footprint_rows
        )
    )

    print()
    print(
        f"Step 2O rows loaded: "
        f"{len(ranking_rows)}"
    )

    print(
        f"AlphaFold confidence rows loaded: "
        f"{len(af_rows)}"
    )

    print(
        f"Compact footprints loaded: "
        f"{len(footprint_rows)}"
    )

    print(
        "AlphaFold model interpretation:"
    )

    print(
        f"  {af_model_decision}"
    )

    output_rows = []

    for row in ranking_rows:

        human_residue = as_int(
            row.get(
                "human_NKG2A_residue"
            )
        )

        rhesus_residue = as_int(
            row.get(
                "rhesus_NKG2A_residue"
            )
        )

        pigtail_residue = as_int(
            row.get(
                "pigtail_NKG2A_residue"
            )
        )

        alignment_column = as_int(
            row.get(
                "alignment_column"
            )
        )

        # ---------------------------------------------------------------------
        # Within-species discrimination
        # ---------------------------------------------------------------------

        human_disc = yes(
            row.get(
                "human_NKG2A_vs_NKG2C_discriminatory"
            )
        )

        rhesus_iso1_disc = yes(
            row.get(
                "rhesus_NKG2A_vs_NKG2C1_discriminatory"
            )
        )

        rhesus_iso2_disc = yes(
            row.get(
                "rhesus_NKG2A_vs_NKG2C2_discriminatory"
            )
        )

        rhesus_both_disc = yes(
            row.get(
                "rhesus_NKG2A_vs_both_NKG2C_isoforms_discriminatory"
            )
        )

        pigtail_disc = yes(
            row.get(
                "pigtail_NKG2A_vs_NKG2C_discriminatory"
            )
        )

        # ---------------------------------------------------------------------
        # Experimental structure
        # ---------------------------------------------------------------------

        structure_class = (
            experimental_structure_class(
                row
            )
        )

        # ---------------------------------------------------------------------
        # AlphaFold confidence for human numbering
        # ---------------------------------------------------------------------

        af_row = None

        if human_residue is not None:

            af_row = af_lookup.get(
                human_residue
            )

        af_plddt = ""

        af_confidence_class = ""

        if af_row is not None:

            af_plddt = clean(
                af_row.get(
                    "plddt"
                )
            )

            af_confidence_class = clean(
                af_row.get(
                    "confidence_class"
                )
            )

        af_evidence = (
            alphafold_evidence_class(
                human_residue,
                af_row,
                af_model_decision,
            )
        )

        # ---------------------------------------------------------------------
        # Glycosylation
        # ---------------------------------------------------------------------

        sequon_status = clean(
            row.get(
                "human_NKG2A_canonical_sequon_status"
            )
        )

        sequon = clean(
            row.get(
                "human_NKG2A_canonical_sequon"
            )
        )

        # ---------------------------------------------------------------------
        # Footprint membership
        # ---------------------------------------------------------------------

        footprint_ids = []

        if human_residue is not None:

            footprint_ids = (
                footprint_lookup.get(
                    human_residue,
                    []
                )
            )

        # ---------------------------------------------------------------------
        # Species-specific candidate tiers
        # ---------------------------------------------------------------------

        human_evidence_tier = (
            classify_species_candidate(
                discriminatory=human_disc,
                structure_class=structure_class,
                sequon_status=sequon_status,
            )
        )

        rhesus_evidence_tier = (
            classify_species_candidate(
                discriminatory=rhesus_both_disc,
                structure_class=structure_class,
                sequon_status=sequon_status,
            )
        )

        pigtail_evidence_tier = (
            classify_species_candidate(
                discriminatory=pigtail_disc,
                structure_class=structure_class,
                sequon_status=sequon_status,
            )
        )

        # ---------------------------------------------------------------------
        # Experimental structural interpretation
        # ---------------------------------------------------------------------

        if (
            structure_class
            == "no_experimental_coordinate"
        ):

            geometry_use = (
                "sequence_candidate_only_"
                "experimental_geometry_unavailable"
            )

        else:

            geometry_use = (
                "experimental_3CDG_human_geometry"
            )

        # ---------------------------------------------------------------------
        # Overall row
        # ---------------------------------------------------------------------

        output_rows.append(
            {
                "alignment_column":
                    (
                        alignment_column
                        if alignment_column
                        is not None
                        else ""
                    ),

                # =============================================================
                # HUMAN
                # =============================================================

                "human_NKG2A_residue":
                    (
                        human_residue
                        if human_residue
                        is not None
                        else ""
                    ),

                "human_NKG2A_aa":
                    clean(
                        row.get(
                            "human_NKG2A_aa"
                        )
                    ),

                "human_NKG2C_aa":
                    clean(
                        row.get(
                            "human_NKG2C_aa"
                        )
                    ),

                "human_within_species_discriminatory":
                    (
                        "yes"
                        if human_disc
                        else "no"
                    ),

                "human_evidence_tier":
                    human_evidence_tier,

                # =============================================================
                # RHESUS
                # =============================================================

                "rhesus_NKG2A_residue":
                    (
                        rhesus_residue
                        if rhesus_residue
                        is not None
                        else ""
                    ),

                "rhesus_NKG2A_aa":
                    clean(
                        row.get(
                            "rhesus_NKG2A_aa"
                        )
                    ),

                "rhesus_NKG2C1_aa":
                    clean(
                        row.get(
                            "rhesus_NKG2C1_aa"
                        )
                    ),

                "rhesus_NKG2C2_aa":
                    clean(
                        row.get(
                            "rhesus_NKG2C2_aa"
                        )
                    ),

                "rhesus_vs_isoform1_discriminatory":
                    (
                        "yes"
                        if rhesus_iso1_disc
                        else "no"
                    ),

                "rhesus_vs_isoform2_discriminatory":
                    (
                        "yes"
                        if rhesus_iso2_disc
                        else "no"
                    ),

                "rhesus_vs_both_isoforms_discriminatory":
                    (
                        "yes"
                        if rhesus_both_disc
                        else "no"
                    ),

                "rhesus_evidence_tier":
                    rhesus_evidence_tier,

                # =============================================================
                # PIGTAIL
                # =============================================================

                "pigtail_NKG2A_residue":
                    (
                        pigtail_residue
                        if pigtail_residue
                        is not None
                        else ""
                    ),

                "pigtail_NKG2A_aa":
                    clean(
                        row.get(
                            "pigtail_NKG2A_aa"
                        )
                    ),

                "pigtail_NKG2C_aa":
                    clean(
                        row.get(
                            "pigtail_NKG2C_aa"
                        )
                    ),

                "pigtail_within_species_discriminatory":
                    (
                        "yes"
                        if pigtail_disc
                        else "no"
                    ),

                "pigtail_evidence_tier":
                    pigtail_evidence_tier,

                # =============================================================
                # EXPERIMENTAL STRUCTURE
                # =============================================================

                "experimental_structure_class":
                    structure_class,

                "human_3CDG_coordinate_status":
                    clean(
                        row.get(
                            "human_3CDG_coordinate_status"
                        )
                    ),

                "complex_sasa_A2":
                    clean(
                        row.get(
                            "complex_sasa_A2"
                        )
                    ),

                "complex_rsa":
                    clean(
                        row.get(
                            "complex_rsa"
                        )
                    ),

                "complex_exposure_class":
                    clean(
                        row.get(
                            "complex_exposure_class"
                        )
                    ),

                "isolated_nkg2a_rsa":
                    clean(
                        row.get(
                            "isolated_nkg2a_rsa"
                        )
                    ),

                "buried_surface_area_A2":
                    clean(
                        row.get(
                            "buried_surface_area_A2"
                        )
                    ),

                "contact_CD94":
                    clean(
                        row.get(
                            "contact_CD94"
                        )
                    ),

                "contact_HLA_E":
                    clean(
                        row.get(
                            "contact_HLA_E"
                        )
                    ),

                "contact_B2M":
                    clean(
                        row.get(
                            "contact_B2M"
                        )
                    ),

                "contact_peptide":
                    clean(
                        row.get(
                            "contact_peptide"
                        )
                    ),

                "any_interface_contact":
                    clean(
                        row.get(
                            "any_interface_contact"
                        )
                    ),

                # =============================================================
                # GLYCOSYLATION
                # =============================================================

                "human_NKG2A_sequon_status":
                    sequon_status,

                "human_NKG2A_sequon":
                    sequon,

                "minimum_modeled_glycan_distance_A":
                    clean(
                        row.get(
                            "minimum_modeled_glycan_distance_A"
                        )
                    ),

                "within_5A_of_modeled_glycan":
                    clean(
                        row.get(
                            "within_5A_of_modeled_glycan"
                        )
                    ),

                # =============================================================
                # ALPHAFOLD
                # =============================================================

                "alphafold_plddt":
                    af_plddt,

                "alphafold_confidence_class":
                    af_confidence_class,

                "alphafold_geometry_evidence":
                    af_evidence,

                "alphafold_model_use_decision":
                    af_model_decision,

                # =============================================================
                # FOOTPRINT / SEQUENCE CONTEXT
                # =============================================================

                "compact_footprint_ids":
                    ",".join(
                        footprint_ids
                    ),

                "sequence_classification":
                    clean(
                        row.get(
                            "sequence_classification"
                        )
                    ),

                "sequence_priority":
                    clean(
                        row.get(
                            "sequence_priority"
                        )
                    ),

                # =============================================================
                # FINAL INTERPRETATION
                # =============================================================

                "geometry_evidence_basis":
                    geometry_use,

                "structural_annotation_note":
                    (
                        "direct_human_3CDG"
                        if structure_class
                        != "no_experimental_coordinate"
                        else
                        "no_experimental_coordinates;"
                        "AlphaFold_94_112_low_confidence_high_PAE"
                    ),
            }
        )

    # =========================================================================
    # WRITE MASTER MATRIX
    # =========================================================================

    fields = [
        "alignment_column",

        "human_NKG2A_residue",
        "human_NKG2A_aa",
        "human_NKG2C_aa",
        "human_within_species_discriminatory",
        "human_evidence_tier",

        "rhesus_NKG2A_residue",
        "rhesus_NKG2A_aa",
        "rhesus_NKG2C1_aa",
        "rhesus_NKG2C2_aa",
        "rhesus_vs_isoform1_discriminatory",
        "rhesus_vs_isoform2_discriminatory",
        "rhesus_vs_both_isoforms_discriminatory",
        "rhesus_evidence_tier",

        "pigtail_NKG2A_residue",
        "pigtail_NKG2A_aa",
        "pigtail_NKG2C_aa",
        "pigtail_within_species_discriminatory",
        "pigtail_evidence_tier",

        "experimental_structure_class",
        "human_3CDG_coordinate_status",
        "complex_sasa_A2",
        "complex_rsa",
        "complex_exposure_class",
        "isolated_nkg2a_rsa",
        "buried_surface_area_A2",

        "contact_CD94",
        "contact_HLA_E",
        "contact_B2M",
        "contact_peptide",
        "any_interface_contact",

        "human_NKG2A_sequon_status",
        "human_NKG2A_sequon",
        "minimum_modeled_glycan_distance_A",
        "within_5A_of_modeled_glycan",

        "alphafold_plddt",
        "alphafold_confidence_class",
        "alphafold_geometry_evidence",
        "alphafold_model_use_decision",

        "compact_footprint_ids",

        "sequence_classification",
        "sequence_priority",

        "geometry_evidence_basis",
        "structural_annotation_note",
    ]

    write_tsv(
        OUTPUT_FILE,
        output_rows,
        fields,
    )

    # =========================================================================
    # SPECIES-SPECIFIC OUTPUTS
    # =========================================================================

    human_rows = [
        row
        for row in output_rows
        if row[
            "human_within_species_discriminatory"
        ] == "yes"
    ]

    rhesus_rows = [
        row
        for row in output_rows
        if row[
            "rhesus_vs_both_isoforms_discriminatory"
        ] == "yes"
    ]

    pigtail_rows = [
        row
        for row in output_rows
        if row[
            "pigtail_within_species_discriminatory"
        ] == "yes"
    ]

    human_rows.sort(
        key=lambda row: (
            tier_sort_key(
                row[
                    "human_evidence_tier"
                ]
            ),
            -(
                as_float(
                    row[
                        "complex_rsa"
                    ]
                )
                or -1
            ),
            as_int(
                row[
                    "human_NKG2A_residue"
                ]
            )
            or 9999,
        )
    )

    rhesus_rows.sort(
        key=lambda row: (
            tier_sort_key(
                row[
                    "rhesus_evidence_tier"
                ]
            ),
            -(
                as_float(
                    row[
                        "complex_rsa"
                    ]
                )
                or -1
            ),
            as_int(
                row[
                    "rhesus_NKG2A_residue"
                ]
            )
            or 9999,
        )
    )

    pigtail_rows.sort(
        key=lambda row: (
            tier_sort_key(
                row[
                    "pigtail_evidence_tier"
                ]
            ),
            -(
                as_float(
                    row[
                        "complex_rsa"
                    ]
                )
                or -1
            ),
            as_int(
                row[
                    "pigtail_NKG2A_residue"
                ]
            )
            or 9999,
        )
    )

    write_tsv(
        HUMAN_OUT,
        human_rows,
        fields,
    )

    write_tsv(
        RHESUS_OUT,
        rhesus_rows,
        fields,
    )

    write_tsv(
        PIGTAIL_OUT,
        pigtail_rows,
        fields,
    )

    # =========================================================================
    # REPORT
    # =========================================================================

    print()
    print("=" * 78)
    print("HUMAN WITHIN-SPECIES SPECIFICITY")
    print("=" * 78)

    for row in human_rows:

        residue = row[
            "human_NKG2A_residue"
        ]

        aa_a = row[
            "human_NKG2A_aa"
        ]

        aa_c = row[
            "human_NKG2C_aa"
        ]

        print(
            f"{str(residue):>4} "
            f"{aa_a}>{aa_c:<3} "
            f"{row['human_evidence_tier']:<42} "
            f"RSA={row['complex_rsa'] or 'NA':<7} "
            f"AF={row['alphafold_plddt'] or 'NA'}"
        )

    print()
    print("=" * 78)
    print(
        "RHESUS WITHIN-SPECIES SPECIFICITY "
        "(BOTH NKG2C ISOFORMS)"
    )
    print("=" * 78)

    for row in rhesus_rows:

        residue = row[
            "rhesus_NKG2A_residue"
        ]

        aa_a = row[
            "rhesus_NKG2A_aa"
        ]

        aa_c1 = row[
            "rhesus_NKG2C1_aa"
        ]

        aa_c2 = row[
            "rhesus_NKG2C2_aa"
        ]

        print(
            f"{str(residue):>4} "
            f"{aa_a}>{aa_c1}/{aa_c2:<5} "
            f"{row['rhesus_evidence_tier']:<42} "
            f"RSA(human3CDG)="
            f"{row['complex_rsa'] or 'NA'}"
        )

    print()
    print("=" * 78)
    print("PIGTAIL WITHIN-SPECIES SPECIFICITY")
    print("=" * 78)

    for row in pigtail_rows:

        residue = row[
            "pigtail_NKG2A_residue"
        ]

        aa_a = row[
            "pigtail_NKG2A_aa"
        ]

        aa_c = row[
            "pigtail_NKG2C_aa"
        ]

        print(
            f"{str(residue):>4} "
            f"{aa_a}>{aa_c:<3} "
            f"{row['pigtail_evidence_tier']:<42} "
            f"RSA(human3CDG)="
            f"{row['complex_rsa'] or 'NA'}"
        )

    # =========================================================================
    # EVIDENCE SUMMARY
    # =========================================================================

    print()
    print("=" * 78)
    print("EVIDENCE SUMMARY")
    print("=" * 78)

    print(
        f"\nHuman discriminatory candidates: "
        f"{len(human_rows)}"
    )

    print(
        f"Rhesus stringent discriminatory candidates: "
        f"{len(rhesus_rows)}"
    )

    print(
        f"Pigtail discriminatory candidates: "
        f"{len(pigtail_rows)}"
    )

    human_resolved = sum(
        row[
            "experimental_structure_class"
        ]
        != "no_experimental_coordinate"
        for row in human_rows
    )

    human_unresolved = (
        len(human_rows)
        - human_resolved
    )

    print()
    print(
        "Human candidates with experimental "
        f"structural coordinates: {human_resolved}"
    )

    print(
        "Human candidates without experimental "
        f"structural coordinates: {human_unresolved}"
    )

    print()
    print(
        "AlphaFold Step 2Q interpretation:"
    )

    print(
        f"  {af_model_decision}"
    )

    print()
    print("=" * 78)
    print("OUTPUTS")
    print("=" * 78)

    print()
    print(OUTPUT_FILE)
    print(HUMAN_OUT)
    print(RHESUS_OUT)
    print(PIGTAIL_OUT)

    print()
    print(
        "NOTE: human, rhesus, and pigtail are ranked "
        "independently for within-species NKG2A-vs-NKG2C specificity."
    )

    print(
        "No cross-species NKG2A conservation or antibody "
        "cross-reactivity criterion is used."
    )

    print(
        "For residues 94-112, AlphaFold confidence is retained "
        "as local prediction evidence only; fixed epitope geometry "
        "must not be inferred because Step 2Q found low confidence "
        "and high relative-position PAE."
    )

    print(
        "For rhesus and pigtail, experimental structural annotations "
        "are based on homologous positions in human 3CDG rather than "
        "direct macaque structures."
    )


if __name__ == "__main__":
    main()
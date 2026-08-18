from pathlib import Path
import csv


# =============================================================================
# STEP 2O - RANK WITHIN-SPECIES NKG2A-vs-NKG2C ANTIBODY TARGETS
# =============================================================================
#
# Purpose
# -------
# Rank NKG2A-vs-NKG2C discriminatory residues independently for:
#
#   1. Human:
#        human NKG2A vs human NKG2C
#
#   2. Rhesus macaque:
#        rhesus NKG2A vs rhesus NKG2C isoform 1
#        rhesus NKG2A vs rhesus NKG2C isoform 2
#
#      A stringent rhesus discriminator must differ from BOTH NKG2C isoforms.
#
#   3. Pigtail macaque:
#        pigtail NKG2A vs pigtail NKG2C
#
# This step does NOT rank cross-species NKG2A conservation or antibody
# cross-reactivity.
#
# Structural evidence from human NKG2A in 3CDG is used as a screening
# annotation only. Sequence-discriminatory residues lacking coordinates
# in 3CDG are retained as structurally unresolved candidates rather than
# being assigned low priority.
# =============================================================================


ROOT = Path(__file__).resolve().parent.parent


# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------

CANDIDATE_FILE = (
    ROOT
    / "results"
    / "tables"
    / "candidate_specificity_positions.tsv"
)

CONTACT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_structural_contacts.tsv"
)

SASA_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_solvent_accessibility.tsv"
)

GLYCOSYLATION_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2_ectodomain_n_glycosylation_sites.tsv"
)

GLYCAN_PROXIMITY_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_candidate_glycan_proximity.tsv"
)

FOOTPRINT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_antibody_footprint_candidates.tsv"
)


# -----------------------------------------------------------------------------
# Output file
# -----------------------------------------------------------------------------

OUTPUT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_within_species_specificity_ranking.tsv"
)


# =============================================================================
# Utility functions
# =============================================================================


def read_tsv(path):
    """Read a TSV file into a list of dictionaries."""

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


def clean(value):
    """Return a stripped string; convert None to an empty string."""

    if value is None:
        return ""

    return str(value).strip()


def is_real_residue(value):
    """
    Return True if a residue field contains an amino acid rather than
    an empty value or alignment gap.
    """

    value = clean(value)

    return value not in {
        "",
        "-",
        "NA",
        "N/A",
        "None",
        "none",
    }


def yes_no(value):
    """
    Normalize common truth-like values to yes/no.
    """

    value = clean(value).lower()

    if value in {
        "yes",
        "y",
        "true",
        "1",
    }:
        return "yes"

    return "no"


def aa_differs(aa_a, aa_c):
    """
    Compare two aligned amino-acid states.

    A gap is allowed as a discriminatory state if the NKG2A residue itself
    exists. This function therefore simply compares the normalized strings.
    """

    aa_a = clean(aa_a)
    aa_c = clean(aa_c)

    if aa_a == "":
        return False

    if aa_c == "":
        return False

    return aa_a != aa_c


def numeric_or_blank(value):
    """
    Convert a numeric-looking value to float, otherwise return None.
    """

    value = clean(value)

    if value == "":
        return None

    try:
        return float(value)
    except ValueError:
        return None


def integer_or_blank(value):
    """
    Convert an integer-looking field to int, otherwise return None.
    """

    value = clean(value)

    if value == "":
        return None

    try:
        return int(float(value))
    except ValueError:
        return None


def first_present(row, names):
    """
    Return the first existing, non-empty value among alternative column names.

    This makes the script somewhat tolerant of small naming changes in
    upstream files.
    """

    for name in names:
        if name in row:
            value = clean(row.get(name))

            if value != "":
                return value

    return ""


# =============================================================================
# Load structural annotations
# =============================================================================


def build_contact_lookup(rows):
    """
    Index structural contact rows by human NKG2A full-length residue number.
    """

    lookup = {}

    for row in rows:

        residue = integer_or_blank(
            row.get("full_length_residue")
        )

        if residue is None:
            continue

        lookup[residue] = row

    return lookup


def build_sasa_lookup(rows):
    """
    Index solvent accessibility rows by human NKG2A full-length residue.
    """

    lookup = {}

    for row in rows:

        residue = integer_or_blank(
            row.get("full_length_residue")
        )

        if residue is None:
            continue

        lookup[residue] = row

    return lookup


def build_glycan_lookup(rows):
    """
    Index modeled-glycan proximity rows by human NKG2A full-length residue.
    """

    lookup = {}

    for row in rows:

        residue = integer_or_blank(
            row.get("full_length_residue")
        )

        if residue is None:
            continue

        lookup[residue] = row

    return lookup


def build_human_sequon_lookup(rows):
    """
    Record canonical human NKG2A N-X-S/T sequon positions.

    We store all three residues belonging to each sequon so that a candidate
    can be annotated as:
        - sequon_asparagine
        - within_sequon
        - no
    """

    sequon_asparagines = {}
    sequon_members = {}

    for row in rows:

        if clean(row.get("record_id")) != "human_NKG2A":
            continue

        start = integer_or_blank(
            row.get("sequon_full_length_residue")
        )

        motif = clean(
            row.get("motif")
        )

        if start is None:
            continue

        sequon_asparagines[start] = motif

        for residue in range(
            start,
            start + 3,
        ):
            sequon_members.setdefault(
                residue,
                [],
            ).append(
                f"{start}:{motif}"
            )

    return (
        sequon_asparagines,
        sequon_members,
    )


def build_footprint_lookup(rows):
    """
    Map human NKG2A residue numbers to compact footprint IDs.

    The footprint table contains residue labels such as:
        171P,213V,214N,225I
    """

    lookup = {}

    for row in rows:

        footprint_id = clean(
            row.get("footprint_id")
        )

        residue_labels = clean(
            row.get("residue_labels")
        )

        if not footprint_id:
            continue

        if not residue_labels:
            continue

        labels = [
            item.strip()
            for item in residue_labels.split(",")
            if item.strip()
        ]

        for label in labels:

            digits = ""

            for character in label:

                if character.isdigit():
                    digits += character

                else:
                    break

            if not digits:
                continue

            residue = int(digits)

            lookup.setdefault(
                residue,
                [],
            ).append(
                footprint_id
            )

    return lookup


# =============================================================================
# Structural tier classification
# =============================================================================


def classify_structural_tier(
    resolved,
    discriminatory,
    exposure_class,
    any_interface_contact,
):
    """
    Assign a transparent structural screening tier.

    IMPORTANT:
    structurally unresolved discriminatory residues are NOT assigned a low
    priority. They remain unresolved candidates requiring another structure,
    model, or experimental accessibility assessment.
    """

    if not discriminatory:
        return "not_discriminatory_for_species"

    if not resolved:
        return "structurally_unresolved_discriminator"

    exposure_class = clean(
        exposure_class
    ).lower()

    interface = (
        clean(any_interface_contact).lower()
        == "yes"
    )

    if (
        exposure_class == "exposed"
        and not interface
    ):
        return "tier_1_exposed_noninterface_discriminator"

    if (
        exposure_class == "exposed"
        and interface
    ):
        return "tier_2_exposed_interface_discriminator"

    if (
        exposure_class == "partially_exposed"
        and not interface
    ):
        return "tier_3_partially_exposed_noninterface_discriminator"

    if (
        exposure_class == "partially_exposed"
        and interface
    ):
        return "tier_4_partially_exposed_interface_discriminator"

    if exposure_class == "buried":
        return "tier_5_buried_discriminator"

    return "resolved_structural_status_uncertain"


def tier_sort_key(tier):
    """
    Sorting order used only for console presentation.
    """

    order = {
        "tier_1_exposed_noninterface_discriminator": 1,
        "tier_2_exposed_interface_discriminator": 2,
        "tier_3_partially_exposed_noninterface_discriminator": 3,
        "tier_4_partially_exposed_interface_discriminator": 4,
        "structurally_unresolved_discriminator": 5,
        "tier_5_buried_discriminator": 6,
        "resolved_structural_status_uncertain": 7,
        "not_discriminatory_for_species": 99,
    }

    return order.get(
        tier,
        50,
    )


# =============================================================================
# Main
# =============================================================================


def main():

    print("=" * 78)
    print(
        "STEP 2O - RANK WITHIN-SPECIES NKG2A-vs-NKG2C ANTIBODY TARGETS"
    )
    print("=" * 78)

    # -------------------------------------------------------------------------
    # Load data
    # -------------------------------------------------------------------------

    candidate_rows = read_tsv(
        CANDIDATE_FILE
    )

    contact_rows = read_tsv(
        CONTACT_FILE
    )

    sasa_rows = read_tsv(
        SASA_FILE
    )

    glycosylation_rows = read_tsv(
        GLYCOSYLATION_FILE
    )

    glycan_rows = read_tsv(
        GLYCAN_PROXIMITY_FILE
    )

    footprint_rows = read_tsv(
        FOOTPRINT_FILE
    )

    print(
        f"\nSequence candidate positions loaded: "
        f"{len(candidate_rows)}"
    )

    # -------------------------------------------------------------------------
    # Build lookups
    # -------------------------------------------------------------------------

    contact_lookup = build_contact_lookup(
        contact_rows
    )

    sasa_lookup = build_sasa_lookup(
        sasa_rows
    )

    glycan_lookup = build_glycan_lookup(
        glycan_rows
    )

    (
        human_sequon_asparagines,
        human_sequon_members,
    ) = build_human_sequon_lookup(
        glycosylation_rows
    )

    footprint_lookup = build_footprint_lookup(
        footprint_rows
    )

    # -------------------------------------------------------------------------
    # Build integrated residue table
    # -------------------------------------------------------------------------

    output_rows = []

    for candidate in candidate_rows:

        alignment_column = integer_or_blank(
            candidate.get(
                "alignment_column"
            )
        )

        # ---------------------------------------------------------------------
        # Human
        # ---------------------------------------------------------------------

        human_a_residue = integer_or_blank(
            candidate.get(
                "human_NKG2A_residue"
            )
        )

        human_a_aa = clean(
            candidate.get(
                "human_NKG2A_aa"
            )
        )

        human_c_residue = integer_or_blank(
            candidate.get(
                "human_NKG2C_residue"
            )
        )

        human_c_aa = clean(
            candidate.get(
                "human_NKG2C_aa"
            )
        )

        human_discriminatory = (
            is_real_residue(human_a_aa)
            and aa_differs(
                human_a_aa,
                human_c_aa,
            )
        )

        # ---------------------------------------------------------------------
        # Rhesus
        # ---------------------------------------------------------------------

        rhesus_a_residue = integer_or_blank(
            candidate.get(
                "rhesus_NKG2A_residue"
            )
        )

        rhesus_a_aa = clean(
            candidate.get(
                "rhesus_NKG2A_aa"
            )
        )

        rhesus_c1_residue = integer_or_blank(
            candidate.get(
                "rhesus_NKG2C1_residue"
            )
        )

        rhesus_c1_aa = clean(
            candidate.get(
                "rhesus_NKG2C1_aa"
            )
        )

        rhesus_c2_residue = integer_or_blank(
            candidate.get(
                "rhesus_NKG2C2_residue"
            )
        )

        rhesus_c2_aa = clean(
            candidate.get(
                "rhesus_NKG2C2_aa"
            )
        )

        rhesus_vs_iso1 = (
            is_real_residue(rhesus_a_aa)
            and aa_differs(
                rhesus_a_aa,
                rhesus_c1_aa,
            )
        )

        rhesus_vs_iso2 = (
            is_real_residue(rhesus_a_aa)
            and aa_differs(
                rhesus_a_aa,
                rhesus_c2_aa,
            )
        )

        rhesus_vs_both = (
            rhesus_vs_iso1
            and rhesus_vs_iso2
        )

        # ---------------------------------------------------------------------
        # Pigtail
        # ---------------------------------------------------------------------

        pigtail_a_residue = integer_or_blank(
            candidate.get(
                "pigtail_NKG2A_residue"
            )
        )

        pigtail_a_aa = clean(
            candidate.get(
                "pigtail_NKG2A_aa"
            )
        )

        pigtail_c_residue = integer_or_blank(
            candidate.get(
                "pigtail_NKG2C_residue"
            )
        )

        pigtail_c_aa = clean(
            candidate.get(
                "pigtail_NKG2C_aa"
            )
        )

        pigtail_discriminatory = (
            is_real_residue(pigtail_a_aa)
            and aa_differs(
                pigtail_a_aa,
                pigtail_c_aa,
            )
        )

        # ---------------------------------------------------------------------
        # Structural annotation is based on the human 3CDG NKG2A reference.
        # Therefore it can only be attached directly when a human NKG2A
        # residue exists at this alignment position.
        # ---------------------------------------------------------------------

        contact = {}

        sasa = {}

        glycan = {}

        if human_a_residue is not None:

            contact = contact_lookup.get(
                human_a_residue,
                {},
            )

            sasa = sasa_lookup.get(
                human_a_residue,
                {},
            )

            glycan = glycan_lookup.get(
                human_a_residue,
                {},
            )

        resolved = (
            human_a_residue is not None
            and (
                human_a_residue in contact_lookup
                or human_a_residue in sasa_lookup
            )
        )

        exposure_class = first_present(
            sasa,
            [
                "complex_exposure_class",
                "exposure_class",
            ],
        )

        complex_rsa = first_present(
            sasa,
            [
                "complex_rsa",
            ],
        )

        complex_sasa = first_present(
            sasa,
            [
                "complex_sasa_A2",
            ],
        )

        isolated_rsa = first_present(
            sasa,
            [
                "isolated_nkg2a_rsa",
                "isolated_NKG2A_rsa",
            ],
        )

        buried_surface_area = first_present(
            sasa,
            [
                "buried_surface_area_A2",
            ],
        )

        contact_cd94 = yes_no(
            contact.get(
                "contact_CD94"
            )
        )

        contact_hla_e = yes_no(
            contact.get(
                "contact_HLA_E"
            )
        )

        contact_b2m = yes_no(
            contact.get(
                "contact_B2M"
            )
        )

        contact_peptide = yes_no(
            contact.get(
                "contact_peptide"
            )
        )

        any_interface_contact = (
            "yes"
            if any(
                value == "yes"
                for value in [
                    contact_cd94,
                    contact_hla_e,
                    contact_b2m,
                    contact_peptide,
                ]
            )
            else "no"
        )

        # ---------------------------------------------------------------------
        # Glycosylation annotation
        # ---------------------------------------------------------------------

        sequon_status = "no"

        sequon_description = ""

        if human_a_residue is not None:

            if human_a_residue in human_sequon_asparagines:

                sequon_status = (
                    "sequon_asparagine"
                )

                sequon_description = (
                    f"{human_a_residue}:"
                    f"{human_sequon_asparagines[human_a_residue]}"
                )

            elif human_a_residue in human_sequon_members:

                sequon_status = (
                    "within_sequon"
                )

                sequon_description = ";".join(
                    human_sequon_members[
                        human_a_residue
                    ]
                )

        minimum_glycan_distance = first_present(
            glycan,
            [
                "minimum_glycan_distance_A",
            ],
        )

        within_modeled_glycan_cutoff = first_present(
            glycan,
            [
                "within_5A_of_modeled_glycan",
            ],
        )

        if not within_modeled_glycan_cutoff:
            within_modeled_glycan_cutoff = "no"

        # ---------------------------------------------------------------------
        # Footprint memberships
        # ---------------------------------------------------------------------

        footprint_ids = ""

        if human_a_residue is not None:

            footprint_ids = ",".join(
                footprint_lookup.get(
                    human_a_residue,
                    [],
                )
            )

        # ---------------------------------------------------------------------
        # Species-specific structural classifications
        #
        # Because the structure is human NKG2A, these structural annotations
        # describe the homologous alignment position on the human reference.
        # They are not direct rhesus/pigtail structural measurements.
        # ---------------------------------------------------------------------

        human_tier = classify_structural_tier(
            resolved=resolved,
            discriminatory=human_discriminatory,
            exposure_class=exposure_class,
            any_interface_contact=any_interface_contact,
        )

        rhesus_tier = classify_structural_tier(
            resolved=resolved,
            discriminatory=rhesus_vs_both,
            exposure_class=exposure_class,
            any_interface_contact=any_interface_contact,
        )

        pigtail_tier = classify_structural_tier(
            resolved=resolved,
            discriminatory=pigtail_discriminatory,
            exposure_class=exposure_class,
            any_interface_contact=any_interface_contact,
        )

        # ---------------------------------------------------------------------
        # Preserve original sequence-analysis annotations
        # ---------------------------------------------------------------------

        original_classification = clean(
            candidate.get(
                "classification"
            )
        )

        original_sequence_priority = clean(
            candidate.get(
                "sequence_priority"
            )
        )

        # ---------------------------------------------------------------------
        # Output row
        # ---------------------------------------------------------------------

        output_rows.append(
            {
                "alignment_column":
                    alignment_column
                    if alignment_column is not None
                    else "",

                # Human comparison
                "human_NKG2A_residue":
                    human_a_residue
                    if human_a_residue is not None
                    else "",
                "human_NKG2A_aa":
                    human_a_aa,
                "human_NKG2C_residue":
                    human_c_residue
                    if human_c_residue is not None
                    else "",
                "human_NKG2C_aa":
                    human_c_aa,
                "human_NKG2A_vs_NKG2C_discriminatory":
                    "yes"
                    if human_discriminatory
                    else "no",

                # Rhesus comparison
                "rhesus_NKG2A_residue":
                    rhesus_a_residue
                    if rhesus_a_residue is not None
                    else "",
                "rhesus_NKG2A_aa":
                    rhesus_a_aa,
                "rhesus_NKG2C1_residue":
                    rhesus_c1_residue
                    if rhesus_c1_residue is not None
                    else "",
                "rhesus_NKG2C1_aa":
                    rhesus_c1_aa,
                "rhesus_NKG2C2_residue":
                    rhesus_c2_residue
                    if rhesus_c2_residue is not None
                    else "",
                "rhesus_NKG2C2_aa":
                    rhesus_c2_aa,
                "rhesus_NKG2A_vs_NKG2C1_discriminatory":
                    "yes"
                    if rhesus_vs_iso1
                    else "no",
                "rhesus_NKG2A_vs_NKG2C2_discriminatory":
                    "yes"
                    if rhesus_vs_iso2
                    else "no",
                "rhesus_NKG2A_vs_both_NKG2C_isoforms_discriminatory":
                    "yes"
                    if rhesus_vs_both
                    else "no",

                # Pigtail comparison
                "pigtail_NKG2A_residue":
                    pigtail_a_residue
                    if pigtail_a_residue is not None
                    else "",
                "pigtail_NKG2A_aa":
                    pigtail_a_aa,
                "pigtail_NKG2C_residue":
                    pigtail_c_residue
                    if pigtail_c_residue is not None
                    else "",
                "pigtail_NKG2C_aa":
                    pigtail_c_aa,
                "pigtail_NKG2A_vs_NKG2C_discriminatory":
                    "yes"
                    if pigtail_discriminatory
                    else "no",

                # Structural annotation
                "human_3CDG_coordinate_status":
                    "resolved"
                    if resolved
                    else "unresolved",
                "complex_sasa_A2":
                    complex_sasa,
                "complex_rsa":
                    complex_rsa,
                "complex_exposure_class":
                    exposure_class,
                "isolated_nkg2a_rsa":
                    isolated_rsa,
                "buried_surface_area_A2":
                    buried_surface_area,

                # Interface annotation
                "contact_CD94":
                    contact_cd94,
                "contact_HLA_E":
                    contact_hla_e,
                "contact_B2M":
                    contact_b2m,
                "contact_peptide":
                    contact_peptide,
                "any_interface_contact":
                    any_interface_contact,

                # Glycosylation annotation
                "human_NKG2A_canonical_sequon_status":
                    sequon_status,
                "human_NKG2A_canonical_sequon":
                    sequon_description,
                "minimum_modeled_glycan_distance_A":
                    minimum_glycan_distance,
                "within_5A_of_modeled_glycan":
                    within_modeled_glycan_cutoff,

                # Existing footprint membership
                "compact_footprint_ids":
                    footprint_ids,

                # Original sequence annotations
                "sequence_classification":
                    original_classification,
                "sequence_priority":
                    original_sequence_priority,

                # Species-specific screening classifications
                "human_specificity_tier":
                    human_tier,
                "rhesus_specificity_tier":
                    rhesus_tier,
                "pigtail_specificity_tier":
                    pigtail_tier,

                # Important interpretation note
                "structural_annotation_basis":
                    (
                        "human_NKG2A_3CDG"
                        if resolved
                        else
                        "no_3CDG_coordinate_for_human_alignment_position"
                    ),
            }
        )

    # -------------------------------------------------------------------------
    # Write output
    # -------------------------------------------------------------------------

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = [
        "alignment_column",

        "human_NKG2A_residue",
        "human_NKG2A_aa",
        "human_NKG2C_residue",
        "human_NKG2C_aa",
        "human_NKG2A_vs_NKG2C_discriminatory",

        "rhesus_NKG2A_residue",
        "rhesus_NKG2A_aa",
        "rhesus_NKG2C1_residue",
        "rhesus_NKG2C1_aa",
        "rhesus_NKG2C2_residue",
        "rhesus_NKG2C2_aa",
        "rhesus_NKG2A_vs_NKG2C1_discriminatory",
        "rhesus_NKG2A_vs_NKG2C2_discriminatory",
        "rhesus_NKG2A_vs_both_NKG2C_isoforms_discriminatory",

        "pigtail_NKG2A_residue",
        "pigtail_NKG2A_aa",
        "pigtail_NKG2C_residue",
        "pigtail_NKG2C_aa",
        "pigtail_NKG2A_vs_NKG2C_discriminatory",

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

        "human_NKG2A_canonical_sequon_status",
        "human_NKG2A_canonical_sequon",
        "minimum_modeled_glycan_distance_A",
        "within_5A_of_modeled_glycan",

        "compact_footprint_ids",

        "sequence_classification",
        "sequence_priority",

        "human_specificity_tier",
        "rhesus_specificity_tier",
        "pigtail_specificity_tier",

        "structural_annotation_basis",
    ]

    with OUTPUT_FILE.open(
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
            output_rows
        )

    # =========================================================================
    # Console summaries
    # =========================================================================

    human_rows = [
        row
        for row in output_rows
        if row[
            "human_NKG2A_vs_NKG2C_discriminatory"
        ] == "yes"
    ]

    rhesus_rows = [
        row
        for row in output_rows
        if row[
            "rhesus_NKG2A_vs_both_NKG2C_isoforms_discriminatory"
        ] == "yes"
    ]

    pigtail_rows = [
        row
        for row in output_rows
        if row[
            "pigtail_NKG2A_vs_NKG2C_discriminatory"
        ] == "yes"
    ]

    unresolved_rows = [
        row
        for row in output_rows
        if (
            row[
                "human_3CDG_coordinate_status"
            ] == "unresolved"
            and (
                row[
                    "human_NKG2A_vs_NKG2C_discriminatory"
                ] == "yes"
                or row[
                    "rhesus_NKG2A_vs_both_NKG2C_isoforms_discriminatory"
                ] == "yes"
                or row[
                    "pigtail_NKG2A_vs_NKG2C_discriminatory"
                ] == "yes"
            )
        )
    ]

    print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)

    print(
        f"\nHuman NKG2A-vs-NKG2C discriminatory positions: "
        f"{len(human_rows)}"
    )

    print(
        f"Rhesus NKG2A-vs-both-NKG2C-isoforms discriminatory positions: "
        f"{len(rhesus_rows)}"
    )

    print(
        f"Pigtail NKG2A-vs-NKG2C discriminatory positions: "
        f"{len(pigtail_rows)}"
    )

    print(
        f"Structurally unresolved discriminatory positions: "
        f"{len(unresolved_rows)}"
    )

    # =========================================================================
    # Human
    # =========================================================================

    print()
    print("=" * 78)
    print("HUMAN NKG2A-vs-NKG2C TARGETS")
    print("=" * 78)

    human_rows_sorted = sorted(
        human_rows,
        key=lambda row: (
            tier_sort_key(
                row[
                    "human_specificity_tier"
                ]
            ),
            -(
                numeric_or_blank(
                    row["complex_rsa"]
                )
                or -1
            ),
            integer_or_blank(
                row[
                    "human_NKG2A_residue"
                ]
            )
            or 9999,
        ),
    )

    for row in human_rows_sorted:

        residue = row[
            "human_NKG2A_residue"
        ]

        aa_a = row[
            "human_NKG2A_aa"
        ]

        aa_c = row[
            "human_NKG2C_aa"
        ]

        rsa = (
            row["complex_rsa"]
            if row["complex_rsa"]
            else "NA"
        )

        print(
            f"{str(residue):>4} "
            f"{aa_a}>{aa_c}  "
            f"RSA={rsa:<7} "
            f"interface={row['any_interface_contact']:<3}  "
            f"{row['human_specificity_tier']}"
        )

    # =========================================================================
    # Rhesus
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "RHESUS NKG2A-vs-NKG2C TARGETS "
        "(DISCRIMINATES BOTH ISOFORMS)"
    )
    print("=" * 78)

    rhesus_rows_sorted = sorted(
        rhesus_rows,
        key=lambda row: (
            tier_sort_key(
                row[
                    "rhesus_specificity_tier"
                ]
            ),
            -(
                numeric_or_blank(
                    row["complex_rsa"]
                )
                or -1
            ),
            integer_or_blank(
                row[
                    "rhesus_NKG2A_residue"
                ]
            )
            or 9999,
        ),
    )

    for row in rhesus_rows_sorted:

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

        rsa = (
            row["complex_rsa"]
            if row["complex_rsa"]
            else "NA"
        )

        print(
            f"{str(residue):>4} "
            f"{aa_a}>{aa_c1}/{aa_c2}  "
            f"RSA(human3CDG)={rsa:<7} "
            f"interface={row['any_interface_contact']:<3}  "
            f"{row['rhesus_specificity_tier']}"
        )

    # =========================================================================
    # Pigtail
    # =========================================================================

    print()
    print("=" * 78)
    print("PIGTAIL NKG2A-vs-NKG2C TARGETS")
    print("=" * 78)

    pigtail_rows_sorted = sorted(
        pigtail_rows,
        key=lambda row: (
            tier_sort_key(
                row[
                    "pigtail_specificity_tier"
                ]
            ),
            -(
                numeric_or_blank(
                    row["complex_rsa"]
                )
                or -1
            ),
            integer_or_blank(
                row[
                    "pigtail_NKG2A_residue"
                ]
            )
            or 9999,
        ),
    )

    for row in pigtail_rows_sorted:

        residue = row[
            "pigtail_NKG2A_residue"
        ]

        aa_a = row[
            "pigtail_NKG2A_aa"
        ]

        aa_c = row[
            "pigtail_NKG2C_aa"
        ]

        rsa = (
            row["complex_rsa"]
            if row["complex_rsa"]
            else "NA"
        )

        print(
            f"{str(residue):>4} "
            f"{aa_a}>{aa_c}  "
            f"RSA(human3CDG)={rsa:<7} "
            f"interface={row['any_interface_contact']:<3}  "
            f"{row['pigtail_specificity_tier']}"
        )

    # =========================================================================
    # Structurally unresolved discriminatory positions
    # =========================================================================

    print()
    print("=" * 78)
    print("STRUCTURALLY UNRESOLVED DISCRIMINATORY POSITIONS")
    print("=" * 78)

    if not unresolved_rows:

        print("None.")

    else:

        for row in sorted(
            unresolved_rows,
            key=lambda row: (
                integer_or_blank(
                    row[
                        "alignment_column"
                    ]
                )
                or 9999
            ),
        ):

            species = []

            if (
                row[
                    "human_NKG2A_vs_NKG2C_discriminatory"
                ]
                == "yes"
            ):
                species.append("human")

            if (
                row[
                    "rhesus_NKG2A_vs_both_NKG2C_isoforms_discriminatory"
                ]
                == "yes"
            ):
                species.append("rhesus")

            if (
                row[
                    "pigtail_NKG2A_vs_NKG2C_discriminatory"
                ]
                == "yes"
            ):
                species.append("pigtail")

            human_label = ""

            if row["human_NKG2A_residue"]:

                human_label = (
                    f"{row['human_NKG2A_residue']}"
                    f"{row['human_NKG2A_aa']}"
                )

            else:

                human_label = (
                    f"alignment_column_"
                    f"{row['alignment_column']}"
                )

            print(
                f"{human_label:<22} "
                f"species={','.join(species):<22} "
                f"{row['sequence_classification']}"
            )

    # =========================================================================
    # Output
    # =========================================================================

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)

    print()
    print(OUTPUT_FILE)

    print()
    print(
        "NOTE: species are evaluated independently for within-species "
        "NKG2A-vs-NKG2C discrimination."
    )

    print(
        "Rhesus stringent specificity requires the rhesus NKG2A state "
        "to differ from BOTH rhesus NKG2C isoforms."
    )

    print(
        "Structural exposure/interface annotations are derived from "
        "human NKG2A in 3CDG and are therefore homologous-position "
        "screening evidence for macaque residues, not direct macaque "
        "structural measurements."
    )

    print(
        "Structurally unresolved discriminatory positions are retained "
        "as candidates and are not interpreted as low-priority residues."
    )

    print(
        "No cross-species NKG2A conservation or antibody cross-reactivity "
        "criterion is used in this ranking."
    )


if __name__ == "__main__":
    main()
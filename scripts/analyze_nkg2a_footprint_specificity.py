from pathlib import Path
import csv


# =============================================================================
# STEP 2N - NKG2A FOOTPRINT TARGET CONSERVATION VS NKG2C SEPARATION
# =============================================================================

ROOT = Path(__file__).resolve().parent.parent

FOOTPRINT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_antibody_footprint_candidates.tsv"
)

INTEGRATION_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_epitope_candidate_integration.tsv"
)

OUTPUT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_footprint_cross_species_specificity.tsv"
)


# =============================================================================
# BASIC HELPERS
# =============================================================================

def read_tsv(path):
    """
    Read a tab-delimited file into a list of dictionaries.
    """

    if not path.exists():
        raise FileNotFoundError(
            f"Required file not found:\n{path}"
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
    """
    Normalize a table value to a stripped string.
    """

    if value is None:
        return ""

    return str(value).strip()


def as_int(value):
    """
    Convert a numeric table value to int.

    Blank values return None.
    """

    text = clean(value)

    if not text:
        return None

    return int(float(text))


def parse_residue_labels(value):
    """
    Convert a footprint residue-label field such as:

        171P,213V,214N,225I

    into:

        [171, 213, 214, 225]

    This deliberately extracts the leading residue number
    rather than assuming the field contains pure integers.
    """

    residues = []

    for label in clean(value).split(","):

        label = label.strip()

        if not label:
            continue

        digits = ""

        for char in label:

            if char.isdigit():
                digits += char

            else:
                break

        if digits:
            residues.append(
                int(digits)
            )

    return sorted(
        set(residues)
    )


def same_valid_aa(*values):
    """
    Return True if every supplied amino-acid state:

    1. is present,
    2. is not a gap,
    3. is identical.
    """

    cleaned = [
        clean(value)
        for value in values
    ]

    if any(
        value in {"", "-"}
        for value in cleaned
    ):
        return False

    return len(
        set(cleaned)
    ) == 1


# =============================================================================
# SPECIES / RECEPTOR STATE HELPERS
# =============================================================================

def target_state(row, species):
    """
    Return the NKG2A amino acid for the requested species.
    """

    columns = {
        "human": "human_NKG2A_aa",
        "rhesus": "rhesus_NKG2A_aa",
        "pigtail": "pigtail_NKG2A_aa",
    }

    if species not in columns:
        raise ValueError(
            f"Unknown species: {species}"
        )

    return clean(
        row.get(
            columns[species]
        )
    )


def negative_states(row, species):
    """
    Return the corresponding NKG2C amino-acid state(s).

    Human:
        one NKG2C sequence

    Rhesus:
        two NKG2C isoforms

    Pigtail:
        one NKG2C sequence
    """

    if species == "human":

        return [
            clean(
                row.get(
                    "human_NKG2C_aa"
                )
            )
        ]

    if species == "rhesus":

        return [
            clean(
                row.get(
                    "rhesus_NKG2C1_aa"
                )
            ),
            clean(
                row.get(
                    "rhesus_NKG2C2_aa"
                )
            ),
        ]

    if species == "pigtail":

        return [
            clean(
                row.get(
                    "pigtail_NKG2C_aa"
                )
            )
        ]

    raise ValueError(
        f"Unknown species: {species}"
    )


def target_differs_from_all_negatives(
    target_aa,
    negative_aas,
):
    """
    Return True when the NKG2A amino acid differs from
    every available corresponding NKG2C state.

    For rhesus, this requires discrimination from BOTH
    NKG2C isoforms.
    """

    target_aa = clean(
        target_aa
    )

    if target_aa in {"", "-"}:
        return False

    valid_negatives = [
        clean(aa)
        for aa in negative_aas
        if clean(aa) not in {"", "-"}
    ]

    if not valid_negatives:
        return False

    return all(
        target_aa != aa
        for aa in valid_negatives
    )


# =============================================================================
# SIGNATURE GENERATION
# =============================================================================

def sequence_signature(
    residues,
    lookup,
    species,
    receptor,
):
    """
    Produce a compact residue signature.

    Examples:

        171P|213V|214N|225I

    For rhesus NKG2C, the two isoforms are represented
    at each position as:

        171A/T
    """

    pieces = []

    for residue in residues:

        row = lookup.get(
            residue
        )

        if row is None:
            pieces.append(
                f"{residue}?"
            )
            continue

        if receptor == "NKG2A":

            aa = target_state(
                row,
                species,
            )

            if not aa:
                aa = "?"

        elif receptor == "NKG2C":

            negatives = negative_states(
                row,
                species,
            )

            if species == "rhesus":

                formatted = [
                    aa if aa else "?"
                    for aa in negatives
                ]

                aa = "/".join(
                    formatted
                )

            else:

                if negatives:
                    aa = (
                        negatives[0]
                        if negatives[0]
                        else "?"
                    )
                else:
                    aa = "?"

        else:

            raise ValueError(
                f"Unknown receptor: {receptor}"
            )

        pieces.append(
            f"{residue}{aa}"
        )

    return "|".join(
        pieces
    )


# =============================================================================
# FOOTPRINT CLASSIFICATION
# =============================================================================

def classify_footprint(
    target_conserved_count,
    human_discrimination,
    rhesus_discrimination,
    pigtail_discrimination,
    residue_count,
):
    """
    Assign a transparent screening classification.

    Important:
    This does NOT claim experimental antibody specificity.

    The classification separates two ideas:

    1. conservation of the NKG2A target surface across species
    2. separation of NKG2A from NKG2C within each species
    """

    all_species_discrimination = (
        human_discrimination > 0
        and rhesus_discrimination > 0
        and pigtail_discrimination > 0
    )

    if (
        all_species_discrimination
        and residue_count > 0
        and target_conserved_count == residue_count
    ):
        return (
            "strong_cross_species_target_conservation_"
            "with_NKG2C_separation"
        )

    if (
        all_species_discrimination
        and target_conserved_count > 0
    ):
        return (
            "mixed_cross_species_target_conservation_"
            "with_NKG2C_separation"
        )

    if all_species_discrimination:
        return (
            "cross_species_discrimination_"
            "but_target_surface_varies"
        )

    if (
        rhesus_discrimination > 0
        and pigtail_discrimination > 0
    ):
        return (
            "macaque_NKG2A_specificity_candidate"
        )

    if human_discrimination > 0:
        return (
            "human_NKG2A_specificity_candidate"
        )

    if pigtail_discrimination > 0:
        return (
            "pigtail_NKG2A_specificity_candidate"
        )

    if rhesus_discrimination > 0:
        return (
            "rhesus_NKG2A_specificity_candidate"
        )

    return (
        "weak_NKG2A_NKG2C_separation"
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "STEP 2N - NKG2A FOOTPRINT TARGET CONSERVATION "
        "VS NKG2C SEPARATION"
    )
    print("=" * 78)

    footprints = read_tsv(
        FOOTPRINT_FILE
    )

    integration = read_tsv(
        INTEGRATION_FILE
    )

    # -------------------------------------------------------------------------
    # Build lookup keyed by human NKG2A full-length residue number.
    #
    # The compact structural footprints are defined using human NKG2A
    # full-length numbering from the 3CDG structural mapping.
    # -------------------------------------------------------------------------

    residue_lookup = {}

    for row in integration:

        residue = as_int(
            row.get(
                "human_NKG2A_residue"
            )
        )

        if residue is None:
            continue

        residue_lookup[
            residue
        ] = row

    print(
        f"\nFootprints loaded: "
        f"{len(footprints)}"
    )

    print(
        "Integrated human NKG2A residue candidates loaded: "
        f"{len(residue_lookup)}"
    )

    output_rows = []

    # -------------------------------------------------------------------------
    # Analyze each compact candidate footprint
    # -------------------------------------------------------------------------

    for footprint in footprints:

        footprint_id = clean(
            footprint.get(
                "footprint_id"
            )
        )

        residue_labels = clean(
            footprint.get(
                "residue_labels"
            )
        )

        residues = parse_residue_labels(
            residue_labels
        )

        if not residues:
            print(
                f"\nWARNING: {footprint_id} contains "
                "no parseable residue labels."
            )
            continue

        human_disc = 0
        rhesus_disc = 0
        pigtail_disc = 0

        conserved_all_three = 0
        conserved_macaques = 0

        varying_target_positions = []

        human_disc_positions = []
        rhesus_disc_positions = []
        pigtail_disc_positions = []

        missing_lookup_positions = []

        # ---------------------------------------------------------------------
        # Examine every discriminatory residue included in the footprint
        # ---------------------------------------------------------------------

        for residue in residues:

            row = residue_lookup.get(
                residue
            )

            if row is None:

                missing_lookup_positions.append(
                    residue
                )

                continue

            human_a = target_state(
                row,
                "human",
            )

            rhesus_a = target_state(
                row,
                "rhesus",
            )

            pigtail_a = target_state(
                row,
                "pigtail",
            )

            # -----------------------------------------------------------------
            # Positive-target conservation:
            # Are the NKG2A amino acids identical across species?
            # -----------------------------------------------------------------

            if same_valid_aa(
                human_a,
                rhesus_a,
                pigtail_a,
            ):

                conserved_all_three += 1

            else:

                varying_target_positions.append(
                    residue
                )

            if same_valid_aa(
                rhesus_a,
                pigtail_a,
            ):

                conserved_macaques += 1

            # -----------------------------------------------------------------
            # Negative-target separation:
            # Does NKG2A differ from corresponding NKG2C?
            # -----------------------------------------------------------------

            if target_differs_from_all_negatives(
                human_a,
                negative_states(
                    row,
                    "human",
                ),
            ):

                human_disc += 1

                human_disc_positions.append(
                    residue
                )

            if target_differs_from_all_negatives(
                rhesus_a,
                negative_states(
                    row,
                    "rhesus",
                ),
            ):

                rhesus_disc += 1

                rhesus_disc_positions.append(
                    residue
                )

            if target_differs_from_all_negatives(
                pigtail_a,
                negative_states(
                    row,
                    "pigtail",
                ),
            ):

                pigtail_disc += 1

                pigtail_disc_positions.append(
                    residue
                )

        residue_count = len(
            residues
        )

        mapped_residue_count = (
            residue_count
            - len(
                missing_lookup_positions
            )
        )

        classification = classify_footprint(
            conserved_all_three,
            human_disc,
            rhesus_disc,
            pigtail_disc,
            mapped_residue_count,
        )

        cross_species_discrimination = (
            "yes"
            if (
                human_disc > 0
                and rhesus_disc > 0
                and pigtail_disc > 0
            )
            else "no"
        )

        complete_target_conservation = (
            "yes"
            if (
                mapped_residue_count > 0
                and conserved_all_three
                == mapped_residue_count
            )
            else "no"
        )

        output_rows.append(
            {
                "footprint_id":
                    footprint_id,

                "residue_count":
                    residue_count,

                "mapped_residue_count":
                    mapped_residue_count,

                "residue_labels":
                    residue_labels,

                "residues":
                    ",".join(
                        str(x)
                        for x in residues
                    ),

                "missing_lookup_positions":
                    ",".join(
                        str(x)
                        for x
                        in missing_lookup_positions
                    ),

                # -------------------------------------------------------------
                # Positive-target signatures
                # -------------------------------------------------------------

                "human_NKG2A_signature":
                    sequence_signature(
                        residues,
                        residue_lookup,
                        "human",
                        "NKG2A",
                    ),

                "rhesus_NKG2A_signature":
                    sequence_signature(
                        residues,
                        residue_lookup,
                        "rhesus",
                        "NKG2A",
                    ),

                "pigtail_NKG2A_signature":
                    sequence_signature(
                        residues,
                        residue_lookup,
                        "pigtail",
                        "NKG2A",
                    ),

                # -------------------------------------------------------------
                # Negative-target signatures
                # -------------------------------------------------------------

                "human_NKG2C_signature":
                    sequence_signature(
                        residues,
                        residue_lookup,
                        "human",
                        "NKG2C",
                    ),

                "rhesus_NKG2C_signature":
                    sequence_signature(
                        residues,
                        residue_lookup,
                        "rhesus",
                        "NKG2C",
                    ),

                "pigtail_NKG2C_signature":
                    sequence_signature(
                        residues,
                        residue_lookup,
                        "pigtail",
                        "NKG2C",
                    ),

                # -------------------------------------------------------------
                # Target conservation
                # -------------------------------------------------------------

                "NKG2A_positions_conserved_all_three":
                    conserved_all_three,

                "NKG2A_positions_conserved_rhesus_pigtail":
                    conserved_macaques,

                "NKG2A_positions_varying_across_species":
                    ",".join(
                        str(x)
                        for x
                        in varying_target_positions
                    ),

                "complete_NKG2A_target_conservation":
                    complete_target_conservation,

                # -------------------------------------------------------------
                # NKG2A vs NKG2C discrimination
                # -------------------------------------------------------------

                "human_NKG2A_vs_NKG2C_discriminators":
                    human_disc,

                "human_discriminatory_positions":
                    ",".join(
                        str(x)
                        for x
                        in human_disc_positions
                    ),

                "rhesus_NKG2A_vs_both_NKG2C_discriminators":
                    rhesus_disc,

                "rhesus_discriminatory_positions":
                    ",".join(
                        str(x)
                        for x
                        in rhesus_disc_positions
                    ),

                "pigtail_NKG2A_vs_NKG2C_discriminators":
                    pigtail_disc,

                "pigtail_discriminatory_positions":
                    ",".join(
                        str(x)
                        for x
                        in pigtail_disc_positions
                    ),

                "footprint_cross_species_NKG2A_vs_NKG2C_discrimination":
                    cross_species_discrimination,

                # -------------------------------------------------------------
                # Overall screening classification
                # -------------------------------------------------------------

                "cross_species_specificity_class":
                    classification,

                "original_footprint_priority":
                    clean(
                        footprint.get(
                            "footprint_priority"
                        )
                    ),

                "exposed_non_interface_residues":
                    clean(
                        footprint.get(
                            "exposed_non_interface_residues"
                        )
                    ),

                "interface_residues":
                    clean(
                        footprint.get(
                            "interface_residues"
                        )
                    ),

                "mean_complex_rsa":
                    clean(
                        footprint.get(
                            "mean_complex_rsa"
                        )
                    ),
            }
        )

    # =========================================================================
    # WRITE OUTPUT
    # =========================================================================

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fields = [
        "footprint_id",
        "residue_count",
        "mapped_residue_count",
        "residue_labels",
        "residues",
        "missing_lookup_positions",

        "human_NKG2A_signature",
        "rhesus_NKG2A_signature",
        "pigtail_NKG2A_signature",

        "human_NKG2C_signature",
        "rhesus_NKG2C_signature",
        "pigtail_NKG2C_signature",

        "NKG2A_positions_conserved_all_three",
        "NKG2A_positions_conserved_rhesus_pigtail",
        "NKG2A_positions_varying_across_species",
        "complete_NKG2A_target_conservation",

        "human_NKG2A_vs_NKG2C_discriminators",
        "human_discriminatory_positions",

        "rhesus_NKG2A_vs_both_NKG2C_discriminators",
        "rhesus_discriminatory_positions",

        "pigtail_NKG2A_vs_NKG2C_discriminators",
        "pigtail_discriminatory_positions",

        "footprint_cross_species_NKG2A_vs_NKG2C_discrimination",

        "cross_species_specificity_class",

        "original_footprint_priority",
        "exposed_non_interface_residues",
        "interface_residues",
        "mean_complex_rsa",
    ]

    with OUTPUT_FILE.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        writer.writerows(
            output_rows
        )

    # =========================================================================
    # TERMINAL SUMMARY
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "FOOTPRINT CROSS-SPECIES SPECIFICITY"
    )
    print("=" * 78)

    for row in output_rows:

        print()

        print(
            f"{row['footprint_id']:<13} "
            f"{row['residue_labels']}"
        )

        print(
            "  NKG2A conservation: "
            f"{row['NKG2A_positions_conserved_all_three']}"
            f"/{row['mapped_residue_count']} "
            "mapped positions identical across "
            "human/rhesus/pigtail"
        )

        print(
            "  NKG2A-vs-NKG2C discrimination: "
            f"H="
            f"{row['human_NKG2A_vs_NKG2C_discriminators']}  "
            f"R="
            f"{row['rhesus_NKG2A_vs_both_NKG2C_discriminators']}  "
            f"P="
            f"{row['pigtail_NKG2A_vs_NKG2C_discriminators']}"
        )

        print(
            "  Cross-species discrimination: "
            f"{row['footprint_cross_species_NKG2A_vs_NKG2C_discrimination']}"
        )

        print(
            "  Class: "
            f"{row['cross_species_specificity_class']}"
        )

    # =========================================================================
    # HIGHLIGHT CROSS-SPECIES CANDIDATES
    # =========================================================================

    cross_species_rows = [
        row
        for row in output_rows
        if (
            row[
                "footprint_cross_species_NKG2A_vs_NKG2C_discrimination"
            ]
            == "yes"
        )
    ]

    print()
    print("=" * 78)
    print(
        "CROSS-SPECIES NKG2A-vs-NKG2C "
        "DISCRIMINATORY FOOTPRINTS"
    )
    print("=" * 78)

    if cross_species_rows:

        for row in cross_species_rows:

            print()

            print(
                f"{row['footprint_id']}: "
                f"{row['residue_labels']}"
            )

            print(
                "  Human NKG2A:   "
                f"{row['human_NKG2A_signature']}"
            )

            print(
                "  Human NKG2C:   "
                f"{row['human_NKG2C_signature']}"
            )

            print(
                "  Rhesus NKG2A:  "
                f"{row['rhesus_NKG2A_signature']}"
            )

            print(
                "  Rhesus NKG2C:  "
                f"{row['rhesus_NKG2C_signature']}"
            )

            print(
                "  Pigtail NKG2A: "
                f"{row['pigtail_NKG2A_signature']}"
            )

            print(
                "  Pigtail NKG2C: "
                f"{row['pigtail_NKG2C_signature']}"
            )

            print(
                "  Classification: "
                f"{row['cross_species_specificity_class']}"
            )

    else:

        print()
        print(
            "No compact footprint currently contains "
            "NKG2A-vs-NKG2C discrimination in all "
            "three species."
        )

    # =========================================================================
    # OUTPUT
    # =========================================================================

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)

    print()
    print(OUTPUT_FILE)

    print()

    print(
        "NOTE: this step separates positive-target NKG2A "
        "conservation from negative-target NKG2C separation."
    )

    print(
        "A cross-species discriminatory footprint means the "
        "local candidate neighborhood contains at least one "
        "NKG2A-vs-NKG2C discriminator for human, rhesus, "
        "and pigtail."
    )

    print(
        "It does NOT require the same individual residue to "
        "provide discrimination in every species."
    )

    print(
        "The analysis is still restricted to discriminatory "
        "candidate residues contained in each compact footprint; "
        "it is not yet a complete physical antibody-contact "
        "surface model."
    )


if __name__ == "__main__":
    main()
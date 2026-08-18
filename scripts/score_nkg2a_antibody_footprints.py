from pathlib import Path
import csv
from collections import defaultdict


ROOT = Path(__file__).resolve().parent.parent

INTEGRATION = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_epitope_candidate_integration.tsv"
)

DISTANCES = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_candidate_spatial_distances.tsv"
)

OUTPUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_antibody_footprint_candidates.tsv"
)

NEIGHBOR_CUTOFF_A = 8.0


def read_tsv(path):

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

    if value is None:
        return ""

    return str(value).strip()


def yes(value):

    return (
        clean(value).lower()
        == "yes"
    )


def as_float(value):

    value = clean(value)

    if not value:
        return None

    try:
        return float(value)

    except ValueError:
        return None


def as_int(value):

    value = clean(value)

    if not value:
        return None

    try:
        return int(value)

    except ValueError:
        return int(float(value))


def nkg2a_label(row, residue):

    aa = clean(
        row.get(
            "human_NKG2A_aa"
        )
    )

    return f"{residue}{aa}"


def species_discrimination(row):
    """
    Read the actual Step 2K column names.

    Returns within-species NKG2A-vs-NKG2C
    discrimination for human, rhesus, and
    pig-tailed macaque.
    """

    human = yes(
        row.get(
            "human_NKG2A_vs_NKG2C_diff"
        )
    )

    rhesus_c1 = yes(
        row.get(
            "rhesus_NKG2A_vs_NKG2C1_diff"
        )
    )

    rhesus_c2 = yes(
        row.get(
            "rhesus_NKG2A_vs_NKG2C2_diff"
        )
    )

    pigtail = yes(
        row.get(
            "pigtail_NKG2A_vs_NKG2C_diff"
        )
    )

    return {
        "human": human,
        "rhesus_c1": rhesus_c1,
        "rhesus_c2": rhesus_c2,
        "rhesus_both": (
            rhesus_c1
            and rhesus_c2
        ),
        "pigtail": pigtail,
    }


def main():

    print("=" * 78)
    print(
        "STEP 2M - COMPACT NKG2A ANTIBODY FOOTPRINT ANALYSIS"
    )
    print("=" * 78)

    integration_rows = read_tsv(
        INTEGRATION
    )

    distance_rows = read_tsv(
        DISTANCES
    )

    # ---------------------------------------------------------------------
    # Structurally resolved discriminatory candidates
    # ---------------------------------------------------------------------

    candidates = {}

    for row in integration_rows:

        residue = as_int(
            row.get(
                "human_NKG2A_residue"
            )
        )

        if residue is None:
            continue

        rsa = as_float(
            row.get(
                "complex_rsa"
            )
        )

        # No RSA means the human NKG2A residue
        # is not resolved in 3CDG.
        if rsa is None:
            continue

        candidates[
            residue
        ] = row

    print(
        f"\nStructurally resolved candidate residues: "
        f"{len(candidates)}"
    )

    # ---------------------------------------------------------------------
    # QC discrimination counts BEFORE footprint construction
    # ---------------------------------------------------------------------

    human_total = 0
    rhesus_total = 0
    pigtail_total = 0

    for row in candidates.values():

        disc = species_discrimination(
            row
        )

        if disc["human"]:
            human_total += 1

        if disc["rhesus_both"]:
            rhesus_total += 1

        if disc["pigtail"]:
            pigtail_total += 1

    print()
    print(
        "Resolved discriminatory residues:"
    )

    print(
        f"  Human NKG2A vs NKG2C: "
        f"{human_total}"
    )

    print(
        f"  Rhesus NKG2A vs both NKG2C isoforms: "
        f"{rhesus_total}"
    )

    print(
        f"  Pigtail NKG2A vs NKG2C: "
        f"{pigtail_total}"
    )

    # ---------------------------------------------------------------------
    # Candidate-candidate distance graph
    # ---------------------------------------------------------------------

    neighbors = defaultdict(
        dict
    )

    for residue in candidates:

        neighbors[
            residue
        ][residue] = 0.0

    for row in distance_rows:

        r1 = as_int(
            row.get(
                "residue_1"
            )
        )

        r2 = as_int(
            row.get(
                "residue_2"
            )
        )

        distance = as_float(
            row.get(
                "minimum_heavy_atom_distance_A"
            )
        )

        if (
            r1 is None
            or r2 is None
            or distance is None
        ):
            continue

        if (
            r1 not in candidates
            or r2 not in candidates
        ):
            continue

        neighbors[
            r1
        ][r2] = distance

        neighbors[
            r2
        ][r1] = distance

    # ---------------------------------------------------------------------
    # Construct compact neighborhoods
    # ---------------------------------------------------------------------

    footprints = {}

    for center in sorted(
        candidates
    ):

        members = {
            center
        }

        for (
            other,
            distance,
        ) in neighbors[
            center
        ].items():

            if (
                distance
                <= NEIGHBOR_CUTOFF_A
            ):
                members.add(
                    other
                )

        key = tuple(
            sorted(members)
        )

        footprints.setdefault(
            key,
            [],
        ).append(
            center
        )

    print(
        f"\nUnique compact candidate neighborhoods: "
        f"{len(footprints)}"
    )

    # ---------------------------------------------------------------------
    # Summarize footprints
    # ---------------------------------------------------------------------

    output_rows = []

    for index, (
        member_tuple,
        centers,
    ) in enumerate(
        sorted(
            footprints.items(),
            key=lambda item: (
                -len(item[0]),
                item[0],
            ),
        ),
        start=1,
    ):

        member_rows = [
            candidates[
                residue
            ]
            for residue
            in member_tuple
        ]

        labels = [
            nkg2a_label(
                candidates[
                    residue
                ],
                residue,
            )
            for residue
            in member_tuple
        ]

        human_disc = 0
        rhesus_both_disc = 0
        pigtail_disc = 0

        pan_species_single_residue = 0

        exposed_non_interface = 0
        interface = 0

        conserved_all_species = 0
        conserved_macaques = 0

        rsas = []

        for row in member_rows:

            disc = (
                species_discrimination(
                    row
                )
            )

            if disc["human"]:
                human_disc += 1

            if disc[
                "rhesus_both"
            ]:
                rhesus_both_disc += 1

            if disc[
                "pigtail"
            ]:
                pigtail_disc += 1

            if (
                disc["human"]
                and disc[
                    "rhesus_both"
                ]
                and disc[
                    "pigtail"
                ]
            ):
                pan_species_single_residue += 1

            if yes(
                row.get(
                    "exposed_non_interface"
                )
            ):
                exposed_non_interface += 1

            if yes(
                row.get(
                    "any_interface_contact"
                )
            ):
                interface += 1

            if yes(
                row.get(
                    "NKG2A_state_conserved_"
                    "human_rhesus_pigtail"
                )
            ):
                conserved_all_species += 1

            if yes(
                row.get(
                    "NKG2A_state_conserved_"
                    "rhesus_pigtail"
                )
            ):
                conserved_macaques += 1

            rsa = as_float(
                row.get(
                    "complex_rsa"
                )
            )

            if rsa is not None:
                rsas.append(
                    rsa
                )

        mean_rsa = (
            sum(rsas)
            / len(rsas)
            if rsas
            else None
        )

        # -----------------------------------------------------------------
        # Does the footprint as a WHOLE contain discrimination
        # for all three species?
        #
        # These discriminatory residues do not have to be the
        # same residue in all species.
        # -----------------------------------------------------------------

        footprint_cross_species = (
            human_disc >= 1
            and rhesus_both_disc >= 1
            and pigtail_disc >= 1
        )

        footprint_macaque = (
            rhesus_both_disc >= 1
            and pigtail_disc >= 1
        )

        # -----------------------------------------------------------------
        # Transparent footprint priority
        # -----------------------------------------------------------------

        if (
            footprint_cross_species
            and exposed_non_interface >= 1
            and interface == 0
        ):

            footprint_priority = (
                "cross_species_NKG2A_"
                "discriminatory_noninterface"
            )

        elif (
            footprint_cross_species
            and exposed_non_interface >= 1
        ):

            footprint_priority = (
                "cross_species_NKG2A_"
                "discriminatory_mixed_interface"
            )

        elif (
            footprint_macaque
            and exposed_non_interface >= 1
            and interface == 0
        ):

            footprint_priority = (
                "macaque_NKG2A_"
                "discriminatory_noninterface"
            )

        elif (
            pigtail_disc >= 1
            and exposed_non_interface >= 1
            and interface == 0
        ):

            footprint_priority = (
                "pigtail_NKG2A_"
                "discriminatory_noninterface"
            )

        elif (
            exposed_non_interface >= 1
            and interface == 0
        ):

            footprint_priority = (
                "accessible_low_cross_species_"
                "discrimination"
            )

        elif interface >= 1:

            footprint_priority = (
                "interface_containing"
            )

        else:

            footprint_priority = (
                "low_structural_priority"
            )

        output_rows.append(
            {
                "footprint_id":
                    f"footprint_{index}",

                "center_residues":
                    ",".join(
                        str(x)
                        for x in centers
                    ),

                "residue_count":
                    len(
                        member_tuple
                    ),

                "residue_labels":
                    ",".join(
                        labels
                    ),

                "human_discriminatory_residues":
                    human_disc,

                "rhesus_both_isoforms_"
                "discriminatory_residues":
                    rhesus_both_disc,

                "pigtail_discriminatory_residues":
                    pigtail_disc,

                "single_residue_pan_species_"
                "discriminators":
                    pan_species_single_residue,

                "footprint_has_human_"
                "discrimination":
                    (
                        "yes"
                        if human_disc >= 1
                        else "no"
                    ),

                "footprint_has_rhesus_"
                "discrimination":
                    (
                        "yes"
                        if rhesus_both_disc >= 1
                        else "no"
                    ),

                "footprint_has_pigtail_"
                "discrimination":
                    (
                        "yes"
                        if pigtail_disc >= 1
                        else "no"
                    ),

                "footprint_cross_species_"
                "discrimination":
                    (
                        "yes"
                        if footprint_cross_species
                        else "no"
                    ),

                "cross_species_conserved_"
                "NKG2A_residues":
                    conserved_all_species,

                "macaque_conserved_"
                "NKG2A_residues":
                    conserved_macaques,

                "exposed_non_interface_residues":
                    exposed_non_interface,

                "interface_residues":
                    interface,

                "mean_complex_rsa":
                    (
                        f"{mean_rsa:.4f}"
                        if mean_rsa
                        is not None
                        else ""
                    ),

                "footprint_priority":
                    footprint_priority,
            }
        )

    # ---------------------------------------------------------------------
    # Write
    # ---------------------------------------------------------------------

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fields = [
        "footprint_id",
        "center_residues",
        "residue_count",
        "residue_labels",

        "human_discriminatory_residues",
        "rhesus_both_isoforms_discriminatory_residues",
        "pigtail_discriminatory_residues",

        "single_residue_pan_species_discriminators",

        "footprint_has_human_discrimination",
        "footprint_has_rhesus_discrimination",
        "footprint_has_pigtail_discrimination",
        "footprint_cross_species_discrimination",

        "cross_species_conserved_NKG2A_residues",
        "macaque_conserved_NKG2A_residues",

        "exposed_non_interface_residues",
        "interface_residues",

        "mean_complex_rsa",

        "footprint_priority",
    ]

    with OUTPUT.open(
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

    # ---------------------------------------------------------------------
    # Console ranking
    # ---------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "COMPACT FOOTPRINT CANDIDATES"
    )
    print("=" * 78)

    priority_order = {
        "cross_species_NKG2A_"
        "discriminatory_noninterface":
            0,

        "cross_species_NKG2A_"
        "discriminatory_mixed_interface":
            1,

        "macaque_NKG2A_"
        "discriminatory_noninterface":
            2,

        "pigtail_NKG2A_"
        "discriminatory_noninterface":
            3,

        "accessible_low_cross_species_"
        "discrimination":
            4,

        "interface_containing":
            5,

        "low_structural_priority":
            6,
    }

    ranked = sorted(
        output_rows,
        key=lambda row: (
            priority_order.get(
                row[
                    "footprint_priority"
                ],
                99,
            ),

            -int(
                row[
                    "exposed_non_interface_residues"
                ]
            ),

            -int(
                row[
                    "human_discriminatory_residues"
                ]
            ),

            -int(
                row[
                    "rhesus_both_isoforms_"
                    "discriminatory_residues"
                ]
            ),

            -int(
                row[
                    "pigtail_discriminatory_residues"
                ]
            ),

            -int(
                row[
                    "residue_count"
                ]
            ),
        ),
    )

    for row in ranked:

        print(
            f"{row['footprint_id']:<14} "
            f"{row['residue_labels']:<45} "
            f"H={row['human_discriminatory_residues']} "
            f"R={row['rhesus_both_isoforms_discriminatory_residues']} "
            f"P={row['pigtail_discriminatory_residues']} "
            f"EXP={row['exposed_non_interface_residues']} "
            f"INT={row['interface_residues']}  "
            f"{row['footprint_priority']}"
        )

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)

    print()
    print(OUTPUT)

    print()
    print(
        "NOTE: footprint discrimination means "
        "the local surface neighborhood contains "
        "at least one NKG2A-vs-NKG2C difference "
        "for the indicated species."
    )

    print(
        "It does not require the same individual "
        "residue to discriminate NKG2A from NKG2C "
        "in every species."
    )

    print(
        "This is still a structural screening "
        "analysis, not experimental antibody "
        "specificity evidence."
    )


if __name__ == "__main__":
    main()
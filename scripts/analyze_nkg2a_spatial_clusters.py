from pathlib import Path
import csv
import math

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa


ROOT = Path(__file__).resolve().parent.parent

PDB_ID = "3CDG"

CIF_FILE = (
    ROOT
    / "structures"
    / "reference"
    / f"{PDB_ID}.cif"
)

SASA_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_solvent_accessibility.tsv"
)

RESIDUE_MAP_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_structure_residue_map.tsv"
)

NEIGHBOR_OUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_candidate_spatial_neighbors.tsv"
)

CLUSTER_OUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_candidate_surface_clusters.tsv"
)


NKG2A_CHAIN = "F"

DISTANCE_THRESHOLDS = [
    6.0,
    8.0,
    10.0,
    12.0,
]


def read_tsv(path):
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
        writer.writerows(rows)


def clean(value):
    if value is None:
        return ""

    return str(value).strip()


def safe_int(value):
    value = clean(value)

    if not value:
        return None

    try:
        return int(value)

    except ValueError:
        return None


def safe_float(value):
    value = clean(value)

    if not value:
        return None

    try:
        return float(value)

    except ValueError:
        return None


def atom_distance(atom1, atom2):
    delta = atom1.coord - atom2.coord

    return math.sqrt(
        float(
            delta.dot(delta)
        )
    )


def heavy_atoms(residue):
    atoms = []

    for atom in residue.get_atoms():

        element = clean(
            getattr(
                atom,
                "element",
                "",
            )
        ).upper()

        if element == "H":
            continue

        atoms.append(atom)

    return atoms


def minimum_residue_distance(
    residue1,
    residue2,
):
    atoms1 = heavy_atoms(
        residue1
    )

    atoms2 = heavy_atoms(
        residue2
    )

    if not atoms1 or not atoms2:
        return None

    best = None

    for atom1 in atoms1:

        for atom2 in atoms2:

            distance = atom_distance(
                atom1,
                atom2,
            )

            if (
                best is None
                or distance < best
            ):
                best = distance

    return best


def build_mapping_lookup(rows):
    lookup = {}

    for row in rows:

        if clean(
            row.get("pdb_id")
        ).upper() != PDB_ID:
            continue

        if clean(
            row.get("chain_id")
        ) != NKG2A_CHAIN:
            continue

        full_residue = safe_int(
            row.get(
                "full_length_residue"
            )
        )

        auth_residue = safe_int(
            row.get(
                "auth_residue_number"
            )
        )

        if (
            full_residue is None
            or auth_residue is None
        ):
            continue

        lookup[
            full_residue
        ] = row

    return lookup


def build_coordinate_lookup(chain):
    lookup = {}

    for residue in chain:

        if not is_aa(
            residue,
            standard=False,
        ):
            continue

        auth_number = residue.id[1]

        insertion_code = clean(
            residue.id[2]
        )

        lookup[
            (
                auth_number,
                insertion_code,
            )
        ] = residue

    return lookup


def get_coordinate_residue(
    mapping_row,
    coordinate_lookup,
):
    auth_number = safe_int(
        mapping_row.get(
            "auth_residue_number"
        )
    )

    insertion_code = clean(
        mapping_row.get(
            "insertion_code"
        )
    )

    key = (
        auth_number,
        insertion_code,
    )

    residue = coordinate_lookup.get(
        key
    )

    if residue is not None:
        return residue

    candidates = [
        residue
        for (
            number,
            ins_code
        ), residue
        in coordinate_lookup.items()
        if number == auth_number
    ]

    if len(candidates) == 1:
        return candidates[0]

    return None


def candidate_score(row):
    """
    Sequence/structure screening score only.

    This is not a binding-energy score.
    """

    score = 0

    classification = clean(
        row.get(
            "classification"
        )
    )

    priority = clean(
        row.get(
            "sequence_priority"
        )
    )

    exposure = clean(
        row.get(
            "complex_exposure_class"
        )
    )

    interface_contact = clean(
        row.get(
            "any_interface_contact"
        )
    )

    if classification == (
        "pan_species_NKG2A_vs_NKG2C_difference"
    ):
        score += 4

    elif classification == (
        "macaque_shared_difference"
    ):
        score += 4

    elif classification == (
        "pigtail_specific_difference"
    ):
        score += 3

    elif classification == (
        "human_pigtail_shared_difference"
    ):
        score += 2

    else:
        score += 1

    if priority.startswith("high"):
        score += 3

    elif priority.startswith("moderate"):
        score += 2

    elif priority:
        score += 1

    if exposure == "exposed":
        score += 3

    elif exposure == "partially_exposed":
        score += 1

    if interface_contact == "no":
        score += 2

    return score


def main():

    print("=" * 78)
    print("STEP 2H - NKG2A SPATIAL CANDIDATE CLUSTER ANALYSIS")
    print("=" * 78)

    parser = MMCIFParser(
        QUIET=True
    )

    structure = parser.get_structure(
        PDB_ID,
        str(CIF_FILE),
    )

    model = structure[0]

    chain = model[
        NKG2A_CHAIN
    ]

    mapping_lookup = (
        build_mapping_lookup(
            read_tsv(
                RESIDUE_MAP_FILE
            )
        )
    )

    coordinate_lookup = (
        build_coordinate_lookup(
            chain
        )
    )

    sasa_rows = read_tsv(
        SASA_FILE
    )

    candidate_rows = [
        row
        for row in sasa_rows
        if clean(
            row.get(
                "candidate_specificity_position"
            )
        ) == "yes"
    ]

    print()
    print(
        f"Candidate positions in SASA table: "
        f"{len(candidate_rows)}"
    )

    candidates = []

    for row in candidate_rows:

        full_residue = safe_int(
            row.get(
                "full_length_residue"
            )
        )

        if full_residue is None:
            continue

        mapping = mapping_lookup.get(
            full_residue
        )

        if mapping is None:
            continue

        residue = get_coordinate_residue(
            mapping,
            coordinate_lookup,
        )

        if residue is None:
            continue

        candidates.append(
            {
                "full_residue":
                    full_residue,

                "aa":
                    clean(
                        row.get(
                            "nkg2a_aa"
                        )
                    ),

                "classification":
                    clean(
                        row.get(
                            "classification"
                        )
                    ),

                "sequence_priority":
                    clean(
                        row.get(
                            "sequence_priority"
                        )
                    ),

                "complex_rsa":
                    safe_float(
                        row.get(
                            "complex_rsa"
                        )
                    ),

                "exposure_class":
                    clean(
                        row.get(
                            "complex_exposure_class"
                        )
                    ),

                "interface_contact":
                    clean(
                        row.get(
                            "any_interface_contact"
                        )
                    ),

                "screening_score":
                    candidate_score(
                        row
                    ),

                "residue":
                    residue,
            }
        )

    print(
        "Candidate positions with coordinates: "
        f"{len(candidates)}"
    )

    # ---------------------------------------------------------------------
    # Pairwise candidate distances
    # ---------------------------------------------------------------------

    neighbor_rows = []

    for i in range(
        len(candidates)
    ):

        for j in range(
            i + 1,
            len(candidates),
        ):

            c1 = candidates[i]
            c2 = candidates[j]

            distance = minimum_residue_distance(
                c1["residue"],
                c2["residue"],
            )

            if distance is None:
                continue

            row = {
                "residue_1":
                    c1["full_residue"],

                "aa_1":
                    c1["aa"],

                "classification_1":
                    c1["classification"],

                "residue_2":
                    c2["full_residue"],

                "aa_2":
                    c2["aa"],

                "classification_2":
                    c2["classification"],

                "minimum_heavy_atom_distance_A":
                    f"{distance:.3f}",
            }

            for cutoff in (
                DISTANCE_THRESHOLDS
            ):

                field = (
                    f"within_{int(cutoff)}A"
                )

                row[field] = (
                    "yes"
                    if distance <= cutoff
                    else "no"
                )

            neighbor_rows.append(
                row
            )

    neighbor_fields = [
        "residue_1",
        "aa_1",
        "classification_1",
        "residue_2",
        "aa_2",
        "classification_2",
        "minimum_heavy_atom_distance_A",
    ]

    for cutoff in (
        DISTANCE_THRESHOLDS
    ):

        neighbor_fields.append(
            f"within_{int(cutoff)}A"
        )

    write_tsv(
        NEIGHBOR_OUT,
        neighbor_rows,
        neighbor_fields,
    )

    # ---------------------------------------------------------------------
    # Per-residue cluster summaries
    # ---------------------------------------------------------------------

    cluster_rows = []

    for candidate in candidates:

        residue_number = (
            candidate[
                "full_residue"
            ]
        )

        neighbor_counts = {}

        neighbor_lists = {}

        for cutoff in (
            DISTANCE_THRESHOLDS
        ):

            neighbors = []

            for other in candidates:

                if (
                    other["full_residue"]
                    == residue_number
                ):
                    continue

                distance = (
                    minimum_residue_distance(
                        candidate[
                            "residue"
                        ],
                        other[
                            "residue"
                        ],
                    )
                )

                if (
                    distance is not None
                    and distance <= cutoff
                ):
                    neighbors.append(
                        other[
                            "full_residue"
                        ]
                    )

            neighbor_counts[
                cutoff
            ] = len(
                neighbors
            )

            neighbor_lists[
                cutoff
            ] = ",".join(
                str(value)
                for value
                in sorted(neighbors)
            )

        # Simple cluster score:
        # candidate quality + density of nearby discriminatory sites.
        cluster_score = (
            candidate[
                "screening_score"
            ]
            + 2
            * neighbor_counts[
                8.0
            ]
            + neighbor_counts[
                12.0
            ]
        )

        cluster_rows.append(
            {
                "full_length_residue":
                    residue_number,

                "nkg2a_aa":
                    candidate[
                        "aa"
                    ],

                "classification":
                    candidate[
                        "classification"
                    ],

                "sequence_priority":
                    candidate[
                        "sequence_priority"
                    ],

                "complex_rsa":
                    (
                        ""
                        if candidate[
                            "complex_rsa"
                        ]
                        is None
                        else (
                            f"{candidate['complex_rsa']:.4f}"
                        )
                    ),

                "exposure_class":
                    candidate[
                        "exposure_class"
                    ],

                "interface_contact":
                    candidate[
                        "interface_contact"
                    ],

                "screening_score":
                    candidate[
                        "screening_score"
                    ],

                "neighbors_within_6A":
                    neighbor_counts[
                        6.0
                    ],

                "neighbors_6A":
                    neighbor_lists[
                        6.0
                    ],

                "neighbors_within_8A":
                    neighbor_counts[
                        8.0
                    ],

                "neighbors_8A":
                    neighbor_lists[
                        8.0
                    ],

                "neighbors_within_10A":
                    neighbor_counts[
                        10.0
                    ],

                "neighbors_10A":
                    neighbor_lists[
                        10.0
                    ],

                "neighbors_within_12A":
                    neighbor_counts[
                        12.0
                    ],

                "neighbors_12A":
                    neighbor_lists[
                        12.0
                    ],

                "cluster_score":
                    cluster_score,
            }
        )

    cluster_rows.sort(
        key=lambda row: (
            -int(
                row[
                    "cluster_score"
                ]
            ),
            -int(
                row[
                    "neighbors_within_8A"
                ]
            ),
            int(
                row[
                    "full_length_residue"
                ]
            ),
        )
    )

    cluster_fields = [
        "full_length_residue",
        "nkg2a_aa",
        "classification",
        "sequence_priority",
        "complex_rsa",
        "exposure_class",
        "interface_contact",
        "screening_score",
        "neighbors_within_6A",
        "neighbors_6A",
        "neighbors_within_8A",
        "neighbors_8A",
        "neighbors_within_10A",
        "neighbors_10A",
        "neighbors_within_12A",
        "neighbors_12A",
        "cluster_score",
    ]

    write_tsv(
        CLUSTER_OUT,
        cluster_rows,
        cluster_fields,
    )

    # ---------------------------------------------------------------------
    # Console summary
    # ---------------------------------------------------------------------

    print()
    print("=" * 78)
    print("TOP SPATIAL CANDIDATE CENTERS")
    print("=" * 78)

    for row in cluster_rows[:10]:

        print(
            f"{row['full_length_residue']:>3} "
            f"{row['nkg2a_aa']:<2} "
            f"cluster={row['cluster_score']:<3} "
            f"8A={row['neighbors_within_8A']:<2} "
            f"12A={row['neighbors_within_12A']:<2} "
            f"RSA={row['complex_rsa']:<7} "
            f"interface={row['interface_contact']:<3} "
            f"{row['classification']}"
        )

    print()
    print("=" * 78)
    print("OUTPUTS")
    print("=" * 78)

    print(
        NEIGHBOR_OUT
    )

    print(
        CLUSTER_OUT
    )

    print()
    print(
        "NOTE: cluster_score is a heuristic screening metric, "
        "not an affinity or binding-energy prediction."
    )


if __name__ == "__main__":
    main()
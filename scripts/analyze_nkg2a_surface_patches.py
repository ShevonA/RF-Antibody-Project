from pathlib import Path
from collections import defaultdict
import csv
import math

from Bio.PDB import MMCIFParser


# =============================================================================
# PATHS
# =============================================================================

ROOT = Path(__file__).resolve().parent.parent

STRUCTURE_FILE = (
    ROOT
    / "structures"
    / "reference"
    / "3CDG.cif"
)

INTEGRATION_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_epitope_candidate_integration.tsv"
)

STRUCTURE_MAP_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_structure_residue_map.tsv"
)

OUTPUT_DIR = (
    ROOT
    / "results"
    / "tables"
    / "structure"
)

PAIRWISE_OUT = (
    OUTPUT_DIR
    / "nkg2a_candidate_spatial_distances.tsv"
)

PATCH_OUT = (
    OUTPUT_DIR
    / "nkg2a_candidate_surface_patches.tsv"
)


# =============================================================================
# ANALYSIS SETTINGS
# =============================================================================

PDB_ID = "3CDG"

# Biological assembly 2
NKG2A_CHAIN = "F"

# Candidate-candidate neighborhood cutoffs.
PATCH_CUTOFFS = [
    6.0,
    8.0,
    10.0,
    12.0,
]

# Primary cutoff used to define connected candidate patches.
PRIMARY_PATCH_CUTOFF = 10.0


# =============================================================================
# HELPERS
# =============================================================================

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


def write_tsv(path, rows, fieldnames):
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


def yes(value):
    return clean(value).lower() == "yes"


def as_int(value):
    value = clean(value)

    if not value:
        return None

    try:
        return int(value)
    except ValueError:
        return int(float(value))


def as_float(value):
    value = clean(value)

    if not value:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def get_first(row, names, default=""):
    """
    Return the first matching column value.

    This makes the script tolerant of small naming differences
    between earlier Step 2 output tables.
    """

    for name in names:
        if name in row:
            return row[name]

    return default


def minimum_residue_distance(residue1, residue2):
    """
    Minimum heavy-atom distance between two residues.
    """

    minimum = None

    for atom1 in residue1.get_atoms():

        if atom1.element == "H":
            continue

        for atom2 in residue2.get_atoms():

            if atom2.element == "H":
                continue

            distance = atom1 - atom2

            if minimum is None or distance < minimum:
                minimum = distance

    return minimum


def residue_one_letter(residue):
    """
    Convert a standard amino-acid residue name to one-letter code.
    """

    mapping = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLN": "Q",
        "GLU": "E",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        "MSE": "M",
    }

    return mapping.get(
        residue.get_resname().upper(),
        "X",
    )


def connected_components(nodes, adjacency):
    """
    Return connected components of an undirected graph.
    """

    unseen = set(nodes)
    components = []

    while unseen:

        start = min(unseen)
        stack = [start]
        component = set()

        while stack:

            node = stack.pop()

            if node in component:
                continue

            component.add(node)
            unseen.discard(node)

            for neighbor in adjacency.get(node, set()):

                if neighbor not in component:
                    stack.append(neighbor)

        components.append(
            sorted(component)
        )

    components.sort(
        key=lambda x: (
            -len(x),
            x[0],
        )
    )

    return components


def species_discrimination(row):
    """
    Extract within-species NKG2A-vs-NKG2C discrimination flags.
    """

    human = yes(
        get_first(
            row,
            [
                "human_A_vs_C_diff",
                "human_NKG2A_vs_NKG2C_diff",
            ],
        )
    )

    rhesus1 = yes(
        get_first(
            row,
            [
                "rhesus_A_vs_C1_diff",
                "rhesus_NKG2A_vs_NKG2C1_diff",
            ],
        )
    )

    rhesus2 = yes(
        get_first(
            row,
            [
                "rhesus_A_vs_C2_diff",
                "rhesus_NKG2A_vs_NKG2C2_diff",
            ],
        )
    )

    pigtail = yes(
        get_first(
            row,
            [
                "pigtail_A_vs_C_diff",
                "pigtail_NKG2A_vs_NKG2C_diff",
            ],
        )
    )

    rhesus_any = rhesus1 or rhesus2
    rhesus_both = rhesus1 and rhesus2

    return {
        "human": human,
        "rhesus1": rhesus1,
        "rhesus2": rhesus2,
        "rhesus_any": rhesus_any,
        "rhesus_both": rhesus_both,
        "pigtail": pigtail,
    }


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("STEP 2L - NKG2A CANDIDATE SURFACE PATCH ANALYSIS")
    print("=" * 78)

    print()
    print(f"PDB: {PDB_ID}")
    print("Biological assembly: 2")
    print(f"NKG2A author chain: {NKG2A_CHAIN}")
    print(
        "Patch distance cutoffs: "
        + ", ".join(
            f"{x:.1f} A"
            for x in PATCH_CUTOFFS
        )
    )
    print(
        f"Primary connected-patch cutoff: "
        f"{PRIMARY_PATCH_CUTOFF:.1f} A"
    )

    # -------------------------------------------------------------------------
    # Load integrated candidate table
    # -------------------------------------------------------------------------

    integration_rows = read_tsv(
        INTEGRATION_FILE
    )

    candidates = {}

    insertion_only = []

    for row in integration_rows:

        residue_number = as_int(
            get_first(
                row,
                [
                    "human_NKG2A_residue",
                    "full_length_residue",
                ],
            )
        )

        if residue_number is None:

            insertion_only.append(row)
            continue

        candidate_flag = get_first(
            row,
            [
                "candidate_specificity_position",
            ],
        )

        # If the integration table contains candidate flag,
        # honor it. Otherwise every row in this table is assumed
        # to represent an integrated candidate position.
        if candidate_flag and not yes(candidate_flag):
            continue

        candidates[residue_number] = row

    print()
    print(
        "Human NKG2A candidate residues loaded: "
        f"{len(candidates)}"
    )

    print(
        "Insertion-only candidate positions: "
        f"{len(insertion_only)}"
    )

    # -------------------------------------------------------------------------
    # Load structure residue map
    # -------------------------------------------------------------------------

    map_rows = read_tsv(
        STRUCTURE_MAP_FILE
    )

    structure_mapping = {}

    for row in map_rows:

        if clean(row.get("pdb_id")) != PDB_ID:
            continue

        if clean(row.get("chain_id")) != NKG2A_CHAIN:
            continue

        full_residue = as_int(
            row.get("full_length_residue")
        )

        auth_residue = as_int(
            row.get("auth_residue_number")
        )

        if (
            full_residue is None
            or auth_residue is None
        ):
            continue

        structure_mapping[
            full_residue
        ] = auth_residue

    print(
        "Mapped NKG2A structural residues loaded: "
        f"{len(structure_mapping)}"
    )

    # -------------------------------------------------------------------------
    # Load 3CDG coordinates
    # -------------------------------------------------------------------------

    parser = MMCIFParser(
        QUIET=True
    )

    structure = parser.get_structure(
        PDB_ID,
        str(STRUCTURE_FILE),
    )

    model = structure[0]

    if NKG2A_CHAIN not in model:
        raise ValueError(
            f"NKG2A chain {NKG2A_CHAIN} "
            "was not found in 3CDG."
        )

    chain = model[NKG2A_CHAIN]

    coordinate_residues = {}

    for residue in chain:

        hetero_flag, residue_number, insertion_code = (
            residue.id
        )

        if hetero_flag.strip():
            continue

        coordinate_residues[
            residue_number
        ] = residue

    # -------------------------------------------------------------------------
    # Match candidates to coordinate residues
    # -------------------------------------------------------------------------

    candidate_coordinates = {}

    unresolved_candidates = []

    for full_residue in sorted(candidates):

        auth_residue = structure_mapping.get(
            full_residue
        )

        if auth_residue is None:

            unresolved_candidates.append(
                full_residue
            )
            continue

        residue = coordinate_residues.get(
            auth_residue
        )

        if residue is None:

            unresolved_candidates.append(
                full_residue
            )
            continue

        candidate_coordinates[
            full_residue
        ] = {
            "auth_residue": auth_residue,
            "residue": residue,
        }

    print(
        "Candidate residues with coordinates: "
        f"{len(candidate_coordinates)}"
    )

    print(
        "Candidate residues unresolved structurally: "
        f"{len(unresolved_candidates)}"
    )

    # -------------------------------------------------------------------------
    # Pairwise candidate distances
    # -------------------------------------------------------------------------

    pairwise_rows = []

    adjacency_by_cutoff = {
        cutoff: defaultdict(set)
        for cutoff in PATCH_CUTOFFS
    }

    coordinate_numbers = sorted(
        candidate_coordinates
    )

    for i, residue1_number in enumerate(
        coordinate_numbers
    ):

        residue1 = (
            candidate_coordinates[
                residue1_number
            ]["residue"]
        )

        for residue2_number in (
            coordinate_numbers[i + 1:]
        ):

            residue2 = (
                candidate_coordinates[
                    residue2_number
                ]["residue"]
            )

            distance = minimum_residue_distance(
                residue1,
                residue2,
            )

            if distance is None:
                continue

            row1 = candidates[
                residue1_number
            ]

            row2 = candidates[
                residue2_number
            ]

            output_row = {
                "residue_1":
                    residue1_number,

                "aa_1":
                    get_first(
                        row1,
                        [
                            "human_NKG2A_aa",
                            "nkg2a_aa",
                        ],
                        residue_one_letter(
                            residue1
                        ),
                    ),

                "classification_1":
                    get_first(
                        row1,
                        ["classification"],
                    ),

                "integrated_priority_1":
                    get_first(
                        row1,
                        [
                            "integrated_priority",
                            "sequence_priority",
                        ],
                    ),

                "residue_2":
                    residue2_number,

                "aa_2":
                    get_first(
                        row2,
                        [
                            "human_NKG2A_aa",
                            "nkg2a_aa",
                        ],
                        residue_one_letter(
                            residue2
                        ),
                    ),

                "classification_2":
                    get_first(
                        row2,
                        ["classification"],
                    ),

                "integrated_priority_2":
                    get_first(
                        row2,
                        [
                            "integrated_priority",
                            "sequence_priority",
                        ],
                    ),

                "minimum_heavy_atom_distance_A":
                    f"{distance:.3f}",
            }

            for cutoff in PATCH_CUTOFFS:

                within = (
                    distance <= cutoff
                )

                output_row[
                    f"within_{int(cutoff)}A"
                ] = (
                    "yes"
                    if within
                    else "no"
                )

                if within:

                    adjacency_by_cutoff[
                        cutoff
                    ][residue1_number].add(
                        residue2_number
                    )

                    adjacency_by_cutoff[
                        cutoff
                    ][residue2_number].add(
                        residue1_number
                    )

            pairwise_rows.append(
                output_row
            )

    pairwise_fields = [
        "residue_1",
        "aa_1",
        "classification_1",
        "integrated_priority_1",
        "residue_2",
        "aa_2",
        "classification_2",
        "integrated_priority_2",
        "minimum_heavy_atom_distance_A",
    ]

    for cutoff in PATCH_CUTOFFS:
        pairwise_fields.append(
            f"within_{int(cutoff)}A"
        )

    write_tsv(
        PAIRWISE_OUT,
        pairwise_rows,
        pairwise_fields,
    )

    # -------------------------------------------------------------------------
    # Connected patches at primary cutoff
    # -------------------------------------------------------------------------

    primary_adjacency = (
        adjacency_by_cutoff[
            PRIMARY_PATCH_CUTOFF
        ]
    )

    components = connected_components(
        coordinate_numbers,
        primary_adjacency,
    )

    patch_rows = []

    print()
    print("=" * 78)
    print(
        f"CANDIDATE PATCHES AT "
        f"{PRIMARY_PATCH_CUTOFF:.1f} A"
    )
    print("=" * 78)

    for patch_index, component in enumerate(
        components,
        start=1,
    ):

        component_rows = [
            candidates[x]
            for x in component
        ]

        human_count = 0
        rhesus_any_count = 0
        rhesus_both_count = 0
        pigtail_count = 0

        pan_species_count = 0

        exposed_count = 0
        non_interface_count = 0
        exposed_non_interface_count = 0

        interface_count = 0

        high_priority_count = 0

        rsa_values = []

        labels = []

        for residue_number, row in zip(
            component,
            component_rows,
        ):

            discrimination = (
                species_discrimination(row)
            )

            if discrimination["human"]:
                human_count += 1

            if discrimination["rhesus_any"]:
                rhesus_any_count += 1

            if discrimination["rhesus_both"]:
                rhesus_both_count += 1

            if discrimination["pigtail"]:
                pigtail_count += 1

            if (
                discrimination["human"]
                and discrimination["rhesus_both"]
                and discrimination["pigtail"]
            ):
                pan_species_count += 1

            exposure = get_first(
                row,
                [
                    "complex_exposure_class",
                ],
            ).lower()

            if exposure == "exposed":
                exposed_count += 1

            interface = yes(
                get_first(
                    row,
                    [
                        "any_interface_contact",
                    ],
                )
            )

            if interface:
                interface_count += 1
            else:
                non_interface_count += 1

            exposed_non_interface = yes(
                get_first(
                    row,
                    [
                        "exposed_non_interface",
                    ],
                )
            )

            if exposed_non_interface:
                exposed_non_interface_count += 1

            rsa = as_float(
                get_first(
                    row,
                    [
                        "complex_rsa",
                    ],
                )
            )

            if rsa is not None:
                rsa_values.append(rsa)

            priority = get_first(
                row,
                [
                    "integrated_priority",
                    "sequence_priority",
                ],
            )

            if priority.startswith(
                (
                    "highest",
                    "high_",
                    "high ",
                )
            ):
                high_priority_count += 1

            aa = get_first(
                row,
                [
                    "human_NKG2A_aa",
                    "nkg2a_aa",
                ],
                "?",
            )

            labels.append(
                f"{residue_number}{aa}"
            )

        mean_rsa = (
            sum(rsa_values)
            / len(rsa_values)
            if rsa_values
            else None
        )

        # ---------------------------------------------------------------------
        # Transparent patch classification
        # ---------------------------------------------------------------------

        if (
            human_count > 0
            and rhesus_both_count > 0
            and pigtail_count > 0
            and exposed_non_interface_count > 0
        ):
            patch_priority = (
                "cross_species_discriminatory_patch"
            )

        elif (
            rhesus_both_count > 0
            and pigtail_count > 0
            and exposed_non_interface_count > 0
        ):
            patch_priority = (
                "macaque_discriminatory_patch"
            )

        elif (
            pigtail_count > 0
            and exposed_non_interface_count > 0
        ):
            patch_priority = (
                "pigtail_discriminatory_patch"
            )

        elif interface_count > 0:
            patch_priority = (
                "interface_associated_patch"
            )

        else:
            patch_priority = (
                "lower_priority_patch"
            )

        patch_row = {
            "patch_id":
                f"patch_{patch_index}",

            "distance_cutoff_A":
                f"{PRIMARY_PATCH_CUTOFF:.1f}",

            "residue_count":
                len(component),

            "residues":
                ",".join(
                    str(x)
                    for x in component
                ),

            "residue_labels":
                ",".join(labels),

            "human_discriminatory_residues":
                human_count,

            "rhesus_any_discriminatory_residues":
                rhesus_any_count,

            "rhesus_both_isoforms_discriminatory_residues":
                rhesus_both_count,

            "pigtail_discriminatory_residues":
                pigtail_count,

            "single_residue_pan_species_discriminators":
                pan_species_count,

            "exposed_residues":
                exposed_count,

            "non_interface_residues":
                non_interface_count,

            "exposed_non_interface_residues":
                exposed_non_interface_count,

            "interface_residues":
                interface_count,

            "high_priority_residues":
                high_priority_count,

            "mean_complex_rsa":
                (
                    f"{mean_rsa:.4f}"
                    if mean_rsa is not None
                    else ""
                ),

            "patch_priority":
                patch_priority,
        }

        patch_rows.append(
            patch_row
        )

        print()
        print(
            f"Patch {patch_index}: "
            + ", ".join(labels)
        )

        print(
            f"  residues: "
            f"{len(component)}"
        )

        print(
            "  discrimination: "
            f"human={human_count}, "
            f"rhesus_both={rhesus_both_count}, "
            f"pigtail={pigtail_count}"
        )

        print(
            "  structural: "
            f"exposed_non_interface="
            f"{exposed_non_interface_count}, "
            f"interface={interface_count}"
        )

        print(
            f"  priority: "
            f"{patch_priority}"
        )

    patch_fields = [
        "patch_id",
        "distance_cutoff_A",
        "residue_count",
        "residues",
        "residue_labels",
        "human_discriminatory_residues",
        "rhesus_any_discriminatory_residues",
        "rhesus_both_isoforms_discriminatory_residues",
        "pigtail_discriminatory_residues",
        "single_residue_pan_species_discriminators",
        "exposed_residues",
        "non_interface_residues",
        "exposed_non_interface_residues",
        "interface_residues",
        "high_priority_residues",
        "mean_complex_rsa",
        "patch_priority",
    ]

    write_tsv(
        PATCH_OUT,
        patch_rows,
        patch_fields,
    )

    # -------------------------------------------------------------------------
    # Neighborhood summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("NEIGHBORHOOD SUMMARY")
    print("=" * 78)

    for cutoff in PATCH_CUTOFFS:

        adjacency = adjacency_by_cutoff[
            cutoff
        ]

        pair_count = sum(
            len(neighbors)
            for neighbors in adjacency.values()
        ) // 2

        components_at_cutoff = (
            connected_components(
                coordinate_numbers,
                adjacency,
            )
        )

        multi_residue_components = [
            component
            for component
            in components_at_cutoff
            if len(component) > 1
        ]

        print(
            f"{cutoff:>4.1f} A: "
            f"{pair_count:>3} candidate pairs, "
            f"{len(multi_residue_components):>2} "
            "multi-residue connected patch(es)"
        )

    # -------------------------------------------------------------------------
    # Structurally unresolved candidates
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("STRUCTURALLY UNRESOLVED CANDIDATES")
    print("=" * 78)

    if unresolved_candidates:

        for residue_number in (
            unresolved_candidates
        ):

            row = candidates[
                residue_number
            ]

            aa = get_first(
                row,
                [
                    "human_NKG2A_aa",
                    "nkg2a_aa",
                ],
                "?",
            )

            classification = get_first(
                row,
                ["classification"],
            )

            priority = get_first(
                row,
                [
                    "integrated_priority",
                    "sequence_priority",
                ],
            )

            print(
                f"{residue_number:>3} "
                f"{aa:<2} "
                f"{classification:<40} "
                f"{priority}"
            )

    else:
        print("None.")

    # -------------------------------------------------------------------------
    # Insertion-only candidates
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("NKG2C INSERTION-ONLY CANDIDATES")
    print("=" * 78)

    if insertion_only:

        for row in insertion_only:

            alignment_column = (
                get_first(
                    row,
                    ["alignment_column"],
                )
            )

            human_c_aa = get_first(
                row,
                ["human_NKG2C_aa"],
            )

            classification = get_first(
                row,
                ["classification"],
            )

            print(
                f"Alignment column "
                f"{alignment_column}: "
                f"NKG2A=- "
                f"NKG2C={human_c_aa} "
                f"{classification}"
            )

    else:
        print("None.")

    # -------------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)

    print()
    print(PAIRWISE_OUT)
    print(PATCH_OUT)

    print()
    print(
        "NOTE: candidate patches are geometric "
        "screening neighborhoods, not experimentally "
        "validated antibody epitopes."
    )

    print(
        "A connected patch means that candidate residues "
        "can be linked through pairwise heavy-atom "
        f"distances <= {PRIMARY_PATCH_CUTOFF:.1f} A."
    )

    print(
        "Structurally unresolved N-terminal candidates "
        "remain candidates and are not assigned low "
        "priority merely because 3CDG lacks coordinates."
    )


if __name__ == "__main__":
    main()
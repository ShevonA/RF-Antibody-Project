from pathlib import Path
import csv
import math

from Bio.PDB import MMCIFParser


# =============================================================================
# PATHS
# =============================================================================

ROOT = Path(__file__).resolve().parent.parent

PDB_ID = "3CDG"

STRUCTURE_FILE = (
    ROOT
    / "structures"
    / "reference"
    / "3CDG.cif"
)

CONTACT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_structural_contacts.tsv"
)

GLYCOSYLATION_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2_ectodomain_n_glycosylation_sites.tsv"
)

GLYCAN_INVENTORY_OUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "3CDG_glycan_inventory.tsv"
)

GLYCAN_PROXIMITY_OUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_candidate_glycan_proximity.tsv"
)


# =============================================================================
# STRUCTURAL SETTINGS
# =============================================================================

# Verified biological assembly 2 author-chain assignments from Step 2E.
ASSEMBLY_ID = "2"

CHAIN_COMPONENTS = {
    "C": "HLA-E",
    "D": "B2M",
    "E": "CD94",
    "F": "NKG2A",
    "Q": "peptide",
}

NKG2A_CHAIN = "F"

# Distance used only as a screening threshold for possible steric proximity.
GLYCAN_PROXIMITY_CUTOFF_A = 5.0


# =============================================================================
# CARBOHYDRATE DEFINITIONS
# =============================================================================

# Common carbohydrate residue names encountered in PDB/mmCIF structures.
#
# We intentionally use a defined list rather than treating every HETATM as a
# carbohydrate. This prevents waters, ions, buffer molecules, etc. from being
# misclassified as glycans.
#
# The inventory also prints any unclassified hetero residues so that unexpected
# carbohydrate residue names can be noticed and added deliberately.

CARBOHYDRATE_RESNAMES = {
    # N-acetylglucosamine
    "NAG",
    "NDG",

    # Mannose
    "MAN",
    "BMA",

    # Fucose
    "FUC",
    "FUL",

    # Galactose
    "GAL",
    "GLA",

    # Glucose
    "GLC",
    "BGC",

    # N-acetylgalactosamine
    "NGA",
    "A2G",

    # Sialic acids
    "SIA",
    "NAN",
    "NGC",

    # Other common monosaccharides
    "XYS",
    "XYL",
    "ARA",
    "RIB",
    "FRU",
}


# =============================================================================
# BASIC HELPERS
# =============================================================================

def is_hydrogen(atom):
    """
    Return True for hydrogen/deuterium atoms.
    """

    element = (atom.element or "").strip().upper()

    if element in {"H", "D"}:
        return True

    atom_name = atom.get_name().strip().upper()

    return atom_name.startswith("H")


def heavy_atoms(residue):
    """
    Return all non-hydrogen atoms in a residue.
    """

    return [
        atom
        for atom in residue.get_atoms()
        if not is_hydrogen(atom)
    ]


def residue_number(residue):
    """
    Return the author residue number used by Bio.PDB.
    """

    return residue.id[1]


def insertion_code(residue):
    """
    Return insertion code or blank string.
    """

    code = residue.id[2]

    if code is None:
        return ""

    return str(code).strip()


def residue_identifier(residue):
    """
    Human-readable residue identifier.
    """

    number = residue_number(residue)
    icode = insertion_code(residue)

    if icode:
        return f"{number}{icode}"

    return str(number)


def euclidean_distance(atom1, atom2):
    """
    Euclidean distance between two Bio.PDB atoms.
    """

    c1 = atom1.coord
    c2 = atom2.coord

    return math.sqrt(
        float(
            (c1[0] - c2[0]) ** 2
            + (c1[1] - c2[1]) ** 2
            + (c1[2] - c2[2]) ** 2
        )
    )


def minimum_residue_distance(residue1, residue2):
    """
    Minimum heavy-atom distance between two residues.

    Returns:
        distance,
        atom_name_residue1,
        atom_name_residue2
    """

    atoms1 = heavy_atoms(residue1)
    atoms2 = heavy_atoms(residue2)

    if not atoms1 or not atoms2:
        return None, "", ""

    best_distance = None
    best_atom1 = ""
    best_atom2 = ""

    for atom1 in atoms1:
        for atom2 in atoms2:

            distance = euclidean_distance(
                atom1,
                atom2,
            )

            if (
                best_distance is None
                or distance < best_distance
            ):
                best_distance = distance
                best_atom1 = atom1.get_name().strip()
                best_atom2 = atom2.get_name().strip()

    return (
        best_distance,
        best_atom1,
        best_atom2,
    )


def format_distance(value):
    """
    Format structural distances consistently.
    """

    if value is None:
        return ""

    return f"{value:.3f}"


# =============================================================================
# INPUT TABLES
# =============================================================================

def read_candidate_positions():
    """
    Read candidate NKG2A residues from the Step 2F structural-contact table.

    This is preferable to reconstructing the candidate list from scratch,
    because Step 2F already established which candidate positions have
    coordinates in 3CDG NKG2A chain F.
    """

    if not CONTACT_FILE.exists():
        raise FileNotFoundError(
            f"Missing structural contact table:\n{CONTACT_FILE}"
        )

    candidates = {}

    with CONTACT_FILE.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        required = {
            "full_length_residue",
            "nkg2a_aa",
            "candidate_specificity_position",
            "classification",
            "sequence_priority",
        }

        missing = required - set(
            reader.fieldnames or []
        )

        if missing:
            raise ValueError(
                "Structural contact table is missing required "
                f"columns: {sorted(missing)}"
            )

        for row in reader:

            if (
                row[
                    "candidate_specificity_position"
                ].strip().lower()
                != "yes"
            ):
                continue

            residue_text = row[
                "full_length_residue"
            ].strip()

            if not residue_text:
                continue

            residue_num = int(residue_text)

            candidates[residue_num] = {
                "full_length_residue":
                    residue_num,
                "nkg2a_aa":
                    row["nkg2a_aa"].strip(),
                "classification":
                    row["classification"].strip(),
                "sequence_priority":
                    row["sequence_priority"].strip(),
            }

    return candidates


def read_human_nkg2a_sequons():
    """
    Read human NKG2A canonical N-X-S/T sequons from Step 2I.

    These are included in the output for interpretation only. A modeled
    carbohydrate can be absent even when a sequon exists.
    """

    sequons = {}

    if not GLYCOSYLATION_FILE.exists():
        return sequons

    with GLYCOSYLATION_FILE.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:

            if (
                row.get(
                    "record_id",
                    "",
                ).strip()
                != "human_NKG2A"
            ):
                continue

            residue_text = row.get(
                "sequon_full_length_residue",
                "",
            ).strip()

            if not residue_text:
                continue

            residue_num = int(residue_text)

            sequons[residue_num] = (
                row.get(
                    "motif",
                    "",
                ).strip()
            )

    return sequons


# =============================================================================
# STRUCTURE HELPERS
# =============================================================================

def load_structure():
    """
    Load 3CDG using Bio.PDB.
    """

    if not STRUCTURE_FILE.exists():
        raise FileNotFoundError(
            f"Missing structure file:\n{STRUCTURE_FILE}"
        )

    parser = MMCIFParser(
        QUIET=True,
    )

    structure = parser.get_structure(
        PDB_ID,
        str(STRUCTURE_FILE),
    )

    model = structure[0]

    return structure, model


def verify_expected_chains(model):
    """
    Verify that all author chains expected for biological assembly 2
    exist in the coordinate model.
    """

    observed = {
        chain.id
        for chain in model
    }

    missing = (
        set(CHAIN_COMPONENTS)
        - observed
    )

    if missing:
        raise ValueError(
            "3CDG is missing expected biological assembly 2 "
            f"author chain(s): {sorted(missing)}"
        )


def classify_hetero_residue(residue):
    """
    Classify a hetero residue.

    Returns:
        carbohydrate
        water
        other_hetero
        polymer
    """

    hetflag = residue.id[0]

    if hetflag == " ":
        return "polymer"

    resname = residue.get_resname().strip().upper()

    if resname in {"HOH", "WAT", "DOD"}:
        return "water"

    if resname in CARBOHYDRATE_RESNAMES:
        return "carbohydrate"

    return "other_hetero"


def collect_glycans(model):
    """
    Collect carbohydrate residues from the complete coordinate model.

    We scan all chains, not only assembly 2 chains, because the same asymmetric
    unit contains both verified biological assemblies. The assembly/component
    column makes the context explicit.
    """

    glycans = []
    other_hetero = []

    for chain in model:

        component = CHAIN_COMPONENTS.get(
            chain.id,
            "outside_selected_assembly",
        )

        for residue in chain:

            classification = (
                classify_hetero_residue(
                    residue
                )
            )

            if classification == "carbohydrate":

                glycans.append(
                    {
                        "chain_id": chain.id,
                        "component": component,
                        "residue": residue,
                    }
                )

            elif classification == "other_hetero":

                other_hetero.append(
                    {
                        "chain_id": chain.id,
                        "component": component,
                        "residue": residue,
                    }
                )

    return glycans, other_hetero


def build_nkg2a_residue_lookup(model):
    """
    Map author residue number -> NKG2A residue object for chain F.

    Step 2D established the relationship between coordinate residues and
    full-length human NKG2A numbering. In 3CDG chain F, Bio.PDB author residue
    numbering corresponds to the structural construct numbering rather than
    necessarily full-length KLRC1 numbering, so we do NOT assume that the
    coordinate residue ID itself equals the full-length residue number.

    Instead, candidate full-length residue numbers are linked to chain-F
    coordinates through the Step 2F contact table below.
    """

    chain = model[NKG2A_CHAIN]

    lookup = {}

    for residue in chain:

        if residue.id[0] != " ":
            continue

        lookup[
            residue_number(residue)
        ] = residue

    return lookup


def read_structural_candidate_coordinate_numbers():
    """
    Read the Step 2F contact table and recover the coordinate residue number
    corresponding to each full-length candidate residue.

    The exact coordinate-number column name is detected from the known names
    used during the structural mapping workflow.
    """

    with CONTACT_FILE.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        fields = set(
            reader.fieldnames or []
        )

        coordinate_column = None

        possible_columns = [
            "auth_residue_number",
            "coordinate_residue_number",
            "structure_residue_number",
        ]

        for column in possible_columns:
            if column in fields:
                coordinate_column = column
                break

        if coordinate_column is None:
            raise ValueError(
                "Could not find a coordinate residue-number column "
                "in nkg2a_structural_contacts.tsv.\n"
                "Expected one of:\n"
                + "\n".join(
                    f"  {x}"
                    for x in possible_columns
                )
            )

        mapping = {}

        for row in reader:

            if (
                row.get(
                    "candidate_specificity_position",
                    "",
                ).strip().lower()
                != "yes"
            ):
                continue

            full_text = row.get(
                "full_length_residue",
                "",
            ).strip()

            coord_text = row.get(
                coordinate_column,
                "",
            ).strip()

            if not full_text or not coord_text:
                continue

            # Permit values such as "113" or "113.0".
            full_number = int(
                float(full_text)
            )

            coordinate_number = int(
                float(coord_text)
            )

            mapping[
                full_number
            ] = coordinate_number

    return mapping, coordinate_column


# =============================================================================
# GLYCAN ATTACHMENT SCREEN
# =============================================================================

def find_nearest_polymer_residue(
    glycan_residue,
    model,
):
    """
    Find the nearest polymer residue to a carbohydrate residue.

    This is a geometric screen, not a chemical-bond parser.

    A very short distance to ASN is consistent with an N-linked attachment,
    but the output deliberately labels this as a nearest-polymer relationship
    unless covalent connectivity is independently established.
    """

    best = None

    for chain in model:

        for residue in chain:

            if residue.id[0] != " ":
                continue

            distance, atom1, atom2 = (
                minimum_residue_distance(
                    glycan_residue,
                    residue,
                )
            )

            if distance is None:
                continue

            if (
                best is None
                or distance < best["distance"]
            ):
                best = {
                    "distance": distance,
                    "glycan_atom": atom1,
                    "polymer_atom": atom2,
                    "chain_id": chain.id,
                    "residue": residue,
                }

    return best


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "STEP 2J - 3CDG MODELED GLYCAN PROXIMITY ANALYSIS"
    )
    print("=" * 78)

    print(f"\nPDB: {PDB_ID}")
    print(
        f"Biological assembly used for NKG2A analysis: "
        f"{ASSEMBLY_ID}"
    )
    print(
        f"NKG2A author chain: {NKG2A_CHAIN}"
    )
    print(
        f"Glycan proximity screening cutoff: "
        f"{GLYCAN_PROXIMITY_CUTOFF_A:.1f} A"
    )

    candidates = read_candidate_positions()

    print(
        f"\nCandidate specificity residues loaded: "
        f"{len(candidates)}"
    )

    human_sequons = (
        read_human_nkg2a_sequons()
    )

    print(
        f"Human NKG2A canonical sequons loaded: "
        f"{len(human_sequons)}"
    )

    structure, model = load_structure()

    verify_expected_chains(model)

    glycans, other_hetero = (
        collect_glycans(model)
    )

    print(
        f"Modeled carbohydrate residues found: "
        f"{len(glycans)}"
    )

    print(
        f"Other non-water hetero residues found: "
        f"{len(other_hetero)}"
    )

    # -------------------------------------------------------------------------
    # Inventory modeled carbohydrates
    # -------------------------------------------------------------------------

    inventory_rows = []

    print("\n" + "=" * 78)
    print("MODELED CARBOHYDRATE INVENTORY")
    print("=" * 78)

    if not glycans:
        print(
            "\nNo recognized carbohydrate residues were "
            "found in the coordinate model."
        )

    for glycan in glycans:

        chain_id = glycan["chain_id"]
        component = glycan["component"]
        residue = glycan["residue"]

        resname = (
            residue.get_resname()
            .strip()
            .upper()
        )

        nearest = (
            find_nearest_polymer_residue(
                residue,
                model,
            )
        )

        if nearest is None:

            nearest_chain = ""
            nearest_component = ""
            nearest_residue = ""
            nearest_resname = ""
            nearest_distance = None
            glycan_atom = ""
            polymer_atom = ""

        else:

            nearest_chain = (
                nearest["chain_id"]
            )

            nearest_component = (
                CHAIN_COMPONENTS.get(
                    nearest_chain,
                    "outside_selected_assembly",
                )
            )

            nearest_residue_obj = (
                nearest["residue"]
            )

            nearest_residue = (
                residue_identifier(
                    nearest_residue_obj
                )
            )

            nearest_resname = (
                nearest_residue_obj
                .get_resname()
                .strip()
                .upper()
            )

            nearest_distance = (
                nearest["distance"]
            )

            glycan_atom = (
                nearest["glycan_atom"]
            )

            polymer_atom = (
                nearest["polymer_atom"]
            )

        inventory_rows.append(
            {
                "pdb_id": PDB_ID,
                "assembly_context":
                    ASSEMBLY_ID,
                "glycan_chain_id":
                    chain_id,
                "glycan_component_context":
                    component,
                "glycan_resname":
                    resname,
                "glycan_residue_number":
                    residue_identifier(
                        residue
                    ),
                "heavy_atom_count":
                    len(
                        heavy_atoms(
                            residue
                        )
                    ),
                "nearest_polymer_chain":
                    nearest_chain,
                "nearest_polymer_component":
                    nearest_component,
                "nearest_polymer_residue":
                    nearest_residue,
                "nearest_polymer_resname":
                    nearest_resname,
                "minimum_distance_A":
                    format_distance(
                        nearest_distance
                    ),
                "glycan_atom":
                    glycan_atom,
                "polymer_atom":
                    polymer_atom,
                "nearest_polymer_is_asn":
                    (
                        "yes"
                        if nearest_resname == "ASN"
                        else "no"
                    ),
            }
        )

        print(
            f"{chain_id:>2} "
            f"{resname:<4} "
            f"{residue_identifier(residue):<6} "
            f"context={component:<25} "
            f"nearest="
            f"{nearest_component}:"
            f"{nearest_resname}"
            f"{nearest_residue} "
            f"{format_distance(nearest_distance)} A"
        )

    GLYCAN_INVENTORY_OUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with GLYCAN_INVENTORY_OUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        fields = [
            "pdb_id",
            "assembly_context",
            "glycan_chain_id",
            "glycan_component_context",
            "glycan_resname",
            "glycan_residue_number",
            "heavy_atom_count",
            "nearest_polymer_chain",
            "nearest_polymer_component",
            "nearest_polymer_residue",
            "nearest_polymer_resname",
            "minimum_distance_A",
            "glycan_atom",
            "polymer_atom",
            "nearest_polymer_is_asn",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            inventory_rows
        )

    # -------------------------------------------------------------------------
    # Candidate -> coordinate mapping
    # -------------------------------------------------------------------------

    coordinate_mapping, coordinate_column = (
        read_structural_candidate_coordinate_numbers()
    )

    nkg2a_lookup = (
        build_nkg2a_residue_lookup(
            model
        )
    )

    print("\n" + "=" * 78)
    print("CANDIDATE NKG2A RESIDUE GLYCAN PROXIMITY")
    print("=" * 78)

    print(
        f"\nCoordinate mapping column: "
        f"{coordinate_column}"
    )

    proximity_rows = []

    candidate_count_with_coordinates = 0
    candidate_count_near_glycan = 0

    for full_residue in sorted(
        candidates
    ):

        candidate = candidates[
            full_residue
        ]

        coordinate_number = (
            coordinate_mapping.get(
                full_residue
            )
        )

        if coordinate_number is None:

            proximity_rows.append(
                {
                    "pdb_id": PDB_ID,
                    "assembly_id":
                        ASSEMBLY_ID,
                    "nkg2a_chain":
                        NKG2A_CHAIN,
                    "full_length_residue":
                        full_residue,
                    "nkg2a_aa":
                        candidate["nkg2a_aa"],
                    "classification":
                        candidate["classification"],
                    "sequence_priority":
                        candidate["sequence_priority"],
                    "coordinate_residue_number":
                        "",
                    "coordinate_present":
                        "no",
                    "nearest_glycan_chain":
                        "",
                    "nearest_glycan_resname":
                        "",
                    "nearest_glycan_residue":
                        "",
                    "minimum_glycan_distance_A":
                        "",
                    "nkg2a_atom":
                        "",
                    "glycan_atom":
                        "",
                    "within_5A_of_modeled_glycan":
                        "unknown",
                    "human_nkg2a_sequon_at_candidate":
                        (
                            human_sequons.get(
                                full_residue,
                                ""
                            )
                        ),
                }
            )

            print(
                f"{full_residue:>3} "
                f"{candidate['nkg2a_aa']}  "
                f"no coordinate mapping"
            )

            continue

        residue = nkg2a_lookup.get(
            coordinate_number
        )

        if residue is None:

            proximity_rows.append(
                {
                    "pdb_id": PDB_ID,
                    "assembly_id":
                        ASSEMBLY_ID,
                    "nkg2a_chain":
                        NKG2A_CHAIN,
                    "full_length_residue":
                        full_residue,
                    "nkg2a_aa":
                        candidate["nkg2a_aa"],
                    "classification":
                        candidate["classification"],
                    "sequence_priority":
                        candidate["sequence_priority"],
                    "coordinate_residue_number":
                        coordinate_number,
                    "coordinate_present":
                        "no",
                    "nearest_glycan_chain":
                        "",
                    "nearest_glycan_resname":
                        "",
                    "nearest_glycan_residue":
                        "",
                    "minimum_glycan_distance_A":
                        "",
                    "nkg2a_atom":
                        "",
                    "glycan_atom":
                        "",
                    "within_5A_of_modeled_glycan":
                        "unknown",
                    "human_nkg2a_sequon_at_candidate":
                        (
                            human_sequons.get(
                                full_residue,
                                ""
                            )
                        ),
                }
            )

            print(
                f"{full_residue:>3} "
                f"{candidate['nkg2a_aa']}  "
                f"coordinate residue "
                f"{coordinate_number} not found"
            )

            continue

        candidate_count_with_coordinates += 1

        best = None

        for glycan in glycans:

            distance, nkg2a_atom, glycan_atom = (
                minimum_residue_distance(
                    residue,
                    glycan["residue"],
                )
            )

            if distance is None:
                continue

            if (
                best is None
                or distance < best["distance"]
            ):
                best = {
                    "distance":
                        distance,
                    "nkg2a_atom":
                        nkg2a_atom,
                    "glycan_atom":
                        glycan_atom,
                    "glycan":
                        glycan,
                }

        if best is None:

            nearest_chain = ""
            nearest_resname = ""
            nearest_residue = ""
            nearest_distance = None
            nkg2a_atom = ""
            glycan_atom = ""
            within_cutoff = "no"

        else:

            nearest_glycan = (
                best["glycan"]
            )

            nearest_glycan_residue = (
                nearest_glycan[
                    "residue"
                ]
            )

            nearest_chain = (
                nearest_glycan[
                    "chain_id"
                ]
            )

            nearest_resname = (
                nearest_glycan_residue
                .get_resname()
                .strip()
                .upper()
            )

            nearest_residue = (
                residue_identifier(
                    nearest_glycan_residue
                )
            )

            nearest_distance = (
                best["distance"]
            )

            nkg2a_atom = (
                best["nkg2a_atom"]
            )

            glycan_atom = (
                best["glycan_atom"]
            )

            within_cutoff = (
                "yes"
                if nearest_distance
                <= GLYCAN_PROXIMITY_CUTOFF_A
                else "no"
            )

        if within_cutoff == "yes":
            candidate_count_near_glycan += 1

        proximity_rows.append(
            {
                "pdb_id": PDB_ID,
                "assembly_id":
                    ASSEMBLY_ID,
                "nkg2a_chain":
                    NKG2A_CHAIN,
                "full_length_residue":
                    full_residue,
                "nkg2a_aa":
                    candidate["nkg2a_aa"],
                "classification":
                    candidate["classification"],
                "sequence_priority":
                    candidate["sequence_priority"],
                "coordinate_residue_number":
                    coordinate_number,
                "coordinate_present":
                    "yes",
                "nearest_glycan_chain":
                    nearest_chain,
                "nearest_glycan_resname":
                    nearest_resname,
                "nearest_glycan_residue":
                    nearest_residue,
                "minimum_glycan_distance_A":
                    format_distance(
                        nearest_distance
                    ),
                "nkg2a_atom":
                    nkg2a_atom,
                "glycan_atom":
                    glycan_atom,
                "within_5A_of_modeled_glycan":
                    within_cutoff,
                "human_nkg2a_sequon_at_candidate":
                    (
                        human_sequons.get(
                            full_residue,
                            ""
                        )
                    ),
            }
        )

        glycan_label = "none"

        if best is not None:
            glycan_label = (
                f"{nearest_chain}:"
                f"{nearest_resname}"
                f"{nearest_residue}"
            )

        print(
            f"{full_residue:>3} "
            f"{candidate['nkg2a_aa']}  "
            f"nearest glycan: "
            f"{glycan_label:<15} "
            f"distance="
            f"{format_distance(nearest_distance):>7} "
            f"A  "
            f"within cutoff={within_cutoff}"
        )

    with GLYCAN_PROXIMITY_OUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        fields = [
            "pdb_id",
            "assembly_id",
            "nkg2a_chain",
            "full_length_residue",
            "nkg2a_aa",
            "classification",
            "sequence_priority",
            "coordinate_residue_number",
            "coordinate_present",
            "nearest_glycan_chain",
            "nearest_glycan_resname",
            "nearest_glycan_residue",
            "minimum_glycan_distance_A",
            "nkg2a_atom",
            "glycan_atom",
            "within_5A_of_modeled_glycan",
            "human_nkg2a_sequon_at_candidate",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            proximity_rows
        )

    # -------------------------------------------------------------------------
    # Other hetero-residue audit
    # -------------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("UNCLASSIFIED NON-WATER HETERO RESIDUES")
    print("=" * 78)

    if not other_hetero:

        print(
            "\nNone."
        )

    else:

        counts = {}

        for item in other_hetero:

            residue = item["residue"]

            key = (
                residue
                .get_resname()
                .strip()
                .upper()
            )

            counts[key] = (
                counts.get(
                    key,
                    0,
                )
                + 1
            )

        for resname in sorted(
            counts
        ):
            print(
                f"{resname:<6} "
                f"{counts[resname]}"
            )

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)

    print(
        f"\nModeled carbohydrate residues: "
        f"{len(glycans)}"
    )

    print(
        f"Candidate specificity residues: "
        f"{len(candidates)}"
    )

    print(
        f"Candidate residues with coordinates: "
        f"{candidate_count_with_coordinates}"
    )

    print(
        f"Candidates within "
        f"{GLYCAN_PROXIMITY_CUTOFF_A:.1f} A "
        f"of a modeled carbohydrate: "
        f"{candidate_count_near_glycan}"
    )

    print("\nOutputs:")
    print(GLYCAN_INVENTORY_OUT)
    print(GLYCAN_PROXIMITY_OUT)

    print(
        "\nNOTE: modeled carbohydrate proximity is structural "
        "evidence only."
    )

    print(
        "Absence of a modeled glycan does not establish absence "
        "of glycosylation in native NKG2A."
    )

    print(
        "Nearest-polymer assignments are geometric screens and "
        "should not automatically be interpreted as covalent "
        "glycosidic attachments."
    )


if __name__ == "__main__":
    main()
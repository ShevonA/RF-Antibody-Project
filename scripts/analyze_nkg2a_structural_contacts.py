from pathlib import Path
import csv
import math

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa


# =============================================================================
# PATHS / SETTINGS
# =============================================================================

ROOT = Path(__file__).resolve().parent.parent

PDB_ID = "3CDG"

CIF_FILE = (
    ROOT
    / "structures"
    / "reference"
    / f"{PDB_ID}.cif"
)

RESIDUE_MAP_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_structure_residue_map.tsv"
)

CANDIDATE_FILE = (
    ROOT
    / "results"
    / "tables"
    / "candidate_specificity_positions.tsv"
)

OUTPUT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_structural_contacts.tsv"
)


# Biological assembly 2, verified in Step 2E.
CHAINS = {
    "HLA-E": "C",
    "B2M": "D",
    "CD94": "E",
    "NKG2A": "F",
    "peptide": "Q",
}

CONTACT_CUTOFF = 4.5


# =============================================================================
# BASIC HELPERS
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


def safe_int(value):
    value = clean(value)

    if not value:
        return None

    try:
        return int(value)
    except ValueError:
        return None


def residue_label(residue):
    """
    Human-readable author residue identifier.
    """

    hetflag, number, insertion_code = residue.id

    insertion_code = clean(
        insertion_code
    )

    if insertion_code:
        return f"{number}{insertion_code}"

    return str(number)


def atom_distance(atom1, atom2):
    """
    Euclidean distance between two Bio.PDB atoms.
    """

    delta = atom1.coord - atom2.coord

    return math.sqrt(
        float(
            delta.dot(delta)
        )
    )


def heavy_atoms(residue):
    """
    Return non-hydrogen atoms for one residue.
    """

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


def chain_heavy_atoms(chain):
    """
    Return all heavy atoms from a chain.

    Protein chains use standard amino-acid residues.
    The peptide is also represented as amino acids.

    Waters are excluded.
    """

    atoms = []

    for residue in chain:

        if not is_aa(
            residue,
            standard=False,
        ):
            continue

        atoms.extend(
            heavy_atoms(residue)
        )

    return atoms


# =============================================================================
# DISTANCE CALCULATION
# =============================================================================

def minimum_distance(
    source_atoms,
    target_atoms,
):
    """
    Minimum heavy-atom distance between two atom collections.
    """

    if not source_atoms:
        return None

    if not target_atoms:
        return None

    best = None

    for atom1 in source_atoms:

        for atom2 in target_atoms:

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


def contact_flag(distance):
    if distance is None:
        return "no_coordinates"

    if distance <= CONTACT_CUTOFF:
        return "yes"

    return "no"


def format_distance(distance):
    if distance is None:
        return ""

    return f"{distance:.3f}"


# =============================================================================
# CANDIDATE TABLE
# =============================================================================

def build_candidate_lookup(rows):
    """
    Build lookup by human NKG2A full-length residue number.

    candidate_specificity_positions.tsv was generated from the
    ectodomain alignment, and human_NKG2A_residue refers to the
    full-length human NKG2A residue number.
    """

    lookup = {}

    for row in rows:

        residue_number = safe_int(
            row.get(
                "human_NKG2A_residue"
            )
        )

        if residue_number is None:
            continue

        lookup[
            residue_number
        ] = row

    return lookup


# =============================================================================
# STRUCTURE RESIDUE MAP
# =============================================================================

def get_chain_f_mapping(rows):
    """
    Select the verified 3CDG NKG2A chain F mapping produced in Step 2D.
    """

    selected = []

    for row in rows:

        if clean(
            row.get("pdb_id")
        ).upper() != PDB_ID:
            continue

        if clean(
            row.get("chain_id")
        ) != CHAINS["NKG2A"]:
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

        if full_residue is None:
            continue

        if auth_residue is None:
            continue

        selected.append(row)

    return selected


def build_coordinate_residue_lookup(chain):
    """
    Lookup coordinate residues by author residue number.

    Step 2D mapping already established which coordinate residue maps
    to which full-length NKG2A residue.
    """

    lookup = {}

    for residue in chain:

        if not is_aa(
            residue,
            standard=False,
        ):
            continue

        number = residue.id[1]

        insertion_code = clean(
            residue.id[2]
        )

        key = (
            number,
            insertion_code,
        )

        lookup[key] = residue

    return lookup


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("STEP 2F - NKG2A STRUCTURAL INTERFACE CONTACT ANALYSIS")
    print("=" * 78)

    print()
    print(f"PDB: {PDB_ID}")
    print("Biological assembly: 2")
    print(
        f"Contact cutoff: "
        f"{CONTACT_CUTOFF:.1f} A"
    )

    print()
    print("Chains:")

    for component, chain_id in CHAINS.items():
        print(
            f"  {component:<8} {chain_id}"
        )

    # -------------------------------------------------------------------------
    # Load structure
    # -------------------------------------------------------------------------

    parser = MMCIFParser(
        QUIET=True
    )

    structure = parser.get_structure(
        PDB_ID,
        str(CIF_FILE),
    )

    model = structure[0]

    for component, chain_id in CHAINS.items():

        if chain_id not in model:
            raise ValueError(
                f"Required chain {chain_id} "
                f"({component}) not found in "
                f"{PDB_ID}."
            )

    nkg2a_chain = model[
        CHAINS["NKG2A"]
    ]

    target_chains = {
        component: model[chain_id]
        for component, chain_id
        in CHAINS.items()
        if component != "NKG2A"
    }

    print()
    print("Loading target-chain atoms...")

    target_atoms = {}

    for component, chain in target_chains.items():

        atoms = chain_heavy_atoms(
            chain
        )

        target_atoms[
            component
        ] = atoms

        print(
            f"  {component:<8} "
            f"{len(atoms):>5} heavy atoms"
        )

    # -------------------------------------------------------------------------
    # Load Step 2D mapping
    # -------------------------------------------------------------------------

    residue_map_rows = read_tsv(
        RESIDUE_MAP_FILE
    )

    mapping_rows = (
        get_chain_f_mapping(
            residue_map_rows
        )
    )

    if not mapping_rows:
        raise ValueError(
            "No 3CDG chain F residue mappings "
            "were found in "
            "nkg2a_structure_residue_map.tsv."
        )

    print()
    print(
        "Mapped NKG2A coordinate residues: "
        f"{len(mapping_rows)}"
    )

    # -------------------------------------------------------------------------
    # Coordinate residue lookup
    # -------------------------------------------------------------------------

    coordinate_lookup = (
        build_coordinate_residue_lookup(
            nkg2a_chain
        )
    )

    # -------------------------------------------------------------------------
    # Candidate specificity table
    # -------------------------------------------------------------------------

    candidate_rows = read_tsv(
        CANDIDATE_FILE
    )

    candidate_lookup = (
        build_candidate_lookup(
            candidate_rows
        )
    )

    print(
        "Sequence candidate positions loaded: "
        f"{len(candidate_lookup)}"
    )

    # -------------------------------------------------------------------------
    # Calculate contacts
    # -------------------------------------------------------------------------

    output_rows = []

    missing_coordinate_lookup = 0
    candidate_count = 0

    contact_counts = {
        "CD94": 0,
        "HLA-E": 0,
        "B2M": 0,
        "peptide": 0,
    }

    for mapping in mapping_rows:

        full_residue = safe_int(
            mapping.get(
                "full_length_residue"
            )
        )

        auth_number = safe_int(
            mapping.get(
                "auth_residue_number"
            )
        )

        insertion_code = clean(
            mapping.get(
                "insertion_code"
            )
        )

        # Some Step 2D versions may not contain an insertion_code
        # column. In that case the PDB residues here use blank insertion
        # codes, which is appropriate for the mapped NKG2A residues.
        key = (
            auth_number,
            insertion_code,
        )

        residue = coordinate_lookup.get(
            key
        )

        if residue is None:

            # Fallback for a mapping table that did not preserve
            # insertion-code information.
            candidates = [
                res
                for (
                    number,
                    ins_code
                ), res
                in coordinate_lookup.items()
                if number == auth_number
            ]

            if len(candidates) == 1:
                residue = candidates[0]

        if residue is None:
            missing_coordinate_lookup += 1
            continue

        source_atoms = heavy_atoms(
            residue
        )

        distances = {}

        for component in [
            "CD94",
            "HLA-E",
            "B2M",
            "peptide",
        ]:

            distances[
                component
            ] = minimum_distance(
                source_atoms,
                target_atoms[
                    component
                ],
            )

        contacts = {
            component:
                contact_flag(
                    distances[
                        component
                    ]
                )
            for component in distances
        }

        for component, flag in contacts.items():

            if flag == "yes":
                contact_counts[
                    component
                ] += 1

        candidate = (
            candidate_lookup.get(
                full_residue
            )
        )

        is_candidate = (
            candidate is not None
        )

        if is_candidate:
            candidate_count += 1

        row = {
            "pdb_id":
                PDB_ID,

            "assembly_id":
                "2",

            "nkg2a_chain":
                CHAINS["NKG2A"],

            "auth_residue_number":
                residue_label(
                    residue
                ),

            "full_length_residue":
                full_residue,

            "nkg2a_aa":
                clean(
                    mapping.get(
                        "reference_aa"
                    )
                )
                or clean(
                    mapping.get(
                        "structure_aa"
                    )
                ),

            "candidate_specificity_position":
                (
                    "yes"
                    if is_candidate
                    else "no"
                ),

            "classification":
                (
                    clean(
                        candidate.get(
                            "classification"
                        )
                    )
                    if candidate
                    else ""
                ),

            "sequence_priority":
                (
                    clean(
                        candidate.get(
                            "sequence_priority"
                        )
                    )
                    if candidate
                    else ""
                ),

            "distance_to_CD94_A":
                format_distance(
                    distances["CD94"]
                ),

            "contact_CD94":
                contacts["CD94"],

            "distance_to_HLA_E_A":
                format_distance(
                    distances["HLA-E"]
                ),

            "contact_HLA_E":
                contacts["HLA-E"],

            "distance_to_B2M_A":
                format_distance(
                    distances["B2M"]
                ),

            "contact_B2M":
                contacts["B2M"],

            "distance_to_peptide_A":
                format_distance(
                    distances["peptide"]
                ),

            "contact_peptide":
                contacts["peptide"],
        }

        output_rows.append(
            row
        )

    # -------------------------------------------------------------------------
    # Sort by full-length residue
    # -------------------------------------------------------------------------

    output_rows.sort(
        key=lambda row:
            int(
                row[
                    "full_length_residue"
                ]
            )
    )

    # -------------------------------------------------------------------------
    # Write output
    # -------------------------------------------------------------------------

    fieldnames = [
        "pdb_id",
        "assembly_id",
        "nkg2a_chain",
        "auth_residue_number",
        "full_length_residue",
        "nkg2a_aa",
        "candidate_specificity_position",
        "classification",
        "sequence_priority",
        "distance_to_CD94_A",
        "contact_CD94",
        "distance_to_HLA_E_A",
        "contact_HLA_E",
        "distance_to_B2M_A",
        "contact_B2M",
        "distance_to_peptide_A",
        "contact_peptide",
    ]

    write_tsv(
        OUTPUT_FILE,
        output_rows,
        fieldnames,
    )

    # -------------------------------------------------------------------------
    # Console summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("CONTACT SUMMARY")
    print("=" * 78)

    print(
        f"NKG2A residues analyzed: "
        f"{len(output_rows)}"
    )

    print(
        f"Candidate specificity residues "
        f"with coordinates: "
        f"{candidate_count}"
    )

    print()

    for component in [
        "CD94",
        "HLA-E",
        "B2M",
        "peptide",
    ]:

        print(
            f"NKG2A residues contacting "
            f"{component:<7}: "
            f"{contact_counts[component]}"
        )

    if missing_coordinate_lookup:
        print()
        print(
            "WARNING: coordinate lookup failed "
            f"for {missing_coordinate_lookup} "
            "mapped residues."
        )

    # -------------------------------------------------------------------------
    # Candidate contacts
    # -------------------------------------------------------------------------

    candidate_output = [
        row
        for row in output_rows
        if row[
            "candidate_specificity_position"
        ] == "yes"
    ]

    print()
    print("=" * 78)
    print("CANDIDATE SPECIFICITY POSITIONS WITH STRUCTURAL COORDINATES")
    print("=" * 78)

    for row in candidate_output:

        contacts_here = [
            component
            for component, field
            in [
                ("CD94", "contact_CD94"),
                ("HLA-E", "contact_HLA_E"),
                ("B2M", "contact_B2M"),
                ("peptide", "contact_peptide"),
            ]
            if row[field] == "yes"
        ]

        contact_text = (
            ", ".join(contacts_here)
            if contacts_here
            else "none"
        )

        print(
            f"{row['full_length_residue']:>3} "
            f"{row['nkg2a_aa']:<2} "
            f"{row['classification']:<40} "
            f"contacts: {contact_text}"
        )

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)
    print(OUTPUT_FILE)

    print()
    print(
        "NOTE: contact status is based on minimum "
        f"heavy-atom distance <= {CONTACT_CUTOFF:.1f} A."
    )

    print(
        "This step measures intermolecular contacts only; "
        "it does not measure solvent accessibility."
    )


if __name__ == "__main__":
    main()
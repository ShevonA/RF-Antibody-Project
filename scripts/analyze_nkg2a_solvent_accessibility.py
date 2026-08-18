from pathlib import Path
import csv

from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
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

CONTACT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_structural_contacts.tsv"
)

OUTPUT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_solvent_accessibility.tsv"
)


# Verified biological assembly 2.
CHAINS = {
    "HLA-E": "C",
    "B2M": "D",
    "CD94": "E",
    "NKG2A": "F",
    "peptide": "Q",
}


# Maximum theoretical residue SASA values in Angstrom^2.
#
# These are used to convert absolute SASA to approximate relative
# solvent accessibility (RSA).
#
# Values are from commonly used extended Gly-X-Gly reference-state
# maximum accessibility estimates.
MAX_SASA = {
    "A": 129.0,
    "R": 274.0,
    "N": 195.0,
    "D": 193.0,
    "C": 167.0,
    "Q": 225.0,
    "E": 223.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "L": 201.0,
    "K": 236.0,
    "M": 224.0,
    "F": 240.0,
    "P": 159.0,
    "S": 155.0,
    "T": 172.0,
    "W": 285.0,
    "Y": 263.0,
    "V": 174.0,
}


# RSA interpretation.
#
# These are screening categories rather than absolute biological
# boundaries.
BURIED_RSA_MAX = 0.10
PARTIAL_RSA_MAX = 0.25


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


def residue_label(residue):

    number = residue.id[1]

    insertion_code = clean(
        residue.id[2]
    )

    if insertion_code:
        return (
            f"{number}"
            f"{insertion_code}"
        )

    return str(number)


def classify_rsa(rsa):

    if rsa is None:
        return ""

    if rsa < BURIED_RSA_MAX:
        return "buried"

    if rsa < PARTIAL_RSA_MAX:
        return "partially_exposed"

    return "exposed"


# =============================================================================
# LOOKUPS
# =============================================================================

def build_contact_lookup(rows):

    lookup = {}

    for row in rows:

        full_residue = safe_int(
            row.get(
                "full_length_residue"
            )
        )

        if full_residue is None:
            continue

        lookup[
            full_residue
        ] = row

    return lookup


def get_chain_f_mapping(rows):

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

        lookup[
            (
                number,
                insertion_code,
            )
        ] = residue

    return lookup


# =============================================================================
# SASA
# =============================================================================

def calculate_complex_sasa(model):

    sr = ShrakeRupley(
        probe_radius=1.4,
        n_points=100,
    )

    sr.compute(
        model,
        level="R",
    )


def calculate_isolated_chain_sasa(chain):

    sr = ShrakeRupley(
        probe_radius=1.4,
        n_points=100,
    )

    sr.compute(
        chain,
        level="R",
    )


def get_residue_sasa(residue):

    sasa = getattr(
        residue,
        "sasa",
        None,
    )

    if sasa is None:
        return None

    return float(sasa)


def calculate_rsa(
    aa,
    sasa,
):

    if sasa is None:
        return None

    maximum = MAX_SASA.get(
        aa
    )

    if maximum is None:
        return None

    return sasa / maximum


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("STEP 2G - NKG2A SOLVENT ACCESSIBILITY ANALYSIS")
    print("=" * 78)

    print()
    print(f"PDB: {PDB_ID}")
    print("Biological assembly: 2")
    print(
        f"NKG2A chain: "
        f"{CHAINS['NKG2A']}"
    )

    # -------------------------------------------------------------------------
    # Parse structure
    # -------------------------------------------------------------------------

    parser = MMCIFParser(
        QUIET=True
    )

    structure = parser.get_structure(
        PDB_ID,
        str(CIF_FILE),
    )

    model = structure[0]

    nkg2a_chain = model[
        CHAINS["NKG2A"]
    ]

    # -------------------------------------------------------------------------
    # Mapping and contact data
    # -------------------------------------------------------------------------

    mapping_rows = get_chain_f_mapping(
        read_tsv(
            RESIDUE_MAP_FILE
        )
    )

    contact_lookup = (
        build_contact_lookup(
            read_tsv(
                CONTACT_FILE
            )
        )
    )

    coordinate_lookup = (
        build_coordinate_residue_lookup(
            nkg2a_chain
        )
    )

    print()
    print(
        "Mapped NKG2A residues: "
        f"{len(mapping_rows)}"
    )

    # -------------------------------------------------------------------------
    # Complex SASA
    # -------------------------------------------------------------------------

    print()
    print(
        "Calculating SASA in complete "
        "3CDG coordinate model..."
    )

    calculate_complex_sasa(
        model
    )

    complex_sasa = {}

    for key, residue in (
        coordinate_lookup.items()
    ):

        complex_sasa[key] = (
            get_residue_sasa(
                residue
            )
        )

    # -------------------------------------------------------------------------
    # Isolated NKG2A SASA
    #
    # Reparse structure so the complex calculation does not interfere
    # with the isolated-chain calculation.
    # -------------------------------------------------------------------------

    print(
        "Calculating SASA for isolated "
        "NKG2A chain F..."
    )

    isolated_structure = (
        parser.get_structure(
            f"{PDB_ID}_isolated",
            str(CIF_FILE),
        )
    )

    isolated_model = (
        isolated_structure[0]
    )

    isolated_chain = (
        isolated_model[
            CHAINS["NKG2A"]
        ]
    )

    calculate_isolated_chain_sasa(
        isolated_chain
    )

    isolated_lookup = (
        build_coordinate_residue_lookup(
            isolated_chain
        )
    )

    isolated_sasa = {}

    for key, residue in (
        isolated_lookup.items()
    ):

        isolated_sasa[key] = (
            get_residue_sasa(
                residue
            )
        )

    # -------------------------------------------------------------------------
    # Assemble output
    # -------------------------------------------------------------------------

    output_rows = []

    missing = 0

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

        key = (
            auth_number,
            insertion_code,
        )

        residue = coordinate_lookup.get(
            key
        )

        # Fallback if Step 2D did not preserve insertion code.
        if residue is None:

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

                key = (
                    residue.id[1],
                    clean(
                        residue.id[2]
                    ),
                )

        if residue is None:

            missing += 1
            continue

        contact = contact_lookup.get(
            full_residue,
            {},
        )

        aa = (
            clean(
                mapping.get(
                    "reference_aa"
                )
            )
            or clean(
                mapping.get(
                    "structure_aa"
                )
            )
        )

        sasa_complex = (
            complex_sasa.get(
                key
            )
        )

        sasa_isolated = (
            isolated_sasa.get(
                key
            )
        )

        rsa_complex = calculate_rsa(
            aa,
            sasa_complex,
        )

        rsa_isolated = calculate_rsa(
            aa,
            sasa_isolated,
        )

        buried_surface_area = None

        if (
            sasa_complex is not None
            and sasa_isolated is not None
        ):

            buried_surface_area = (
                sasa_isolated
                - sasa_complex
            )

        any_interface_contact = any(
            contact.get(field) == "yes"
            for field in [
                "contact_CD94",
                "contact_HLA_E",
                "contact_B2M",
                "contact_peptide",
            ]
        )

        output_rows.append(
            {
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
                    aa,

                "candidate_specificity_position":
                    clean(
                        contact.get(
                            "candidate_specificity_position"
                        )
                    ),

                "classification":
                    clean(
                        contact.get(
                            "classification"
                        )
                    ),

                "sequence_priority":
                    clean(
                        contact.get(
                            "sequence_priority"
                        )
                    ),

                "complex_sasa_A2":
                    (
                        ""
                        if sasa_complex is None
                        else f"{sasa_complex:.3f}"
                    ),

                "complex_rsa":
                    (
                        ""
                        if rsa_complex is None
                        else f"{rsa_complex:.4f}"
                    ),

                "complex_exposure_class":
                    classify_rsa(
                        rsa_complex
                    ),

                "isolated_nkg2a_sasa_A2":
                    (
                        ""
                        if sasa_isolated is None
                        else f"{sasa_isolated:.3f}"
                    ),

                "isolated_nkg2a_rsa":
                    (
                        ""
                        if rsa_isolated is None
                        else f"{rsa_isolated:.4f}"
                    ),

                "isolated_exposure_class":
                    classify_rsa(
                        rsa_isolated
                    ),

                "buried_surface_area_A2":
                    (
                        ""
                        if buried_surface_area is None
                        else (
                            f"{buried_surface_area:.3f}"
                        )
                    ),

                "contact_CD94":
                    clean(
                        contact.get(
                            "contact_CD94"
                        )
                    ),

                "contact_HLA_E":
                    clean(
                        contact.get(
                            "contact_HLA_E"
                        )
                    ),

                "contact_B2M":
                    clean(
                        contact.get(
                            "contact_B2M"
                        )
                    ),

                "contact_peptide":
                    clean(
                        contact.get(
                            "contact_peptide"
                        )
                    ),

                "any_interface_contact":
                    (
                        "yes"
                        if any_interface_contact
                        else "no"
                    ),
            }
        )

    output_rows.sort(
        key=lambda row:
            int(
                row[
                    "full_length_residue"
                ]
            )
    )

    # -------------------------------------------------------------------------
    # Write table
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
        "complex_sasa_A2",
        "complex_rsa",
        "complex_exposure_class",
        "isolated_nkg2a_sasa_A2",
        "isolated_nkg2a_rsa",
        "isolated_exposure_class",
        "buried_surface_area_A2",
        "contact_CD94",
        "contact_HLA_E",
        "contact_B2M",
        "contact_peptide",
        "any_interface_contact",
    ]

    write_tsv(
        OUTPUT_FILE,
        output_rows,
        fieldnames,
    )

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------

    candidates = [
        row
        for row in output_rows
        if row[
            "candidate_specificity_position"
        ] == "yes"
    ]

    exposed_candidates = [
        row
        for row in candidates
        if row[
            "complex_exposure_class"
        ] == "exposed"
    ]

    exposed_noninterface = [
        row
        for row in exposed_candidates
        if row[
            "any_interface_contact"
        ] == "no"
    ]

    print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)

    print(
        f"Residues analyzed: "
        f"{len(output_rows)}"
    )

    print(
        f"Candidate specificity residues: "
        f"{len(candidates)}"
    )

    print(
        f"Exposed candidate residues: "
        f"{len(exposed_candidates)}"
    )

    print(
        "Exposed candidates without "
        "4.5 A interface contact: "
        f"{len(exposed_noninterface)}"
    )

    if missing:

        print(
            f"WARNING: {missing} mapped residues "
            "could not be found in chain F."
        )

    print()
    print("=" * 78)
    print("EXPOSED NON-INTERFACE CANDIDATES")
    print("=" * 78)

    for row in exposed_noninterface:

        print(
            f"{row['full_length_residue']:>3} "
            f"{row['nkg2a_aa']:<2} "
            f"RSA={row['complex_rsa']:<7} "
            f"{row['classification']:<40} "
            f"{row['sequence_priority']}"
        )

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)
    print(OUTPUT_FILE)

    print()
    print(
        "NOTE: RSA categories are screening classifications:"
    )

    print(
        f"  buried:            RSA < {BURIED_RSA_MAX:.2f}"
    )

    print(
        f"  partially exposed: {BURIED_RSA_MAX:.2f} <= "
        f"RSA < {PARTIAL_RSA_MAX:.2f}"
    )

    print(
        f"  exposed:           RSA >= {PARTIAL_RSA_MAX:.2f}"
    )

    print()
    print(
        "Complex SASA measures accessibility in the complete "
        "3CDG coordinate model."
    )

    print(
        "Isolated-chain SASA measures accessibility after removing "
        "CD94/HLA-E/B2M/peptide neighbors."
    )


if __name__ == "__main__":
    main()
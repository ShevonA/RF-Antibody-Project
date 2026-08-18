from pathlib import Path
import csv

from Bio.PDB.MMCIF2Dict import MMCIF2Dict


ROOT = Path(__file__).resolve().parent.parent

STRUCTURE_DIR = (
    ROOT
    / "structures"
    / "reference"
)

OUTPUT_DIR = (
    ROOT
    / "results"
    / "tables"
    / "structure"
)

OUTPUT_FILE = (
    OUTPUT_DIR
    / "structural_chain_inventory.tsv"
)

PDB_IDS = [
    "3BDW",
    "3CDG",
]


def as_list(value):
    """
    MMCIF2Dict sometimes returns a single string and sometimes
    a list. Normalize everything to a list.
    """

    if value is None:
        return []

    if isinstance(value, list):
        return value

    return [value]


def clean(value):
    """
    Convert mmCIF missing-value markers to empty strings.
    """

    if value in {
        None,
        ".",
        "?",
    }:
        return ""

    return str(value)


def classify_component(
    description,
    polymer_type,
):
    """
    Assign a useful biological component label from the
    mmCIF entity description.
    """

    text = description.lower()

    if (
        "nkg2a" in text
        or "nkg2-a" in text
        or "klrc1" in text
        or "killer cell lectin-like receptor subfamily c member 1"
        in text
    ):
        return "NKG2A"

    if (
        "cd94" in text
        or "klrd1" in text
        or "killer cell lectin-like receptor subfamily d member 1"
        in text
    ):
        return "CD94"

    if (
        "hla-e" in text
        or "hla e" in text
        or "histocompatibility antigen, alpha chain e"
        in text
        or "histocompatibility antigen alpha chain e"
        in text
    ):
        return "HLA-E"

    if (
        "beta-2-microglobulin" in text
        or "beta-2 microglobulin" in text
        or "beta-2-microglobulin" in text
    ):
        return "B2M"

    if (
        "leader peptide" in text
        or polymer_type == "polypeptide(L)"
    ):
        return "peptide"

    return "other"


def get_entity_descriptions(cif):
    """
    Return:
        entity_id -> description
    """

    entity_ids = as_list(
        cif.get("_entity.id")
    )

    descriptions = as_list(
        cif.get(
            "_entity.pdbx_description"
        )
    )

    result = {}

    for index, entity_id in enumerate(
        entity_ids
    ):

        description = ""

        if index < len(descriptions):
            description = clean(
                descriptions[index]
            )

        result[
            clean(entity_id)
        ] = description

    return result


def inspect_structure(pdb_id):

    path = (
        STRUCTURE_DIR
        / f"{pdb_id}.cif"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Missing structure: {path}"
        )

    print()
    print("=" * 78)
    print(f"PDB {pdb_id}")
    print("=" * 78)

    cif = MMCIF2Dict(
        str(path)
    )

    descriptions = (
        get_entity_descriptions(cif)
    )

    entity_ids = as_list(
        cif.get(
            "_entity_poly.entity_id"
        )
    )

    polymer_types = as_list(
        cif.get(
            "_entity_poly.type"
        )
    )

    strand_ids = as_list(
        cif.get(
            "_entity_poly.pdbx_strand_id"
        )
    )

    sequences = as_list(
        cif.get(
            "_entity_poly.pdbx_seq_one_letter_code_can"
        )
    )

    rows = []

    print("\nPOLYMER ENTITIES")

    for index, entity_id in enumerate(
        entity_ids
    ):

        entity_id = clean(
            entity_id
        )

        polymer_type = ""

        if index < len(polymer_types):
            polymer_type = clean(
                polymer_types[index]
            )

        chains = ""

        if index < len(strand_ids):
            chains = clean(
                strand_ids[index]
            )

        sequence = ""

        if index < len(sequences):
            sequence = clean(
                sequences[index]
            )

            sequence = "".join(
                sequence.split()
            )

        description = (
            descriptions.get(
                entity_id,
                "",
            )
        )

        component = (
            classify_component(
                description,
                polymer_type,
            )
        )

        sequence_length = len(
            sequence
        )

        print()
        print(
            f"Entity:      {entity_id}"
        )
        print(
            f"Component:   {component}"
        )
        print(
            f"Description: {description}"
        )
        print(
            f"Polymer:     {polymer_type}"
        )
        print(
            f"Chains:      {chains}"
        )
        print(
            f"Seq length:  {sequence_length}"
        )

        rows.append(
            {
                "pdb_id": pdb_id,
                "entity_id": entity_id,
                "component": component,
                "description": description,
                "polymer_type": polymer_type,
                "chain_ids": chains,
                "sequence_length": sequence_length,
                "canonical_sequence": sequence,
            }
        )

    return rows


def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 78)
    print(
        "STEP 2C - STRUCTURAL CHAIN / ENTITY INSPECTION"
    )
    print("=" * 78)

    all_rows = []

    for pdb_id in PDB_IDS:

        rows = inspect_structure(
            pdb_id
        )

        all_rows.extend(rows)

    fields = [
        "pdb_id",
        "entity_id",
        "component",
        "description",
        "polymer_type",
        "chain_ids",
        "sequence_length",
        "canonical_sequence",
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
            all_rows
        )

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)

    print(
        OUTPUT_FILE
    )

    print(
        f"\nPolymer entities written: "
        f"{len(all_rows)}"
    )


if __name__ == "__main__":
    main()
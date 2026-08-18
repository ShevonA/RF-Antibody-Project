from pathlib import Path
import csv

from Bio.PDB.MMCIF2Dict import MMCIF2Dict


ROOT = Path(__file__).resolve().parent.parent

CIF_FILE = (
    ROOT
    / "structures"
    / "reference"
    / "3CDG.cif"
)

OUTPUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "3CDG_biological_assembly_chains.tsv"
)


def as_list(value):

    if value is None:
        return []

    if isinstance(value, list):
        return value

    return [value]


def clean(value):

    if value in {
        None,
        ".",
        "?",
    }:
        return ""

    return str(value)


def classify_component(description):

    text = description.lower()

    if (
        "nkg2a" in text
        or "nkg2-a" in text
    ):
        return "NKG2A"

    if "cd94" in text:
        return "CD94"

    if (
        "alpha chain e" in text
        or "hla-e" in text
    ):
        return "HLA-E"

    if "beta-2-microglobulin" in text:
        return "B2M"

    if "leader peptide" in text:
        return "peptide"

    return "other"


def main():

    print("=" * 78)
    print("STEP 2E - VERIFY 3CDG BIOLOGICAL ASSEMBLIES")
    print("=" * 78)

    cif = MMCIF2Dict(
        str(CIF_FILE)
    )

    # ---------------------------------------------------------
    # Entity descriptions
    # ---------------------------------------------------------

    entity_ids = as_list(
        cif.get("_entity.id")
    )

    descriptions = as_list(
        cif.get(
            "_entity.pdbx_description"
        )
    )

    entity_description = {}

    for index, entity_id in enumerate(
        entity_ids
    ):

        description = ""

        if index < len(descriptions):
            description = clean(
                descriptions[index]
            )

        entity_description[
            clean(entity_id)
        ] = description

    # ---------------------------------------------------------
    # label_asym_id -> entity
    # ---------------------------------------------------------

    struct_asym_ids = as_list(
        cif.get(
            "_struct_asym.id"
        )
    )

    struct_asym_entities = as_list(
        cif.get(
            "_struct_asym.entity_id"
        )
    )

    label_to_entity = {}

    for label_id, entity_id in zip(
        struct_asym_ids,
        struct_asym_entities,
    ):

        label_to_entity[
            clean(label_id)
        ] = clean(entity_id)

    # ---------------------------------------------------------
    # label_asym_id -> author chain ID
    #
    # Extract this from atom_site records.
    # ---------------------------------------------------------

    atom_label_asym = as_list(
        cif.get(
            "_atom_site.label_asym_id"
        )
    )

    atom_auth_asym = as_list(
        cif.get(
            "_atom_site.auth_asym_id"
        )
    )

    label_to_auth = {}

    for label_id, auth_id in zip(
        atom_label_asym,
        atom_auth_asym,
    ):

        label_id = clean(
            label_id
        )

        auth_id = clean(
            auth_id
        )

        if label_id and auth_id:

            if label_id not in label_to_auth:
                label_to_auth[
                    label_id
                ] = auth_id

            elif (
                label_to_auth[label_id]
                != auth_id
            ):
                raise ValueError(
                    f"label asym ID {label_id} "
                    "maps to multiple author chains."
                )

    # ---------------------------------------------------------
    # Biological assemblies
    # ---------------------------------------------------------

    assembly_ids = as_list(
        cif.get(
            "_pdbx_struct_assembly_gen.assembly_id"
        )
    )

    asym_lists = as_list(
        cif.get(
            "_pdbx_struct_assembly_gen.asym_id_list"
        )
    )

    oper_expressions = as_list(
        cif.get(
            "_pdbx_struct_assembly_gen.oper_expression"
        )
    )

    rows = []

    for index, assembly_id in enumerate(
        assembly_ids
    ):

        asym_string = clean(
            asym_lists[index]
        )

        operation = ""

        if index < len(
            oper_expressions
        ):
            operation = clean(
                oper_expressions[index]
            )

        label_ids = [
            item.strip()
            for item in asym_string.split(",")
            if item.strip()
        ]

        print()
        print(
            f"Assembly {assembly_id}"
        )
        print("-" * 78)

        for label_id in label_ids:

            entity_id = (
                label_to_entity.get(
                    label_id,
                    "",
                )
            )

            auth_chain = (
                label_to_auth.get(
                    label_id,
                    "",
                )
            )

            description = (
                entity_description.get(
                    entity_id,
                    "",
                )
            )

            component = (
                classify_component(
                    description
                )
            )

            print(
                f"label={label_id:<3} "
                f"auth={auth_chain:<3} "
                f"entity={entity_id:<3} "
                f"{component:<8} "
                f"{description}"
            )

            rows.append(
                {
                    "assembly_id":
                        clean(assembly_id),
                    "operation":
                        operation,
                    "label_asym_id":
                        label_id,
                    "author_chain_id":
                        auth_chain,
                    "entity_id":
                        entity_id,
                    "component":
                        component,
                    "description":
                        description,
                }
            )

    # ---------------------------------------------------------
    # Write output
    # ---------------------------------------------------------

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fields = [
        "assembly_id",
        "operation",
        "label_asym_id",
        "author_chain_id",
        "entity_id",
        "component",
        "description",
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
            rows
        )

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)

    print(OUTPUT)


if __name__ == "__main__":
    main()
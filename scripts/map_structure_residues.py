from pathlib import Path
import csv

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import PairwiseAligner


ROOT = Path(__file__).resolve().parent.parent

STRUCTURE_DIR = ROOT / "structures" / "reference"

OUTPUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_structure_residue_map.tsv"
)


# ---------------------------------------------------------
# Reference sequence
# ---------------------------------------------------------

HUMAN_NKG2A_ECTODOMAIN = (
    "PSTLIQRHNNSSLNTRTQK"
    "ARHCGHCPEEWITYSNSCYYIGKERRTWEESLLACTSKNSSLLSIDNEEEM"
    "KFLSIISPSSWIGVFRNSSHHPWVTMNGLAFKHEIKDSDNAELNCAVLQVN"
    "RLKSAQCGSSIIYHCKHKL"
)

ECTODOMAIN_START = 94


STRUCTURES = {
    "3BDW": {
        "file": STRUCTURE_DIR / "3BDW.cif",
        "nkg2a_chains": ["B", "D"],
    },
    "3CDG": {
        "file": STRUCTURE_DIR / "3CDG.cif",
        "nkg2a_chains": ["K", "F"],
    },
}


# ---------------------------------------------------------
# Structural entity sequence
# ---------------------------------------------------------

STRUCTURE_SEQUENCE = (
    "ARHCGHCPEEWITYSNSCYYIGKERRTWEESLLACTSKNSSLLSIDNEEEM"
    "KFLSIISPSSWIGVFRNSSHHPWVTMNGLAFKHEIKDSDNAELNCAVLQVN"
    "RLKSAQCGSSIIYHCKHKL"
)


# ---------------------------------------------------------
# Amino-acid conversion
# ---------------------------------------------------------

THREE_TO_ONE = {
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
}


def residue_to_one_letter(residue):
    return THREE_TO_ONE.get(
        residue.get_resname().upper(),
        "X",
    )


# ---------------------------------------------------------
# Locate crystallographic sequence in reference ectodomain
# ---------------------------------------------------------

def locate_structure_sequence():

    index = HUMAN_NKG2A_ECTODOMAIN.find(
        STRUCTURE_SEQUENCE
    )

    if index == -1:
        raise ValueError(
            "Structural NKG2A sequence was not found "
            "inside the human NKG2A ectodomain."
        )

    full_start = (
        ECTODOMAIN_START
        + index
    )

    full_end = (
        full_start
        + len(STRUCTURE_SEQUENCE)
        - 1
    )

    return index, full_start, full_end


# ---------------------------------------------------------
# Build coordinate-to-entity sequence mapping
# ---------------------------------------------------------

def align_coordinate_sequence(
    pdb_id,
    chain_id,
    coordinate_sequence,
):

    aligner = PairwiseAligner()

    aligner.mode = "global"

    # Favor exact amino-acid matches.
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0

    # Permit missing coordinate residues without shifting
    # all downstream residue numbers.
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(
        STRUCTURE_SEQUENCE,
        coordinate_sequence,
    )

    if len(alignments) == 0:
        raise ValueError(
            f"{pdb_id} chain {chain_id}: "
            "no sequence alignment was produced."
        )

    alignment = alignments[0]

    coordinate_to_structure = {}

    target_blocks = alignment.aligned[0]
    query_blocks = alignment.aligned[1]

    for target_block, query_block in zip(
        target_blocks,
        query_blocks,
    ):

        target_start, target_end = target_block
        query_start, query_end = query_block

        target_length = (
            target_end
            - target_start
        )

        query_length = (
            query_end
            - query_start
        )

        if target_length != query_length:
            raise ValueError(
                f"{pdb_id} chain {chain_id}: "
                "unexpected unequal aligned block lengths."
            )

        for offset in range(target_length):

            structure_index = (
                target_start
                + offset
            )

            coordinate_index = (
                query_start
                + offset
            )

            coordinate_to_structure[
                coordinate_index
            ] = structure_index

    return alignment, coordinate_to_structure


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    print("=" * 78)
    print("STEP 2D - MAP STRUCTURAL NKG2A RESIDUES")
    print("=" * 78)

    (
        structure_ecto_index,
        structure_full_start,
        structure_full_end,
    ) = locate_structure_sequence()

    ectodomain_end = (
        ECTODOMAIN_START
        + len(HUMAN_NKG2A_ECTODOMAIN)
        - 1
    )

    print()
    print(
        "Human NKG2A ectodomain: "
        f"{ECTODOMAIN_START}-{ectodomain_end}"
    )

    print(
        "Structural NKG2A sequence maps to: "
        f"{structure_full_start}-{structure_full_end}"
    )

    print(
        "Structural sequence length: "
        f"{len(STRUCTURE_SEQUENCE)} aa"
    )

    print(
        "Structural sequence begins at ectodomain "
        f"residue {structure_ecto_index + 1}"
    )

    rows = []

    parser = MMCIFParser(
        QUIET=True
    )

    # -----------------------------------------------------
    # Process each structure
    # -----------------------------------------------------

    for pdb_id, info in STRUCTURES.items():

        print()
        print("-" * 78)
        print(f"PDB {pdb_id}")
        print("-" * 78)

        structure_file = info["file"]

        if not structure_file.exists():
            raise FileNotFoundError(
                f"Structure file not found: "
                f"{structure_file}"
            )

        structure = parser.get_structure(
            pdb_id,
            structure_file,
        )

        model = structure[0]

        # -------------------------------------------------
        # Process each NKG2A chain
        # -------------------------------------------------

        for chain_id in info["nkg2a_chains"]:

            if chain_id not in model:
                raise KeyError(
                    f"{pdb_id}: expected NKG2A "
                    f"chain {chain_id} was not found."
                )

            chain = model[chain_id]

            residues = [
                residue
                for residue in chain
                if is_aa(
                    residue,
                    standard=False,
                )
            ]

            coordinate_sequence = "".join(
                residue_to_one_letter(residue)
                for residue in residues
            )

            print()
            print(f"Chain {chain_id}")

            print(
                "Coordinate amino-acid residues: "
                f"{len(residues)}"
            )

            if not residues:
                raise ValueError(
                    f"{pdb_id} chain {chain_id}: "
                    "no amino-acid coordinate residues found."
                )

            # -------------------------------------------------
            # Align coordinate sequence to entity sequence
            # -------------------------------------------------

            (
                alignment,
                coordinate_to_structure,
            ) = align_coordinate_sequence(
                pdb_id,
                chain_id,
                coordinate_sequence,
            )

            print(
                "Best alignment score: "
                f"{alignment.score:.1f}"
            )

            # Every coordinate residue should map to an
            # entity-sequence residue.
            unmapped_coordinate_indices = [
                i
                for i in range(len(residues))
                if i not in coordinate_to_structure
            ]

            if unmapped_coordinate_indices:
                raise ValueError(
                    f"{pdb_id} chain {chain_id}: "
                    f"{len(unmapped_coordinate_indices)} "
                    "coordinate residue(s) could not be mapped "
                    "to the structural entity sequence."
                )

            mapped_structure_positions = [
                coordinate_to_structure[i] + 1
                for i in range(len(residues))
            ]

            print(
                "Mapped structural sequence positions: "
                f"{min(mapped_structure_positions)}-"
                f"{max(mapped_structure_positions)}"
            )

            mapped_full_positions = [
                structure_full_start
                + position
                - 1
                for position
                in mapped_structure_positions
            ]

            print(
                "Mapped full-length residues: "
                f"{min(mapped_full_positions)}-"
                f"{max(mapped_full_positions)}"
            )

            # -------------------------------------------------
            # Entity residues lacking coordinates
            # -------------------------------------------------

            mapped_structure_indices = set(
                coordinate_to_structure.values()
            )

            missing_structure_positions = [
                i + 1
                for i in range(
                    len(STRUCTURE_SEQUENCE)
                )
                if i not in mapped_structure_indices
            ]

            if missing_structure_positions:
                print(
                    "Entity residues without coordinates: "
                    + ", ".join(
                        str(position)
                        for position
                        in missing_structure_positions
                    )
                )
            else:
                print(
                    "Entity residues without coordinates: none"
                )

            # -------------------------------------------------
            # Generate individual residue mappings
            # -------------------------------------------------

            for coordinate_index, residue in enumerate(
                residues
            ):

                structure_index = (
                    coordinate_to_structure[
                        coordinate_index
                    ]
                )

                structure_sequence_position = (
                    structure_index
                    + 1
                )

                full_length_residue = (
                    structure_full_start
                    + structure_index
                )

                ectodomain_residue = (
                    full_length_residue
                    - ECTODOMAIN_START
                    + 1
                )

                structure_aa = (
                    residue_to_one_letter(
                        residue
                    )
                )

                entity_aa = (
                    STRUCTURE_SEQUENCE[
                        structure_index
                    ]
                )

                reference_aa = (
                    HUMAN_NKG2A_ECTODOMAIN[
                        ectodomain_residue
                        - 1
                    ]
                )

                (
                    hetflag,
                    auth_residue_number,
                    insertion_code,
                ) = residue.id

                insertion_code = (
                    insertion_code.strip()
                )

                rows.append(
                    {
                        "pdb_id": pdb_id,
                        "chain_id": chain_id,
                        "coordinate_index":
                            coordinate_index + 1,
                        "auth_residue_number":
                            auth_residue_number,
                        "insertion_code":
                            insertion_code,
                        "structure_sequence_position":
                            structure_sequence_position,
                        "ectodomain_residue":
                            ectodomain_residue,
                        "full_length_residue":
                            full_length_residue,
                        "structure_aa":
                            structure_aa,
                        "entity_aa":
                            entity_aa,
                        "reference_aa":
                            reference_aa,
                        "sequence_match": (
                            "yes"
                            if (
                                structure_aa
                                == entity_aa
                                == reference_aa
                            )
                            else "no"
                        ),
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
        "pdb_id",
        "chain_id",
        "coordinate_index",
        "auth_residue_number",
        "insertion_code",
        "structure_sequence_position",
        "ectodomain_residue",
        "full_length_residue",
        "structure_aa",
        "entity_aa",
        "reference_aa",
        "sequence_match",
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
        writer.writerows(rows)

    # ---------------------------------------------------------
    # QC summary
    # ---------------------------------------------------------

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)

    print(OUTPUT)

    print()
    print(
        f"Residue mapping rows: {len(rows)}"
    )

    mismatches = [
        row
        for row in rows
        if row["sequence_match"] != "yes"
    ]

    print(
        f"Sequence mismatches: {len(mismatches)}"
    )

    if mismatches:

        print()
        print(
            "WARNING: sequence mismatches detected."
        )

        for row in mismatches:

            print(
                f"{row['pdb_id']} "
                f"chain {row['chain_id']} "
                f"auth {row['auth_residue_number']}: "
                f"coordinate={row['structure_aa']} "
                f"entity={row['entity_aa']} "
                f"reference={row['reference_aa']} "
                f"(full-length "
                f"{row['full_length_residue']})"
            )


if __name__ == "__main__":
    main()
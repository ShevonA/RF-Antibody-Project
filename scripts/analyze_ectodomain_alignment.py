from pathlib import Path
from itertools import combinations
import csv


ROOT = Path(__file__).resolve().parent.parent

ALIGNMENT = (
    ROOT
    / "alignments"
    / "primary_ectodomain_aligned.fasta"
)

IDENTITY_OUT = (
    ROOT
    / "results"
    / "tables"
    / "ectodomain_pairwise_identity.tsv"
)

DIFF_OUT = (
    ROOT
    / "results"
    / "tables"
    / "ectodomain_within_species_differences.tsv"
)

MAP_OUT = (
    ROOT
    / "results"
    / "tables"
    / "ectodomain_alignment_position_map.tsv"
)


WITHIN_SPECIES_PAIRS = [
    (
        "human",
        "human_NKG2A|P26715|ecto_94-233",
        "human_NKG2C|NP_002251.2|ecto_94-231",
    ),
    (
        "rhesus_iso1",
        "rhesus_NKG2A|NP_001028001.3|ecto_94-233",
        "rhesus_NKG2C_isoform1|NP_001305287.3|ecto_94-231",
    ),
    (
        "rhesus_iso2",
        "rhesus_NKG2A|NP_001028001.3|ecto_94-233",
        "rhesus_NKG2C_isoform2|NP_001098647.3|ecto_94-231",
    ),
    (
        "pigtail",
        "pigtail_NKG2A|XP_070928357.1|ecto_94-233",
        "pigtail_NKG2C|XP_070928345.1|ecto_94-231",
    ),
]


def read_alignment(path):
    """
    Read a FASTA alignment into a dictionary.

    utf-8-sig is used so the script works even if
    PowerShell added a UTF-8 byte-order mark.
    """

    records = {}

    current = None
    pieces = []

    for line in path.read_text(
        encoding="utf-8-sig"
    ).splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if current is not None:
                records[current] = "".join(pieces)

            current = line[1:]
            pieces = []

        else:
            pieces.append(line)

    if current is not None:
        records[current] = "".join(pieces)

    return records


def aligned_identity(seq1, seq2):
    """
    Calculate percent identity over alignment columns
    in which both sequences contain an amino acid.

    Columns containing a gap in either sequence are
    excluded from the denominator.
    """

    matches = 0
    compared = 0

    for aa1, aa2 in zip(seq1, seq2):

        if aa1 == "-" or aa2 == "-":
            continue

        compared += 1

        if aa1 == aa2:
            matches += 1

    identity = (
        100.0 * matches / compared
        if compared
        else 0.0
    )

    return matches, compared, identity


def residue_positions(
    aligned_sequence,
    full_length_start=94,
):
    """
    For every alignment column, calculate:

    1. Ectodomain-relative residue number
       Example: first ectodomain residue = 1

    2. Original full-length protein residue number
       Example: first ectodomain residue = 94

    Gap columns receive None.
    """

    ecto_positions = []
    full_positions = []

    ecto_number = 0

    for aa in aligned_sequence:

        if aa == "-":

            ecto_positions.append(None)
            full_positions.append(None)

        else:

            ecto_number += 1

            ecto_positions.append(
                ecto_number
            )

            full_positions.append(
                full_length_start
                + ecto_number
                - 1
            )

    return ecto_positions, full_positions


def main():

    print("=" * 78)
    print("STEP 1D - ECTODOMAIN ALIGNMENT ANALYSIS")
    print("=" * 78)

    if not ALIGNMENT.exists():
        raise FileNotFoundError(
            f"Alignment file not found:\n{ALIGNMENT}"
        )

    records = read_alignment(
        ALIGNMENT
    )

    print(
        f"\nSequences loaded: {len(records)}"
    )

    print("\nSequence headers:")

    for name in records:
        print(f"  - {name}")

    # -----------------------------------------------------
    # Validate alignment
    # -----------------------------------------------------

    alignment_lengths = {
        len(sequence)
        for sequence in records.values()
    }

    if len(alignment_lengths) != 1:
        raise ValueError(
            "Aligned sequences have unequal lengths."
        )

    alignment_length = next(
        iter(alignment_lengths)
    )

    print(
        f"\nAlignment length: "
        f"{alignment_length} columns"
    )

    expected_names = set()

    for _, nkg2a_name, nkg2c_name in (
        WITHIN_SPECIES_PAIRS
    ):

        expected_names.add(
            nkg2a_name
        )

        expected_names.add(
            nkg2c_name
        )

    missing_names = (
        expected_names
        - set(records.keys())
    )

    if missing_names:

        print(
            "\nERROR: Expected sequence headers "
            "are missing:"
        )

        for name in sorted(
            missing_names
        ):
            print(f"  - {name}")

        raise KeyError(
            "One or more required sequences "
            "are missing from the alignment."
        )

    # -----------------------------------------------------
    # Build ectodomain and full-length residue maps
    # -----------------------------------------------------

    ecto_position_maps = {}
    full_position_maps = {}

    for name, sequence in records.items():

        ecto_positions, full_positions = (
            residue_positions(
                sequence,
                full_length_start=94,
            )
        )

        ecto_position_maps[
            name
        ] = ecto_positions

        full_position_maps[
            name
        ] = full_positions

    # -----------------------------------------------------
    # Alignment position map
    # -----------------------------------------------------

    map_fields = [
        "alignment_column"
    ]

    for name in records:

        map_fields.extend(
            [
                f"{name}_aa",
                f"{name}_ecto_residue",
                f"{name}_full_residue",
            ]
        )

    map_rows = []

    for column_index in range(
        alignment_length
    ):

        row = {
            "alignment_column":
                column_index + 1
        }

        for name, sequence in (
            records.items()
        ):

            aa = sequence[
                column_index
            ]

            ecto_residue = (
                ecto_position_maps[
                    name
                ][column_index]
            )

            full_residue = (
                full_position_maps[
                    name
                ][column_index]
            )

            row[
                f"{name}_aa"
            ] = aa

            row[
                f"{name}_ecto_residue"
            ] = (
                ""
                if ecto_residue is None
                else ecto_residue
            )

            row[
                f"{name}_full_residue"
            ] = (
                ""
                if full_residue is None
                else full_residue
            )

        map_rows.append(
            row
        )

    MAP_OUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with MAP_OUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=map_fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            map_rows
        )

    # -----------------------------------------------------
    # Pairwise identity
    # -----------------------------------------------------

    identity_rows = []

    print(
        "\nPAIRWISE ALIGNMENT-BASED IDENTITY"
    )

    for name1, name2 in combinations(
        records.keys(),
        2,
    ):

        matches, compared, identity = (
            aligned_identity(
                records[name1],
                records[name2],
            )
        )

        identity_rows.append(
            {
                "sequence_1":
                    name1,
                "sequence_2":
                    name2,
                "matches":
                    matches,
                "positions_compared":
                    compared,
                "identity_percent":
                    f"{identity:.3f}",
            }
        )

    with IDENTITY_OUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sequence_1",
                "sequence_2",
                "matches",
                "positions_compared",
                "identity_percent",
            ],
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            identity_rows
        )

    # -----------------------------------------------------
    # Within-species NKG2A vs NKG2C differences
    # -----------------------------------------------------

    diff_rows = []

    print(
        "\nWITHIN-SPECIES "
        "NKG2A vs NKG2C"
    )

    for (
        species,
        nkg2a_name,
        nkg2c_name,
    ) in WITHIN_SPECIES_PAIRS:

        seq_a = records[
            nkg2a_name
        ]

        seq_c = records[
            nkg2c_name
        ]

        ecto_pos_a = (
            ecto_position_maps[
                nkg2a_name
            ]
        )

        ecto_pos_c = (
            ecto_position_maps[
                nkg2c_name
            ]
        )

        full_pos_a = (
            full_position_maps[
                nkg2a_name
            ]
        )

        full_pos_c = (
            full_position_maps[
                nkg2c_name
            ]
        )

        (
            matches,
            compared,
            identity,
        ) = aligned_identity(
            seq_a,
            seq_c,
        )

        print(
            f"{species:<15} "
            f"{identity:>7.3f}% identity "
            f"({matches}/{compared})"
        )

        for i, (
            aa_a,
            aa_c,
        ) in enumerate(
            zip(
                seq_a,
                seq_c,
            )
        ):

            if aa_a == aa_c:
                continue

            diff_rows.append(
                {
                    "comparison":
                        species,

                    "alignment_column":
                        i + 1,

                    "NKG2A_ecto_residue":
                        (
                            ""
                            if ecto_pos_a[i]
                            is None
                            else ecto_pos_a[i]
                        ),

                    "NKG2A_full_residue":
                        (
                            ""
                            if full_pos_a[i]
                            is None
                            else full_pos_a[i]
                        ),

                    "NKG2A_aa":
                        aa_a,

                    "NKG2C_ecto_residue":
                        (
                            ""
                            if ecto_pos_c[i]
                            is None
                            else ecto_pos_c[i]
                        ),

                    "NKG2C_full_residue":
                        (
                            ""
                            if full_pos_c[i]
                            is None
                            else full_pos_c[i]
                        ),

                    "NKG2C_aa":
                        aa_c,

                    "difference_type":
                        (
                            "gap"
                            if (
                                aa_a == "-"
                                or aa_c == "-"
                            )
                            else
                            "substitution"
                        ),
                }
            )

    with DIFF_OUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "comparison",
                "alignment_column",
                "NKG2A_ecto_residue",
                "NKG2A_full_residue",
                "NKG2A_aa",
                "NKG2C_ecto_residue",
                "NKG2C_full_residue",
                "NKG2C_aa",
                "difference_type",
            ],
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            diff_rows
        )

    # -----------------------------------------------------
    # Summary
    # -----------------------------------------------------

    print("\nOutputs:")

    print(
        IDENTITY_OUT
    )

    print(
        DIFF_OUT
    )

    print(
        MAP_OUT
    )

    print(
        f"\nTotal within-species "
        f"difference rows: "
        f"{len(diff_rows)}"
    )


if __name__ == "__main__":
    main()
from pathlib import Path
from itertools import combinations
import csv


ROOT = Path(__file__).resolve().parent.parent

ALIGNMENT = (
    ROOT
    / "alignments"
    / "primary_full_length_aligned.fasta"
)

IDENTITY_OUT = (
    ROOT
    / "results"
    / "tables"
    / "primary_alignment_pairwise_identity.tsv"
)

DIFF_OUT = (
    ROOT
    / "results"
    / "tables"
    / "within_species_NKG2A_vs_NKG2C_differences.tsv"
)

MAP_OUT = (
    ROOT
    / "results"
    / "tables"
    / "alignment_position_map.tsv"
)


WITHIN_SPECIES_PAIRS = [
    (
        "human",
        "human_NKG2A_P26715",
        "human_NKG2C_NP_002251.2",
    ),
    (
        "rhesus_iso1",
        "rhesus_NKG2A_NP_001028001.3",
        "rhesus_NKG2C_iso1_NP_001305287.3",
    ),
    (
        "rhesus_iso2",
        "rhesus_NKG2A_NP_001028001.3",
        "rhesus_NKG2C_iso2_NP_001098647.3",
    ),
    (
        "pigtail",
        "pigtail_NKG2A_XP_070928357.1",
        "pigtail_NKG2C_XP_070928345.1",
    ),
]


def read_alignment(path):
    records = {}

    current = None
    pieces = []

    for line in path.read_text(
        encoding="utf-8"
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
    Identity calculated only over columns where
    both sequences contain an amino acid.
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


def residue_positions(aligned_sequence):
    """
    For every alignment column, return the original
    residue number. Gap columns receive None.
    """

    positions = []

    residue_number = 0

    for aa in aligned_sequence:

        if aa == "-":
            positions.append(None)

        else:
            residue_number += 1
            positions.append(residue_number)

    return positions


def main():

    records = read_alignment(ALIGNMENT)

    print("=" * 78)
    print("STEP 1D - PRIMARY ALIGNMENT ANALYSIS")
    print("=" * 78)

    print(
        f"\nSequences loaded: {len(records)}"
    )

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
        f"Alignment length: {alignment_length} columns"
    )

    # -----------------------------------------------------
    # Alignment position map
    # -----------------------------------------------------

    position_maps = {
        name: residue_positions(sequence)
        for name, sequence in records.items()
    }

    map_fields = [
        "alignment_column"
    ]

    for name in records:
        map_fields.extend(
            [
                f"{name}_aa",
                f"{name}_residue",
            ]
        )

    map_rows = []

    for column_index in range(alignment_length):

        row = {
            "alignment_column":
                column_index + 1
        }

        for name, sequence in records.items():

            aa = sequence[column_index]
            residue_number = (
                position_maps[name][column_index]
            )

            row[f"{name}_aa"] = aa

            row[f"{name}_residue"] = (
                ""
                if residue_number is None
                else residue_number
            )

        map_rows.append(row)

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
        writer.writerows(map_rows)

    # -----------------------------------------------------
    # Pairwise identity
    # -----------------------------------------------------

    identity_rows = []

    print("\nPAIRWISE ALIGNMENT-BASED IDENTITY")

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
                "sequence_1": name1,
                "sequence_2": name2,
                "matches": matches,
                "positions_compared": compared,
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
        writer.writerows(identity_rows)

    # -----------------------------------------------------
    # Within-species NKG2A vs NKG2C differences
    # -----------------------------------------------------

    diff_rows = []

    print(
        "\nWITHIN-SPECIES NKG2A vs NKG2C"
    )

    for species, nkg2a_name, nkg2c_name in (
        WITHIN_SPECIES_PAIRS
    ):

        seq_a = records[nkg2a_name]
        seq_c = records[nkg2c_name]

        pos_a = position_maps[nkg2a_name]
        pos_c = position_maps[nkg2c_name]

        matches, compared, identity = (
            aligned_identity(
                seq_a,
                seq_c,
            )
        )

        print(
            f"{species:<15} "
            f"{identity:>7.3f}% identity "
            f"({matches}/{compared})"
        )

        for i, (aa_a, aa_c) in enumerate(
            zip(seq_a, seq_c)
        ):

            if aa_a == aa_c:
                continue

            diff_rows.append(
                {
                    "comparison": species,
                    "alignment_column": i + 1,
                    "NKG2A_residue": (
                        ""
                        if pos_a[i] is None
                        else pos_a[i]
                    ),
                    "NKG2A_aa": aa_a,
                    "NKG2C_residue": (
                        ""
                        if pos_c[i] is None
                        else pos_c[i]
                    ),
                    "NKG2C_aa": aa_c,
                    "difference_type": (
                        "gap"
                        if aa_a == "-"
                        or aa_c == "-"
                        else "substitution"
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
                "NKG2A_residue",
                "NKG2A_aa",
                "NKG2C_residue",
                "NKG2C_aa",
                "difference_type",
            ],
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(diff_rows)

    print("\nOutputs:")
    print(IDENTITY_OUT)
    print(DIFF_OUT)
    print(MAP_OUT)


if __name__ == "__main__":
    main()
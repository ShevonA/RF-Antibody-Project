from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parent.parent

ECTODOMAIN_FILE = (
    ROOT
    / "results"
    / "tables"
    / "primary_ectodomain_sequences.tsv"
)

CLUSTER_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_candidate_surface_clusters.tsv"
)

OUT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2_ectodomain_n_glycosylation_sites.tsv"
)


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


def write_tsv(path, rows, fields):
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
            fieldnames=fields,
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

    return int(value)


def find_n_linked_sequons(sequence):
    """
    Find canonical N-linked glycosylation sequons:

        N-X-S/T

    where X is any amino acid except proline.

    Returns zero-based sequence indexes.
    """

    hits = []

    for i in range(len(sequence) - 2):

        aa1 = sequence[i]
        aa2 = sequence[i + 1]
        aa3 = sequence[i + 2]

        if (
            aa1 == "N"
            and aa2 != "P"
            and aa3 in {"S", "T"}
        ):
            hits.append(i)

    return hits


def main():

    print("=" * 78)
    print("STEP 2I - NKG2 ECTODOMAIN N-LINKED GLYCOSYLATION ANALYSIS")
    print("=" * 78)

    ectodomain_rows = read_tsv(
        ECTODOMAIN_FILE
    )

    cluster_rows = read_tsv(
        CLUSTER_FILE
    )

    candidate_positions = set()

    for row in cluster_rows:

        residue = safe_int(
            row.get(
                "full_length_residue"
            )
        )

        if residue is not None:
            candidate_positions.add(
                residue
            )

    output_rows = []

    print()

    for row in ectodomain_rows:

        record_id = clean(
            row.get("record_id")
        )

        accession = clean(
            row.get("accession")
        )

        sequence = clean(
            row.get("sequence")
        ).upper()

        ectodomain_start = safe_int(
            row.get(
                "ectodomain_start"
            )
        )

        ectodomain_end = safe_int(
            row.get(
                "ectodomain_end"
            )
        )

        if (
            not sequence
            or ectodomain_start is None
        ):
            continue

        sequons = find_n_linked_sequons(
            sequence
        )

        print(
            f"{record_id:<25} "
            f"{len(sequons)} sequon(s)"
        )

        for index in sequons:

            ectodomain_residue = (
                index + 1
            )

            full_length_residue = (
                ectodomain_start
                + index
            )

            motif = sequence[
                index:index + 3
            ]

            # Flag whether the sequon itself overlaps
            # one of the structurally analyzed human
            # NKG2A candidate positions.
            motif_positions = {
                full_length_residue,
                full_length_residue + 1,
                full_length_residue + 2,
            }

            candidate_overlap = (
                "yes"
                if (
                    record_id
                    == "human_NKG2A"
                    and motif_positions
                    & candidate_positions
                )
                else "no"
            )

            output_rows.append(
                {
                    "record_id":
                        record_id,

                    "accession":
                        accession,

                    "ectodomain_start":
                        ectodomain_start,

                    "ectodomain_end":
                        ectodomain_end,

                    "sequon_ecto_residue":
                        ectodomain_residue,

                    "sequon_full_length_residue":
                        full_length_residue,

                    "motif":
                        motif,

                    "candidate_position_overlap":
                        candidate_overlap,
                }
            )

            print(
                f"    {full_length_residue:>3} "
                f"{motif}"
            )

    fields = [
        "record_id",
        "accession",
        "ectodomain_start",
        "ectodomain_end",
        "sequon_ecto_residue",
        "sequon_full_length_residue",
        "motif",
        "candidate_position_overlap",
    ]

    write_tsv(
        OUT_FILE,
        output_rows,
        fields,
    )

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)

    print(
        OUT_FILE
    )

    print()
    print(
        "NOTE: N-X-S/T identifies canonical sequence "
        "sequons only. It does not establish that a site "
        "is experimentally glycosylated."
    )


if __name__ == "__main__":
    main()
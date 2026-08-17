from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
BOUNDARIES = ROOT / "data" / "curated" / "ectodomain_boundaries.tsv"

OUT_FASTA = (
    ROOT
    / "alignments"
    / "primary_ectodomain_input.fasta"
)

OUT_TABLE = (
    ROOT
    / "results"
    / "tables"
    / "primary_ectodomain_sequences.tsv"
)


FILES = {
    "human_NKG2A":
        RAW / "human_NKG2A_uniprot_P26715.fasta",

    "human_NKG2C":
        RAW / "human_NKG2C_refseq_NP_002251_2.fasta",

    "rhesus_NKG2A":
        RAW / "rhesus_NKG2A_refseq_NP_001028001_3.fasta",

    "rhesus_NKG2C_isoform1":
        RAW / "rhesus_NKG2C_isoform1_NP_001305287_3.fasta",

    "rhesus_NKG2C_isoform2":
        RAW / "rhesus_NKG2C_isoform2_NP_001098647_3.fasta",

    "pigtail_NKG2A":
        RAW / "pigtail_NKG2A_X1_XP_070928357_1.fasta",

    "pigtail_NKG2C":
        RAW / "pigtail_NKG2C_XP_070928345_1.fasta",
}


def read_fasta(path):
    return "".join(
        line.strip()
        for line in path.read_text(
            encoding="utf-8"
        ).splitlines()
        if line.strip()
        and not line.startswith(">")
    )


def main():

    with BOUNDARIES.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        boundaries = {
            row["record_id"]: row
            for row in csv.DictReader(
                handle,
                delimiter="\t",
            )
        }

    fasta_records = []
    table_rows = []

    print("=" * 78)
    print("STEP 1D - EXTRACT PRIMARY ECTODOMAINS")
    print("=" * 78)

    for record_id, path in FILES.items():

        if not path.exists():
            raise FileNotFoundError(path)

        if record_id not in boundaries:
            raise KeyError(
                f"No boundary for {record_id}"
            )

        sequence = read_fasta(path)
        info = boundaries[record_id]

        start = int(
            info["ectodomain_start"]
        )

        end = int(
            info["ectodomain_end"]
        )

        ectodomain = sequence[
            start - 1:end
        ]

        expected_length = (
            end - start + 1
        )

        if len(ectodomain) != expected_length:
            raise ValueError(
                f"{record_id}: expected "
                f"{expected_length} aa, got "
                f"{len(ectodomain)}"
            )

        label = (
            f"{record_id}|{info['accession']}"
            f"|ecto_{start}-{end}"
        )

        fasta_records.append(
            f">{label}\n{ectodomain}\n"
        )

        table_rows.append({
            "record_id": record_id,
            "accession": info["accession"],
            "full_length": len(sequence),
            "ectodomain_start": start,
            "ectodomain_end": end,
            "ectodomain_length":
                len(ectodomain),
            "boundary_status":
                info["boundary_status"],
            "sequence": ectodomain,
        })

        print(
            f"{record_id:<25} "
            f"{start:>3}-{end:<3} "
            f"{len(ectodomain):>3} aa  "
            f"{info['boundary_status']}"
        )

    OUT_FASTA.write_text(
        "".join(fasta_records),
        encoding="utf-8",
    )

    with OUT_TABLE.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        fields = [
            "record_id",
            "accession",
            "full_length",
            "ectodomain_start",
            "ectodomain_end",
            "ectodomain_length",
            "boundary_status",
            "sequence",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(table_rows)

    print("\nFASTA:")
    print(OUT_FASTA)

    print("\nTable:")
    print(OUT_TABLE)


if __name__ == "__main__":
    main()
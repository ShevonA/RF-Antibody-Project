from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "alignments" / "primary_full_length_input.fasta"


SEQUENCES = [
    (
        "human_NKG2A_P26715",
        RAW / "human_NKG2A_uniprot_P26715.fasta",
    ),
    (
        "human_NKG2C_NP_002251.2",
        RAW / "human_NKG2C_refseq_NP_002251_2.fasta",
    ),
    (
        "rhesus_NKG2A_NP_001028001.3",
        RAW / "rhesus_NKG2A_refseq_NP_001028001_3.fasta",
    ),
    (
        "rhesus_NKG2C_iso1_NP_001305287.3",
        RAW / "rhesus_NKG2C_isoform1_NP_001305287_3.fasta",
    ),
    (
        "rhesus_NKG2C_iso2_NP_001098647.3",
        RAW / "rhesus_NKG2C_isoform2_NP_001098647_3.fasta",
    ),
    (
        "pigtail_NKG2A_XP_070928357.1",
        RAW / "pigtail_NKG2A_X1_XP_070928357_1.fasta",
    ),
    (
        "pigtail_NKG2C_XP_070928345.1",
        RAW / "pigtail_NKG2C_XP_070928345_1.fasta",
    ),
]


def read_fasta(path):
    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    return "".join(
        line.strip()
        for line in lines
        if line.strip() and not line.startswith(">")
    )


def main():
    OUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    records = []

    print("=" * 72)
    print("BUILD PRIMARY FULL-LENGTH ALIGNMENT INPUT")
    print("=" * 72)

    for label, path in SEQUENCES:

        if not path.exists():
            raise FileNotFoundError(
                f"Missing required FASTA: {path}"
            )

        sequence = read_fasta(path)

        print(
            f"{label:<42} {len(sequence):>3} aa"
        )

        records.append(
            f">{label}\n{sequence}\n"
        )

    OUT.write_text(
        "".join(records),
        encoding="utf-8",
    )

    print(f"\nWrote {len(records)} sequences to:")
    print(OUT)


if __name__ == "__main__":
    main()
from pathlib import Path
from itertools import combinations
import csv


ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"

OUT_SUMMARY = ROOT / "results" / "tables" / "sequence_comparison_summary.tsv"
OUT_DIFFS = ROOT / "results" / "tables" / "sequence_pairwise_differences.tsv"


SEQUENCES = {
    # Human
    "human_NKG2A_UniProt_P26715":
        RAW / "human_NKG2A_uniprot_P26715.fasta",

    "human_NKG2A_RefSeq_NP_002250_2":
        RAW / "human_NKG2A_refseq_NP_002250_2.fasta",

    "human_NKG2A_isoformC_JC":
        RAW / "human_NKG2A_isoform_C_JC.fasta",

    "human_NKG2C_NP_002251_2":
        RAW / "human_NKG2C_refseq_NP_002251_2.fasta",

    # Rhesus
    "rhesus_NKG2A_old_NP_001028001_2":
        RAW / "rhesus_NKG2A_refseq_NP_001028001_2.fasta",

    "rhesus_NKG2C_old_NP_001098647_2":
        RAW / "rhesus_NKG2C_refseq_NP_001098647_2.fasta",

    "rhesus_NKG2C_current_iso1_NP_001305287_3":
        RAW / "rhesus_NKG2C_isoform1_NP_001305287_3.fasta",

    "rhesus_NKG2C_current_iso2_NP_001098647_3":
        RAW / "rhesus_NKG2C_isoform2_NP_001098647_3.fasta",

    "rhesus_NKG2E_legacy_NP_001038194_2":
        RAW / "rhesus_NKG2C2_legacy_NP_001038194_2.fasta",

    # Pig-tailed macaque
    "pigtail_NKG2A_XP_070928357_1":
        RAW / "pigtail_NKG2A_X1_XP_070928357_1.fasta",

    "pigtail_NKG2A_XP_070928358_1":
        RAW / "pigtail_NKG2A_X1_XP_070928358_1.fasta",

    "pigtail_NKG2C_XP_070928345_1":
        RAW / "pigtail_NKG2C_XP_070928345_1.fasta",

        "rhesus_NKG2A_current_NP_001028001_3":
    RAW / "rhesus_NKG2A_refseq_NP_001028001_3.fasta",
}


def read_fasta(path):
    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    sequence = "".join(
        line.strip()
        for line in lines
        if line.strip()
        and not line.startswith(">")
    ).upper()

    return sequence


def compare_equal_length(seq1, seq2):
    if len(seq1) != len(seq2):
        return None

    differences = []

    matches = 0

    for position, (aa1, aa2) in enumerate(
        zip(seq1, seq2),
        start=1,
    ):
        if aa1 == aa2:
            matches += 1
        else:
            differences.append(
                (position, aa1, aa2)
            )

    identity = (
        matches / len(seq1) * 100
        if seq1
        else 0
    )

    return identity, differences


def main():
    OUT_SUMMARY.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    loaded = {}

    print("=" * 80)
    print("STEP 1C - SEQUENCE QC AND COMPARISON")
    print("=" * 80)

    print("\nLoading sequences...")

    for name, path in SEQUENCES.items():
        if not path.exists():
            print(
                f"MISSING  {name}: {path.name}"
            )
            continue

        seq = read_fasta(path)
        loaded[name] = seq

        print(
            f"OK  {name:<48} "
            f"{len(seq):>3} aa"
        )

    summary_rows = []
    difference_rows = []

    print("\n" + "=" * 80)
    print("PAIRWISE COMPARISON")
    print("=" * 80)

    for name1, name2 in combinations(
        loaded.keys(),
        2,
    ):
        seq1 = loaded[name1]
        seq2 = loaded[name2]

        result = compare_equal_length(
            seq1,
            seq2,
        )

        if result is None:
            summary_rows.append({
                "sequence_1": name1,
                "length_1": len(seq1),
                "sequence_2": name2,
                "length_2": len(seq2),
                "same_length": "no",
                "identical": "",
                "identity_percent": "",
                "n_differences": "",
            })

            continue

        identity, differences = result

        identical = (
            "yes"
            if len(differences) == 0
            else "no"
        )

        summary_rows.append({
            "sequence_1": name1,
            "length_1": len(seq1),
            "sequence_2": name2,
            "length_2": len(seq2),
            "same_length": "yes",
            "identical": identical,
            "identity_percent": f"{identity:.3f}",
            "n_differences": len(differences),
        })

        for position, aa1, aa2 in differences:
            difference_rows.append({
                "sequence_1": name1,
                "sequence_2": name2,
                "position": position,
                "aa_sequence_1": aa1,
                "aa_sequence_2": aa2,
            })

    with OUT_SUMMARY.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sequence_1",
                "length_1",
                "sequence_2",
                "length_2",
                "same_length",
                "identical",
                "identity_percent",
                "n_differences",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    with OUT_DIFFS.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sequence_1",
                "sequence_2",
                "position",
                "aa_sequence_1",
                "aa_sequence_2",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(difference_rows)

    print(
        f"\nSummary written to:\n{OUT_SUMMARY}"
    )

    print(
        f"\nDifferences written to:\n{OUT_DIFFS}"
    )

    # Key comparisons we specifically care about
    KEY_PAIRS = [
        (
            "pigtail_NKG2A_XP_070928357_1",
            "pigtail_NKG2A_XP_070928358_1",
        ),
        (
            "rhesus_NKG2C_old_NP_001098647_2",
            "rhesus_NKG2C_current_iso2_NP_001098647_3",
        ),
        (
            "human_NKG2A_UniProt_P26715",
            "human_NKG2A_RefSeq_NP_002250_2",
        ),
        (
            "human_NKG2A_UniProt_P26715",
            "rhesus_NKG2A_old_NP_001028001_2",
        ),
        (
            "human_NKG2A_UniProt_P26715",
            "pigtail_NKG2A_XP_070928357_1",
        ),
        (
            "rhesus_NKG2A_old_NP_001028001_2",
            "pigtail_NKG2A_XP_070928357_1",
        ),

        (
    "rhesus_NKG2A_old_NP_001028001_2",
    "rhesus_NKG2A_current_NP_001028001_3",
),
(
    "human_NKG2A_UniProt_P26715",
    "rhesus_NKG2A_current_NP_001028001_3",
),
(
    "rhesus_NKG2A_current_NP_001028001_3",
    "pigtail_NKG2A_XP_070928357_1",
),
    ]

    print("\n" + "=" * 80)
    print("KEY COMPARISONS")
    print("=" * 80)

    for name1, name2 in KEY_PAIRS:
        if (
            name1 not in loaded
            or name2 not in loaded
        ):
            continue

        seq1 = loaded[name1]
        seq2 = loaded[name2]

        result = compare_equal_length(
            seq1,
            seq2,
        )

        print(f"\n{name1}")
        print(f"vs")
        print(f"{name2}")

        if result is None:
            print(
                f"Different lengths: "
                f"{len(seq1)} vs {len(seq2)}"
            )
            continue

        identity, differences = result

        print(
            f"Identity: {identity:.3f}%"
        )
        print(
            f"Differences: {len(differences)}"
        )

        if len(differences) <= 20:
            for position, aa1, aa2 in differences:
                print(
                    f"  {position:>3}: "
                    f"{aa1} -> {aa2}"
                )


if __name__ == "__main__":
    main()
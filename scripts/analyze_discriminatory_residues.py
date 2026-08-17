from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parent.parent

ALIGNMENT = (
    ROOT
    / "alignments"
    / "primary_ectodomain_aligned.fasta"
)

OUT_TABLE = (
    ROOT
    / "results"
    / "tables"
    / "discriminatory_residue_analysis.tsv"
)

OUT_CANDIDATES = (
    ROOT
    / "results"
    / "tables"
    / "candidate_specificity_positions.tsv"
)


SEQUENCE_NAMES = {
    "human_A":
        "human_NKG2A|P26715|ecto_94-233",

    "human_C":
        "human_NKG2C|NP_002251.2|ecto_94-231",

    "rhesus_A":
        "rhesus_NKG2A|NP_001028001.3|ecto_94-233",

    "rhesus_C1":
        "rhesus_NKG2C_isoform1|NP_001305287.3|ecto_94-231",

    "rhesus_C2":
        "rhesus_NKG2C_isoform2|NP_001098647.3|ecto_94-231",

    "pigtail_A":
        "pigtail_NKG2A|XP_070928357.1|ecto_94-233",

    "pigtail_C":
        "pigtail_NKG2C|XP_070928345.1|ecto_94-231",
}


def read_alignment(path):
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


def residue_number_map(
    aligned_sequence,
    full_length_start=94,
):
    """
    Return the original full-length residue number
    corresponding to each alignment column.

    Gaps receive None.
    """

    numbers = []

    residue_number = full_length_start - 1

    for aa in aligned_sequence:

        if aa == "-":
            numbers.append(None)

        else:
            residue_number += 1
            numbers.append(residue_number)

    return numbers


def differs(aa_a, aa_c):
    """
    True if NKG2A and NKG2C differ at the aligned site.

    Gap/amino-acid differences count as differences.
    Double gaps do not.
    """

    if aa_a == "-" and aa_c == "-":
        return False

    return aa_a != aa_c


def classify_position(
    human_diff,
    rhesus_c1_diff,
    rhesus_c2_diff,
    pigtail_diff,
):
    """
    Assign each alignment position to a useful
    discrimination category.
    """

    rhesus_both = (
        rhesus_c1_diff
        and rhesus_c2_diff
    )

    rhesus_any = (
        rhesus_c1_diff
        or rhesus_c2_diff
    )

    if (
        human_diff
        and rhesus_both
        and pigtail_diff
    ):
        return "pan_species_NKG2A_vs_NKG2C_difference"

    if (
        not human_diff
        and rhesus_both
        and pigtail_diff
    ):
        return "macaque_shared_difference"

    if (
        not human_diff
        and not rhesus_any
        and pigtail_diff
    ):
        return "pigtail_specific_difference"

    if (
        human_diff
        and not rhesus_any
        and pigtail_diff
    ):
        return "human_pigtail_shared_difference"

    if (
        human_diff
        and rhesus_both
        and not pigtail_diff
    ):
        return "human_rhesus_shared_difference"

    if (
        not human_diff
        and rhesus_both
        and not pigtail_diff
    ):
        return "rhesus_specific_difference"

    if (
        rhesus_c1_diff
        != rhesus_c2_diff
    ):
        return "rhesus_isoform_dependent_difference"

    if (
        human_diff
        or rhesus_any
        or pigtail_diff
    ):
        return "other_difference_pattern"

    return "fully_conserved_A_C"


def same_nkg2a_across_macaques(
    rhesus_aa,
    pigtail_aa,
):
    return (
        rhesus_aa != "-"
        and pigtail_aa != "-"
        and rhesus_aa == pigtail_aa
    )


def candidate_priority(
    classification,
    rhesus_a,
    pigtail_a,
    human_a,
):
    """
    Simple sequence-only priority category.

    This does NOT imply surface accessibility or
    structural suitability.
    """

    if classification == "pan_species_NKG2A_vs_NKG2C_difference":
        if (
            rhesus_a == pigtail_a
            and rhesus_a != "-"
        ):
            return "high_sequence_interest"

        return "moderate_sequence_interest"

    if classification == "macaque_shared_difference":
        if (
            rhesus_a == pigtail_a
            and rhesus_a != "-"
        ):
            return "high_macaque_interest"

        return "moderate_macaque_interest"

    if classification == "pigtail_specific_difference":
        return "pigtail_specific_interest"

    if classification == "human_pigtail_shared_difference":
        return "pigtail_specific_interest"

    if classification == "rhesus_isoform_dependent_difference":
        return "caution_rhesus_isoform_variability"

    return ""


def main():

    print("=" * 80)
    print("STEP 1E - DISCRIMINATORY RESIDUE ANALYSIS")
    print("=" * 80)

    if not ALIGNMENT.exists():
        raise FileNotFoundError(
            f"Alignment not found:\n{ALIGNMENT}"
        )

    records = read_alignment(
        ALIGNMENT
    )

    print(
        f"\nSequences loaded: {len(records)}"
    )

    for short_name, full_name in (
        SEQUENCE_NAMES.items()
    ):

        if full_name not in records:
            raise KeyError(
                f"Missing expected sequence:\n"
                f"{full_name}"
            )

    selected = {
        short_name: records[full_name]
        for short_name, full_name
        in SEQUENCE_NAMES.items()
    }

    lengths = {
        len(seq)
        for seq in selected.values()
    }

    if len(lengths) != 1:
        raise ValueError(
            "Alignment lengths are inconsistent."
        )

    alignment_length = next(
        iter(lengths)
    )

    print(
        f"Alignment length: "
        f"{alignment_length} columns"
    )

    residue_maps = {
        name: residue_number_map(
            sequence,
            full_length_start=94,
        )
        for name, sequence
        in selected.items()
    }

    rows = []

    class_counts = {}

    for i in range(
        alignment_length
    ):

        human_a = selected[
            "human_A"
        ][i]

        human_c = selected[
            "human_C"
        ][i]

        rhesus_a = selected[
            "rhesus_A"
        ][i]

        rhesus_c1 = selected[
            "rhesus_C1"
        ][i]

        rhesus_c2 = selected[
            "rhesus_C2"
        ][i]

        pigtail_a = selected[
            "pigtail_A"
        ][i]

        pigtail_c = selected[
            "pigtail_C"
        ][i]

        human_diff = differs(
            human_a,
            human_c,
        )

        rhesus_c1_diff = differs(
            rhesus_a,
            rhesus_c1,
        )

        rhesus_c2_diff = differs(
            rhesus_a,
            rhesus_c2,
        )

        pigtail_diff = differs(
            pigtail_a,
            pigtail_c,
        )

        classification = classify_position(
            human_diff,
            rhesus_c1_diff,
            rhesus_c2_diff,
            pigtail_diff,
        )

        priority = candidate_priority(
            classification,
            rhesus_a,
            pigtail_a,
            human_a,
        )

        class_counts[
            classification
        ] = (
            class_counts.get(
                classification,
                0,
            )
            + 1
        )

        row = {
            "alignment_column":
                i + 1,

            "human_NKG2A_residue":
                residue_maps[
                    "human_A"
                ][i]
                or "",

            "human_NKG2A_aa":
                human_a,

            "human_NKG2C_residue":
                residue_maps[
                    "human_C"
                ][i]
                or "",

            "human_NKG2C_aa":
                human_c,

            "human_A_vs_C_diff":
                "yes"
                if human_diff
                else "no",

            "rhesus_NKG2A_residue":
                residue_maps[
                    "rhesus_A"
                ][i]
                or "",

            "rhesus_NKG2A_aa":
                rhesus_a,

            "rhesus_NKG2C1_residue":
                residue_maps[
                    "rhesus_C1"
                ][i]
                or "",

            "rhesus_NKG2C1_aa":
                rhesus_c1,

            "rhesus_A_vs_C1_diff":
                "yes"
                if rhesus_c1_diff
                else "no",

            "rhesus_NKG2C2_residue":
                residue_maps[
                    "rhesus_C2"
                ][i]
                or "",

            "rhesus_NKG2C2_aa":
                rhesus_c2,

            "rhesus_A_vs_C2_diff":
                "yes"
                if rhesus_c2_diff
                else "no",

            "pigtail_NKG2A_residue":
                residue_maps[
                    "pigtail_A"
                ][i]
                or "",

            "pigtail_NKG2A_aa":
                pigtail_a,

            "pigtail_NKG2C_residue":
                residue_maps[
                    "pigtail_C"
                ][i]
                or "",

            "pigtail_NKG2C_aa":
                pigtail_c,

            "pigtail_A_vs_C_diff":
                "yes"
                if pigtail_diff
                else "no",

            "macaque_NKG2A_same":
                (
                    "yes"
                    if same_nkg2a_across_macaques(
                        rhesus_a,
                        pigtail_a,
                    )
                    else "no"
                ),

            "classification":
                classification,

            "sequence_priority":
                priority,
        }

        rows.append(
            row
        )

    fields = list(
        rows[0].keys()
    )

    OUT_TABLE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with OUT_TABLE.open(
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

    candidate_rows = [
        row
        for row in rows
        if row[
            "sequence_priority"
        ]
    ]

    with OUT_CANDIDATES.open(
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
            candidate_rows
        )

    print("\nClassification counts:")

    for category, count in sorted(
        class_counts.items(),
        key=lambda item: (
            -item[1],
            item[0],
        ),
    ):

        print(
            f"  {category:<45} "
            f"{count:>3}"
        )

    print(
        f"\nCandidate sequence positions: "
        f"{len(candidate_rows)}"
    )

    print("\nOutputs:")

    print(
        OUT_TABLE
    )

    print(
        OUT_CANDIDATES
    )

    print(
        "\nNOTE: sequence_priority reflects "
        "sequence discrimination only."
    )

    print(
        "Surface exposure, CD94 contacts, "
        "HLA-E contacts, glycosylation, and "
        "known antibody epitopes still need "
        "structural annotation."
    )


if __name__ == "__main__":
    main()
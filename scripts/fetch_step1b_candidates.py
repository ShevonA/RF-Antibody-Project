from pathlib import Path
import csv
import time
import urllib.parse
import urllib.request


ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "data" / "raw"
LOGFILE = ROOT / "results" / "qc" / "step1b_download_log.tsv"

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

RECORDS = [
    # Rhesus CD94
    ("rhesus_CD94_NP_001028000_1", "NP_001028000.1"),

    # Current rhesus KLRC2 products
    ("rhesus_NKG2C_isoform1_NP_001305287_3", "NP_001305287.3"),
    ("rhesus_NKG2C_isoform2_NP_001098647_3", "NP_001098647.3"),

    # Pig-tailed macaque KLRC1 / NKG2A candidates
    ("pigtail_NKG2A_X1_XP_070928357_1", "XP_070928357.1"),
    ("pigtail_NKG2A_X1_XP_070928358_1", "XP_070928358.1"),

    # Pig-tailed macaque KLRC2 / NKG2C
    ("pigtail_NKG2C_XP_070928345_1", "XP_070928345.1"),

    # Pig-tailed macaque KLRD1 / CD94
    ("pigtail_CD94_X1_XP_011771121_1", "XP_011771121.1"),
    ("pigtail_CD94_X2_XP_011771122_1", "XP_011771122.1"),
]


def fetch_fasta(accession):
    params = {
        "db": "protein",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }

    url = BASE + "?" + urllib.parse.urlencode(params)

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent":
                "NKG2A-ectodomain-comparison/1.0"
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=30,
    ) as response:
        return response.read().decode("utf-8")


def sequence_length(fasta):
    sequence = "".join(
        line.strip()
        for line in fasta.splitlines()
        if line and not line.startswith(">")
    )

    return len(sequence)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    LOGFILE.parent.mkdir(parents=True, exist_ok=True)

    rows = []

    print("=" * 72)
    print("STEP 1B - FETCH RESOLVED CANDIDATES")
    print("=" * 72)

    for i, (record_id, accession) in enumerate(RECORDS, start=1):

        outfile = OUTDIR / f"{record_id}.fasta"

        print(f"\n[{i}/{len(RECORDS)}] {record_id}")
        print(f"Accession: {accession}")

        try:
            fasta = fetch_fasta(accession)
            length = sequence_length(fasta)

            outfile.write_text(
                fasta,
                encoding="utf-8",
            )

            status = "downloaded"

            print(f"OK - {length} aa")

        except Exception as error:
            length = ""
            status = f"ERROR: {error}"

            print(status)

        rows.append({
            "record_id": record_id,
            "accession": accession,
            "sequence_length": length,
            "status": status,
            "output_file": str(outfile),
        })

        time.sleep(0.4)

    with LOGFILE.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "record_id",
                "accession",
                "sequence_length",
                "status",
                "output_file",
            ],
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)

    print("\n" + "=" * 72)
    print("COMPLETE")
    print("=" * 72)

    print(f"\nLog:\n{LOGFILE}")


if __name__ == "__main__":
    main()
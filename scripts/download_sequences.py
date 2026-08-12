from pathlib import Path
import csv
import sys


PROJECT_ROOT = Path(__file__).resolve().parent.parent

MANIFEST = PROJECT_ROOT / "sequence_manifest.tsv"
RAW_DIR = PROJECT_ROOT / "data" / "raw"
CURATED_DIR = PROJECT_ROOT / "data" / "curated"
REFERENCE_DIR = PROJECT_ROOT / "data" / "reference"


def read_manifest(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames

    if not rows:
        raise ValueError("sequence_manifest.tsv is empty.")

    return rows, fieldnames


def clean(value):
    if value is None:
        return ""
    return value.strip()


def get_best_accession(row):
    """
    Return the best available accession for retrieval.

    Priority:
    1. fetch_id
    2. refseq_protein
    3. uniprot_id
    """

    fetch_id = clean(row.get("fetch_id"))
    refseq = clean(row.get("refseq_protein"))
    uniprot = clean(row.get("uniprot_id"))

    if fetch_id:
        return fetch_id, clean(row.get("source_database")) or "unspecified"

    if refseq:
        return refseq, "NCBI RefSeq"

    if uniprot:
        return uniprot, "UniProt"

    return "", ""


def main():
    print("=" * 72)
    print("NKG2A/NKG2C Sequence Download Pipeline")
    print("=" * 72)

    print(f"\nPython version:\n{sys.version}")
    print(f"\nProject root:\n{PROJECT_ROOT}")

    required_paths = {
        "sequence_manifest.tsv": MANIFEST,
        "data/raw": RAW_DIR,
        "data/curated": CURATED_DIR,
        "data/reference": REFERENCE_DIR,
    }

    print("\nChecking required project files and folders...")

    for label, path in required_paths.items():
        if path.exists():
            print(f"OK       {label:25} -> {path}")
        else:
            print(f"MISSING  {label:25} -> {path}")
            sys.exit(1)

    rows, fieldnames = read_manifest(MANIFEST)

    print("\nManifest loaded successfully.")
    print(f"Rows: {len(rows)}")

    print("\nColumns found:")
    for field in fieldnames:
        print(f"  - {field}")

    resolved = []
    unresolved = []

    for row_number, row in enumerate(rows, start=1):

        record_id = clean(row.get("record_id"))
        species_common = clean(row.get("species_common"))
        species_scientific = clean(row.get("species_scientific"))
        receptor = clean(row.get("receptor"))
        gene_symbol = clean(row.get("gene_symbol"))
        locus = clean(row.get("locus_or_paralog"))
        curation_status = clean(row.get("curation_status"))
        fetch_now = clean(row.get("fetch_now"))

        accession, database = get_best_accession(row)

        record = {
            "row": row_number,
            "record_id": record_id,
            "species_common": species_common,
            "species_scientific": species_scientific,
            "receptor": receptor,
            "gene_symbol": gene_symbol,
            "locus": locus,
            "accession": accession,
            "database": database,
            "curation_status": curation_status,
            "fetch_now": fetch_now,
        }

        if accession:
            resolved.append(record)
        else:
            unresolved.append(record)

    print("\n" + "=" * 72)
    print("ACCESSION SUMMARY")
    print("=" * 72)

    print(f"\nResolved records:   {len(resolved)}")
    print(f"Unresolved records: {len(unresolved)}")

    if resolved:
        print("\nResolved records:")
        for record in resolved:
            print(
                f"  Row {record['row']:>2} | "
                f"{record['record_id']:<18} | "
                f"{record['species_common']:<20} | "
                f"{record['receptor']:<6} | "
                f"{record['accession']:<15} | "
                f"{record['database']}"
            )

    if unresolved:
        print("\nUnresolved records:")
        for record in unresolved:
            print(
                f"  Row {record['row']:>2} | "
                f"{record['record_id']:<18} | "
                f"{record['species_common']:<20} | "
                f"{record['receptor']:<6} | "
                f"{record['locus'] or '[no locus/paralog]'}"
            )

    print("\n" + "=" * 72)
    print("RECORDS MARKED FOR FETCHING")
    print("=" * 72)

    fetch_records = []

    for record in resolved:
        flag = record["fetch_now"].lower()

        if flag in {"yes", "true", "1", "y"}:
            fetch_records.append(record)

    if fetch_records:
        for record in fetch_records:
            print(
                f"  {record['record_id']:<18} "
                f"{record['species_common']:<20} "
                f"{record['receptor']:<6} "
                f"{record['accession']}"
            )
    else:
        print("  No resolved records are currently marked fetch_now=yes.")

    print("\nNo network requests were made.")
    print("Next stage: download only validated records marked for retrieval.")


if __name__ == "__main__":
    main()
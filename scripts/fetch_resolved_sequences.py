from pathlib import Path
import argparse
import csv
import sys
import time
import urllib.error
import urllib.request


PROJECT_ROOT = Path(__file__).resolve().parent.parent

MANIFEST = PROJECT_ROOT / "sequence_manifest.tsv"
RAW_DIR = PROJECT_ROOT / "data" / "raw"
LOG_FILE = PROJECT_ROOT / "results" / "qc" / "download_log.tsv"

USER_AGENT = "NKG2A-NKG2C-sequence-curation/1.0"


def clean(value):
    if value is None:
        return ""
    return value.strip()


def truthy(value):
    return clean(value).lower() in {
        "yes",
        "y",
        "true",
        "1",
    }


def read_manifest():
    with MANIFEST.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )
        return list(reader)


def get_accession(row):
    fetch_id = clean(row.get("fetch_id"))

    if fetch_id:
        return fetch_id

    refseq = clean(row.get("refseq_protein"))

    if refseq:
        return refseq

    uniprot = clean(row.get("uniprot_id"))

    if uniprot:
        return uniprot

    return ""


def infer_route(row, accession):
    source = clean(
        row.get("source_database")
    ).lower()

    # Explicit database assignments first.
    if "uniprot" in source:
        return "uniprot"

    if (
        "ncbi" in source
        or "refseq" in source
    ):
        return "ncbi"

    # A user-provided sequence with an NCBI accession
    # can still be independently verified against NCBI.
    if accession.startswith(
        (
            "NP_",
            "XP_",
            "YP_",
            "WP_",
        )
    ):
        return "ncbi"

    # UniProt accessions in this project do not use
    # RefSeq-style prefixes.
    if accession:
        return "uniprot"

    return ""


def build_url(route, accession):
    if route == "ncbi":
        return (
            "https://eutils.ncbi.nlm.nih.gov/"
            "entrez/eutils/efetch.fcgi"
            f"?db=protein&id={accession}"
            "&rettype=fasta&retmode=text"
        )

    if route == "uniprot":
        return (
            "https://rest.uniprot.org/"
            f"uniprotkb/{accession}.fasta"
        )

    raise ValueError(
        f"Unknown retrieval route: {route}"
    )


def download_text(url, retries=3):
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
        },
    )

    last_error = None

    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(
                request,
                timeout=30,
            ) as response:
                return response.read().decode(
                    "utf-8"
                )

        except (
            urllib.error.URLError,
            urllib.error.HTTPError,
            TimeoutError,
        ) as error:

            last_error = error

            if attempt < retries:
                wait_seconds = attempt * 2
                print(
                    f"    Retry {attempt}/{retries} "
                    f"after error: {error}"
                )
                time.sleep(wait_seconds)

    raise RuntimeError(
        f"Download failed after {retries} "
        f"attempts: {last_error}"
    )


def validate_fasta(text):
    lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip()
    ]

    if not lines:
        return False, "Empty response"

    if not lines[0].startswith(">"):
        return False, "Response is not FASTA"

    sequence = "".join(
        line
        for line in lines[1:]
        if not line.startswith(">")
    ).upper()

    if not sequence:
        return False, "FASTA contains no sequence"

    allowed = set(
        "ACDEFGHIKLMNPQRSTVWY"
        "BXZJUO"
        "*"
    )

    bad = sorted(
        set(sequence) - allowed
    )

    if bad:
        return (
            False,
            "Unexpected characters: "
            + ",".join(bad),
        )

    return True, str(len(sequence))


def write_log(records):
    LOG_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fields = [
        "record_id",
        "species_common",
        "receptor",
        "accession",
        "route",
        "status",
        "sequence_length",
        "output_file",
        "message",
    ]

    with LOG_FILE.open(
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
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download resolved protein sequences "
            "listed in sequence_manifest.tsv."
        )
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Show what would be downloaded "
            "without making network requests."
        ),
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Overwrite FASTA files that "
            "already exist."
        ),
    )

    args = parser.parse_args()

    print("=" * 72)
    print(
        "NKG2A/NKG2C Resolved Sequence Fetcher"
    )
    print("=" * 72)

    if not MANIFEST.exists():
        sys.exit(
            f"Manifest not found: {MANIFEST}"
        )

    RAW_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    rows = read_manifest()

    selected = []

    for row in rows:
        accession = get_accession(row)

        if (
            accession
            and truthy(row.get("fetch_now"))
        ):
            selected.append(
                (row, accession)
            )

    print(
        f"\nRecords selected for retrieval: "
        f"{len(selected)}"
    )

    log_rows = []

    for number, (
        row,
        accession,
    ) in enumerate(
        selected,
        start=1,
    ):

        record_id = clean(
            row.get("record_id")
        )

        species = clean(
            row.get("species_common")
        )

        receptor = clean(
            row.get("receptor")
        )

        route = infer_route(
            row,
            accession,
        )

        output_file = (
            RAW_DIR
            / f"{record_id}.fasta"
        )

        print(
            f"\n[{number}/{len(selected)}] "
            f"{record_id}"
        )

        print(
            f"    Species:   {species}"
        )
        print(
            f"    Receptor:  {receptor}"
        )
        print(
            f"    Accession: {accession}"
        )
        print(
            f"    Route:     {route}"
        )
        print(
            f"    Output:    {output_file}"
        )

        if args.dry_run:
            print(
                "    DRY RUN - not downloaded"
            )
            continue

        if (
            output_file.exists()
            and not args.force
        ):
            print(
                "    SKIP - file already exists"
            )

            log_rows.append(
                {
                    "record_id": record_id,
                    "species_common": species,
                    "receptor": receptor,
                    "accession": accession,
                    "route": route,
                    "status": "skipped_existing",
                    "sequence_length": "",
                    "output_file": str(
                        output_file
                    ),
                    "message": (
                        "Existing file preserved"
                    ),
                }
            )

            continue

        try:
            url = build_url(
                route,
                accession,
            )

            fasta = download_text(url)

            valid, detail = validate_fasta(
                fasta
            )

            if not valid:
                raise ValueError(detail)

            output_file.write_text(
                fasta.rstrip() + "\n",
                encoding="utf-8",
            )

            print(
                f"    OK - {detail} aa"
            )

            log_rows.append(
                {
                    "record_id": record_id,
                    "species_common": species,
                    "receptor": receptor,
                    "accession": accession,
                    "route": route,
                    "status": "downloaded",
                    "sequence_length": detail,
                    "output_file": str(
                        output_file
                    ),
                    "message": "",
                }
            )

        except Exception as error:
            print(
                f"    ERROR - {error}"
            )

            log_rows.append(
                {
                    "record_id": record_id,
                    "species_common": species,
                    "receptor": receptor,
                    "accession": accession,
                    "route": route,
                    "status": "error",
                    "sequence_length": "",
                    "output_file": str(
                        output_file
                    ),
                    "message": str(error),
                }
            )

        # Be polite to public databases.
        time.sleep(0.4)

    if not args.dry_run:
        write_log(log_rows)

        print("\n" + "=" * 72)
        print("FETCH COMPLETE")
        print("=" * 72)

        print(
            f"\nLog written to:\n{LOG_FILE}"
        )

        downloaded = sum(
            row["status"] == "downloaded"
            for row in log_rows
        )

        errors = sum(
            row["status"] == "error"
            for row in log_rows
        )

        skipped = sum(
            row["status"]
            == "skipped_existing"
            for row in log_rows
        )

        print(
            f"\nDownloaded: {downloaded}"
        )
        print(
            f"Skipped:    {skipped}"
        )
        print(
            f"Errors:     {errors}"
        )


if __name__ == "__main__":
    main()
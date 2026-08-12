from pathlib import Path
import csv
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parent.parent

QUERY_FILE = PROJECT_ROOT / "ncbi_search_queries.tsv"
OUTPUT_FILE = PROJECT_ROOT / "results" / "tables" / "ncbi_candidate_hits.tsv"

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

USER_AGENT = "NKG2A-NKG2C-sequence-curation/1.0"

MAX_HITS_PER_QUERY = 20


# These are deliberately broad.
# We want candidate discovery, not automatic gene assignment.
SEARCHES = [
    {
        "search_id": "rhesus_NKG2C",
        "species_common": "rhesus_macaque",
        "species_scientific": "Macaca mulatta",
        "target": "NKG2C",
        "query": (
            'Macaca mulatta[Organism] AND '
            '("NKG2-C"[Protein Name] OR '
            '"NKG2C"[All Fields] OR '
            '"KLRC2"[All Fields] OR '
            '"killer cell lectin like receptor C2"[All Fields])'
        ),
    },
    {
        "search_id": "rhesus_KLRD1",
        "species_common": "rhesus_macaque",
        "species_scientific": "Macaca mulatta",
        "target": "CD94",
        "query": (
            'Macaca mulatta[Organism] AND '
            '("CD94"[All Fields] OR '
            '"KLRD1"[All Fields] OR '
            '"killer cell lectin like receptor D1"[All Fields])'
        ),
    },
    {
        "search_id": "pigtail_NKG2A",
        "species_common": "pig_tailed_macaque",
        "species_scientific": "Macaca nemestrina",
        "target": "NKG2A",
        "query": (
            'Macaca nemestrina[Organism] AND '
            '("NKG2-A"[Protein Name] OR '
            '"NKG2A"[All Fields] OR '
            '"KLRC1"[All Fields] OR '
            '"killer cell lectin like receptor C1"[All Fields])'
        ),
    },
    {
        "search_id": "pigtail_NKG2C",
        "species_common": "pig_tailed_macaque",
        "species_scientific": "Macaca nemestrina",
        "target": "NKG2C",
        "query": (
            'Macaca nemestrina[Organism] AND '
            '("NKG2-C"[Protein Name] OR '
            '"NKG2C"[All Fields] OR '
            '"KLRC2"[All Fields] OR '
            '"killer cell lectin like receptor C2"[All Fields])'
        ),
    },
    {
        "search_id": "pigtail_KLRD1",
        "species_common": "pig_tailed_macaque",
        "species_scientific": "Macaca nemestrina",
        "target": "CD94",
        "query": (
            'Macaca nemestrina[Organism] AND '
            '("CD94"[All Fields] OR '
            '"KLRD1"[All Fields] OR '
            '"killer cell lectin like receptor D1"[All Fields])'
        ),
    },
]


def request_text(url):
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT},
    )

    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8")


def ncbi_esearch(query):
    params = {
        "db": "protein",
        "term": query,
        "retmode": "json",
        "retmax": str(MAX_HITS_PER_QUERY),
    }

    url = (
        f"{NCBI_EUTILS}/esearch.fcgi?"
        + urllib.parse.urlencode(params)
    )

    text = request_text(url)
    data = json.loads(text)

    return data["esearchresult"].get("idlist", [])


def ncbi_esummary(ids):
    if not ids:
        return {}

    params = {
        "db": "protein",
        "id": ",".join(ids),
        "retmode": "json",
    }

    url = (
        f"{NCBI_EUTILS}/esummary.fcgi?"
        + urllib.parse.urlencode(params)
    )

    text = request_text(url)

    return json.loads(text)


def get_fasta_sequence_length(uid):
    params = {
        "db": "protein",
        "id": uid,
        "rettype": "fasta",
        "retmode": "text",
    }

    url = (
        f"{NCBI_EUTILS}/efetch.fcgi?"
        + urllib.parse.urlencode(params)
    )

    fasta = request_text(url)

    sequence = "".join(
        line.strip()
        for line in fasta.splitlines()
        if line and not line.startswith(">")
    )

    return len(sequence)


def extract_accession(record):
    accession = record.get("accessionversion", "")

    if accession:
        return accession

    caption = record.get("caption", "")
    return caption


def main():
    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_rows = []

    print("=" * 78)
    print("STEP 1B - NCBI CANDIDATE DISCOVERY")
    print("=" * 78)

    for search_number, search in enumerate(
        SEARCHES,
        start=1,
    ):
        print(
            f"\n[{search_number}/{len(SEARCHES)}] "
            f"{search['search_id']}"
        )

        print(f"Target: {search['target']}")
        print(f"Species: {search['species_scientific']}")
        print(f"Query: {search['query']}")

        try:
            ids = ncbi_esearch(search["query"])

        except Exception as error:
            print(f"SEARCH ERROR: {error}")
            continue

        print(f"Hits: {len(ids)}")

        if not ids:
            continue

        time.sleep(0.4)

        try:
            summary_data = ncbi_esummary(ids)

        except Exception as error:
            print(f"SUMMARY ERROR: {error}")
            continue

        result = summary_data.get("result", {})

        for rank, uid in enumerate(ids, start=1):
            record = result.get(uid, {})

            accession = extract_accession(record)
            title = record.get("title", "")
            taxid = record.get("taxid", "")
            update_date = record.get("updatedate", "")

            length = ""

            try:
                length = get_fasta_sequence_length(uid)

            except Exception:
                length = ""

            row = {
                "search_id": search["search_id"],
                "rank": rank,
                "species_common": search["species_common"],
                "species_scientific": search["species_scientific"],
                "target": search["target"],
                "protein_uid": uid,
                "accession": accession,
                "protein_title": title,
                "taxid": taxid,
                "sequence_length": length,
                "update_date": update_date,
                "query": search["query"],
                "curation_decision": "",
                "curation_notes": "",
            }

            output_rows.append(row)

            print(
                f"  {rank:>2}. "
                f"{accession:<18} "
                f"{str(length):>4} aa  "
                f"{title}"
            )

            time.sleep(0.35)

    fields = [
        "search_id",
        "rank",
        "species_common",
        "species_scientific",
        "target",
        "protein_uid",
        "accession",
        "protein_title",
        "taxid",
        "sequence_length",
        "update_date",
        "query",
        "curation_decision",
        "curation_notes",
    ]

    with OUTPUT_FILE.open(
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
        writer.writerows(output_rows)

    print("\n" + "=" * 78)
    print("SEARCH COMPLETE")
    print("=" * 78)

    print(f"\nCandidate records saved to:\n{OUTPUT_FILE}")

    print(f"\nTotal candidate hits saved: {len(output_rows)}")

    print(
        "\nIMPORTANT: Nothing has been added to "
        "sequence_manifest.tsv automatically."
    )

    print(
        "Candidate identity must be reviewed before "
        "assigning NKG2C paralog names."
    )


if __name__ == "__main__":
    main()
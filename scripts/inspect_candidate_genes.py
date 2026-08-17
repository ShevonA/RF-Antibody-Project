from pathlib import Path
import json
import time
import urllib.parse
import urllib.request


PROJECT_ROOT = Path(__file__).resolve().parent.parent

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "NKG2A-NKG2C-sequence-curation/1.0"


GENE_IDS = {
    "rhesus_KLRC2_current": "114670780",
    "rhesus_KLRC2_old": "709885",
    "pigtail_KLRC1": "105499932",
}


def request_json(url):
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT},
    )

    with urllib.request.urlopen(request, timeout=30) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_gene_summary(gene_id):
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json",
    }

    url = (
        f"{NCBI_EUTILS}/esummary.fcgi?"
        + urllib.parse.urlencode(params)
    )

    return request_json(url)


def main():
    print("=" * 78)
    print("STEP 1B - NCBI GENE RECORD INSPECTION")
    print("=" * 78)

    for label, gene_id in GENE_IDS.items():

        print(f"\n{label}")
        print("-" * 78)
        print(f"Gene ID: {gene_id}")

        try:
            data = fetch_gene_summary(gene_id)

            result = data.get("result", {})
            record = result.get(gene_id, {})

            print(f"Name:        {record.get('name', '')}")
            print(f"Description: {record.get('description', '')}")
            print(f"Status:      {record.get('status', '')}")
            print(f"Current ID:  {record.get('currentid', '')}")

            chromosome = record.get("chromosome", "")
            print(f"Chromosome:  {chromosome}")

            map_location = record.get("maplocation", "")
            print(f"Map loc:     {map_location}")

            organism = record.get("organism", {})
            print(
                "Organism:    "
                f"{organism.get('scientificname', '')}"
            )

            aliases = record.get("otheraliases", "")
            print(f"Aliases:     {aliases}")

            other_designations = record.get(
                "otherdesignations",
                "",
            )
            print(
                f"Other names: {other_designations}"
            )

        except Exception as error:
            print(f"ERROR: {error}")

        time.sleep(0.4)


if __name__ == "__main__":
    main()
from pathlib import Path
import time
import urllib.parse
import urllib.request


ROOT = Path(__file__).resolve().parent.parent

OUTDIR = (
    ROOT
    / "data"
    / "reference"
    / "topology_annotations"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


RECORDS = [
    # Human reviewed UniProt
    {
        "name": "human_NKG2A_P26715",
        "database": "uniprot",
        "accession": "P26715",
    },

    # NCBI protein records
    {
        "name": "human_NKG2C_NP_002251_2",
        "database": "ncbi",
        "accession": "NP_002251.2",
    },
    {
        "name": "rhesus_NKG2A_NP_001028001_3",
        "database": "ncbi",
        "accession": "NP_001028001.3",
    },
    {
        "name": "rhesus_NKG2C_iso1_NP_001305287_3",
        "database": "ncbi",
        "accession": "NP_001305287.3",
    },
    {
        "name": "rhesus_NKG2C_iso2_NP_001098647_3",
        "database": "ncbi",
        "accession": "NP_001098647.3",
    },
    {
        "name": "pigtail_NKG2A_XP_070928357_1",
        "database": "ncbi",
        "accession": "XP_070928357.1",
    },
    {
        "name": "pigtail_NKG2C_XP_070928345_1",
        "database": "ncbi",
        "accession": "XP_070928345.1",
    },
]


def fetch(url):
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

        return response.read().decode(
            "utf-8"
        )


def fetch_uniprot(accession):
    return fetch(
        "https://rest.uniprot.org/"
        f"uniprotkb/{accession}.txt"
    )


def fetch_ncbi(accession):

    params = {
        "db": "protein",
        "id": accession,
        "rettype": "gb",
        "retmode": "text",
    }

    url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/efetch.fcgi?"
        + urllib.parse.urlencode(params)
    )

    return fetch(url)


def print_relevant_lines(text):

    keywords = [
        "TRANSMEM",
        "TOPO_DOM",
        "transmembrane",
        "extracellular",
        "cytoplasmic",
        "Region",
        "/region_name",
    ]

    found = False

    for line in text.splitlines():

        if any(
            keyword.lower() in line.lower()
            for keyword in keywords
        ):
            print("   ", line.rstrip())
            found = True

    if not found:
        print(
            "    No obvious topology feature "
            "found in annotation."
        )


def main():

    print("=" * 76)
    print("STEP 1D - TOPOLOGY ANNOTATION RETRIEVAL")
    print("=" * 76)

    for index, record in enumerate(
        RECORDS,
        start=1,
    ):

        print(
            f"\n[{index}/{len(RECORDS)}] "
            f"{record['name']}"
        )

        print(
            f"Accession: {record['accession']}"
        )

        try:

            if record["database"] == "uniprot":

                text = fetch_uniprot(
                    record["accession"]
                )

                suffix = ".uniprot.txt"

            else:

                text = fetch_ncbi(
                    record["accession"]
                )

                suffix = ".genpept.txt"

            outfile = (
                OUTDIR
                / f"{record['name']}{suffix}"
            )

            outfile.write_text(
                text,
                encoding="utf-8",
            )

            print_relevant_lines(text)

            print(
                f"    Saved: {outfile.name}"
            )

        except Exception as error:

            print(
                f"    ERROR: {error}"
            )

        time.sleep(0.4)


if __name__ == "__main__":
    main()
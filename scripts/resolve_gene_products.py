from pathlib import Path
import csv
import json
import time
import urllib.parse
import urllib.request


PROJECT_ROOT = Path(__file__).resolve().parent.parent

OUTPUT_FILE = (
    PROJECT_ROOT
    / "results"
    / "tables"
    / "gene_product_resolution.tsv"
)

NCBI_EUTILS = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
)

USER_AGENT = "NKG2A-NKG2C-sequence-curation/1.0"


SEARCH_TARGETS = [
    {
        "record_id": "rhesus_KLRC2",
        "species": "Macaca mulatta",
        "gene_term": "KLRC2",
    },
    {
        "record_id": "rhesus_KLRD1",
        "species": "Macaca mulatta",
        "gene_term": "KLRD1",
    },
    {
        "record_id": "pigtail_KLRC1",
        "species": "Macaca nemestrina",
        "gene_term": "KLRC1",
    },
    {
        "record_id": "pigtail_KLRC2",
        "species": "Macaca nemestrina",
        "gene_term": "KLRC2",
    },
    {
        "record_id": "pigtail_KLRD1",
        "species": "Macaca nemestrina",
        "gene_term": "KLRD1",
    },
]


def request_json(endpoint, params):
    url = (
        f"{NCBI_EUTILS}/{endpoint}?"
        + urllib.parse.urlencode(params)
    )

    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT},
    )

    with urllib.request.urlopen(
        request,
        timeout=30,
    ) as response:
        return json.loads(
            response.read().decode("utf-8")
        )


def search_gene(species, gene_term):
    query = (
        f'{species}[Organism] AND '
        f'{gene_term}[All Fields]'
    )

    data = request_json(
        "esearch.fcgi",
        {
            "db": "gene",
            "term": query,
            "retmode": "json",
            "retmax": "20",
        },
    )

    return (
        data["esearchresult"]["idlist"],
        query,
    )


def gene_summary(gene_id):
    data = request_json(
        "esummary.fcgi",
        {
            "db": "gene",
            "id": gene_id,
            "retmode": "json",
        },
    )

    return data.get(
        "result",
        {},
    ).get(
        gene_id,
        {},
    )


def linked_protein_ids(gene_id):
    data = request_json(
        "elink.fcgi",
        {
            "dbfrom": "gene",
            "db": "protein",
            "id": gene_id,
            "retmode": "json",
        },
    )

    protein_ids = []

    for linkset in data.get(
        "linksets",
        [],
    ):
        for dbset in linkset.get(
            "linksetdbs",
            [],
        ):
            if dbset.get("dbto") == "protein":
                protein_ids.extend(
                    dbset.get("links", [])
                )

    return protein_ids


def protein_summaries(protein_ids):
    if not protein_ids:
        return {}

    data = request_json(
        "esummary.fcgi",
        {
            "db": "protein",
            "id": ",".join(protein_ids),
            "retmode": "json",
        },
    )

    return data.get(
        "result",
        {},
    )


def main():

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 80)
    print("STEP 1B - GENE TO PROTEIN RESOLUTION")
    print("=" * 80)

    rows = []

    for target in SEARCH_TARGETS:

        print(
            f"\n### {target['record_id']}"
        )

        print(
            f"Species: {target['species']}"
        )

        print(
            f"Gene term: {target['gene_term']}"
        )

        gene_ids, query = search_gene(
            target["species"],
            target["gene_term"],
        )

        print(
            f"Gene hits: {len(gene_ids)}"
        )

        for gene_id in gene_ids:

            time.sleep(0.35)

            summary = gene_summary(gene_id)

            gene_name = summary.get(
                "name",
                "",
            )

            description = summary.get(
                "description",
                "",
            )

            chromosome = summary.get(
                "chromosome",
                "",
            )

            aliases = summary.get(
                "otheraliases",
                "",
            )

            status = summary.get(
                "status",
                "",
            )

            current_id = summary.get(
                "currentid",
                "",
            )

            print(
                f"\n  Gene ID: {gene_id}"
            )
            print(
                f"  Name: {gene_name}"
            )
            print(
                f"  Description: {description}"
            )
            print(
                f"  Chromosome: {chromosome}"
            )

            protein_ids = linked_protein_ids(
                gene_id
            )

            print(
                f"  Linked proteins: "
                f"{len(protein_ids)}"
            )

            time.sleep(0.35)

            protein_data = protein_summaries(
                protein_ids
            )

            if not protein_ids:

                rows.append(
                    {
                        "search_record_id":
                            target["record_id"],
                        "species":
                            target["species"],
                        "gene_search_term":
                            target["gene_term"],
                        "gene_id":
                            gene_id,
                        "gene_name":
                            gene_name,
                        "gene_description":
                            description,
                        "chromosome":
                            chromosome,
                        "aliases":
                            aliases,
                        "gene_status":
                            status,
                        "current_gene_id":
                            current_id,
                        "protein_uid":
                            "",
                        "protein_accession":
                            "",
                        "protein_title":
                            "",
                        "protein_length":
                            "",
                        "query":
                            query,
                    }
                )

            for protein_uid in protein_ids:

                protein = protein_data.get(
                    str(protein_uid),
                    {},
                )

                accession = protein.get(
                    "accessionversion",
                    protein.get(
                        "caption",
                        "",
                    ),
                )

                title = protein.get(
                    "title",
                    "",
                )

                length = protein.get(
                    "slen",
                    "",
                )

                print(
                    f"    {accession:<18} "
                    f"{str(length):>4} aa  "
                    f"{title}"
                )

                rows.append(
                    {
                        "search_record_id":
                            target["record_id"],
                        "species":
                            target["species"],
                        "gene_search_term":
                            target["gene_term"],
                        "gene_id":
                            gene_id,
                        "gene_name":
                            gene_name,
                        "gene_description":
                            description,
                        "chromosome":
                            chromosome,
                        "aliases":
                            aliases,
                        "gene_status":
                            status,
                        "current_gene_id":
                            current_id,
                        "protein_uid":
                            protein_uid,
                        "protein_accession":
                            accession,
                        "protein_title":
                            title,
                        "protein_length":
                            length,
                        "query":
                            query,
                    }
                )

        time.sleep(0.5)

    fields = [
        "search_record_id",
        "species",
        "gene_search_term",
        "gene_id",
        "gene_name",
        "gene_description",
        "chromosome",
        "aliases",
        "gene_status",
        "current_gene_id",
        "protein_uid",
        "protein_accession",
        "protein_title",
        "protein_length",
        "query",
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
        writer.writerows(rows)

    print("\n" + "=" * 80)
    print("COMPLETE")
    print("=" * 80)

    print(
        f"\nResults written to:\n"
        f"{OUTPUT_FILE}"
    )


if __name__ == "__main__":
    main()
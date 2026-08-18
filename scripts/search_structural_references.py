from pathlib import Path
import csv
import json
import urllib.request
import urllib.error


ROOT = Path(__file__).resolve().parent.parent

OUTPUT_DIR = (
    ROOT
    / "results"
    / "tables"
    / "structure"
)

OUTPUT_FILE = (
    OUTPUT_DIR
    / "rcsb_nkg2a_structural_references.tsv"
)


# ---------------------------------------------------------------------
# RCSB PDB Search API
# ---------------------------------------------------------------------

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

RCSB_ENTRY_URL = (
    "https://data.rcsb.org/rest/v1/core/entry/{}"
)


# ---------------------------------------------------------------------
# Search terms
#
# We deliberately use several searches rather than assuming that every
# PDB record uses exactly the same receptor terminology.
# ---------------------------------------------------------------------

SEARCH_TERMS = [
    ("NKG2A", "NKG2A"),
    ("KLRC1", "KLRC1"),
    ("CD94 NKG2A", "CD94 NKG2A"),
    ("HLA-E NKG2A", "HLA-E NKG2A"),
]


def rcsb_text_search(term):
    """
    Search RCSB PDB entry text for a term.

    Returns a list of PDB entry IDs.
    """

    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": term
            },
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True
        },
    }

    data = json.dumps(query).encode("utf-8")

    request = urllib.request.Request(
        RCSB_SEARCH_URL,
        data=data,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "NKG2A-structural-reference-search/1.0",
        },
        method="POST",
    )

    try:
        with urllib.request.urlopen(
            request,
            timeout=60,
        ) as response:

            result = json.loads(
                response.read().decode("utf-8")
            )

    except urllib.error.HTTPError as exc:
        print(
            f"HTTP error searching '{term}': "
            f"{exc.code} {exc.reason}"
        )
        return []

    except Exception as exc:
        print(
            f"Error searching '{term}': {exc}"
        )
        return []

    pdb_ids = []

    for item in result.get(
        "result_set",
        [],
    ):
        identifier = item.get("identifier")

        if identifier:
            pdb_ids.append(identifier)

    return pdb_ids


def fetch_entry_metadata(pdb_id):
    """
    Retrieve core metadata for one PDB entry.
    """

    url = RCSB_ENTRY_URL.format(pdb_id)

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "NKG2A-structural-reference-search/1.0"
        },
    )

    try:
        with urllib.request.urlopen(
            request,
            timeout=60,
        ) as response:

            return json.loads(
                response.read().decode("utf-8")
            )

    except Exception as exc:
        print(
            f"Could not retrieve metadata for "
            f"{pdb_id}: {exc}"
        )

        return {}


def get_title(metadata):
    struct = metadata.get(
        "struct",
        {},
    )

    return struct.get(
        "title",
        "",
    )


def get_method(metadata):
    exptl = metadata.get(
        "exptl",
        [],
    )

    methods = []

    for item in exptl:
        method = item.get("method")

        if method:
            methods.append(method)

    return "; ".join(methods)


def get_resolution(metadata):
    info = metadata.get(
        "rcsb_entry_info",
        {},
    )

    resolution = info.get(
        "resolution_combined",
        [],
    )

    if not resolution:
        return ""

    return "; ".join(
        str(value)
        for value in resolution
    )


def get_release_date(metadata):
    accession = metadata.get(
        "rcsb_accession_info",
        {},
    )

    return accession.get(
        "initial_release_date",
        "",
    )


def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 78)
    print("STEP 2A - SEARCH STRUCTURAL REFERENCES")
    print("=" * 78)

    hits = {}

    # -------------------------------------------------------------
    # Search RCSB
    # -------------------------------------------------------------

    for label, term in SEARCH_TERMS:

        print(
            f"\nSearching RCSB PDB for: {term}"
        )

        pdb_ids = rcsb_text_search(term)

        print(
            f"Hits: {len(pdb_ids)}"
        )

        for pdb_id in pdb_ids:

            if pdb_id not in hits:
                hits[pdb_id] = set()

            hits[pdb_id].add(label)

    print(
        f"\nUnique PDB entries found: "
        f"{len(hits)}"
    )

    # -------------------------------------------------------------
    # Retrieve metadata
    # -------------------------------------------------------------

    rows = []

    sorted_ids = sorted(hits)

    for index, pdb_id in enumerate(
        sorted_ids,
        start=1,
    ):

        print(
            f"[{index}/{len(sorted_ids)}] "
            f"{pdb_id}"
        )

        metadata = fetch_entry_metadata(
            pdb_id
        )

        title = get_title(metadata)
        method = get_method(metadata)
        resolution = get_resolution(metadata)
        release_date = get_release_date(metadata)

        search_matches = ";".join(
            sorted(hits[pdb_id])
        )

        title_lower = title.lower()

        # ---------------------------------------------------------
        # Simple preliminary relevance flags.
        #
        # These are intentionally conservative. We will inspect
        # actual polymer/entity composition in the next stage.
        # ---------------------------------------------------------

        mentions_nkg2a = (
            "nkg2a" in title_lower
            or "nkg2-a" in title_lower
        )

        mentions_cd94 = (
            "cd94" in title_lower
        )

        mentions_hlae = (
            "hla-e" in title_lower
            or "hla e" in title_lower
        )

        rows.append(
            {
                "pdb_id": pdb_id,
                "search_matches": search_matches,
                "title": title,
                "experimental_method": method,
                "resolution_angstrom": resolution,
                "initial_release_date": release_date,
                "title_mentions_NKG2A":
                    "yes" if mentions_nkg2a else "no",
                "title_mentions_CD94":
                    "yes" if mentions_cd94 else "no",
                "title_mentions_HLA_E":
                    "yes" if mentions_hlae else "no",
            }
        )

    # -------------------------------------------------------------
    # Write table
    # -------------------------------------------------------------

    fields = [
        "pdb_id",
        "search_matches",
        "title",
        "experimental_method",
        "resolution_angstrom",
        "initial_release_date",
        "title_mentions_NKG2A",
        "title_mentions_CD94",
        "title_mentions_HLA_E",
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

    print("\nSearch complete.")

    print(
        f"\nResults written to:\n"
        f"{OUTPUT_FILE}"
    )

    # -------------------------------------------------------------
    # Display likely relevant entries
    # -------------------------------------------------------------

    print(
        "\nLIKELY NKG2A STRUCTURAL REFERENCES"
    )

    relevant = []

    for row in rows:

        if (
            row["title_mentions_NKG2A"] == "yes"
            or (
                row["title_mentions_CD94"] == "yes"
                and row["title_mentions_HLA_E"] == "yes"
            )
        ):
            relevant.append(row)

    if not relevant:

        print(
            "No entries passed the preliminary "
            "title filter."
        )

    else:

        for row in relevant:

            print()
            print(
                f"PDB:        {row['pdb_id']}"
            )
            print(
                f"Method:     "
                f"{row['experimental_method']}"
            )
            print(
                f"Resolution: "
                f"{row['resolution_angstrom']}"
            )
            print(
                f"Released:   "
                f"{row['initial_release_date']}"
            )
            print(
                f"Title:      {row['title']}"
            )

    print(
        "\nNOTE: These are search candidates only."
    )
    print(
        "Do not select the primary structural reference "
        "until polymer/entity composition and chain "
        "assignments have been verified."
    )


if __name__ == "__main__":
    main()
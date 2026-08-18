from pathlib import Path
import urllib.request
import hashlib
import csv


ROOT = Path(__file__).resolve().parent.parent

STRUCTURE_DIR = (
    ROOT
    / "structures"
    / "reference"
)

TABLE_DIR = (
    ROOT
    / "results"
    / "tables"
    / "structure"
)

MANIFEST = (
    TABLE_DIR
    / "structural_reference_manifest.tsv"
)


STRUCTURES = [
    {
        "pdb_id": "3BDW",
        "role": "unliganded_receptor_reference",
        "description":
            "Human CD94/NKG2A receptor; 2.5 A X-ray structure",
    },
    {
        "pdb_id": "3CDG",
        "role": "ligand_bound_primary_reference",
        "description":
            "Human CD94/NKG2A bound to HLA-E and leader peptide; "
            "3.4 A X-ray structure",
    },
]


def sha256_file(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def download_file(url, output_path):

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent":
                "NKG2A-structural-analysis/1.0"
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=60,
    ) as response:

        data = response.read()

    output_path.write_bytes(data)

    return len(data)


def main():

    STRUCTURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    TABLE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 78)
    print("STEP 2B - DOWNLOAD STRUCTURAL REFERENCES")
    print("=" * 78)

    rows = []

    for index, entry in enumerate(
        STRUCTURES,
        start=1,
    ):

        pdb_id = entry["pdb_id"]

        print(
            f"\n[{index}/{len(STRUCTURES)}] "
            f"{pdb_id}"
        )

        # ---------------------------------------------------------
        # Download PDBx/mmCIF.
        #
        # mmCIF is preferred over legacy PDB because it preserves
        # modern chain/entity identifiers and metadata.
        # ---------------------------------------------------------

        url = (
            f"https://files.rcsb.org/download/"
            f"{pdb_id}.cif"
        )

        output_path = (
            STRUCTURE_DIR
            / f"{pdb_id}.cif"
        )

        print(f"URL:    {url}")
        print(f"Output: {output_path}")

        try:
            size = download_file(
                url,
                output_path,
            )

            checksum = sha256_file(
                output_path
            )

            status = "downloaded"

            print(
                f"OK - {size:,} bytes"
            )

        except Exception as exc:

            size = ""
            checksum = ""
            status = f"error: {exc}"

            print(
                f"ERROR - {exc}"
            )

        rows.append(
            {
                "pdb_id": pdb_id,
                "role": entry["role"],
                "description":
                    entry["description"],
                "format": "mmCIF",
                "source_url": url,
                "local_file":
                    str(output_path.relative_to(ROOT)),
                "file_size_bytes": size,
                "sha256": checksum,
                "status": status,
            }
        )

    fields = [
        "pdb_id",
        "role",
        "description",
        "format",
        "source_url",
        "local_file",
        "file_size_bytes",
        "sha256",
        "status",
    ]

    with MANIFEST.open(
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

    print(
        "\nManifest written to:"
    )
    print(MANIFEST)


if __name__ == "__main__":
    main()
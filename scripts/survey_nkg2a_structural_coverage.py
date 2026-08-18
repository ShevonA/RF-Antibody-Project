from pathlib import Path
import csv
import re
import urllib.request

from Bio import Align
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa


# =============================================================================
# STEP 2P - SURVEY EXPERIMENTAL NKG2A STRUCTURAL COVERAGE
# =============================================================================
#
# Goal
# ----
# Determine whether available human NKG2A experimental structures resolve
# the N-terminal extracellular region that is missing from 3BDW/3CDG.
#
# Of particular interest:
#
#     human NKG2A full-length residues 94-112
#
# These positions contain several strong human NKG2A-vs-NKG2C sequence
# differences identified during Step 2O.
#
# This step surveys experimental structures before considering predicted
# structural models.
# =============================================================================


ROOT = Path(__file__).resolve().parent.parent

STRUCTURE_DIR = (
    ROOT
    / "structures"
    / "reference"
    / "coverage_survey"
)

ECTODOMAIN_FILE = (
    ROOT
    / "results"
    / "tables"
    / "primary_ectodomain_sequences.tsv"
)

OUTPUT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_experimental_structure_coverage.tsv"
)


# =============================================================================
# STRUCTURES TO SURVEY
# =============================================================================

PDB_IDS = [
    "3BDW",
    "3CDG",
    "3CII",
]


TARGET_REGION_START = 94
TARGET_REGION_END = 112


# =============================================================================
# AMINO-ACID HELPERS
# =============================================================================

THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
}


def clean(value):
    if value is None:
        return ""

    return str(value).strip()


def clean_sequence(value):
    """
    Normalize mmCIF polymer sequence text.
    """

    text = clean(value)

    text = re.sub(
        r"\s+",
        "",
        text,
    )

    return text.upper()


# =============================================================================
# TABLE HELPERS
# =============================================================================

def read_tsv(path):

    if not path.exists():
        raise FileNotFoundError(
            f"Required input file not found:\n{path}"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        return list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )


def write_tsv(path, rows, fields):

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
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


# =============================================================================
# DOWNLOAD
# =============================================================================

def download_cif(pdb_id):

    STRUCTURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    path = (
        STRUCTURE_DIR
        / f"{pdb_id}.cif"
    )

    if path.exists():
        return path

    url = (
        "https://files.rcsb.org/download/"
        f"{pdb_id}.cif"
    )

    print(
        f"Downloading {pdb_id}..."
    )

    urllib.request.urlretrieve(
        url,
        path,
    )

    return path


# =============================================================================
# HUMAN NKG2A REFERENCE
# =============================================================================

def load_human_nkg2a_ectodomain():

    rows = read_tsv(
        ECTODOMAIN_FILE
    )

    for row in rows:

        if clean(
            row.get("record_id")
        ) != "human_NKG2A":
            continue

        sequence = clean(
            row.get("sequence")
        ).upper()

        start = int(
            row.get(
                "ectodomain_start"
            )
        )

        end = int(
            row.get(
                "ectodomain_end"
            )
        )

        return {
            "sequence": sequence,
            "start": start,
            "end": end,
        }

    raise ValueError(
        "human_NKG2A was not found in "
        "primary_ectodomain_sequences.tsv"
    )


# =============================================================================
# MMCIF HELPERS
# =============================================================================

def as_list(value):

    if isinstance(
        value,
        list,
    ):
        return value

    return [value]


def entity_table(cif_dict):

    ids = as_list(
        cif_dict.get(
            "_entity_poly.entity_id",
            [],
        )
    )

    sequences = as_list(
        cif_dict.get(
            "_entity_poly.pdbx_seq_one_letter_code_can",
            [],
        )
    )

    descriptions = as_list(
        cif_dict.get(
            "_entity.pdbx_description",
            [],
        )
    )

    entity_ids_for_descriptions = as_list(
        cif_dict.get(
            "_entity.id",
            [],
        )
    )

    description_lookup = {}

    for entity_id, description in zip(
        entity_ids_for_descriptions,
        descriptions,
    ):

        description_lookup[
            clean(entity_id)
        ] = clean(
            description
        )

    rows = []

    for entity_id, sequence in zip(
        ids,
        sequences,
    ):

        entity_id = clean(
            entity_id
        )

        rows.append(
            {
                "entity_id":
                    entity_id,

                "sequence":
                    clean_sequence(
                        sequence
                    ),

                "description":
                    description_lookup.get(
                        entity_id,
                        "",
                    ),
            }
        )

    return rows


def entity_chain_lookup(cif_dict):
    """
    Return:
        entity_id -> author chain IDs

    Uses _atom_site label/auth chain mapping.
    """

    label_asym = as_list(
        cif_dict.get(
            "_atom_site.label_asym_id",
            [],
        )
    )

    auth_asym = as_list(
        cif_dict.get(
            "_atom_site.auth_asym_id",
            [],
        )
    )

    label_entity = as_list(
        cif_dict.get(
            "_atom_site.label_entity_id",
            [],
        )
    )

    lookup = {}

    for (
        label_chain,
        auth_chain,
        entity_id,
    ) in zip(
        label_asym,
        auth_asym,
        label_entity,
    ):

        entity_id = clean(
            entity_id
        )

        auth_chain = clean(
            auth_chain
        )

        if not entity_id:
            continue

        if not auth_chain:
            continue

        lookup.setdefault(
            entity_id,
            set(),
        ).add(
            auth_chain
        )

    return lookup


def find_nkg2a_entity(
    entities,
    reference_sequence,
):
    """
    Identify the polymer entity most consistent with NKG2A.

    Primary cue:
        description contains NKG2-A / NKG2A.

    Fallback:
        local sequence alignment to human NKG2A ectodomain.
    """

    direct = []

    for entity in entities:

        description = (
            entity[
                "description"
            ].upper()
        )

        if (
            "NKG2-A" in description
            or "NKG2A" in description
            or "NKG2-A/NKG2-B" in description
        ):

            direct.append(
                entity
            )

    if len(direct) == 1:
        return direct[0]

    # -------------------------------------------------------------------------
    # Sequence fallback
    # -------------------------------------------------------------------------

    aligner = Align.PairwiseAligner()

    aligner.mode = "local"

    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5

    best = None

    for entity in entities:

        sequence = entity[
            "sequence"
        ]

        if not sequence:
            continue

        score = aligner.score(
            reference_sequence,
            sequence,
        )

        if (
            best is None
            or score > best[
                "score"
            ]
        ):

            best = {
                "entity":
                    entity,
                "score":
                    score,
            }

    if best is None:
        return None

    return best[
        "entity"
    ]


# =============================================================================
# ALIGN ENTITY SEQUENCE TO HUMAN NKG2A
# =============================================================================

def map_entity_to_reference(
    reference_sequence,
    reference_start,
    entity_sequence,
):
    """
    Align structural entity sequence against full human NKG2A ectodomain.

    Returns:
        entity position -> full-length human NKG2A position
    """

    aligner = Align.PairwiseAligner()

    aligner.mode = "local"

    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(
        reference_sequence,
        entity_sequence,
    )

    if len(alignments) == 0:
        raise ValueError(
            "No sequence alignment was produced."
        )

    alignment = alignments[0]

    mapping = {}

    reference_blocks = (
        alignment.aligned[0]
    )

    entity_blocks = (
        alignment.aligned[1]
    )

    for (
        ref_block,
        entity_block,
    ) in zip(
        reference_blocks,
        entity_blocks,
    ):

        ref_start_index = int(
            ref_block[0]
        )

        ref_end_index = int(
            ref_block[1]
        )

        entity_start_index = int(
            entity_block[0]
        )

        entity_end_index = int(
            entity_block[1]
        )

        ref_length = (
            ref_end_index
            - ref_start_index
        )

        entity_length = (
            entity_end_index
            - entity_start_index
        )

        if ref_length != entity_length:
            continue

        for offset in range(
            ref_length
        ):

            entity_position = (
                entity_start_index
                + offset
                + 1
            )

            full_length_position = (
                reference_start
                + ref_start_index
                + offset
            )

            mapping[
                entity_position
            ] = (
                full_length_position
            )

    return mapping, alignment.score


# =============================================================================
# COORDINATE SEQUENCE / COVERAGE
# =============================================================================

def coordinate_residues(chain):
    """
    Return amino-acid coordinate residues from one chain.
    """

    residues = []

    for residue in chain:

        if not is_aa(
            residue,
            standard=False,
        ):
            continue

        aa = THREE_TO_ONE.get(
            residue.get_resname().upper(),
            "X",
        )

        residues.append(
            {
                "residue":
                    residue,

                "aa":
                    aa,

                "auth_number":
                    residue.id[1],

                "insertion_code":
                    clean(
                        residue.id[2]
                    ),
            }
        )

    return residues


def map_coordinate_sequence_to_entity(
    coordinate_rows,
    entity_sequence,
):
    """
    Align the coordinate-derived amino-acid sequence to the polymer entity.

    Returns:
        coordinate row index -> entity sequence position
    """

    coordinate_sequence = "".join(
        row["aa"]
        for row in coordinate_rows
    )

    aligner = Align.PairwiseAligner()

    aligner.mode = "global"

    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(
        entity_sequence,
        coordinate_sequence,
    )[0]

    mapping = {}

    entity_blocks = (
        alignment.aligned[0]
    )

    coordinate_blocks = (
        alignment.aligned[1]
    )

    for (
        entity_block,
        coord_block,
    ) in zip(
        entity_blocks,
        coordinate_blocks,
    ):

        entity_start = int(
            entity_block[0]
        )

        entity_end = int(
            entity_block[1]
        )

        coord_start = int(
            coord_block[0]
        )

        coord_end = int(
            coord_block[1]
        )

        entity_length = (
            entity_end
            - entity_start
        )

        coord_length = (
            coord_end
            - coord_start
        )

        if entity_length != coord_length:
            continue

        for offset in range(
            entity_length
        ):

            coordinate_index = (
                coord_start
                + offset
            )

            entity_position = (
                entity_start
                + offset
                + 1
            )

            mapping[
                coordinate_index
            ] = entity_position

    return mapping, alignment.score


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "STEP 2P - SURVEY EXPERIMENTAL HUMAN NKG2A STRUCTURAL COVERAGE"
    )
    print("=" * 78)

    reference = (
        load_human_nkg2a_ectodomain()
    )

    reference_sequence = (
        reference[
            "sequence"
        ]
    )

    reference_start = (
        reference[
            "start"
        ]
    )

    reference_end = (
        reference[
            "end"
        ]
    )

    print()
    print(
        "Human NKG2A ectodomain: "
        f"{reference_start}-{reference_end}"
    )

    print(
        "Critical unresolved region: "
        f"{TARGET_REGION_START}-{TARGET_REGION_END}"
    )

    output_rows = []

    for pdb_id in PDB_IDS:

        print()
        print("=" * 78)
        print(
            f"PDB {pdb_id}"
        )
        print("=" * 78)

        cif_path = download_cif(
            pdb_id
        )

        cif_dict = MMCIF2Dict(
            str(cif_path)
        )

        entities = entity_table(
            cif_dict
        )

        chains_by_entity = (
            entity_chain_lookup(
                cif_dict
            )
        )

        nkg2a_entity = (
            find_nkg2a_entity(
                entities,
                reference_sequence,
            )
        )

        if nkg2a_entity is None:

            print(
                "No NKG2A entity identified."
            )

            continue

        entity_id = (
            nkg2a_entity[
                "entity_id"
            ]
        )

        entity_sequence = (
            nkg2a_entity[
                "sequence"
            ]
        )

        description = (
            nkg2a_entity[
                "description"
            ]
        )

        chains = sorted(
            chains_by_entity.get(
                entity_id,
                set(),
            )
        )

        print(
            f"Entity:      {entity_id}"
        )

        print(
            f"Description: {description}"
        )

        print(
            f"Entity length: "
            f"{len(entity_sequence)} aa"
        )

        print(
            "Author chains: "
            + (
                ",".join(chains)
                if chains
                else "none"
            )
        )

        (
            entity_to_reference,
            entity_alignment_score,
        ) = map_entity_to_reference(
            reference_sequence,
            reference_start,
            entity_sequence,
        )

        mapped_reference_positions = (
            sorted(
                entity_to_reference.values()
            )
        )

        if mapped_reference_positions:

            entity_ref_start = min(
                mapped_reference_positions
            )

            entity_ref_end = max(
                mapped_reference_positions
            )

        else:

            entity_ref_start = None
            entity_ref_end = None

        print(
            "Entity maps to human NKG2A: "
            f"{entity_ref_start}-{entity_ref_end}"
        )

        parser = MMCIFParser(
            QUIET=True
        )

        structure = parser.get_structure(
            pdb_id,
            str(cif_path),
        )

        model = structure[0]

        for chain_id in chains:

            if chain_id not in model:

                print(
                    f"WARNING: author chain "
                    f"{chain_id} not found "
                    "by Bio.PDB."
                )

                continue

            chain = model[
                chain_id
            ]

            coordinate_rows = (
                coordinate_residues(
                    chain
                )
            )

            if not coordinate_rows:

                continue

            (
                coord_to_entity,
                coordinate_alignment_score,
            ) = map_coordinate_sequence_to_entity(
                coordinate_rows,
                entity_sequence,
            )

            coordinate_full_positions = []

            for coordinate_index in sorted(
                coord_to_entity
            ):

                entity_position = (
                    coord_to_entity[
                        coordinate_index
                    ]
                )

                full_position = (
                    entity_to_reference.get(
                        entity_position
                    )
                )

                if full_position is None:
                    continue

                coordinate_full_positions.append(
                    full_position
                )

            coordinate_full_positions = (
                sorted(
                    set(
                        coordinate_full_positions
                    )
                )
            )

            if coordinate_full_positions:

                coord_start = min(
                    coordinate_full_positions
                )

                coord_end = max(
                    coordinate_full_positions
                )

            else:

                coord_start = None
                coord_end = None

            target_positions_present = [
                position
                for position
                in range(
                    TARGET_REGION_START,
                    TARGET_REGION_END + 1,
                )
                if position
                in coordinate_full_positions
            ]

            target_positions_missing = [
                position
                for position
                in range(
                    TARGET_REGION_START,
                    TARGET_REGION_END + 1,
                )
                if position
                not in coordinate_full_positions
            ]

            target_fraction = (
                len(
                    target_positions_present
                )
                / (
                    TARGET_REGION_END
                    - TARGET_REGION_START
                    + 1
                )
            )

            print()
            print(
                f"Chain {chain_id}"
            )

            print(
                "Coordinate amino-acid residues: "
                f"{len(coordinate_rows)}"
            )

            print(
                "Mapped full-length coverage: "
                f"{coord_start}-{coord_end}"
            )

            print(
                f"94-112 residues resolved: "
                f"{len(target_positions_present)}/"
                f"{TARGET_REGION_END - TARGET_REGION_START + 1}"
            )

            if target_positions_present:

                print(
                    "Resolved target positions: "
                    + ",".join(
                        str(x)
                        for x
                        in target_positions_present
                    )
                )

            if target_positions_missing:

                print(
                    "Missing target positions:  "
                    + ",".join(
                        str(x)
                        for x
                        in target_positions_missing
                    )
                )

            output_rows.append(
                {
                    "pdb_id":
                        pdb_id,

                    "entity_id":
                        entity_id,

                    "description":
                        description,

                    "author_chain_id":
                        chain_id,

                    "entity_sequence_length":
                        len(
                            entity_sequence
                        ),

                    "entity_reference_start":
                        (
                            entity_ref_start
                            if entity_ref_start
                            is not None
                            else ""
                        ),

                    "entity_reference_end":
                        (
                            entity_ref_end
                            if entity_ref_end
                            is not None
                            else ""
                        ),

                    "coordinate_residue_count":
                        len(
                            coordinate_rows
                        ),

                    "coordinate_reference_start":
                        (
                            coord_start
                            if coord_start
                            is not None
                            else ""
                        ),

                    "coordinate_reference_end":
                        (
                            coord_end
                            if coord_end
                            is not None
                            else ""
                        ),

                    "target_region_start":
                        TARGET_REGION_START,

                    "target_region_end":
                        TARGET_REGION_END,

                    "target_positions_resolved_count":
                        len(
                            target_positions_present
                        ),

                    "target_positions_total":
                        (
                            TARGET_REGION_END
                            - TARGET_REGION_START
                            + 1
                        ),

                    "target_fraction_resolved":
                        f"{target_fraction:.3f}",

                    "target_positions_resolved":
                        ",".join(
                            str(x)
                            for x
                            in target_positions_present
                        ),

                    "target_positions_missing":
                        ",".join(
                            str(x)
                            for x
                            in target_positions_missing
                        ),

                    "entity_alignment_score":
                        f"{entity_alignment_score:.3f}",

                    "coordinate_alignment_score":
                        f"{coordinate_alignment_score:.3f}",
                }
            )

    # =========================================================================
    # WRITE OUTPUT
    # =========================================================================

    fields = [
        "pdb_id",
        "entity_id",
        "description",
        "author_chain_id",
        "entity_sequence_length",
        "entity_reference_start",
        "entity_reference_end",
        "coordinate_residue_count",
        "coordinate_reference_start",
        "coordinate_reference_end",
        "target_region_start",
        "target_region_end",
        "target_positions_resolved_count",
        "target_positions_total",
        "target_fraction_resolved",
        "target_positions_resolved",
        "target_positions_missing",
        "entity_alignment_score",
        "coordinate_alignment_score",
    ]

    write_tsv(
        OUTPUT_FILE,
        output_rows,
        fields,
    )

    # =========================================================================
    # SUMMARY
    # =========================================================================

    print()
    print("=" * 78)
    print("EXPERIMENTAL COVERAGE SUMMARY")
    print("=" * 78)

    ranked = sorted(
        output_rows,
        key=lambda row: (
            -int(
                row[
                    "target_positions_resolved_count"
                ]
            ),
            row[
                "pdb_id"
            ],
            row[
                "author_chain_id"
            ],
        )
    )

    for row in ranked:

        print(
            f"{row['pdb_id']} "
            f"chain {row['author_chain_id']}: "
            f"{row['target_positions_resolved_count']}/"
            f"{row['target_positions_total']} "
            f"positions from "
            f"{TARGET_REGION_START}-{TARGET_REGION_END} "
            "resolved"
        )

    best_count = max(
        (
            int(
                row[
                    "target_positions_resolved_count"
                ]
            )
            for row in output_rows
        ),
        default=0,
    )

    print()
    print("=" * 78)
    print("STEP 2P DECISION")
    print("=" * 78)

    if best_count == 0:

        print()
        print(
            "No surveyed experimental human NKG2A "
            f"structure resolves residues "
            f"{TARGET_REGION_START}-{TARGET_REGION_END}."
        )

        print(
            "The next step should evaluate a predicted "
            "full ectodomain model for this region."
        )

    elif best_count < (
        TARGET_REGION_END
        - TARGET_REGION_START
        + 1
    ):

        print()
        print(
            "At least one experimental structure provides "
            "partial coverage of the unresolved region."
        )

        print(
            "Inspect that structure before deciding whether "
            "predicted modeling is still required."
        )

    else:

        print()
        print(
            "At least one experimental structure resolves "
            "the complete 94-112 region."
        )

        print(
            "Use the experimental structure before considering "
            "predicted modeling."
        )

    print()
    print("=" * 78)
    print("OUTPUT")
    print("=" * 78)

    print()
    print(OUTPUT_FILE)

    print()
    print(
        "NOTE: coverage is mapped to canonical human "
        "NKG2A full-length numbering using sequence alignment."
    )

    print(
        "Experimental coordinate coverage is distinguished "
        "from polymer entity sequence coverage; residues present "
        "in the construct but lacking coordinates are treated "
        "as unresolved."
    )


if __name__ == "__main__":
    main()
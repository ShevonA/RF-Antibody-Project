from pathlib import Path
import csv
import json
import math
import urllib.request

from Bio import Align
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa


# =============================================================================
# STEP 2Q - ALPHAFOLD HUMAN NKG2A MODEL ACQUISITION AND VALIDATION
# =============================================================================
#
# Purpose
# -------
# 1. Download the current AlphaFold DB prediction for human NKG2A P26715.
# 2. Verify model sequence/coverage against the curated human NKG2A sequence.
# 3. Extract per-residue pLDDT, especially residues 94-112.
# 4. Compare the AlphaFold model to experimental 3CDG over the region that
#    both structures actually contain.
# 5. Decide whether residues 94-112 are confident enough for subsequent
#    surface-accessibility / epitope-geometry analysis.
#
# Important
# ---------
# AlphaFold coordinates exist even for low-confidence residues. Coordinates
# in a low-pLDDT region must not automatically be interpreted as a stable
# antibody epitope geometry.
# =============================================================================


ROOT = Path(__file__).resolve().parent.parent

UNIPROT_ID = "P26715"

MODEL_DIR = (
    ROOT
    / "structures"
    / "models"
    / "alphafold"
)

MODEL_CIF = (
    MODEL_DIR
    / f"{UNIPROT_ID}_alphafold.cif"
)

MODEL_PAE = (
    MODEL_DIR
    / f"{UNIPROT_ID}_alphafold_pae.json"
)

EXPERIMENTAL_CIF = (
    ROOT
    / "structures"
    / "reference"
    / "3CDG.cif"
)

ECTODOMAIN_FILE = (
    ROOT
    / "results"
    / "tables"
    / "primary_ectodomain_sequences.tsv"
)

RESIDUE_OUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "alphafold_nkg2a_residue_confidence.tsv"
)

SUMMARY_OUT = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "alphafold_nkg2a_model_validation.tsv"
)


TARGET_START = 94
TARGET_END = 112

EXPERIMENTAL_CHAIN = "F"

# AlphaFold DB API endpoint.
API_URL = (
    "https://alphafold.ebi.ac.uk/api/prediction/"
    f"{UNIPROT_ID}"
)


# =============================================================================
# BASIC HELPERS
# =============================================================================


def clean(value):

    if value is None:
        return ""

    return str(value).strip()


def mean(values):

    values = [
        value
        for value in values
        if value is not None
    ]

    if not values:
        return None

    return sum(values) / len(values)


def median(values):

    values = sorted(
        value
        for value in values
        if value is not None
    )

    if not values:
        return None

    n = len(values)

    middle = n // 2

    if n % 2:
        return values[middle]

    return (
        values[middle - 1]
        + values[middle]
    ) / 2.0


def confidence_class(plddt):

    if plddt is None:
        return "unknown"

    if plddt > 90:
        return "very_high"

    if plddt >= 70:
        return "confident"

    if plddt >= 50:
        return "low"

    return "very_low"


def write_tsv(
    path,
    rows,
    fieldnames,
):

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
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            rows
        )


# =============================================================================
# LOAD CURATED HUMAN NKG2A
# =============================================================================


def load_human_nkg2a_ectodomain():

    with ECTODOMAIN_FILE.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        rows = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    for row in rows:

        if clean(
            row.get(
                "record_id"
            )
        ) != "human_NKG2A":
            continue

        return {
            "sequence":
                clean(
                    row.get(
                        "sequence"
                    )
                ).upper(),

            "start":
                int(
                    row[
                        "ectodomain_start"
                    ]
                ),

            "end":
                int(
                    row[
                        "ectodomain_end"
                    ]
                ),
        }

    raise ValueError(
        "human_NKG2A not found in "
        "primary_ectodomain_sequences.tsv"
    )


# =============================================================================
# ALPHAFOLD DOWNLOAD
# =============================================================================


def load_alphafold_metadata():

    print(
        "Querying AlphaFold DB metadata..."
    )

    with urllib.request.urlopen(
        API_URL
    ) as response:

        payload = json.loads(
            response.read().decode(
                "utf-8"
            )
        )

    if not payload:
        raise ValueError(
            f"No AlphaFold DB entry returned for {UNIPROT_ID}"
        )

    # Canonical entry should normally be first.
    entry = payload[0]

    return entry


def download_file(
    url,
    destination,
):

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if destination.exists():

        print(
            f"Using existing file: "
            f"{destination}"
        )

        return

    print(
        f"Downloading:\n  {url}"
    )

    print(
        f"To:\n  {destination}"
    )

    urllib.request.urlretrieve(
        url,
        destination,
    )


# =============================================================================
# STRUCTURE / SEQUENCE EXTRACTION
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


def residue_one_letter(
    residue,
):

    return THREE_TO_ONE.get(
        residue.get_resname().upper(),
        "X",
    )


def amino_acid_residues(
    chain,
):

    return [
        residue
        for residue in chain
        if is_aa(
            residue,
            standard=False,
        )
    ]


def select_alphafold_chain(
    model,
):

    chains = list(
        model.get_chains()
    )

    if not chains:
        raise ValueError(
            "AlphaFold model contains no chains."
        )

    # AFDB monomer predictions normally contain one chain.
    if len(chains) > 1:

        print(
            "WARNING: AlphaFold coordinate file contains "
            f"{len(chains)} chains; using the first."
        )

    return chains[0]


def residue_plddt(
    residue,
):
    """
    AlphaFold DB stores pLDDT in B-factor fields.

    Use the mean over residue atoms. These values should normally
    be identical or nearly identical for all atoms of a residue.
    """

    values = [
        float(
            atom.get_bfactor()
        )
        for atom in residue.get_atoms()
    ]

    return mean(values)


# =============================================================================
# MAP MODEL SEQUENCE TO CANONICAL FULL-LENGTH NUMBERING
# =============================================================================


def map_model_to_reference(
    model_sequence,
    reference_sequence,
    reference_start,
):
    """
    Map AlphaFold model sequence positions to canonical human NKG2A
    full-length numbering.

    The reference sequence supplied here is the curated ectodomain
    (94-233), so the resulting coordinates are full-length positions.
    """

    aligner = Align.PairwiseAligner()

    aligner.mode = "local"

    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(
        reference_sequence,
        model_sequence,
    )[0]

    mapping = {}

    ref_blocks = (
        alignment.aligned[0]
    )

    model_blocks = (
        alignment.aligned[1]
    )

    for (
        ref_block,
        model_block,
    ) in zip(
        ref_blocks,
        model_blocks,
    ):

        ref_start_index = int(
            ref_block[0]
        )

        ref_end_index = int(
            ref_block[1]
        )

        model_start_index = int(
            model_block[0]
        )

        model_end_index = int(
            model_block[1]
        )

        ref_len = (
            ref_end_index
            - ref_start_index
        )

        model_len = (
            model_end_index
            - model_start_index
        )

        if ref_len != model_len:
            continue

        for offset in range(
            ref_len
        ):

            model_index = (
                model_start_index
                + offset
            )

            full_length_position = (
                reference_start
                + ref_start_index
                + offset
            )

            mapping[
                model_index
            ] = (
                full_length_position
            )

    return mapping, alignment.score


# =============================================================================
# EXPERIMENTAL COMPARISON
# =============================================================================


def ca_lookup_by_sequence_position(
    residues,
):
    """
    Map zero-based sequence position to CA coordinate.
    """

    lookup = {}

    for index, residue in enumerate(
        residues
    ):

        if "CA" not in residue:
            continue

        lookup[
            index
        ] = residue["CA"].coord

    return lookup


def kabsch_rmsd(
    mobile_points,
    target_points,
):
    """
    Return optimal C-alpha RMSD after Kabsch superposition.

    Uses NumPy through BioPython's dependency.
    """

    import numpy as np

    mobile = np.asarray(
        mobile_points,
        dtype=float,
    )

    target = np.asarray(
        target_points,
        dtype=float,
    )

    if len(mobile) < 3:
        return None

    mobile_center = mobile.mean(
        axis=0
    )

    target_center = target.mean(
        axis=0
    )

    mobile_centered = (
        mobile
        - mobile_center
    )

    target_centered = (
        target
        - target_center
    )

    covariance = (
        mobile_centered.T
        @ target_centered
    )

    u, s, vt = np.linalg.svd(
        covariance
    )

    rotation = (
        vt.T
        @ u.T
    )

    if np.linalg.det(
        rotation
    ) < 0:

        vt[-1, :] *= -1

        rotation = (
            vt.T
            @ u.T
        )

    fitted = (
        mobile_centered
        @ rotation
    )

    diff = (
        fitted
        - target_centered
    )

    rmsd = math.sqrt(
        float(
            (
                diff * diff
            ).sum()
            / len(mobile)
        )
    )

    return rmsd


def compare_to_3cdg(
    af_residue_rows,
):
    """
    Compare AlphaFold residues to 3CDG chain F over common
    full-length positions.

    3CDG NKG2A chain F corresponds to human residues 113-232.
    """

    parser = MMCIFParser(
        QUIET=True
    )

    experimental = (
        parser.get_structure(
            "3CDG",
            str(
                EXPERIMENTAL_CIF
            ),
        )
    )

    chain = (
        experimental[0][
            EXPERIMENTAL_CHAIN
        ]
    )

    exp_residues = (
        amino_acid_residues(
            chain
        )
    )

    # From validated Step 2D:
    # structural entity position 1 -> full residue 113.
    exp_lookup = {}

    for index, residue in enumerate(
        exp_residues
    ):

        full_position = (
            113
            + index
        )

        if "CA" not in residue:
            continue

        exp_lookup[
            full_position
        ] = residue["CA"].coord

    af_lookup = {}

    for row in af_residue_rows:

        full_position = row[
            "full_length_residue"
        ]

        residue = row[
            "residue_object"
        ]

        if (
            full_position is None
            or "CA" not in residue
        ):
            continue

        af_lookup[
            full_position
        ] = residue["CA"].coord

    shared = sorted(
        set(
            exp_lookup
        )
        & set(
            af_lookup
        )
    )

    mobile = [
        af_lookup[position]
        for position in shared
    ]

    target = [
        exp_lookup[position]
        for position in shared
    ]

    rmsd = kabsch_rmsd(
        mobile,
        target,
    )

    return shared, rmsd


# =============================================================================
# OPTIONAL PAE SUMMARY
# =============================================================================


def load_pae_summary(
    path,
):
    """
    Return the PAE matrix if the downloaded AlphaFold JSON follows
    the standard AFDB monomer format.

    The script does not require PAE to succeed.
    """

    if not path.exists():
        return None

    try:

        payload = json.loads(
            path.read_text(
                encoding="utf-8"
            )
        )

    except Exception:

        return None

    if isinstance(
        payload,
        list,
    ) and payload:

        matrix = payload[0].get(
            "predicted_aligned_error"
        )

        return matrix

    if isinstance(
        payload,
        dict,
    ):

        return payload.get(
            "predicted_aligned_error"
        )

    return None


def target_region_pae(
    matrix,
    target_indices,
    domain_indices,
):
    """
    Summarize PAE between target region 94-112 and the resolved
    C-terminal domain.

    Lower PAE means greater confidence in relative positioning.
    """

    if matrix is None:
        return None

    values = []

    n = len(matrix)

    for target_index in target_indices:

        if target_index >= n:
            continue

        for domain_index in domain_indices:

            if domain_index >= len(
                matrix[
                    target_index
                ]
            ):
                continue

            try:
                values.append(
                    float(
                        matrix[
                            target_index
                        ][
                            domain_index
                        ]
                    )
                )

            except (
                TypeError,
                ValueError,
            ):
                continue

    return mean(values)


# =============================================================================
# MAIN
# =============================================================================


def main():

    print("=" * 78)
    print(
        "STEP 2Q - ALPHAFOLD HUMAN NKG2A MODEL VALIDATION"
    )
    print("=" * 78)

    print()
    print(
        f"UniProt: {UNIPROT_ID}"
    )

    print(
        f"Target region: "
        f"{TARGET_START}-{TARGET_END}"
    )

    reference = (
        load_human_nkg2a_ectodomain()
    )

    # -------------------------------------------------------------------------
    # AlphaFold metadata and download URLs
    # -------------------------------------------------------------------------

    metadata = (
        load_alphafold_metadata()
    )

    print()
    print("AlphaFold DB entry:")

    print(
        f"  Entry ID: "
        f"{clean(metadata.get('entryId'))}"
    )

    print(
        f"  UniProt:  "
        f"{clean(metadata.get('uniprotAccession'))}"
    )

    print(
        f"  Gene:     "
        f"{clean(metadata.get('gene'))}"
    )

    print(
        f"  Organism: "
        f"{clean(metadata.get('organismScientificName'))}"
    )

    cif_url = (
        clean(
            metadata.get(
                "cifUrl"
            )
        )
    )

    pdb_url = (
        clean(
            metadata.get(
                "pdbUrl"
            )
        )
    )

    pae_url = (
        clean(
            metadata.get(
                "paeDocUrl"
            )
        )
    )

    if not cif_url:

        raise ValueError(
            "AlphaFold DB API response did not include cifUrl."
        )

    download_file(
        cif_url,
        MODEL_CIF,
    )

    if pae_url:

        try:

            download_file(
                pae_url,
                MODEL_PAE,
            )

        except Exception as error:

            print()
            print(
                "WARNING: PAE download failed."
            )

            print(
                f"Reason: {error}"
            )

    # -------------------------------------------------------------------------
    # Parse AlphaFold model
    # -------------------------------------------------------------------------

    parser = MMCIFParser(
        QUIET=True
    )

    structure = (
        parser.get_structure(
            "AF_P26715",
            str(
                MODEL_CIF
            ),
        )
    )

    model = structure[0]

    chain = (
        select_alphafold_chain(
            model
        )
    )

    residues = (
        amino_acid_residues(
            chain
        )
    )

    model_sequence = "".join(
        residue_one_letter(
            residue
        )
        for residue in residues
    )

    print()
    print(
        "AlphaFold coordinate residues: "
        f"{len(residues)}"
    )

    (
        model_to_reference,
        sequence_alignment_score,
    ) = map_model_to_reference(
        model_sequence,
        reference[
            "sequence"
        ],
        reference[
            "start"
        ],
    )

    # -------------------------------------------------------------------------
    # Build residue confidence table
    # -------------------------------------------------------------------------

    residue_rows = []

    mapped_positions = []

    for model_index, residue in enumerate(
        residues
    ):

        full_position = (
            model_to_reference.get(
                model_index
            )
        )

        if full_position is not None:

            mapped_positions.append(
                full_position
            )

        plddt = residue_plddt(
            residue
        )

        residue_rows.append(
            {
                "model_index":
                    model_index,

                "model_residue_number":
                    residue.id[1],

                "full_length_residue":
                    full_position,

                "aa":
                    residue_one_letter(
                        residue
                    ),

                "plddt":
                    plddt,

                "confidence_class":
                    confidence_class(
                        plddt
                    ),

                "in_target_94_112":
                    (
                        "yes"
                        if (
                            full_position
                            is not None
                            and TARGET_START
                            <= full_position
                            <= TARGET_END
                        )
                        else "no"
                    ),

                "residue_object":
                    residue,
            }
        )

    if mapped_positions:

        mapped_start = min(
            mapped_positions
        )

        mapped_end = max(
            mapped_positions
        )

    else:

        mapped_start = None
        mapped_end = None

    print(
        "Mapped model coverage: "
        f"{mapped_start}-{mapped_end}"
    )

    # -------------------------------------------------------------------------
    # Target confidence
    # -------------------------------------------------------------------------

    target_rows = [
        row
        for row in residue_rows
        if (
            row[
                "full_length_residue"
            ]
            is not None
            and TARGET_START
            <= row[
                "full_length_residue"
            ]
            <= TARGET_END
        )
    ]

    target_plddts = [
        row[
            "plddt"
        ]
        for row in target_rows
    ]

    target_mean = mean(
        target_plddts
    )

    target_median = median(
        target_plddts
    )

    target_min = (
        min(
            target_plddts
        )
        if target_plddts
        else None
    )

    target_max = (
        max(
            target_plddts
        )
        if target_plddts
        else None
    )

    print()
    print("=" * 78)
    print(
        "TARGET REGION 94-112 pLDDT"
    )
    print("=" * 78)

    for row in target_rows:

        print(
            f"{row['full_length_residue']:>3} "
            f"{row['aa']}  "
            f"pLDDT={row['plddt']:>6.2f}  "
            f"{row['confidence_class']}"
        )

    print()
    print(
        f"Target mean pLDDT:   "
        f"{target_mean:.2f}"
        if target_mean is not None
        else
        "Target mean pLDDT:   NA"
    )

    print(
        f"Target median pLDDT: "
        f"{target_median:.2f}"
        if target_median is not None
        else
        "Target median pLDDT: NA"
    )

    print(
        f"Target minimum pLDDT:"
        f" {target_min:.2f}"
        if target_min is not None
        else
        "Target minimum pLDDT: NA"
    )

    print(
        f"Target maximum pLDDT:"
        f" {target_max:.2f}"
        if target_max is not None
        else
        "Target maximum pLDDT: NA"
    )

    # -------------------------------------------------------------------------
    # Compare model to 3CDG
    # -------------------------------------------------------------------------

    (
        shared_positions,
        rmsd,
    ) = compare_to_3cdg(
        residue_rows
    )

    print()
    print("=" * 78)
    print(
        "ALPHAFOLD vs 3CDG"
    )
    print("=" * 78)

    print(
        "Shared C-alpha positions: "
        f"{len(shared_positions)}"
    )

    if shared_positions:

        print(
            "Shared full-length range: "
            f"{min(shared_positions)}-"
            f"{max(shared_positions)}"
        )

    if rmsd is None:

        print(
            "C-alpha RMSD: NA"
        )

    else:

        print(
            f"C-alpha RMSD after superposition: "
            f"{rmsd:.3f} A"
        )

    # -------------------------------------------------------------------------
    # PAE
    # -------------------------------------------------------------------------

    pae_matrix = (
        load_pae_summary(
            MODEL_PAE
        )
    )

    target_model_indices = [
        row[
            "model_index"
        ]
        for row in target_rows
    ]

    domain_model_indices = [
        row[
            "model_index"
        ]
        for row in residue_rows
        if (
            row[
                "full_length_residue"
            ]
            is not None
            and 113
            <= row[
                "full_length_residue"
            ]
            <= 232
        )
    ]

    mean_target_domain_pae = (
        target_region_pae(
            pae_matrix,
            target_model_indices,
            domain_model_indices,
        )
    )

    print()
    print("=" * 78)
    print(
        "RELATIVE-POSITION CONFIDENCE"
    )
    print("=" * 78)

    if (
        mean_target_domain_pae
        is None
    ):

        print(
            "Mean 94-112 vs 113-232 PAE: "
            "not available"
        )

    else:

        print(
            "Mean 94-112 vs 113-232 PAE: "
            f"{mean_target_domain_pae:.2f} A"
        )

    # -------------------------------------------------------------------------
    # Write residue-level output
    # -------------------------------------------------------------------------

    residue_output = []

    for row in residue_rows:

        residue_output.append(
            {
                "model_residue_number":
                    row[
                        "model_residue_number"
                    ],

                "full_length_residue":
                    (
                        row[
                            "full_length_residue"
                        ]
                        if row[
                            "full_length_residue"
                        ]
                        is not None
                        else ""
                    ),

                "aa":
                    row[
                        "aa"
                    ],

                "plddt":
                    f"{row['plddt']:.3f}",

                "confidence_class":
                    row[
                        "confidence_class"
                    ],

                "in_target_94_112":
                    row[
                        "in_target_94_112"
                    ],
            }
        )

    write_tsv(
        RESIDUE_OUT,
        residue_output,
        [
            "model_residue_number",
            "full_length_residue",
            "aa",
            "plddt",
            "confidence_class",
            "in_target_94_112",
        ],
    )

    # -------------------------------------------------------------------------
    # Decision
    # -------------------------------------------------------------------------

    if target_mean is None:

        decision = (
            "target_region_not_mapped"
        )

    elif target_mean >= 70:

        decision = (
            "target_region_locally_confident_"
            "proceed_with_caution"
        )

    elif target_mean >= 50:

        decision = (
            "target_region_low_confidence_"
            "geometry_should_not_be_overinterpreted"
        )

    else:

        decision = (
            "target_region_very_low_confidence_"
            "do_not_use_as_fixed_epitope_geometry"
        )

    if (
        target_mean is not None
        and target_mean >= 70
        and mean_target_domain_pae
        is not None
        and mean_target_domain_pae > 10
    ):

        decision = (
            "local_structure_confident_but_relative_"
            "placement_uncertain"
        )

    summary_row = {
        "uniprot_id":
            UNIPROT_ID,

        "alphafold_entry_id":
            clean(
                metadata.get(
                    "entryId"
                )
            ),

        "model_coordinate_residues":
            len(
                residues
            ),

        "mapped_reference_start":
            (
                mapped_start
                if mapped_start
                is not None
                else ""
            ),

        "mapped_reference_end":
            (
                mapped_end
                if mapped_end
                is not None
                else ""
            ),

        "sequence_alignment_score":
            f"{sequence_alignment_score:.3f}",

        "target_start":
            TARGET_START,

        "target_end":
            TARGET_END,

        "target_residues_mapped":
            len(
                target_rows
            ),

        "target_mean_plddt":
            (
                f"{target_mean:.3f}"
                if target_mean
                is not None
                else ""
            ),

        "target_median_plddt":
            (
                f"{target_median:.3f}"
                if target_median
                is not None
                else ""
            ),

        "target_min_plddt":
            (
                f"{target_min:.3f}"
                if target_min
                is not None
                else ""
            ),

        "target_max_plddt":
            (
                f"{target_max:.3f}"
                if target_max
                is not None
                else ""
            ),

        "shared_3CDG_CA_positions":
            len(
                shared_positions
            ),

        "alphafold_vs_3CDG_CA_rmsd_A":
            (
                f"{rmsd:.3f}"
                if rmsd
                is not None
                else ""
            ),

        "mean_target_94_112_vs_domain_113_232_PAE_A":
            (
                f"{mean_target_domain_pae:.3f}"
                if mean_target_domain_pae
                is not None
                else ""
            ),

        "model_use_decision":
            decision,
    }

    write_tsv(
        SUMMARY_OUT,
        [
            summary_row
        ],
        list(
            summary_row.keys()
        ),
    )

    # -------------------------------------------------------------------------
    # Final report
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "STEP 2Q DECISION"
    )
    print("=" * 78)

    print()
    print(
        decision
    )

    print()
    print(
        "Interpretation guidance:"
    )

    print(
        "  pLDDT >90: high local confidence"
    )

    print(
        "  pLDDT 70-90: generally reliable backbone"
    )

    print(
        "  pLDDT 50-70: low confidence"
    )

    print(
        "  pLDDT <50: do not interpret as fixed geometry"
    )

    print()
    print(
        "Even with acceptable local pLDDT, high PAE between "
        "94-112 and the lectin domain would mean that the "
        "orientation of this region relative to the domain "
        "is uncertain."
    )

    print()
    print("=" * 78)
    print("OUTPUTS")
    print("=" * 78)

    print()
    print(MODEL_CIF)
    print(MODEL_PAE)
    print(RESIDUE_OUT)
    print(SUMMARY_OUT)


if __name__ == "__main__":
    main()
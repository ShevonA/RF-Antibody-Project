from pathlib import Path
import csv


# =============================================================================
# PATHS
# =============================================================================

ROOT = Path(__file__).resolve().parent.parent

SEQUENCE_FILE = (
    ROOT
    / "results"
    / "tables"
    / "candidate_specificity_positions.tsv"
)

CONTACT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_structural_contacts.tsv"
)

SASA_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_solvent_accessibility.tsv"
)

GLYCOSYLATION_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2_ectodomain_n_glycosylation_sites.tsv"
)

GLYCAN_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_candidate_glycan_proximity.tsv"
)

OUTPUT_FILE = (
    ROOT
    / "results"
    / "tables"
    / "structure"
    / "nkg2a_epitope_candidate_integration.tsv"
)


# =============================================================================
# SETTINGS
# =============================================================================

EXPOSED_RSA_THRESHOLD = 0.25

# Residues at an intermolecular interface are not automatically "bad".
# However, for a first-pass specificity epitope screen, exposed residues away
# from CD94/HLA-E/peptide interfaces are especially attractive because their
# accessibility is less dependent on receptor/ligand engagement.

INTERFACE_COLUMNS = [
    "contact_CD94",
    "contact_HLA_E",
    "contact_B2M",
    "contact_peptide",
]


# =============================================================================
# GENERIC HELPERS
# =============================================================================

def read_tsv(path):
    """
    Read a TSV file into a list of dictionaries.
    """

    if not path.exists():
        raise FileNotFoundError(
            f"Required input file not found:\n{path}"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        return list(reader)


def require_columns(rows, required, source_name):
    """
    Confirm that a table contains required columns.
    """

    if not rows:
        raise ValueError(
            f"{source_name} contains no data rows."
        )

    observed = set(rows[0].keys())

    missing = set(required) - observed

    if missing:
        raise ValueError(
            f"{source_name} is missing required columns:\n"
            + "\n".join(
                f"  {column}"
                for column in sorted(missing)
            )
        )


def clean(value):
    """
    Normalize a table value to stripped text.
    """

    if value is None:
        return ""

    return str(value).strip()


def yes(value):
    """
    Interpret a yes/no field.
    """

    return clean(value).lower() == "yes"


def integer_or_none(value):
    """
    Convert a numeric text field to int.
    """

    text = clean(value)

    if not text:
        return None

    return int(float(text))


def float_or_none(value):
    """
    Convert a numeric text field to float.
    """

    text = clean(value)

    if not text:
        return None

    return float(text)


def format_float(value, decimals=4):
    """
    Format a floating-point value.
    """

    if value is None:
        return ""

    return f"{value:.{decimals}f}"


# =============================================================================
# LOAD SEQUENCE CANDIDATES
# =============================================================================

def load_sequence_candidates():
    """
    Load the Step 1 discriminatory candidate table.

    The alignment is NKG2A-centered, so positions lacking a human NKG2A
    residue (for example an NKG2C insertion relative to NKG2A) cannot be
    mapped to a human NKG2A structural residue and are retained separately
    from residue-based structural candidates.
    """

    rows = read_tsv(
        SEQUENCE_FILE
    )

    required = [
        "alignment_column",
        "human_NKG2A_residue",
        "human_NKG2A_aa",
        "human_NKG2C_residue",
        "human_NKG2C_aa",
        "human_A_vs_C_diff",
        "rhesus_NKG2A_residue",
        "rhesus_NKG2A_aa",
        "rhesus_NKG2C1_residue",
        "rhesus_NKG2C1_aa",
        "rhesus_A_vs_C1_diff",
        "rhesus_NKG2C2_residue",
        "rhesus_NKG2C2_aa",
        "rhesus_A_vs_C2_diff",
        "pigtail_NKG2A_residue",
        "pigtail_NKG2A_aa",
        "pigtail_NKG2C_residue",
        "pigtail_NKG2C_aa",
        "pigtail_A_vs_C_diff",
        "classification",
        "sequence_priority",
    ]

    require_columns(
        rows,
        required,
        "candidate_specificity_positions.tsv",
    )

    candidates = {}
    insertion_only = []

    for row in rows:

        human_residue = integer_or_none(
            row["human_NKG2A_residue"]
        )

        if human_residue is None:
            insertion_only.append(row)
            continue

        if human_residue in candidates:
            raise ValueError(
                "Duplicate human NKG2A candidate residue "
                f"{human_residue} in sequence table."
            )

        candidates[human_residue] = row

    return candidates, insertion_only


# =============================================================================
# LOAD STRUCTURAL TABLES
# =============================================================================

def load_contacts():
    """
    Index Step 2F contact rows by human NKG2A full-length residue.
    """

    rows = read_tsv(
        CONTACT_FILE
    )

    required = [
        "full_length_residue",
        "nkg2a_aa",
        "candidate_specificity_position",
        "contact_CD94",
        "contact_HLA_E",
        "contact_B2M",
        "contact_peptide",
    ]

    require_columns(
        rows,
        required,
        "nkg2a_structural_contacts.tsv",
    )

    result = {}

    for row in rows:

        residue = integer_or_none(
            row["full_length_residue"]
        )

        if residue is None:
            continue

        result[residue] = row

    return result


def load_sasa():
    """
    Index Step 2G solvent-accessibility rows by full-length residue.
    """

    rows = read_tsv(
        SASA_FILE
    )

    required = [
        "full_length_residue",
        "complex_rsa",
        "complex_exposure_class",
        "isolated_nkg2a_rsa",
        "buried_surface_area_A2",
        "any_interface_contact",
    ]

    require_columns(
        rows,
        required,
        "nkg2a_solvent_accessibility.tsv",
    )

    result = {}

    for row in rows:

        residue = integer_or_none(
            row["full_length_residue"]
        )

        if residue is None:
            continue

        result[residue] = row

    return result


def load_glycan_status():
    """
    Index Step 2J candidate rows by full-length residue.

    If 3CDG contains no modeled glycans, blank glycan distances are treated
    as "not_assessable_no_modeled_glycans", not as evidence of >5 A distance.
    """

    rows = read_tsv(
        GLYCAN_FILE
    )

    required = [
        "full_length_residue",
        "coordinate_present",
        "nearest_glycan_chain",
        "nearest_glycan_resname",
        "nearest_glycan_residue",
        "minimum_glycan_distance_A",
        "within_5A_of_modeled_glycan",
    ]

    require_columns(
        rows,
        required,
        "nkg2a_candidate_glycan_proximity.tsv",
    )

    result = {}

    for row in rows:

        residue = integer_or_none(
            row["full_length_residue"]
        )

        if residue is None:
            continue

        result[residue] = row

    return result


def load_human_sequons():
    """
    Load canonical human NKG2A N-X-S/T sequon starts from Step 2I.
    """

    rows = read_tsv(
        GLYCOSYLATION_FILE
    )

    required = [
        "record_id",
        "sequon_full_length_residue",
        "motif",
    ]

    require_columns(
        rows,
        required,
        "nkg2_ectodomain_n_glycosylation_sites.tsv",
    )

    result = {}

    for row in rows:

        if clean(row["record_id"]) != "human_NKG2A":
            continue

        residue = integer_or_none(
            row["sequon_full_length_residue"]
        )

        if residue is None:
            continue

        result[residue] = clean(
            row["motif"]
        )

    return result


# =============================================================================
# SPECIFICITY LOGIC
# =============================================================================

def evaluate_species_specificity(row):
    """
    Evaluate NKG2A-vs-NKG2C discrimination separately in each species.

    For rhesus, both current NKG2C isoforms are considered.
    """

    human_diff = yes(
        row["human_A_vs_C_diff"]
    )

    rhesus_c1_diff = yes(
        row["rhesus_A_vs_C1_diff"]
    )

    rhesus_c2_diff = yes(
        row["rhesus_A_vs_C2_diff"]
    )

    pigtail_diff = yes(
        row["pigtail_A_vs_C_diff"]
    )

    rhesus_both_diff = (
        rhesus_c1_diff
        and rhesus_c2_diff
    )

    all_species_diff = (
        human_diff
        and rhesus_both_diff
        and pigtail_diff
    )

    macaque_diff = (
        rhesus_both_diff
        and pigtail_diff
    )

    return {
        "human_diff":
            human_diff,
        "rhesus_c1_diff":
            rhesus_c1_diff,
        "rhesus_c2_diff":
            rhesus_c2_diff,
        "rhesus_both_diff":
            rhesus_both_diff,
        "pigtail_diff":
            pigtail_diff,
        "all_species_diff":
            all_species_diff,
        "macaque_diff":
            macaque_diff,
    }


def evaluate_nkg2a_state_conservation(row):
    """
    Determine whether the NKG2A amino-acid state itself is conserved across
    human, rhesus, and pig-tailed macaque.

    This is distinct from asking whether NKG2A differs from NKG2C.
    """

    human_aa = clean(
        row["human_NKG2A_aa"]
    )

    rhesus_aa = clean(
        row["rhesus_NKG2A_aa"]
    )

    pigtail_aa = clean(
        row["pigtail_NKG2A_aa"]
    )

    all_present = all(
        aa not in {"", "-"}
        for aa in [
            human_aa,
            rhesus_aa,
            pigtail_aa,
        ]
    )

    pan_species_same = (
        all_present
        and human_aa == rhesus_aa
        and human_aa == pigtail_aa
    )

    macaque_same = (
        rhesus_aa not in {"", "-"}
        and pigtail_aa not in {"", "-"}
        and rhesus_aa == pigtail_aa
    )

    return {
        "pan_species_same":
            pan_species_same,
        "macaque_same":
            macaque_same,
    }


# =============================================================================
# STRUCTURAL LOGIC
# =============================================================================

def evaluate_structure(
    residue,
    contact_row,
    sasa_row,
):
    """
    Summarize structural evidence for one human NKG2A residue.
    """

    coordinate_present = (
        contact_row is not None
        and clean(
            contact_row.get(
                "full_length_residue",
                "",
            )
        )
        != ""
    )

    contacts = {
        "CD94": False,
        "HLA_E": False,
        "B2M": False,
        "peptide": False,
    }

    if contact_row is not None:

        contacts["CD94"] = yes(
            contact_row.get(
                "contact_CD94",
                "",
            )
        )

        contacts["HLA_E"] = yes(
            contact_row.get(
                "contact_HLA_E",
                "",
            )
        )

        contacts["B2M"] = yes(
            contact_row.get(
                "contact_B2M",
                "",
            )
        )

        contacts["peptide"] = yes(
            contact_row.get(
                "contact_peptide",
                "",
            )
        )

    any_interface = any(
        contacts.values()
    )

    complex_rsa = None
    isolated_rsa = None
    exposure_class = ""
    buried_surface_area = None

    if sasa_row is not None:

        complex_rsa = float_or_none(
            sasa_row.get(
                "complex_rsa",
                "",
            )
        )

        isolated_rsa = float_or_none(
            sasa_row.get(
                "isolated_nkg2a_rsa",
                "",
            )
        )

        exposure_class = clean(
            sasa_row.get(
                "complex_exposure_class",
                "",
            )
        )

        buried_surface_area = (
            float_or_none(
                sasa_row.get(
                    "buried_surface_area_A2",
                    "",
                )
            )
        )

    exposed = (
        complex_rsa is not None
        and complex_rsa
        >= EXPOSED_RSA_THRESHOLD
    )

    exposed_non_interface = (
        exposed
        and not any_interface
    )

    return {
        "coordinate_present":
            coordinate_present,
        "contact_CD94":
            contacts["CD94"],
        "contact_HLA_E":
            contacts["HLA_E"],
        "contact_B2M":
            contacts["B2M"],
        "contact_peptide":
            contacts["peptide"],
        "any_interface":
            any_interface,
        "complex_rsa":
            complex_rsa,
        "isolated_rsa":
            isolated_rsa,
        "exposure_class":
            exposure_class,
        "buried_surface_area":
            buried_surface_area,
        "exposed":
            exposed,
        "exposed_non_interface":
            exposed_non_interface,
    }


# =============================================================================
# PRIORITIZATION
# =============================================================================

def assign_priority(
    specificity,
    conservation,
    structure,
):
    """
    Assign a transparent qualitative priority.

    This is deliberately rule-based rather than a numerical score.

    Highest priority:
      - NKG2A differs from NKG2C in human, both rhesus isoform comparisons,
        and pig-tailed macaque
      - the NKG2A amino-acid state itself is conserved across species
      - residue is exposed
      - residue is outside the measured CD94/HLA-E/B2M/peptide interfaces

    Other categories preserve useful candidates without implying that they
    meet the full pan-species antibody-specificity objective.
    """

    if not structure["coordinate_present"]:
        return (
            "unresolved_structurally",
            "candidate lacks mapped coordinates in 3CDG",
        )

    if not structure["exposed"]:
        return (
            "low_structural_priority",
            "NKG2A-vs-NKG2C candidate is not exposed in 3CDG",
        )

    if (
        specificity["all_species_diff"]
        and conservation["pan_species_same"]
        and structure["exposed_non_interface"]
    ):
        return (
            "highest_pan_species_candidate",
            "conserved NKG2A state discriminates NKG2C across human, "
            "rhesus, and pigtail and is exposed/non-interface",
        )

    if (
        specificity["all_species_diff"]
        and structure["exposed_non_interface"]
    ):
        return (
            "high_pan_species_discrimination",
            "NKG2A differs from NKG2C across all species and residue is "
            "exposed/non-interface, but NKG2A amino-acid state varies",
        )

    if (
        specificity["macaque_diff"]
        and conservation["macaque_same"]
        and structure["exposed_non_interface"]
    ):
        return (
            "high_macaque_candidate",
            "shared macaque NKG2A state discriminates rhesus and pigtail "
            "NKG2C and is exposed/non-interface",
        )

    if (
        specificity["pigtail_diff"]
        and structure["exposed_non_interface"]
    ):
        return (
            "pigtail_specific_candidate",
            "pigtail NKG2A differs from pigtail NKG2C and is "
            "exposed/non-interface",
        )

    if (
        specificity["human_diff"]
        and structure["exposed_non_interface"]
    ):
        return (
            "human_specific_candidate",
            "human NKG2A differs from human NKG2C and is "
            "exposed/non-interface",
        )

    if (
        structure["exposed"]
        and structure["any_interface"]
    ):
        return (
            "interface_candidate",
            "NKG2A-vs-NKG2C candidate is exposed but participates in a "
            "measured receptor/ligand interface",
        )

    return (
        "secondary_candidate",
        "candidate retains sequence discrimination but does not satisfy "
        "a higher-priority combined rule",
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "STEP 2K - INTEGRATE NKG2A-SPECIFIC EPITOPE CANDIDATES"
    )
    print("=" * 78)

    sequence_candidates, insertion_only = (
        load_sequence_candidates()
    )

    contacts = load_contacts()
    sasa = load_sasa()
    glycan_status = load_glycan_status()
    human_sequons = load_human_sequons()

    print(
        f"\nSequence candidate positions with human "
        f"NKG2A residues: {len(sequence_candidates)}"
    )

    print(
        f"NKG2C-insertion candidate positions without "
        f"human NKG2A residue: {len(insertion_only)}"
    )

    print(
        f"Structural contact residues loaded: "
        f"{len(contacts)}"
    )

    print(
        f"SASA residues loaded: "
        f"{len(sasa)}"
    )

    print(
        f"Human NKG2A canonical sequons: "
        f"{len(human_sequons)}"
    )

    output_rows = []

    priority_counts = {}

    for residue in sorted(
        sequence_candidates
    ):

        seq = sequence_candidates[
            residue
        ]

        specificity = (
            evaluate_species_specificity(
                seq
            )
        )

        conservation = (
            evaluate_nkg2a_state_conservation(
                seq
            )
        )

        structure = (
            evaluate_structure(
                residue,
                contacts.get(residue),
                sasa.get(residue),
            )
        )

        glycan = glycan_status.get(
            residue
        )

        direct_sequon = (
            human_sequons.get(
                residue,
                ""
            )
        )

        if glycan is None:

            modeled_glycan_status = (
                "not_evaluated"
            )

            nearest_glycan_distance = ""

        else:

            nearest_glycan_distance = clean(
                glycan.get(
                    "minimum_glycan_distance_A",
                    "",
                )
            )

            if not nearest_glycan_distance:
                modeled_glycan_status = (
                    "not_assessable_no_modeled_glycans"
                )

            elif yes(
                glycan.get(
                    "within_5A_of_modeled_glycan",
                    "",
                )
            ):
                modeled_glycan_status = (
                    "within_5A_of_modeled_glycan"
                )

            else:
                modeled_glycan_status = (
                    "not_within_5A_of_modeled_glycan"
                )

        priority, rationale = (
            assign_priority(
                specificity,
                conservation,
                structure,
            )
        )

        priority_counts[priority] = (
            priority_counts.get(
                priority,
                0,
            )
            + 1
        )

        row = {
            # -------------------------------------------------------------
            # Position
            # -------------------------------------------------------------
            "alignment_column":
                clean(
                    seq["alignment_column"]
                ),

            "human_NKG2A_residue":
                residue,

            "human_NKG2A_aa":
                clean(
                    seq["human_NKG2A_aa"]
                ),

            # -------------------------------------------------------------
            # Human NKG2C comparison
            # -------------------------------------------------------------
            "human_NKG2C_residue":
                clean(
                    seq["human_NKG2C_residue"]
                ),

            "human_NKG2C_aa":
                clean(
                    seq["human_NKG2C_aa"]
                ),

            "human_NKG2A_vs_NKG2C_diff":
                (
                    "yes"
                    if specificity[
                        "human_diff"
                    ]
                    else "no"
                ),

            # -------------------------------------------------------------
            # Rhesus comparison
            # -------------------------------------------------------------
            "rhesus_NKG2A_residue":
                clean(
                    seq["rhesus_NKG2A_residue"]
                ),

            "rhesus_NKG2A_aa":
                clean(
                    seq["rhesus_NKG2A_aa"]
                ),

            "rhesus_NKG2C1_residue":
                clean(
                    seq["rhesus_NKG2C1_residue"]
                ),

            "rhesus_NKG2C1_aa":
                clean(
                    seq["rhesus_NKG2C1_aa"]
                ),

            "rhesus_NKG2C2_residue":
                clean(
                    seq["rhesus_NKG2C2_residue"]
                ),

            "rhesus_NKG2C2_aa":
                clean(
                    seq["rhesus_NKG2C2_aa"]
                ),

            "rhesus_NKG2A_vs_NKG2C1_diff":
                (
                    "yes"
                    if specificity[
                        "rhesus_c1_diff"
                    ]
                    else "no"
                ),

            "rhesus_NKG2A_vs_NKG2C2_diff":
                (
                    "yes"
                    if specificity[
                        "rhesus_c2_diff"
                    ]
                    else "no"
                ),

            "rhesus_both_NKG2C_isoforms_discriminated":
                (
                    "yes"
                    if specificity[
                        "rhesus_both_diff"
                    ]
                    else "no"
                ),

            # -------------------------------------------------------------
            # Pig-tailed macaque comparison
            # -------------------------------------------------------------
            "pigtail_NKG2A_residue":
                clean(
                    seq["pigtail_NKG2A_residue"]
                ),

            "pigtail_NKG2A_aa":
                clean(
                    seq["pigtail_NKG2A_aa"]
                ),

            "pigtail_NKG2C_residue":
                clean(
                    seq["pigtail_NKG2C_residue"]
                ),

            "pigtail_NKG2C_aa":
                clean(
                    seq["pigtail_NKG2C_aa"]
                ),

            "pigtail_NKG2A_vs_NKG2C_diff":
                (
                    "yes"
                    if specificity[
                        "pigtail_diff"
                    ]
                    else "no"
                ),

            # -------------------------------------------------------------
            # Cross-species interpretation
            # -------------------------------------------------------------
            "pan_species_NKG2A_vs_NKG2C_discrimination":
                (
                    "yes"
                    if specificity[
                        "all_species_diff"
                    ]
                    else "no"
                ),

            "macaque_NKG2A_vs_NKG2C_discrimination":
                (
                    "yes"
                    if specificity[
                        "macaque_diff"
                    ]
                    else "no"
                ),

            "NKG2A_state_conserved_human_rhesus_pigtail":
                (
                    "yes"
                    if conservation[
                        "pan_species_same"
                    ]
                    else "no"
                ),

            "NKG2A_state_conserved_rhesus_pigtail":
                (
                    "yes"
                    if conservation[
                        "macaque_same"
                    ]
                    else "no"
                ),

            "original_classification":
                clean(
                    seq["classification"]
                ),

            "original_sequence_priority":
                clean(
                    seq["sequence_priority"]
                ),

            # -------------------------------------------------------------
            # Structure
            # -------------------------------------------------------------
            "coordinate_present":
                (
                    "yes"
                    if structure[
                        "coordinate_present"
                    ]
                    else "no"
                ),

            "complex_rsa":
                format_float(
                    structure[
                        "complex_rsa"
                    ]
                ),

            "complex_exposure_class":
                structure[
                    "exposure_class"
                ],

            "isolated_nkg2a_rsa":
                format_float(
                    structure[
                        "isolated_rsa"
                    ]
                ),

            "buried_surface_area_A2":
                format_float(
                    structure[
                        "buried_surface_area"
                    ],
                    decimals=3,
                ),

            "contact_CD94":
                (
                    "yes"
                    if structure[
                        "contact_CD94"
                    ]
                    else "no"
                ),

            "contact_HLA_E":
                (
                    "yes"
                    if structure[
                        "contact_HLA_E"
                    ]
                    else "no"
                ),

            "contact_B2M":
                (
                    "yes"
                    if structure[
                        "contact_B2M"
                    ]
                    else "no"
                ),

            "contact_peptide":
                (
                    "yes"
                    if structure[
                        "contact_peptide"
                    ]
                    else "no"
                ),

            "any_interface_contact":
                (
                    "yes"
                    if structure[
                        "any_interface"
                    ]
                    else "no"
                ),

            "exposed_non_interface":
                (
                    "yes"
                    if structure[
                        "exposed_non_interface"
                    ]
                    else "no"
                ),

            # -------------------------------------------------------------
            # Glycosylation
            # -------------------------------------------------------------
            "human_NKG2A_direct_NXS_T_sequon":
                (
                    direct_sequon
                    if direct_sequon
                    else "no"
                ),

            "modeled_glycan_status":
                modeled_glycan_status,

            "nearest_modeled_glycan_distance_A":
                nearest_glycan_distance,

            # -------------------------------------------------------------
            # Integrated interpretation
            # -------------------------------------------------------------
            "integrated_priority":
                priority,

            "priority_rationale":
                rationale,
        }

        output_rows.append(
            row
        )

    # =========================================================================
    # WRITE OUTPUT
    # =========================================================================

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fields = [
        "alignment_column",
        "human_NKG2A_residue",
        "human_NKG2A_aa",

        "human_NKG2C_residue",
        "human_NKG2C_aa",
        "human_NKG2A_vs_NKG2C_diff",

        "rhesus_NKG2A_residue",
        "rhesus_NKG2A_aa",
        "rhesus_NKG2C1_residue",
        "rhesus_NKG2C1_aa",
        "rhesus_NKG2C2_residue",
        "rhesus_NKG2C2_aa",
        "rhesus_NKG2A_vs_NKG2C1_diff",
        "rhesus_NKG2A_vs_NKG2C2_diff",
        "rhesus_both_NKG2C_isoforms_discriminated",

        "pigtail_NKG2A_residue",
        "pigtail_NKG2A_aa",
        "pigtail_NKG2C_residue",
        "pigtail_NKG2C_aa",
        "pigtail_NKG2A_vs_NKG2C_diff",

        "pan_species_NKG2A_vs_NKG2C_discrimination",
        "macaque_NKG2A_vs_NKG2C_discrimination",
        "NKG2A_state_conserved_human_rhesus_pigtail",
        "NKG2A_state_conserved_rhesus_pigtail",

        "original_classification",
        "original_sequence_priority",

        "coordinate_present",
        "complex_rsa",
        "complex_exposure_class",
        "isolated_nkg2a_rsa",
        "buried_surface_area_A2",

        "contact_CD94",
        "contact_HLA_E",
        "contact_B2M",
        "contact_peptide",
        "any_interface_contact",
        "exposed_non_interface",

        "human_NKG2A_direct_NXS_T_sequon",
        "modeled_glycan_status",
        "nearest_modeled_glycan_distance_A",

        "integrated_priority",
        "priority_rationale",
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
        writer.writerows(
            output_rows
        )

    # =========================================================================
    # REPORT
    # =========================================================================

    print("\n" + "=" * 78)
    print("INTEGRATED PRIORITY COUNTS")
    print("=" * 78)

    priority_order = [
        "highest_pan_species_candidate",
        "high_pan_species_discrimination",
        "high_macaque_candidate",
        "pigtail_specific_candidate",
        "human_specific_candidate",
        "interface_candidate",
        "low_structural_priority",
        "secondary_candidate",
        "unresolved_structurally",
    ]

    for priority in priority_order:

        count = priority_counts.get(
            priority,
            0,
        )

        if count:
            print(
                f"{priority:<38} "
                f"{count}"
            )

    print("\n" + "=" * 78)
    print("TOP EXPOSED NON-INTERFACE CANDIDATES")
    print("=" * 78)

    rank_order = {
        "highest_pan_species_candidate": 1,
        "high_pan_species_discrimination": 2,
        "high_macaque_candidate": 3,
        "pigtail_specific_candidate": 4,
        "human_specific_candidate": 5,
        "secondary_candidate": 6,
        "interface_candidate": 7,
        "low_structural_priority": 8,
        "unresolved_structurally": 9,
    }

    top_rows = [
        row
        for row in output_rows
        if row["exposed_non_interface"] == "yes"
    ]

    top_rows.sort(
        key=lambda row: (
            rank_order.get(
                row["integrated_priority"],
                99,
            ),
            -(
                float(row["complex_rsa"])
                if row["complex_rsa"]
                else -1.0
            ),
            int(
                row["human_NKG2A_residue"]
            ),
        )
    )

    for row in top_rows:

        print(
            f"{int(row['human_NKG2A_residue']):>3} "
            f"{row['human_NKG2A_aa']}  "
            f"{row['integrated_priority']:<34} "
            f"RSA={row['complex_rsa']:<7} "
            f"H:{row['human_NKG2A_vs_NKG2C_diff']} "
            f"R1:{row['rhesus_NKG2A_vs_NKG2C1_diff']} "
            f"R2:{row['rhesus_NKG2A_vs_NKG2C2_diff']} "
            f"P:{row['pigtail_NKG2A_vs_NKG2C_diff']}"
        )

    print("\n" + "=" * 78)
    print("NKG2C INSERTION-ONLY CANDIDATE POSITIONS")
    print("=" * 78)

    if not insertion_only:

        print("\nNone.")

    else:

        for row in insertion_only:

            print(
                f"Alignment column "
                f"{clean(row['alignment_column'])}: "
                f"human NKG2A={clean(row['human_NKG2A_aa'])} "
                f"human NKG2C={clean(row['human_NKG2C_aa'])}  "
                f"{clean(row['classification'])}"
            )

    print("\n" + "=" * 78)
    print("OUTPUT")
    print("=" * 78)

    print(f"\n{OUTPUT_FILE}")

    print(
        "\nNOTE: integrated_priority is a transparent screening "
        "classification, not an experimentally validated antibody epitope."
    )

    print(
        "Residues are prioritized for NKG2A-vs-NKG2C discrimination "
        "first, then cross-species conservation and structural accessibility."
    )

    print(
        "Interface residues are retained rather than discarded because "
        "their value depends on the intended antibody mechanism."
    )


if __name__ == "__main__":
    main()
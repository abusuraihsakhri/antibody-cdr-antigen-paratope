#!/usr/bin/env python3
"""
Antibody CDR-H3 Loop Modeling, Paratope-Epitope Contact & Developability Engine.

Domain: Computational Structural Biology & Immunoinformatics
Standards: IMGT / Kabat / Chothia Schemes, Clackson-Wells Hotspot Thermodynamics,
           Kyte-Doolittle Hydropathy, Henderson-Hasselbalch pI / Charge Analytics.

Pure Python standard library implementation.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Any, Set


# =============================================================================
# CONSTANTS & SCALES
# =============================================================================

# Empirical binding free energy constant: ~22 cal/(mol * A^2) = 0.022 kcal/(mol * A^2)
CAL_PER_MOL_PER_A2 = 22.0
HOTSPOT_DDG_THRESHOLD = 1.5  # kcal/mol (Clackson & Wells benchmark)

# Canonical loop boundaries (Kabat numbering conventions)
CDR_RANGES_KABAT: Dict[str, Tuple[int, int]] = {
    "CDR-H1": (31, 35),
    "CDR-H2": (50, 65),
    "CDR-H3": (95, 102),
    "CDR-L1": (24, 34),
    "CDR-L2": (50, 56),
    "CDR-L3": (89, 97),
}

# Canonical loop boundaries (IMGT numbering conventions)
CDR_RANGES_IMGT: Dict[str, Tuple[int, int]] = {
    "CDR-H1": (27, 38),
    "CDR-H2": (56, 65),
    "CDR-H3": (105, 117),
    "CDR-L1": (27, 38),
    "CDR-L2": (56, 65),
    "CDR-L3": (105, 117),
}

# Kyte-Doolittle hydropathy scale (Kyte & Doolittle, J. Mol. Biol. 1982)
KYTE_DOOLITTLE: Dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.4,
}

# Standard Bjellqvist / EMBOSS pKa values for pI and charge calculations
PKA_VALUES: Dict[str, float] = {
    "N_TERM": 8.6,
    "C_TERM": 3.6,
    "K": 10.8,
    "R": 12.5,
    "H": 6.5,
    "D": 3.9,
    "E": 4.1,
    "C": 8.5,
    "Y": 10.1,
}

# Human germline variable domain consensus references for humanness scoring
# VH: IGHV3-23*01 / IGHJ4 consensus
HUMAN_VH_GERMLINE_CONSENSUS = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKG"
    "RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDYSGMDYWGQGTLVTVSS"
)

# VL: IGKV1-39*01 / IGKJ1 consensus
HUMAN_VL_GERMLINE_CONSENSUS = (
    "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQS"
    "GVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPYTFGQGTKVEIK"
)


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass(frozen=True)
class HotSpotResult:
    chain: str
    residue_id: int
    residue_name: str
    bsa_a2: float
    ddg_kcal_mol: float
    classification: str  # "Critical Hotspot", "Hotspot", "Moderate", "Neutral"


@dataclass(frozen=True)
class LoopDecomposition:
    loop: str
    total_bsa_a2: float
    percentage_of_paratope: float
    hotspots_count: int
    critical_hotspots_count: int


@dataclass(frozen=True)
class InterfaceContactSummary:
    hydrogen_bonds: int
    salt_bridges: int
    pi_stacking: int
    hydrophobic_contacts: int
    total_interactions: int
    estimated_interaction_energy_kcal: float


@dataclass(frozen=True)
class ChemicalLiability:
    liability_type: str  # "Deamidation", "Isomerization", "Oxidation", "Glycosylation", "Unpaired Cys", "Acid Labile"
    location: str        # e.g. "CDR-H3" or "VH:45"
    motif: str
    start_pos: int
    risk_level: str      # "High", "Medium", "Low"


@dataclass(frozen=True)
class CDRH3StructureProfile:
    sequence: str
    length: int
    length_category: str  # "Short", "Medium", "Long", "Ultra-long"
    base_conformation: str  # "Kinked" vs "Extended"
    torso_salt_bridge: bool
    aromatic_fraction: float
    gly_ser_fraction: float
    net_charge_ph74: float


@dataclass(frozen=True)
class DevelopabilityReport:
    antibody_id: str
    isoelectric_point: float
    net_charge_ph74: float
    net_charge_ph60: float
    hydrophobic_patches: Dict[str, List[Tuple[int, float]]]
    aggregation_risk: str  # "Low", "Moderate", "High"
    chemical_liabilities: List[ChemicalLiability]
    liabilities_count: int
    framework_humanness_vh: float
    framework_humanness_vl: float
    overall_humanness: float
    immunogenicity_tier: str  # "Low Risk", "Moderate Risk (Screening Recommended)", "High Risk (Humanization Required)"
    estimated_tm_celsius: float
    composite_developability_score: float  # 0-100
    developability_tier: str  # "High Developability", "Moderate Developability", "High Liability Risk"


# =============================================================================
# DOMAIN FUNCTIONS & ALGORITHMS
# =============================================================================

def calculate_ddg_from_bsa(bsa_a2: float) -> float:
    """
    Calculate estimated Delta-Delta-G of binding from Buried Surface Area (BSA) in A^2.
    Empirical approximation: ~22 cal/mol/A^2 = 0.022 kcal/(mol * A^2).
    """
    if bsa_a2 < 0:
        raise ValueError("BSA cannot be negative")
    return round(CAL_PER_MOL_PER_A2 * bsa_a2 / 1000.0, 3)


def classify_hotspot(ddg_kcal_mol: float) -> str:
    """Classify residue energetic contribution based on ddG."""
    if ddg_kcal_mol >= 2.5:
        return "Critical Hotspot"
    elif ddg_kcal_mol >= HOTSPOT_DDG_THRESHOLD:
        return "Hotspot"
    elif ddg_kcal_mol >= 0.5:
        return "Moderate"
    else:
        return "Neutral"


def perform_alanine_scan(
    residue_bsa: Dict[Tuple[str, int, str], float]
) -> List[HotSpotResult]:
    """
    Perform computational alanine scanning analysis from residue buried surface areas.
    residue_bsa: Dict mapping (chain, residue_number, residue_1letter) -> BSA in A^2.
    """
    if not residue_bsa:
        raise ValueError("No residue BSA data provided for alanine scan")

    results: List[HotSpotResult] = []
    for (chain, rid, rname), bsa in residue_bsa.items():
        if bsa < 0:
            raise ValueError(f"Negative BSA for residue {chain}:{rid}{rname}")
        ddg = calculate_ddg_from_bsa(bsa)
        classification = classify_hotspot(ddg)
        results.append(HotSpotResult(
            chain=chain.upper(),
            residue_id=rid,
            residue_name=rname.upper(),
            bsa_a2=round(bsa, 2),
            ddg_kcal_mol=ddg,
            classification=classification,
        ))

    # Sort descending by ddG
    results.sort(key=lambda x: x.ddg_kcal_mol, reverse=True)
    return results


def decompose_cdr_loops(
    residue_bsa: Dict[Tuple[str, int, str], float],
    numbering_scheme: str = "IMGT",
) -> List[LoopDecomposition]:
    """
    Decompose paratope buried surface area into CDR loop contributions.
    """
    ranges = CDR_RANGES_IMGT if numbering_scheme.upper() == "IMGT" else CDR_RANGES_KABAT
    totals: Dict[str, float] = {k: 0.0 for k in ranges}
    hotspots: Dict[str, int] = {k: 0 for k in ranges}
    crit_hotspots: Dict[str, int] = {k: 0 for k in ranges}

    total_paratope_bsa = 0.0

    for (chain, rid, _), bsa in residue_bsa.items():
        chain_upper = chain.upper()
        prefix = "CDR-H" if chain_upper == "H" else "CDR-L"
        total_paratope_bsa += bsa
        ddg = calculate_ddg_from_bsa(bsa)

        for loop_name, (start, end) in ranges.items():
            if loop_name.startswith(prefix) and start <= rid <= end:
                totals[loop_name] += bsa
                if ddg >= 2.5:
                    crit_hotspots[loop_name] += 1
                    hotspots[loop_name] += 1
                elif ddg >= HOTSPOT_DDG_THRESHOLD:
                    hotspots[loop_name] += 1
                break

    decompositions = []
    for loop_name in ranges:
        bsa_val = totals[loop_name]
        pct = (bsa_val / total_paratope_bsa * 100.0) if total_paratope_bsa > 0 else 0.0
        decompositions.append(LoopDecomposition(
            loop=loop_name,
            total_bsa_a2=round(bsa_val, 2),
            percentage_of_paratope=round(pct, 2),
            hotspots_count=hotspots[loop_name],
            critical_hotspots_count=crit_hotspots[loop_name],
        ))

    decompositions.sort(key=lambda x: x.total_bsa_a2, reverse=True)
    return decompositions


def analyze_interface_contacts(
    hbonds: Sequence[Tuple[str, str]],
    salt_bridges: Sequence[Tuple[str, str]],
    pi_stacking: Sequence[Tuple[str, str]],
    hydrophobic_contacts: Sequence[Tuple[str, str]],
) -> InterfaceContactSummary:
    """
    Summarize interface contact types and estimate total interaction energy.
    Approximations:
    - H-bond: ~1.5 kcal/mol
    - Salt bridge: ~3.0 kcal/mol
    - Pi-stacking: ~2.0 kcal/mol
    - Hydrophobic contact: ~0.8 kcal/mol
    """
    n_hb = len(hbonds)
    n_sb = len(salt_bridges)
    n_pi = len(pi_stacking)
    n_hyd = len(hydrophobic_contacts)
    total = n_hb + n_sb + n_pi + n_hyd

    energy = (n_hb * 1.5) + (n_sb * 3.0) + (n_pi * 2.0) + (n_hyd * 0.8)

    return InterfaceContactSummary(
        hydrogen_bonds=n_hb,
        salt_bridges=n_sb,
        pi_stacking=n_pi,
        hydrophobic_contacts=n_hyd,
        total_interactions=total,
        estimated_interaction_energy_kcal=round(energy, 2),
    )


def profile_cdrh3_loop(
    cdrh3_sequence: str,
    vh_residue_94: str = "R",
    vh_residue_101: str = "D",
    vh_residue_103: str = "W",
) -> CDRH3StructureProfile:
    """
    Profile CDR-H3 loop characteristics and predict base/torso conformation.
    Kinked vs Extended base rules (Shirai & Morea rules):
    - Kinked base is formed when a salt bridge or hydrogen bond occurs between
      Arg/Lys at H94 and Asp/Glu at H101 (or H102), with Trp at H103 packing against it.
    - Extended base occurs when position 101 lacks acidic residue or position 94 lacks basic residue.
    """
    seq = cdrh3_sequence.strip().upper()
    if not seq or not all(c in KYTE_DOOLITTLE for c in seq):
        raise ValueError(f"Invalid CDR-H3 amino acid sequence: {cdrh3_sequence}")

    length = len(seq)
    if length <= 8:
        category = "Short"
    elif length <= 15:
        category = "Medium"
    elif length <= 21:
        category = "Long"
    else:
        category = "Ultra-long"

    r94 = vh_residue_94.upper()
    r101 = vh_residue_101.upper()
    r103 = vh_residue_103.upper()

    salt_bridge = (r94 in ("R", "K")) and (r101 in ("D", "E"))
    is_kinked = salt_bridge and (r103 in ("W", "F", "Y"))
    conformation = "Kinked" if is_kinked else "Extended"

    aromatics = sum(1 for c in seq if c in ("F", "Y", "W", "H"))
    gly_ser = sum(1 for c in seq if c in ("G", "S"))

    arom_frac = round(aromatics / length, 3)
    gs_frac = round(gly_ser / length, 3)

    # Net charge at pH 7.4
    net_charge = compute_net_charge(seq, ph=7.4)

    return CDRH3StructureProfile(
        sequence=seq,
        length=length,
        length_category=category,
        base_conformation=conformation,
        torso_salt_bridge=salt_bridge,
        aromatic_fraction=arom_frac,
        gly_ser_fraction=gs_frac,
        net_charge_ph74=net_charge,
    )


def scan_hydrophobic_patches(
    sequence: str,
    window: int = 7,
    threshold: float = 1.0,
) -> List[Tuple[int, float]]:
    """
    Scan sequence for hydrophobic patches using Kyte-Doolittle sliding window.
    Returns list of (center_position_1_indexed, mean_hydropathy).
    """
    seq = sequence.strip().upper()
    if len(seq) < window:
        # If sequence shorter than window, compute full average
        mean_val = sum(KYTE_DOOLITTLE.get(c, 0.0) for c in seq) / len(seq)
        if mean_val >= threshold:
            return [(len(seq) // 2 + 1, round(mean_val, 3))]
        return []

    patches: List[Tuple[int, float]] = []
    half = window // 2
    for center in range(half, len(seq) - half):
        window_seq = seq[center - half : center + half + 1]
        mean_val = sum(KYTE_DOOLITTLE.get(c, 0.0) for c in window_seq) / window
        if mean_val >= threshold:
            patches.append((center + 1, round(mean_val, 3)))

    return patches


def find_chemical_liabilities(
    sequences: Dict[str, str]
) -> List[ChemicalLiability]:
    """
    Screen antibody variable domains and CDRs for sequence liabilities:
    - Deamidation: NG, NS, NN, NH (High in CDRs, Low/Medium in FR)
    - Isomerization: DG, DS, DD, DH, DT (High in CDRs, Low/Medium in FR)
    - Oxidation: M, W (High for M in CDRs, Medium for W in CDRs)
    - N-Glycosylation: N[^P][ST] (High across all regions)
    - Acid-labile cleavage: DP (Medium across all regions)
    - Unpaired Cysteine: C outside canonical disulfide (High)
    """
    liabilities: List[ChemicalLiability] = []

    for region, seq in sequences.items():
        seq_up = seq.strip().upper()
        is_cdr = "CDR" in region.upper()

        # Deamidation
        for m in re.finditer(r"N[GSNH]", seq_up):
            motif = m.group()
            if is_cdr:
                risk = "High" if motif in ("NG", "NS") else "Medium"
            else:
                risk = "Medium" if motif == "NG" else "Low"
            liabilities.append(ChemicalLiability(
                liability_type="Deamidation",
                location=region,
                motif=motif,
                start_pos=m.start() + 1,
                risk_level=risk,
            ))

        # Isomerization
        for m in re.finditer(r"D[GSDHT]", seq_up):
            motif = m.group()
            if is_cdr:
                risk = "High" if motif in ("DG", "DS") else "Medium"
            else:
                risk = "Medium" if motif == "DG" else "Low"
            liabilities.append(ChemicalLiability(
                liability_type="Isomerization",
                location=region,
                motif=motif,
                start_pos=m.start() + 1,
                risk_level=risk,
            ))

        # N-linked Glycosylation motif N-X-S/T where X != P
        for m in re.finditer(r"N[^P][ST]", seq_up):
            liabilities.append(ChemicalLiability(
                liability_type="Glycosylation",
                location=region,
                motif=m.group(),
                start_pos=m.start() + 1,
                risk_level="High",
            ))

        # Acid-labile cleavage DP
        for m in re.finditer(r"DP", seq_up):
            liabilities.append(ChemicalLiability(
                liability_type="Acid Labile",
                location=region,
                motif="DP",
                start_pos=m.start() + 1,
                risk_level="Medium",
            ))

        # Oxidation
        if is_cdr:
            for m in re.finditer(r"[MW]", seq_up):
                motif = m.group()
                risk = "High" if motif == "M" else "Medium"
                liabilities.append(ChemicalLiability(
                    liability_type="Oxidation",
                    location=region,
                    motif=motif,
                    start_pos=m.start() + 1,
                    risk_level=risk,
                ))
        else:
            for m in re.finditer(r"M", seq_up):
                liabilities.append(ChemicalLiability(
                    liability_type="Oxidation",
                    location=region,
                    motif="M",
                    start_pos=m.start() + 1,
                    risk_level="Low",
                ))

        # Unpaired Cysteines
        cys_count = seq_up.count("C")
        if cys_count % 2 != 0:
            for m in re.finditer(r"C", seq_up):
                liabilities.append(ChemicalLiability(
                    liability_type="Unpaired Cys",
                    location=region,
                    motif="C",
                    start_pos=m.start() + 1,
                    risk_level="High",
                ))

    return liabilities


def compute_net_charge(sequence: str, ph: float = 7.4) -> float:
    """
    Calculate protein net charge at a given pH using Henderson-Hasselbalch equation.
    """
    seq = sequence.strip().upper()
    if not seq:
        return 0.0

    # N-terminus positive charge
    charge = 1.0 / (1.0 + 10.0 ** (ph - PKA_VALUES["N_TERM"]))
    # C-terminus negative charge
    charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["C_TERM"] - ph))

    for aa in seq:
        if aa == "K":
            charge += 1.0 / (1.0 + 10.0 ** (ph - PKA_VALUES["K"]))
        elif aa == "R":
            charge += 1.0 / (1.0 + 10.0 ** (ph - PKA_VALUES["R"]))
        elif aa == "H":
            charge += 1.0 / (1.0 + 10.0 ** (ph - PKA_VALUES["H"]))
        elif aa == "D":
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["D"] - ph))
        elif aa == "E":
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["E"] - ph))
        elif aa == "C":
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["C"] - ph))
        elif aa == "Y":
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["Y"] - ph))

    return round(charge, 2)


def compute_isoelectric_point(sequence: str, precision: float = 0.01) -> float:
    """
    Compute isoelectric point (pI) using bisection method across pH range 0-14.
    """
    seq = sequence.strip().upper()
    if not seq:
        return 7.0

    low_ph = 0.0
    high_ph = 14.0

    while (high_ph - low_ph) > precision:
        mid_ph = (low_ph + high_ph) / 2.0
        charge = compute_net_charge(seq, mid_ph)
        if charge > 0:
            low_ph = mid_ph
        else:
            high_ph = mid_ph

    return round((low_ph + high_ph) / 2.0, 2)


def compute_sequence_identity(seq1: str, seq2: str) -> float:
    """Compute % sequence identity between two aligned or trimmed sequences."""
    s1 = seq1.strip().upper()
    s2 = seq2.strip().upper()
    span = min(len(s1), len(s2))
    if span == 0:
        return 0.0
    matches = sum(1 for a, b in zip(s1[:span], s2[:span]) if a == b)
    return round(100.0 * matches / span, 2)


def generate_developability_report(
    antibody_id: str,
    vh_sequence: str,
    vl_sequence: str,
    cdr_sequences: Optional[Dict[str, str]] = None,
    vh_germline_consensus: str = HUMAN_VH_GERMLINE_CONSENSUS,
    vl_germline_consensus: str = HUMAN_VL_GERMLINE_CONSENSUS,
) -> DevelopabilityReport:
    """
    Generate comprehensive developability analysis and risk profiling.
    """
    vh_clean = vh_sequence.strip().upper()
    vl_clean = vl_sequence.strip().upper()
    full_fv = vh_clean + vl_clean

    if not full_fv:
        raise ValueError("Empty VH and VL sequences provided")

    # Isoelectric point and charges
    pi = compute_isoelectric_point(full_fv)
    charge_74 = compute_net_charge(full_fv, 7.4)
    charge_60 = compute_net_charge(full_fv, 6.0)

    # Hydrophobic patch detection
    patch_dict: Dict[str, List[Tuple[int, float]]] = {}
    if cdr_sequences:
        for cdr_name, cdr_seq in cdr_sequences.items():
            cdr_p = scan_hydrophobic_patches(cdr_seq, window=5, threshold=1.0)
            if cdr_p:
                patch_dict[cdr_name] = cdr_p
    else:
        vh_patches = scan_hydrophobic_patches(vh_clean, window=7, threshold=1.8)
        if vh_patches:
            patch_dict["VH"] = vh_patches
        vl_patches = scan_hydrophobic_patches(vl_clean, window=7, threshold=1.8)
        if vl_patches:
            patch_dict["VL"] = vl_patches

    total_patch_hits = sum(len(p) for p in patch_dict.values())
    if total_patch_hits == 0:
        agg_risk = "Low"
    elif total_patch_hits <= 2:
        agg_risk = "Moderate"
    else:
        agg_risk = "High"

    # Chemical liabilities
    if cdr_sequences:
        seq_dict = dict(cdr_sequences)
    else:
        seq_dict = {"VH": vh_clean, "VL": vl_clean}
    liabilities = find_chemical_liabilities(seq_dict)

    # Humanness scoring against representative human germline variable domains
    hum_vh = compute_sequence_identity(vh_clean, vh_germline_consensus)
    hum_vl = compute_sequence_identity(vl_clean, vl_germline_consensus)
    overall_hum = round((hum_vh + hum_vl) / 2.0, 2)

    if overall_hum >= 85.0:
        immuno_tier = "Low Risk"
    elif overall_hum >= 70.0:
        immuno_tier = "Moderate Risk (Screening Recommended)"
    else:
        immuno_tier = "High Risk (Humanization Required)"

    # Estimated Tm
    # Baseline 75 C, penalized by hydrophobic patches, extreme net charge, high liability count
    tm_est = 75.0 - (total_patch_hits * 1.5) - (abs(charge_74) * 0.2) - (len(liabilities) * 0.3)
    tm_est = max(min(round(tm_est, 1), 85.0), 50.0)

    # Composite developability score (0-100)
    score = 100.0
    # Deductions
    score -= min(total_patch_hits * 6.0, 25.0)
    high_liabs = len([l for l in liabilities if l.risk_level == "High"])
    med_liabs = len([l for l in liabilities if l.risk_level == "Medium"])
    score -= min(high_liabs * 8.0, 30.0)
    score -= min(med_liabs * 3.0, 15.0)
    if immuno_tier.startswith("High"):
        score -= 15.0
    elif immuno_tier.startswith("Moderate"):
        score -= 5.0
    if abs(charge_74) > 10.0:
        score -= 10.0
    if pi < 5.5 or pi > 10.0:
        score -= 5.0

    score = max(min(round(score, 1), 100.0), 0.0)

    if score >= 75.0:
        dev_tier = "High Developability"
    elif score >= 50.0:
        dev_tier = "Moderate Developability"
    else:
        dev_tier = "High Liability Risk"

    return DevelopabilityReport(
        antibody_id=antibody_id,
        isoelectric_point=pi,
        net_charge_ph74=charge_74,
        net_charge_ph60=charge_60,
        hydrophobic_patches=patch_dict,
        aggregation_risk=agg_risk,
        chemical_liabilities=liabilities,
        liabilities_count=len(liabilities),
        framework_humanness_vh=hum_vh,
        framework_humanness_vl=hum_vl,
        overall_humanness=overall_hum,
        immunogenicity_tier=immuno_tier,
        estimated_tm_celsius=tm_est,
        composite_developability_score=score,
        developability_tier=dev_tier,
    )

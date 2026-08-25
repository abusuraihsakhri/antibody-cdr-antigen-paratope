# Antibody CDR-H3 Loop Modeling, Paratope Contacts & Developability

A production-grade Python immunoinformatics engine for paratope-epitope interface characterization, computational alanine scanning, CDR-H3 conformational modeling, and therapeutic antibody developability liability screening.

---

## Domain Overview & Scientific Foundation

Monoclonal antibodies (mAbs) and engineered antibody fragments (Fab, scFv, VHH) bind target antigens primarily via six Complementarity-Determining Region (CDR) loops. Among these, **CDR-H3** exhibits the greatest structural variability, loop length diversity, and energetic contribution to binding affinity.

This package provides a standalone, pure-standard-library computational framework implementing established biophysical and immunoinformatics algorithms:

1. **Computational Alanine Scanning & Hotspot Prediction**:
   - Binding free energy estimation ($\Delta\Delta G$) derived from per-residue Buried Surface Area ($\text{BSA}$ in $\text{Å}^2$) using the empirical surface area solvation parameter:
     $$\Delta\Delta G \approx 0.022 \times \text{BSA} \quad (\text{kcal/mol})$$
   - Hotspot identification using the classic **Clackson & Wells** thermodynamic threshold ($\Delta\Delta G \ge 1.5 \text{ kcal/mol}$).
   - Classification into *Critical Hotspot* ($\ge 2.5 \text{ kcal/mol}$), *Hotspot* ($1.5-2.5 \text{ kcal/mol}$), *Moderate* ($0.5-1.5 \text{ kcal/mol}$), and *Neutral* ($< 0.5 \text{ kcal/mol}$).

2. **CDR Loop Energetic Decomposition**:
   - Canonical loop partition under **IMGT** and **Kabat** schemes.
   - Paratope fractional contribution and CDR dominance rankings.

3. **CDR-H3 Torso & Base Conformation Modeling**:
   - **Shirai-Morea-Chothia rules**: Classification of CDR-H3 base into **Kinked** vs. **Extended** geometries based on the existence of the canonical H94 (Arg/Lys) to H101/H102 (Asp/Glu) torso salt bridge and H103 (Trp/Phe/Tyr) hydrophobic packing.
   - Loop length profiling and aromatic/flexibility fractions ($\text{Tyr}, \text{Gly}, \text{Ser}$).

4. **Antibody Developability & Liability Profiling**:
   - **Kyte-Doolittle Hydropathy Sliding-Window Scan** for hydrophobic aggregation patches.
   - **Chemical Sequence Liabilities**:
     - *Deamidation*: Asn motifs (`NG`, `NS`, `NN`, `NH`).
     - *Isomerization*: Asp motifs (`DG`, `DS`, `DD`, `DH`, `DT`).
     - *Oxidation*: Met/Trp residues in hypervariable loops.
     - *N-linked Glycosylation*: Sequon `N-X-[S/T]` ($X \neq P$).
     - *Acid-Labile Cleavage*: `DP` peptide bond susceptibility.
     - *Unpaired Cysteines*: Free thiol oxidation/scrambling risks.
   - **Physicochemical Properties**:
     - Net charge at physiological pH (7.4) and formulation pH (6.0) via Henderson-Hasselbalch exact equilibrium.
     - Isoelectric point ($pI$) computation via bisection search.
     - Thermal stability ($T_m$) heuristic approximation.
   - **Germline Humanness & Immunogenicity Risk**:
     - Sequence identity against human germline variable consensus (VH3-23/VK1-39).
     - Immunogenicity risk tiering (*Low Risk*, *Moderate Risk*, *High Risk*).

---

## Installation & Requirements

- Python 3.9+ (Pure standard library; no external dependencies required).

```bash
git clone https://github.com/abusuraihsakhri/antibody-cdr-antigen-paratope.git
cd antibody-cdr-antigen-paratope
```

---

## CLI Usage & Examples

### 1. Paratope Alanine Scanning
```bash
python cli.py scan
```
Output:
```
===========================================================================
  PARATOPE ALANINE SCANNING & HOT-SPOT PREDICTION
  Threshold: ddG >= 1.5 kcal/mol (Clackson & Wells)
===========================================================================
Chain | Pos | Res | BSA (A^2) | ddG (kcal/mol) | Classification 
------+-----+-----+-----------+----------------+-----------------
H     | 100 | D   | 115.0     | 2.530          | Critical Hotspot
H     | 33  | Y   | 95.2      | 2.094          | Hotspot        
H     | 101 | Y   | 88.9      | 1.956          | Hotspot        
H     | 31  | S   | 78.5      | 1.727          | Hotspot        
L     | 96  | W   | 70.4      | 1.549          | Hotspot        
H     | 102 | F   | 66.3      | 1.459          | Moderate       
L     | 32  | Y   | 55.0      | 1.210          | Moderate       
H     | 52  | W   | 30.4      | 0.669          | Moderate       
L     | 91  | S   | 20.0      | 0.440          | Neutral        
```

### 2. CDR-H3 Structural Base Profiling
```bash
python cli.py cdrh3 --sequence "ARDGYYYGMDV" --r94 R --r101 D --r103 W
```

### 3. Full Developability Report
```bash
python cli.py developability \
  --antibody-id "mAb-Trastuzumab-analog" \
  --vh "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDYSGMDYWGQGTLVTVSS" \
  --vl "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPYTFGQGTKVEIK"
```

---

## Python API Usage

```python
import antibody_cdr_antigen_paratope as acap

# Alanine scan
bsa_data = {("H", 33, "Y"): 95.2, ("H", 100, "D"): 115.0, ("L", 96, "W"): 70.4}
hotspots = acap.perform_alanine_scan(bsa_data)

# CDR loop decomposition
loops = acap.decompose_cdr_loops(bsa_data, numbering_scheme="IMGT")

# Developability report
report = acap.generate_developability_report(
    antibody_id="mAb-Lead-01",
    vh_sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTF...",
    vl_sequence="DIQMTQSPSSLSASVGDRVTITCRASQSI...",
)
print(f"Developability Score: {report.composite_developability_score}/100 [{report.developability_tier}]")
print(f"Isoelectric Point: {report.isoelectric_point}, Net Charge (pH 7.4): {report.net_charge_ph74}")
```

---

## Test Suite

Run the unit test suite:
```bash
python -m unittest test_antibody_cdr_antigen_paratope.py
```

All 28 test cases validate biophysical calculations, loop classification, boundary checks, and developability rules.

---

## License

MIT License. Copyright (c) 2026 Dr. Abu Suraih Sakhri.

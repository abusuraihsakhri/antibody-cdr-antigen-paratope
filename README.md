# Antibody CDR-H3 Modeling & Paratope-Epitope Contact Analyzer

A pure Python structural immunoinformatics, antibody developability, and paratope-epitope interface contact analysis engine implementing:
- **IMGT / Kabat CDR Loop Decomposition & Buried Surface Area (BSA):**
  - Quantifies total complex buried surface area ($\Delta\text{BSA} = \text{SASA}_{\text{unbound}} - \text{SASA}_{\text{complex}}$) partitioned across CDR-H1, CDR-H2, CDR-H3, CDR-L1, CDR-L2, CDR-L3, and framework regions.
  - Implements statistical free energy mapping from interface buried surface area:
    $$\Delta\Delta G_{\text{bind}} \approx 0.035 \times \Delta\text{BSA}_{\text{res}}\text{ kcal/mol}$$
- **Alanine Scanning & Binding Hotspot Identification:**
  - Classifies paratope residues into Hotspots ($\Delta\Delta G \ge 2.0\text{ kcal/mol}$), Warmspots ($1.0 \le \Delta\Delta G < 2.0\text{ kcal/mol}$), and Null/Minor contributors.
- **Shirai / Chothia CDR-H3 Torso & Base Conformation Profiling:**
  - Classifies CDR-H3 base into "Kinked" vs "Extended" conformations based on the presence of the conserved Arg94-Asp101 salt bridge and Trp103 packing interactions.
  - Analyzes loop length, aromatic fraction, and net formal charge at physiological pH 7.4.
- **Sequence Chemical Liabilities & Composite Developability Index:**
  - Scans for post-translational modification (PTM) risk motifs: Asp isomerization (DG, DS), Asn deamidation (NG, NS), Met oxidation, and unpaired Cys residues.
  - Generates composite developability scores (0–100) and tiers (High Developability, Moderate, High Liability).
- **High-Throughput Batch Sequence Screening:** Evaluates antibody VH/VL sequence libraries from CSV files.

Requires Python standard library only (zero external runtime dependencies).

---

## Paratope-Epitope Interaction Metrics

| Metric | Classification Rule | Biophysical Rationale |
|:-------|:-------------------|:----------------------|
| **Binding Hotspot** | $\Delta\Delta G \ge 2.0\text{ kcal/mol}$ ($\text{BSA} \ge 57\text{ \AA}^2$) | Primary energetic driver of complex stabilization |
| **Binding Warmspot** | $1.0 \le \Delta\Delta G < 2.0\text{ kcal/mol}$ | Secondary peripheral interface contact |
| **CDR-H3 Base Kink** | Arg94 + Asp101 salt bridge & Trp103 | Canonical $\beta$-hairpin turn conformation prevalent in $\approx 85\%$ of human antibodies |
| **Deamidation Liability** | NG, NS motifs in CDR loops | High spontaneous deamidation rate during storage |
| **Isomerization Liability** | DG, DS motifs in CDR loops | Isoaspartate formation impairing antigen binding |

---

## Features

- **Structural Biophysics Compliance:** Rigorous alanine scanning energy approximations and Shirai loop base classification.
- **PTM Liability Surveillance:** Identifies chemical degradation risks before expensive cell-line development.
- **Batch CSV Processing:** High-throughput batch triage for candidate antibody leads.

---

## Installation & Requirements

- Python 3.10+ (tested on 3.10, 3.11, 3.12)
- Zero external runtime dependencies. `pytest` is optional for running tests.

```bash
git clone https://github.com/abusuraihsakhri/antibody-cdr-antigen-paratope.git
cd antibody-cdr-antigen-paratope
```

---

## CLI Usage

### 1. CDR-H3 Conformation & Property Profiling
```bash
python cli.py cdrh3 --sequence "CARDGYYYFDYW" --json
```

### 2. Paratope Alanine Scanning Hotspot Analysis
```bash
python cli.py scan --json
```

### 3. Antibody Developability Evaluation
```bash
python cli.py developability \
  --vh "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDYSGMDYWGQGTLVTVSS" \
  --vl "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPYTFGQGTKVEIK" \
  --json
```

### 4. Batch CSV Processing
```bash
python cli.py batch --input sample.csv --output results.csv
```

---

## Python API Quickstart

```python
import antibody_cdr_antigen_paratope as acap

# Profile CDR-H3
h3_profile = acap.profile_cdrh3_conformation("CARDGYYYFDYW", r94="R", r101="D", r103="W")
print(f"Base Conformation: {h3_profile.base_conformation}")
print(f"Torso Salt Bridge: {h3_profile.torso_salt_bridge}")

# Full developability assessment
rep = acap.generate_developability_report(
    antibody_id="mAb-101",
    vh_sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDYSGMDYWGQGTLVTVSS",
    vl_sequence="DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPYTFGQGTKVEIK",
)
print(f"Score: {rep.composite_developability_score} ({rep.developability_tier})")
print(f"pI: {rep.isoelectric_point}, Humanness: {rep.overall_humanness}%")
```

---

## Running Tests

Run the test suite using standard `unittest` or `pytest`:

```bash
pytest -v
```


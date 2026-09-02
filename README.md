# Antibody Cdr Antigen Paratope

> **Domain:** Computational Biology & AI Drug Discovery  
> **Reference Guidelines & Standards:** `wwPDB, IUPAC & CLSI Computational Guidelines`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

**Antibody Cdr Antigen Paratope** is an advanced analytical and computational platform implementing IMGT antibody CDR-H3 loop modeling & paratope-epitope contact area calculator.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`HotSpotResult`** — dedicated module for hot spot result evaluation and state verification.
- **`LoopDecomposition`** — dedicated module for loop decomposition evaluation and state verification.
- **`InterfaceContactSummary`** — dedicated module for interface contact summary evaluation and state verification.
- **`ChemicalLiability`** — dedicated module for chemical liability evaluation and state verification.
- **`CDRH3StructureProfile`** — dedicated module for c d r h3 structure profile evaluation and state verification.
- **`DevelopabilityReport`** — dedicated module for developability report evaluation and state verification.

---

## 📐 Mathematical Formulation & Logic

```text
  ddg = calculate_ddg_from_bsa(bsa)
  risk = "High" if motif in ("NG", "NS") else "Medium"
  risk = "Medium" if motif == "NG" else "Low"
  risk = "High" if motif in ("DG", "DS") else "Medium"
  risk = "Medium" if motif == "DG" else "Low"
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --- <value> --json-data <value> --file <value> --json <value>
```

### Parameter Reference
- `---`: Specifies input measurement or parameter value.
- `--json-data`: Specifies input measurement or parameter value.
- `--file`: Specifies input measurement or parameter value.
- `--json`: Specifies input measurement or parameter value.
- `--scheme`: Specifies input measurement or parameter value.
- `--sequence`: Specifies input measurement or parameter value.
- `--r94`: Specifies input measurement or parameter value.
- `--r101`: Specifies input measurement or parameter value.
- `--r103`: Specifies input measurement or parameter value.
- `--antibody-id`: Specifies input measurement or parameter value.

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t antibody-cdr-antigen-paratope .
docker run -p 8000:8000 antibody-cdr-antigen-paratope
```

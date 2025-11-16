# 🔬 Strep A Inclusivity / Exclusivity Pipeline

## 🧭 Quick Start

**Main wrapper script:**  
➡️ `scripts/run_all_spa.sh`  
(This is the *only* wrapper used for full S. pyogenes all-genomes analysis.)

**Runbooks**  
- 🧪 **Demo subset workflow:**  
  [docs/runbook_panaroo.md](docs/runbook_panaroo.md)
- 🧬 **Full *S. pyogenes* all-genomes workflow:**  
  [docs/runbook_panaroo_Strep_A_All_Genomes_Analysis.md](docs/runbook_panaroo_Strep_A_All_Genomes_Analysis.md)

**Cheat Sheets**  
- 📘 [docs/Strep_A_Pipeline_cheat_sheet.md](docs/Strep_A_Pipeline_cheat_sheet.md)  
- 📗 [docs/Strep_A_Pipeline_cheat_sheet_2.md](docs/Strep_A_Pipeline_cheat_sheet_2.md)

---

# Strep A in-silico inclusivity/exclusivity pipeline

This repository contains a lightweight, reproducible workflow to:

- identify conserved (inclusive) targets across *Streptococcus pyogenes* genomes  
- verify exclusivity vs near-neighbor streptococci and common oropharyngeal flora  
- propose QC’d primers/probes for qPCR

## Workflow (overview)

A[Download target genomes<br/>S. pyogenes] --> B[Inclusivity analysis<br/>Identify core genes >98%]
A2[Download near-neighbor genomes<br/>S. dysgalactiae, S. agalactiae, etc.] --> C
B --> C[Exclusivity analysis<br/>BLAST vs exclusion panel]
C --> D[Primer/probe design<br/>Tm, GC, hairpin checks]
D --> E[Off-target testing<br/>BLAST primers vs near-neighbor DB]
E --> F[Reporting<br/>Tables + plots + summary]

markdown
Copy code

## Repository Layout

- `pipelines/panaroo/`
  - `1_annotate_prokka/` — genome annotation  
  - `2_panaroo/` — pan-genome construction  
  - `3_inclusivity/` — identity filtering  
  - `4_consensus/` — consensus building, SNP mapping  
  - `5_exclusivity/` — BLAST searches vs exclusion panel  
  - `6_reports/` — outputs, plots, summaries  
- `scripts/` — wrappers + helper utilities  
  - **`run_all_spa.sh` — main full-pipeline wrapper**  
- `docs/` — runbooks, troubleshooting notes  
- `data/` — raw genome files or downloaded datasets  
- `examples/genes/` — example per-gene analyses

---

## Documentation

For full detailed workflows, see:

- **Demo (small test set):**  
  [docs/runbook_panaroo.md](docs/runbook_panaroo.md)

- **Full S. pyogenes analysis (all genomes):**  
  [docs/runbook_panaroo_Strep_A_All_Genomes_Analysis.md](docs/runbook_panaroo_Strep_A_All_Genomes_Analysis.md)

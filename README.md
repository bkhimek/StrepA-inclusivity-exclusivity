# From the repo root
cat > /tmp/links.md <<'MD'
> **Quick links**
> - [Runbook: Panaroo (demo)](docs/runbook_panaroo.md)
> - [Runbook: All genomes](docs/runbook_panaroo_Strep_A_All_Genomes_Analysis.md)

MD

# Prepend the links block to README.md
( cat /tmp/links.md; echo; cat README.md ) > README.new && mv README.new README.md

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

## Repo layout

- `pipelines/panaroo/`
  - `1_annotate_prokka/` — genome annotation
  - `2_panaroo/` — pan-genome construction
  - `3_inclusivity/` — identity filtering
  - `4_consensus/` — consensus building, SNP mapping
  - `5_exclusivity/` — BLAST searches vs exclusion panel
  - `6_reports/` — outputs, plots, summaries
- `scripts/` — CLI wrappers and shared utilities
- `data/`
  - `demo/` — tiny test dataset
  - `results/` — generated outputs (ignored by Git)
- `docs/` — runbooks, troubleshooting notes
- `config/` — parameters (e.g., `params.yaml`)
- `examples/genes/` — example per-gene runs (ftsA, pbp2a, etc.)

---

## Documentation

📄 For the full step-by-step runbook, see:  
[docs/runbook_panaroo.md](docs/runbook_panaroo.md)

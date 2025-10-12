# Strep A in-silico inclusivity/exclusivity pipeline

This repository contains a lightweight, reproducible workflow to:
- identify conserved (inclusive) targets across *Streptococcus pyogenes* genomes
- verify exclusivity vs near-neighbor streptococci and common oropharyngeal flora
- propose QC’d primers/probes for qPCR

## Workflow (overview)

```mermaid
flowchart TD
    A[Download target genomes<br/>S. pyogenes] --> B[Inclusivity analysis<br/>Identify core genes >98%]
    A2[Download near-neighbor genomes<br/>S. dysgalactiae, S. agalactiae, etc.] --> C
    B --> C[Exclusivity analysis<br/>BLAST vs exclusion panel]
    C --> D[Primer/probe design<br/>Tm, GC, hairpin checks]
    D --> E[Off-target testing<br/>BLAST primers vs near-neighbor DB]
    E --> F[Reporting<br/>Tables + plots + summary]




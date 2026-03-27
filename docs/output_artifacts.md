---
layout: page
title: Output Artifacts Reference
permalink: /output-artifacts/
---

# Output artifacts reference (canonical)

This page is the single source of truth for **what MicroSeq writes**, **where it appears**, and **why it matters**.

Use this table first, then follow the linked docs for mode-specific details.

## Artifact table (seeded from paired single-run layout)

| Artifact path/pattern | Produced by (stage/command) | What you get | Why it matters |
| --- | --- | --- | --- |
| `raw_ab1/*.ab1` | Ingest (AB1 input retention) | Copied/symlinked original chromatograms. | Provenance and re-analysis from original traces. |
| `raw_fastq/*.fastq` | AB1 -> FASTQ conversion | Base-called FASTQ before quality trimming. | Auditable conversion boundary between vendor traces and QC. |
| `passed_qc_fastq/*.fastq` | Quality trimming/QC | Reads passing configured QC thresholds (for example avg Q). | Defines which reads are eligible for downstream assembly/BLAST. |
| `failed_qc_fastq/*.fastq` | Quality trimming/QC | Reads failing QC thresholds. | Explains data loss and supports troubleshooting threshold choices. |
| `qc/*_avg_qual.txt` | Per-read QC stats | Per-file average quality/length summaries. | Quick spot-check of individual read quality. |
| `qc/trim_summary.tsv` | QC aggregation | Combined quality/length metrics and primer-trim telemetry. | Primary run-level QC report used for auditing and filtering decisions. |
| `qc/trimmed.fasta` | QC output FASTA | FASTA of QC-passed reads inside `qc/`. | QC provenance copy used to verify what left trimming. |
| `reads.fasta` | QC handoff | Canonical merged FASTA for downstream assembly/BLAST. | Stable handoff artifact so later stages use one expected input name. |
| `qc/pairing_report.tsv` | Paired detection | Forward/reverse pairing decisions, wells, detector signals, missing mates. | Deterministic pairing audit trail in paired mode. |
| `qc/paired_fasta/*.fasta(.qual)` | Paired preprocessing | Direction-resolved FASTA/QUAL after trimming. | Inputs used to build per-sample paired assembly inputs. |
| `qc/paired_fasta_pretrim/*.fasta(.qual)` | Optional paired pretrim capture | Paired FASTA/QUAL snapshots before downstream pairing/assembly steps. | Debug aid for diagnosing trim-vs-pair effects. |
| `qc/overlap_audit.tsv` | Overlap diagnostics (enabled paths) | Candidate overlap decisions and gating outcomes. | Shows why merges succeed/fail for specific samples. |
| `qc/overlap_audit_engines.tsv` | Compare assemblers/engine audit | Per-engine overlap diagnostics (ungapped/biopython/edlib, etc.). | Explains backend differences in structural outcomes. |
| `qc/review_queue.tsv` | Resolution funnel | Samples queued for manual review + reasons/warnings. | Focuses analyst effort on actionable ambiguity/safety cases. |
| `asm/<sample>/<sample>_paired.fasta` | Paired assembly prep | Per-sample concatenated paired FASTA fed to assembler. | Exact per-sample sequence input entering assembly. |
| `asm/<sample>/<sample>_paired.fasta.qual` | Paired assembly prep | QUAL companion file for CAP3 scoring. | Preserves quality-aware assembly behavior. |
| `asm/<sample>/*.cap.contigs` | CAP3 assembly | Assembled contigs. | Main structural consensus output when overlaps are accepted. |
| `asm/<sample>/*.cap.singlets` | CAP3 assembly | Unassembled singlet sequences. | Captures residual reads when contig assembly is partial/absent. |
| `asm/<sample>/*.cap.info` | CAP3 assembly metadata | CAP3 summary details (clipping/overlap stats). | Primary diagnostic for CAP3 acceptance/removal behavior. |
| `asm/assembly_summary.tsv` | Assembly aggregation | Per-sample assembly status plus synchronized resolution fields. | High-level per-sample status table for reporting and troubleshooting. |
| `asm/blast_inputs.fasta` | Assembly -> BLAST handoff | Final sequence sent to BLAST (contigs/singlets by policy). | Defines exactly what was queried downstream. |
| `asm/blast_inputs.tsv` | Assembly -> BLAST manifest | Manifest of sequence-output selection, reasons, IDs, hypothesis/source mapping. | End-to-end provenance from structures to BLAST input sequence IDs. |
| `asm/compare_assemblers.tsv` | Compare mode | Cross-backend comparison summary for selected assemblers/profiles. | Supports backend/profile selection with explicit evidence. |
| `asm/compare/<backend_or_profile>/<sample>/*` | Compare mode per-backend outputs | Backend/profile-specific sequence outputs, merge reports, captured stdout/stderr. | Deep debugging of why one assembler/profile was selected or rejected. |
| `asm/selection_trace/*.selection_trace.tsv` | Compare + selector trace | Structured decision trace for backend/profile selection. | Reproducible explanation of automated selection decisions. |
| `hits.tsv` | `microseq blast` | Raw BLAST hits table. | Core homology evidence before taxonomy joins. |
| `hits_tax.tsv` | `microseq add-taxonomy` | BLAST hits with taxonomy columns appended. | Interpretable taxonomic output for downstream summarization. |
| `results.biom` (or chosen BIOM output) | `microseq postblast` | BIOM table aligned with metadata/taxonomy. | Direct handoff to QIIME2/phyloseq/scikit-bio ecosystems. |
| `<output>_blast_provenance.csv` | `microseq postblast` sidecar | BLAST provenance table alongside BIOM output. | Transparency for post-BIOM record lineage and joins. |
| `<output>_taxonomy_only.csv` | `microseq postblast` sidecar | Taxonomy-focused export. | Lightweight table for taxonomy-only checks and reporting. |
| `<output>_metadata.csv` | `microseq postblast` sidecar | Metadata-aligned output table. | Verifies metadata join behavior and sample-ID normalization effects. |

## Abbreviated paired single-run tree (reference layout)

```text
paired_single_pair_ab1_demo_run_microseq/
├── raw_ab1/
├── raw_fastq/
├── passed_qc_fastq/
├── failed_qc_fastq/
├── qc/
│   ├── trim_summary.tsv
│   ├── pairing_report.tsv
│   ├── overlap_audit.tsv
│   ├── review_queue.tsv
│   └── trimmed.fasta
├── reads.fasta
├── asm/
│   ├── assembly_summary.tsv
│   ├── blast_inputs.fasta
│   ├── blast_inputs.tsv
│   ├── compare_assemblers.tsv
│   ├── compare/
│   └── selection_trace/
├── hits.tsv
└── hits_tax.tsv
```

## Related docs

* GUI run behavior and examples: [GUI Walkthrough](gui_usage.md)
* Paired mode details and CAP3 outputs: [Paired CAP3 Assembly](paired_assembly.md)
* Resolution contracts (`review_action`, `warning_flags`, hypothesis counts): [Workflow Resolution Funnel](workflow_resolution.md)

---
layout: page
title: MicroSeq Docs
---

MicroSeq is a microbiology-focused Sanger workflow for turning raw AB1 or FASTQ reads into QC-passed sequences, paired assemblies when appropriate, BLAST results, taxonomy-aware tables, and optional downstream BIOM outputs.

If you are using MicroSeq for 16S isolate identification or related microbiology workflows, this page is the best place to start.

## Where should I start?

### Most microbiology users
* [GUI Walkthrough](gui_usage.md) | how to run MicroSeq through the graphical interface
* [Output Artifacts Reference](output_artifacts.md) | what files MicroSeq writes and what they mean

### If you are working with paired forward/reverse Sanger reads
* [Paired CAP3 Assembly Tutorial](paired_assembly.md) | what paired mode means biologically and how MicroSeq decides what goes to BLAST
* [GUI Table Reference](gui_table_reference.md) | detailed meanings of GUI table columns

### If you are running from a terminal, remote server, or HPC
* [CLI Workflows](cli_workflow_tutorial.md) | current headless workflow and wrapper usage

### Advanced and reference pages
* [Workflow Resolution Funnel](workflow_resolution.md) | advanced status routing, ambiguity handling, and sample-level evidence collapse
* [MicroSeq's Design: Questions & Answers](microseq_design_questions.md) | practical design explanations
* [VSEARCH Chimera & Replicate Filters](vsearch_filters.md) | optional VSEARCH-related workflows

## Tiny glossary

| Term | Meaning |
| --- | --- |
| `Q-cov` | Query coverage. This is the percentage of your query sequence covered by the alignment to a database hit. |
| `contig` | A merged consensus sequence built from overlapping read evidence, such as forward and reverse reads from the same sample. |
| `singlet` | A standalone sequence record that remains usable when a defended merged contig was not produced. |
| `pair_missing` | A paired-mode outcome where only one direction survived pairing, so the sample never had a valid forward/reverse assembly opportunity. |
| `filename primer labels` | The labels in the filename, such as `27F` or `1492R`, that MicroSeq uses to recognize forward and reverse files during pairing. |

## Quick launch

For most local microbiology users, the GUI is the easiest starting point:

```bash
conda activate MicroSeq
microseq-gui

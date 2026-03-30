---
layout: page
title: CLI Workflows
permalink: /cli-workflows/
---

## What this page is for

Use this page when you want to run MicroSeq from a terminal, inside a bash wrapper, over SSH, or on an HPC cluster. It documents the current headless paired-AB1 workflow and shows how the wrapper maps onto the underlying MicroSeq CLI stages.

This page is operational. For pairing rules, assembly semantics, ambiguity handling, and BLAST-input behavior, see the [Paired CAP3 Assembly Tutorial](paired_assembly.md). For the meaning of output files, see [Output Artifacts Reference](output_artifacts.md).

## When to use this page instead of the GUI

Use the CLI/HPC workflow when:

* you want reproducible batch runs
* you are working on a remote server or cluster
* you want to call MicroSeq from a shell script or scheduler
* you want one run-scoped log across multiple stages
* you want a stable operational workflow that mirrors the paired pipeline without relying on interactive dialogs

# MicroSeq CLI / HPC Quickstart

## Install MicroSeq

For installation, environment setup, and database setup, follow the main README:

- GitHub repo: https://github.com/jose-cantu/MicroSeq
- README: https://github.com/jose-cantu/MicroSeq/blob/master/README.MD

This includes:

- creating and activating the `MicroSeq` conda environment
- installing the package
- running `microseq-setup`
- configuring databases and logging

Once installation is complete, activate the environment:

```bash
conda activate MicroSeq
```

## Current paired AB1 HPC workflow

For using the bash wrappers I have setup refer to this readme I made here:
- AB1 CLI README: https://github.com/jose-cantu/MicroSeq/tree/master/workflows 

The current bash wrapper for paired AB1 runs is:

```bash
bash workflows/paired_ab1_pipeline.sh <input_dir> <db_key> <blast_threads> [out_dir]
```

Arguments:

- `input_dir` = directory containing raw paired `.ab1` files
- `db_key` = `gg2`, `silva`, or `ncbi`
- `blast_threads` = number of threads to pass to BLAST
- `out_dir` = optional output directory; default is `<input_dir>_microseq`

## Minimal example

```bash
conda activate MicroSeq

bash workflows/paired_ab1_pipeline.sh \
  tests/paired_single_pair_ab1_demo_run \
  gg2 \
  4
```

## Optional primer-label override

If your forward and reverse primer labels differ from the demo defaults, override them before running the wrapper:

```bash
conda activate MicroSeq

export FWD_PATTERN="27F"
export REV_PATTERN="1492R"

bash workflows/paired_ab1_pipeline.sh \
  /path/to/paired_ab1_inputs \
  gg2 \
  4
```

## Optional run-scoped logging

If you want every `microseq` stage launched by the wrapper to write into one run-scoped log context:

```bash
export MICROSEQ_SESSION_ID="paired-run-001"
```

Then run the wrapper normally.

## What the wrapper currently does

The wrapper currently runs the paired AB1 workflow in this order:

1. trim AB1 -> QC-passed FASTQ
2. stage paired FASTA/QUAL
3. write `qc/pairing_report.tsv`
4. run paired assembly
5. write `asm/assembly_summary.tsv`
6. write `qc/overlap_audit.tsv`
7. write `asm/blast_inputs.fasta` and `asm/blast_inputs.tsv`
8. run BLAST
9. join taxonomy

## Key outputs

After a successful run, the main outputs are:

- `qc/trim_summary.tsv`
- `qc/pairing_report.tsv`
- `qc/overlap_audit.tsv`
- `asm/assembly_summary.tsv`
- `asm/blast_inputs.fasta`
- `asm/blast_inputs.tsv`
- `hits.tsv`
- `hits_tax.tsv`

## Notes

- Use the README for installation and setup.
- Use the 2nd README in `/workflow` for the preview pairs wrapper and this wrapper for the current headless paired AB1 workflow.
- For pairing/assembly behavior and output semantics, see the paired assembly documentation in `/docs`.



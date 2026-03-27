---
layout: page
title: vsearch chimera checking & replicate collapse

MicroSeq can call vsearch to do the following currently.

Collapse technical replicates per sample (`--replicate-id-th 1.0` = strict exact-match collapse; values `< 1.0` add a high-identity clustering pass).  
Run chimera based reference checks (UCHIME-ref) is what is used.  
Orient sequences against a reference (vsearch --orient) to fix reverse-complemented Sanger sequences.

If used these stages are run after FASTA generation and, when applicable, after paired assembly, before BLAST, so they reduce redundant sequences and filter chimeras early without changing the core MicroSeq workflow utilizing single/paired files and CAP3 for assembly.

## When To Use These Stages

* **Collapse replicates** may be used when you have technical repeats of the same isolate or when you want to reduce duplicate BLAST work. In the standard one-isolate, one-sequence Sanger workflow, this is an advanced utility rather than a core step.
* **Chimera check** is when you want an extra QC pass to remove obvious PCR chimeras before classification.
* **Orient sequences** is when your input FASTA can contain reverse-complemented reads (common in Sanger when primers differ). It normalizes orientation before collapse/chimera so technical replicates match on the same strand.

All are off by default and do not affect the standard pipeline I have setup. When collapse is used, `qc/replicate_weights.tsv` records technical-support counts for identical collapsed sequences. These are support values, not biological abundance.

In the GUI this is exposed as an explicit **Orient reads** toggle (default off) so users can keep runs reproducible and transparent. I recommend turning it on for mixed-orientation Sanger datasets (forward + reverse primer runs).

## CLI Usage

> Note: I made a bash workflow the user can run. The current vsearch wrapper is a minimal pre-BLAST classification workflow and does not include post-BLAST / BIOM generation.

```bash
# Collapse replicates + reference chimera checks
workflows/microseq_full_with_vsearch.sh reads/ gg2 8
```

More commands that can be used.

* `--replicate-id-th`: if set < 1.0, run a high-identity clustering pass after exact-match collapse.
* `--min-replicate-size`: minimum unique size for dereplication.
* `--chimera-db`: override the default chimera reference FASTA.
* `--orient-db`: override the default orient reference FASTA.

## Chimera reference selection

When `--chimera-mode reference` is enabled, MicroSeq uses *the same DB family you selected for BLAST*:

* `--db gg2` -> uses `databases.gg2.chimera_ref`
* `--db silva` -> uses `databases.silva.chimera_ref`

By default, MicroSeq uses the chimera reference associated with the selected database family. For `--db gg2`, this is the Greengenes2 backbone full-length reference; for `--db silva`, this is SILVA 138.1 SSU Ref NR99. In the current MicroSeq reference set, SILVA is the broader reference for UCHIME-ref screening, so SILVA is the recommended choice when maximal parent coverage is more important than keeping the chimera reference matched to the BLAST database. Use the db-matched reference when consistency with the downstream taxonomy database is the higher priority.

> **Note:** `--replicate-id-th = 1.0` performs strict exact-match collapse after FASTA generation and, when applicable, after paired assembly. Values `< 1.0` shift the behavior to within a sample that is in essence centroiding and may merge true variants; keep `1.0` unless you explicitly accept reduced variant resolution. I’ll leave that up to you as the user which you prefer given the nature of your project.

> **Also Note:** When `--chimera-mode reference` is enabled, MicroSeq can preserve `;size=` annotations when replicate sizes are present. This keeps technical-support annotations intact for downstream reporting. In reference chimera mode, `uchime_ref` does not use `--sizein` for scoring, so these tags are preserved for provenance rather than chimera scoring.

## Orient reference selection

When orientation is enabled, MicroSeq uses the configured orientation reference:

* `--db gg2` -> uses `databases.gg2.orient_ref` (falls back to `databases.gg2.chimera_ref`)
* `--db silva` -> uses `databases.silva.orient_ref` (falls back to `databases.silva.chimera_ref`)

The oriented reads are written to `qc/oriented.fasta`, the reads with undetermined orientation are written to `qc/orient_notmatched.fasta`, and a vsearch tabular report is written to `qc/orient_report.tsv`. In the full pipeline path, collapse/chimera runs on the oriented sequences and then merges the notmatched sequences back before BLAST so they are still classified (but skipped for chimera filtering). This is recommended for reverse-primer Sanger datasets where strand flips are common.

## Workflow Scripts

Below are workflow scripts that I have created in bash featuring the pipeline MicroSeq offers to users.

The workflow scripts are thin minimal wrappers for the same datasets:

* `workflows/microseq_full_pipeline.sh` steps: trim -> FASTQ to FASTA if needed -> BLAST -> taxonomy join (no vsearch, no postblast).
* `workflows/microseq_full_with_vsearch.sh` steps: trim -> FASTQ to FASTA if needed -> vsearch collapse -> reference chimera -> BLAST -> taxonomy join (no postblast / BIOM).

## Configuration

`config.yaml` includes a `chimera_ref` key per database and optional `orient_ref`:

```yaml
databases:
  gg2:
    chimera_ref: ${MICROSEQ_DB_HOME}/gg2/dna-sequences.fasta
    orient_ref: ${MICROSEQ_DB_HOME}/gg2/dna-sequences.fasta
  silva:
    chimera_ref: ${MICROSEQ_DB_HOME}/silva/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta
    orient_ref: ${MICROSEQ_DB_HOME}/silva/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta
```

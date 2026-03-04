---
layout: page
title: vsearch chimera checking & replicate collapse 

MicroSeq can call vsearch to do the following currently.

Collapse technical replicates (so dereplication / clustering ) per sample. 
Run chimera based reference checks (UCHIME-ref) is what is used. 
Orient reads against a reference (vsearch --orient) to fix reverse-complemented Sanger reads. 

If used these stages are run after trimming and assembly and before BLAST, so they reduce redundant sequences and filter chimeras early without changing the core MicroSeq workflow utilizing single/paired files and CAP3 for assembly. 

## When To Use These Stages 

* **Collapse replicates** should be used when you have technical repeats of the same isolate or when you want to reduce duplicate BLAST work. 
* **Chimera check** is when you want an extra QC pass to remove obvious PCR chimeras before classification. 
* **Orient reads** is when your input FASTA can contain reverse-complemented reads (common in Sanger when primers differ). It normalizes orientation before collapse/chimera so technical replicates match on the same strand.

All are off by default and do not affect the standard pipeline I have setup. Replicate sizes are recorded in `qc/replicate_weights.tsv` and used to weight BIOM counts. 

In the GUI this is exposed as an explicit **Orient reads** toggle (default off) so users can keep runs reproducible and transparent. I recommend turning it on for mixed-orientation Sanger datasets (forward + reverse primer runs).

## CLI Usage

> Note: I made a bash workflow the user can run. The vsearch variation also includes post-blast hence the call to having a metadata.tsv file. 

```bash
# Collapse replicates + reference chimera checks 
workflows/full_pipeline_with_vsearch.sh reads/ gg2 metadata.tsv 
```

More commands that can be used.

*`--replicate-id-th`: if set < 1.0, run a high‑identity clustering pass after dereplication.
* `--min-replicate-size`: minimum unique size for dereplication.
* `--chimera-db`: override the default chimera reference FASTA.
* `--orient-db`: override the default orient reference FASTA.

## Chimera reference selection

When `--chimera-mode reference` is enabled, MicroSeq uses *the same DB family you selected for BLAST*:

* `--db gg2` -> uses `databases.gg2.chimera_ref`
* `--db silva` -> uses `databases.silva.chimera_ref`

The reference is only as good as the size of the database and the curation of the database as well. Because of this I would recommend using the classical SILVA as the default unless you are using GG2 as your blast default for your project needs. 

> **Note:** `--replicate-id-th = 1.0` performs strict exact‑match collapse after assembly. Values `< 1.0` shift the behavior to within a sample that is in essence centroiding and may merge true variants; keep `1.0` unless you explicitly accept reduced variant resolution. I'll leave that up to you as the user which you perfer given the nature of your project. 

> **Also Note:** When `--chimera-mode reference` is enabled, MicroSeq forwards `--sizein` to vsearch replicate sizes are present. This perserves abundance anntotations without asserting how `uchime_ref` scores them. This keeps the `;size=` annotations intact as vsearch reads the input FASTA, so downstream stages (like replicate weighting in post‑BLAST) can still rely on those abundance tags. MicroSeq does **not** assume that `uchime_ref` uses abundance in its scoring; if you need abundance‑driven chimera scoring, I would consider instead using a de novo chimera step instead which I plan to implement in MicroSeq in a future update soon. 

## Orient reference selection

When orientation is enabled, MicroSeq uses the configured orientation reference:

* `--db gg2` -> uses `databases.gg2.orient_ref` (falls back to `databases.gg2.chimera_ref`)
* `--db silva` -> uses `databases.silva.orient_ref` (falls back to `databases.silva.chimera_ref`)

The oriented reads are written to `qc/oriented.fasta`, the reads with undetermined orientation are written to `qc/orient_notmatched.fasta`, and a vsearch tabular report is written to `qc/orient_report.tsv`. The pipeline runs collapse/chimera on the oriented reads and then merges the notmatched reads back before BLAST so they are still classified (but skipped for chimera filtering). This is recommended for reverse‑primer Sanger datasets where strand flips are common. 

## Workflow Scripts

Below are workflow scripts that I have created in bash featuring the pipeline MicroSeq offers to users. 

The workflow scripts mirror the GUI pipeline for the same datasets:

* `workflows/full_pipeline.sh` steps: trim -> skip assembly/assemble -> BLAST -> taxonomy join (no vsearch).
* `workflows/full_pipeline_with_vsearch.sh` steps: trim -> FASTA -> assembly/skip assembly -> vsearch collapse -> reference chimera -> BLAST -> taxonomy join (no postblast).
* `workflow/full_pipeline_with_vsearch_and_post_blast.sh` steps: ab1/fastq -> trim -> fasta -> skip aseembly/assembly -> vsearch collapse -> reference chimera -> BLAST -> taxonomy join -> postblast

Both `tests/paired_ab1_demo_run/...` and `tests/reverse_orientation_only_ab1_demo_run/...` are AB1 directories, so the same scripts applies to either layout.

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






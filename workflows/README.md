# MicroSeq Paired AB1 Workflows

These workflow scripts are thin wrappers around MicroSeq’s atomic CLI commands for **paired AB1 / paired Sanger** inputs.

They are intentionally minimal: no monolithic pipeline subcommand is required. Each workflow composes the same lower-level MicroSeq commands used by the GUI and CLI.

## Recommended order of use

For paired AB1 data, use this order:

1. `preview_pair_names.sh`
2. authoritative MicroSeq paired preview
3. `paired_ab1_pipeline.sh`

This separates:

- **filename cleanup**
- **pairing validation**
- **full analysis**

which makes pairing failures easier to diagnose before you commit to the full run.

---

## Scripts

### 1) `preview_pair_names.sh`

A filename-hygiene helper for paired inputs.

What it does:

- recursively scans sequence-like files under an input directory
- checks whether forward / reverse primer labels look canonically delimited in filenames
- reports:
  - `CURRENT`
  - `STATUS`
  - `PRIMER_LABEL`
  - `ISSUE`
  - `SUGGESTED`
- optionally renames `FIXABLE` files in place
- prints the recommended authoritative MicroSeq paired preview command afterward

Current default primer labels checked by the helper:

- Forward: `27F`, `8F`, `515F`
- Reverse: `1492R`, `806R`, `926R`

Important notes:

- this script is a **filename helper**, not the authoritative pairing engine
- it ignores well matching
- after using it, run the actual MicroSeq paired preview step

Status meanings:

- `OK` = primer label found and clearly delimited
- `FIXABLE` = primer label found but jammed into the sample name; a rename suggestion is generated
- `REVIEW` = missing or ambiguous primer label; inspect manually

Usage:

```bash
bash workflows/preview_pair_names.sh <input_dir>

Example:
bash workflows/preview_pair_names.sh tests/paired_single_pair_ab1_demo_run

Usage:

paired_ab1_pipeline.sh <input_dir> <db_key> <blast_threads> [out_dir] \
  [--identity N] \
  [--qcov N] \
  [--max-hits N]

Arguments:

<input_dir>: directory containing raw paired AB1 files
<db_key>: gg2 | silva | ncbi
<blast_threads>: positive integer used for BLAST
[out_dir]: optional output directory
default: <input_dir>_microseq

Optional BLAST overrides:

--identity N: override the canonical CLI default BLAST identity threshold
--qcov N: override the canonical CLI default BLAST query coverage threshold
--max-hits N: override the canonical CLI default BLAST max target hits

If these optional BLAST flags are omitted, the wrapper inherits the canonical MicroSeq CLI defaults. (Identity: 97%, Qcov: 80%, max-hits: 5, threads: 4) 

Example of runs of how I would use it:
bash workflows/paired_ab1_pipeline.sh "tests/paired_single_pair_ab1_demo_run" gg2 4

or 

bash workflows/paired_ab1_pipeline.sh \
  "tests/paired_single_pair_ab1_demo_run" \
  gg2 \
  4 \
  "tests/paired_single_pair_ab1_demo_run_microseq"
or 

bash workflows/paired_ab1_pipeline.sh \
  "tests/paired_single_pair_ab1_demo_run" \
  gg2 \
  4 \
  --identity 95 \
  --qcov 85 \
  --max-hits 7

Example stage progression

A typical successful run looks like:

[1/6] trim AB1 -> QC-passed reads
[2/6] stage paired FASTA
[2.5/6] pairing report
[3/6] paired assembly
[3.3/6] assembly summary
[3.6/6] overlap audit
[4/6] build BLAST input FASTA
[5/6] BLAST
[6/6] taxonomy join

Key outputs from a successful run

trim summary      : <out_dir>/qc/trim_summary.tsv
pairing report    : <out_dir>/qc/pairing_report.tsv
assembly summary  : <out_dir>/asm/assembly_summary.tsv
blast input sequences : <out_dir>/asm/blast_inputs.tsv
blast hits        : <out_dir>/hits.tsv
hits + taxonomy   : <out_dir>/hits_tax.tsv

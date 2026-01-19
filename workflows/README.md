# MicroSeq CLI Workflows

These scripts are thin wrappers that run an end-to-end MicroSeq classification workflow from input reads to a taxonomy-annotated hits table.

They are intentionally minimal: no “pipeline” subcommand is required. Each workflow composes MicroSeq’s atomic CLI commands the same way the GUI does.

## Scripts

### 1) `microseq_full.sh`
Runs:
1) trim/QC (if input is FASTQ/AB1)
2) FASTQ → FASTA (if needed)
3) BLAST
4) add_taxonomy

Outputs:
- `reads.fasta`
- `hits.tsv`
- `hits_tax.tsv`

### 2) `microseq_full_with_vsearch.sh`
Runs the same steps as `microseq_full.sh`, plus:
- per-sample replicate collapse via `microseq vsearch-collapse`
- reference-based chimera filtering via `microseq vsearch-chimera`
- removes `;size=` tokens from FASTA headers before BLAST so BLAST IDs remain clean

Outputs (in addition to base outputs):
- `qc/replicates_collapsed.fasta`
- `qc/nonchimeras.fasta`
- `qc/nonchimeras_clean.fasta`

## Dependencies

These workflows assume the following are available in your PATH (typically from the MicroSeq conda environment):
- `microseq`
- `blastn` (BLAST+)
- `cap3` (only if you use `microseq assembly`)
- `vsearch` (required only for the vsearch workflow)

## Usage

Example (no vsearch):
```bash
workflows/microseq_full.sh reads/ gg2 8


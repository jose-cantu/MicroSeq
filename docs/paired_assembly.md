---
layout: page 
title: Paired CAP3 assembly 
permalink: /paired-assembly/
---

MicroSeq now ships a forward/reverse pairing workflow that works in both the CLI and GUI.
It can auto-detect primer tokens in filenames, enforce plate well consistency, and run CAP3 per sample with duplicate-handling policies.

## What changed

* **Primer token registry + suggestions.** Filenames are scanned with mid/prefix/suffix detectors by default; you can also provide your own regex via `--fwd-pattern/--rev-pattern`.
  When pairing fails, MicroSeq will summarise which detectors matched and suggest concrete regexes to try.
* **Well-aware pairing.** Plate positions such as `A1`–`H12` can be required with `--enforce-well`; mismatched or missing wells are skipped with an explicit report.
* **Duplicate policies.** Choose whether duplicate forward/reverse files error, keep-first, keep-last, merge, or **keep-separate** (runs CAP3 per duplicate and preserves separate outputs).
* **Generalised inputs.** Mixed FASTA, FASTQ, or AB1 inputs are normalised to per-read FASTA before pairing so you can point at a folder or a merged FASTA file.
* **Contig traceability.** CAP3 headers are rewritten with sample-aware IDs and a `contig_map.tsv` manifest records which input reads built each contig.
* **GUI support.** The GUI now stages multi-file selections, exposes single vs paired assembly, primer presets, advanced regex toggles, duplicate policy selector, well enforcement, and a pairing preview dialog.

## CLI tutorial: guaranteed forward/reverse pairing

1) **Prepare inputs.** Place your forward/reverse FASTA (or FASTQ/AB1) files under one directory. Default detectors recognise tokens like `_27F_` / `-1492R-` anywhere in the filename.
2) **Run paired assembly with defaults** (auto-detection + CAP3):

```bash
microseq assembly --mode paired \
  -i data/paired_reads/ \
  -o results/paired_cap3 \
  --dup-policy merge
```

* I would recommend just defaulting to `--dup-policy error` unless otherwise. 
* `--dup-policy merge` merges duplicate orientations before CAP3; swap to `keep-first`, `keep-last`, or `keep-separate` as needed.
* Outputs land in `results/paired_cap3/<sample>/` with `*.cap.contigs`, `*.cap.singlets`, `contig_map.tsv`, and a pairing report.

3) **Require wells + custom primer regex** (when filenames include plate codes and unusual primer names):

```bash
microseq assembly --mode paired \
  -i plates/run1/ \
  -o results/well_enforced \
  --fwd-pattern "(27F|8F|customF)" \
  --rev-pattern "(1492R|806R|customR)" \
  --enforce-well --well-pattern "(?i)[A-H](?:0?[1-9]|1[0-2])" \
  --dup-policy keep-separate
```

* Custom regex detectors are tried before the built-in registry so edge-case primer strings can still be paired.
* `--enforce-well` buckets pairs by plate position; files without wells or with mismatched wells are dropped and called out in the report.
* `--dup-policy keep-separate` emits separate CAP3 runs per duplicate orientation instead of merging them.

4) **Diagnose failures quickly.** When no pairs are found, MicroSeq prints a summary of how many files matched each detector, which orientations are missing, and concrete `--fwd-pattern/--rev-pattern` suggestions extracted from filenames. Apply those regexes and re-run.

## GUI tutorial: pairing with presets, wells, and duplicate handling

0) In your local terminal you will need to first activate MicroSeq in the conda environemnt. So open the termimal and Conda activate MicroSeq by first doing `conda activate MicroSeq`.
1) **Launch and pick inputs.** Next, you will type in`microseq-gui` and that will open the gui (wait a few moments it will open), choose **Assembly** and toggle **Paired**. Browse to a folder or select multiple files; the GUI now stages your selection into a temporary folder so CAP3 sees a consistent path.
2) **Primer detection.** Use the primer-set dropdown to pick common token/primers or type your own forward/reverse regex. Toggle **Advanced regex** if you need case-insensitive or grouped patterns. The GUI will auto-detect tokens from filenames and surface them in the pairing preview dialog.
3) **Well enforcement.** Check **Enforce well codes** to require plate positions (A1–H12 by default). You can edit the regex if your lab uses a different format.
4) **Duplicate policy.** Choose whether duplicates **Error**, **Keep first**, **Keep last**, **Merge**, or **Keep separate**. The choice persists between sessions via QSettings.
5) **Preview pairs.** Click **Preview pairs** to see which files were detected as forward/reverse, which detector matched, any missing mates, and suggested regexes. Fix naming if needed, then proceed.
6) **Run assembly.** The GUI launches CAP3 per paired sample, rewrites contig headers, and writes a `contig_map.tsv` manifest alongside contigs and singlets for traceability. Logs include the pairing report so you can audit wells, detectors, and duplicate handling.

## Tips

* Keep forward/reverse tokens adjacent to underscores or dashes so the detectors can split the sample ID cleanly.
* When supplying regexes, wrap alternatives in parentheses (e.g., `(27F|8F)`), and include trailing `F`/`R` so orientations are unambiguous.
* For plate runs, enable well enforcement to avoid mixing reads across wells; the manifest reflects the enforced well key.
* Use **keep-separate** when you intentionally want per-orientation CAP3 outputs (e.g., troubleshooting primer performance) without merging reads.

## Mock Example workflows to get an idea how to use it 

### 96-well plate with standard 16S primers (paired mode)
*Folder layout:* `plates/2025-05-plateA/` contains files like `A01__27F.ab1`, `A01__1492R.ab1`, … `H12__1492R.ab1`.

1) Preview the pairs to validate wells and primer detection:

```bash
microseq assembly --mode paired \
  -i plates/2025-05-plateA/ \
  --dup-policy merge \
  --enforce-well \
  --fwd-pattern "(27F|8F|515F)" \
  --rev-pattern "(1492R|806R)" \
  --dry-run
```

The dry-run prints how many files matched each detector, which wells are missing mates, and suggests regex tweaks if a primer name is unexpected.

2) Run for real and collect per-sample contigs:

```bash
microseq assembly --mode paired \
  -i plates/2025-05-plateA/ \
  -o results/plateA_cap3 \
  --dup-policy merge \
  --enforce-well
```

Outputs include `results/plateA_cap3/<sample>/contig_map.tsv` plus CAP3 contigs and singlets; the manifest records which input reads built each contig with the enforced well code.

Single mode either does forward/reverse and simply runs CAP3 on the merged FASTA; use this when you know each file is already single-direction only.

## Real World Example to Pilot and test run yourself 


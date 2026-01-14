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

### Paired assembly example in tests directory

So using the CLI you can do the same thing that was done in the demo example of the GUI of the Full Pipeline run which will be done in four steps below. 
This keeps outputs in the same `tests/paired_ab1_demo_run/10292025_1080497_microseq/` tree the GUI uses to keep it consistent. Note you can easily just put this in bash script if you wish just changed it for your purposes. I have preferences for using dup policy being error given the settings i have for it. 

```bash
# set paths first 
# Here we will be input and output 
INPUT=tests/paired_ab1_demo_run/10292025_1080497 
OUTPUT=tests/paired_ab1_demo_run/10292025_1080497_microseq 

# step 1) QC + Trim Sanger traces 
microseq trim -i "$INPUT" --sanger --workdir "$OUT" 

# Step 2) Pair + assembly with CAP3 
microseq assembly --mode paired \
   -i "$OUTPUT/qc/paired_fasta" \
   -o "$OUTPUT/asm" 
   --dup-policy error 

# Step 3) Blast assembled contigs against Greengenes2 
microseq blast -i "$OUTPUT/asm/paired_contigs.fasta" -d gg2 -o "$OUTPUT/hits.tsv"

# Step 4) Join taxonomy 
microseq add_taxonomy -i "$OUTPUT/hits.tsv" -d gg2 -o "$OUTPUT/hits_tax.tsv"
```

Notes:

* Step 1 writes `qc/pairing_report.tsv`, per-read quality stats, length, avg Phred, MEE, and trimmed FASTAS. `--sanger` triggers AB1-> FASTQ path. 
* Step 2 reuses the trimmed forward/reverse FASTAS under `qc/paired_fastas/` so CAP3 sees quality-filtered reads; from there adjust `dup-policy` if you need anther option based on your situation. 
* Step 3/4 match the GUI's **Full Pipeline** path CAP3 -> BLAST -> Taxonony. I will reference using `microseq postblast` in a separate tutorial doc on how to use it and its uses cases. This is if you want a BIOM/CSV table from `hits_tax.tsv`.  

For now I will also keep the single file mode example under here. I will reference it in a separate file just for the CLI at a certain point should I feel it necessary. 

### Single mode ab1 example using the test directory dataset 

For the single direction ab1 demo run using the CLI so cd to `tests/reverse_orientation_only_ab1_demo_run/`, the CLI sequence below mirors what is being run in the GUI when single mode full pipeline example is being used from the GUI tutorial. 

```bash
# set paths here 
INPUT=tests/reverse_orientation_only_ab1_demo_run
OUTPUT=tests/reverse_orientation_only_ab1_demo_run_microseq 

# step 1) QC + trim Sanger traces 
microseq trim -i "$INPUT" --sanger --workdir "$OUTPUT" 

# step 2) Single Mode treat this as an optional step the best practice is to not use CAP3 for non paired reads..... I may remove this this the GUI version when single is used skips CAP3 altogether and is there as a highlighter only to understand they are using single orientation only. 
microseq assembly --mode single \
   -i "$OUTPUT/reads.fasta" \
   -o "$OUTPUT/asm" 

# step 3) Blast QC-passed reads 
microseq blast -i "$OUTPUT/reads.fasta" -d gg2 -o "$OUTPUT/hits.tsv" 

# Step 4) Join Taxonomy 
microseq add_taxonomy -i "$OUTPUT/hits.tsv" -d gg2 -o "$OUTPUT/hits_tax.tsv"
```

Some Notes to consider:

* Step 1 produces the same layout the GUI writes so you have `raw_ab1/`, `raw_fastq/`, `passed_qc_fastq/`, `failed_qc_fastq/`, `qc/trim_summary.tsv`, which now includes avg MEE, avg MEE per kb, and qeff/label summaries, and `reads.fasta` (the canonical merged FASTA used by later steps). 
* Step 2 is honestly optional at this point I don't bother using it in the GUI and go straight to readings my reads.fasta that's already QC-passed. So feel free to skips it.
* Step 3 matches the GUI **Full Pipeline** for single mode (BLAST + Taxonomy on the QC passed reads) unless you chose CAP3 for whatever reason. 

## Frequently Asked Questions 

### What happens to forward/reverse reads in single mode? 

`microseq assembly --mode single` skips the pairing detector entirely. Any reads you give it including mixed forward/reverse orientations are treated as independent sequences and passed stright to CAP3 as as single pool; so no `pairing_report.tsv` is emitted for example. If you need pairing orientation awareness matching then use the well and paired mode commands instead so MicroSeq can build explicit F/R pairs before assembly. 

### How can I make sure that the files are auto-paired correctly? 

MicroSeq buckets files by sample ID after stripping the primer tokens. With the default dectectors I have setup (or you can use the custom `--fwd-pattern/--rev-pattern`), names like `KD001_27F.*` and `KD001_1492R.*` form a pair because of the shared prefix `KD001` remains once the tokens are removed. In contrast, `KD001_27F*` and `KD004_1492R` stay separate because their prefixes differ. MicroSeq does not attempt to cross-numnber pairs its unreliable.

Lets consider a mock scenario below:

**Mock-Scenario For My Reader**

* Forward reads: `KD001_27F`, `KD002_27F`, `KD003_27F`
* Reverse reads: `KD004_1492R`, `KD005_1492R`, `KD006_1492R`

These do not auto-pair because the prefixes (`KD001` vs `KD004`) don’t match. Rename to shared sample IDs (e.g., `KD001_27F` with `KD001_1492R`) so MicroSeq can align them. If your primer labels differ from the defaults, point `--fwd-pattern/--rev-pattern` at your tokens (for example, `--fwd-pattern "27F" --rev-pattern "1492R"`).

**Primer token position examples (all auto-pair because the sample ID matches)** 

* Primer in the middle: `KD001_27F_A01.ab1` ↔ `KD001_1492R_B01.ab1`
* Primer in the front: `27F_KD001_A01.ab1` ↔ `1492R_KD001_B01.ab1` 
* Primer in the end:   `KD001_A01_27F.ab1` ↔ `KD001_B01_1492R.ab1`

The default detectors cover the front/middle/end placements; use custom regex if your primer token labeling differs. Remember all of this hinges on the shared sample ID which in this case is `KD001`.

This means MicroSeq also strips well tokens (`A01-H12`) from sample IDs even when primer tokens sit at the front or end of the filename, so 27F_KD001_A01.ab1` and `KD001_B01_1492R.ab1` still resolve to the shared sample ID `KD001` with well enforcement turned off. 

Well codes only influence pairing when `--enforce-well` is enabled: files must then share both the sample ID **and** the same plate position (e.g., `A01`). Leave well enforcement off when sample IDs are already unique and wells are irrelevant; turn it on for plate exports where cross-well swaps would be problematic. 





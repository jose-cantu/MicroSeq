---
layout: page 
title: Gui Walkthrough # <title> element shown in browser tab 
permalink: /gui/
--- 

## Launch 

* Note in order to launch the GUI you need to have conda activated and from there you type in the  terminal `microseq-gui`. 

* The window will show the defaults:
Control	Default	Meaning
DB	gg2	Reference database
ID %	97	BLAST identity cut-off
Q-cov %	80	Query coverage cut-off
Threads	4	CPU threads
Mode	Fast (megablast)	       Alignment algorithm
Assembly Single/Paired (Cap3)      Choose Single or Paired Assembly {Single assumes you will not assemble with a Forward & Reverse primers and will only choose one}. 
Dup Policy Error                   How duplicate orientations are handles between files in a run. 
Well enformcement off              Requries A1-H12 plate well labeling so when enabled it will match the same well between forward and reverse files. 




## Linux/Wayland stability mode

If you are on Linux Wayland and see maximize crashes (`xdg_surface buffer ... does not match configured state`), MicroSeq now defaults to `QT_QPA_PLATFORM=xcb` for stability.

Override options:

* `MICROSEQ_QT_BACKEND=wayland microseq-gui` (opt into native Wayland)
* `MICROSEQ_QT_BACKEND=xcb microseq-gui` (force XWayland)
* `MICROSEQ_QT_BACKEND=offscreen microseq-gui` (headless diagnostics)

The GUI also persists normal geometry and restores maximized state only after first show to reduce Wayland configure/resize timing hazards.

To make this automatic in your MicroSeq conda environment (so you do not export every run), add activate/deactivate hooks:

```bash
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d" "$CONDA_PREFIX/etc/conda/deactivate.d"

cat > "$CONDA_PREFIX/etc/conda/activate.d/microseq_qt_backend.sh" <<'EOF'
export _MICROSEQ_OLD_QT_QPA_PLATFORM="${QT_QPA_PLATFORM-}"
export QT_QPA_PLATFORM=xcb
EOF

cat > "$CONDA_PREFIX/etc/conda/deactivate.d/microseq_qt_backend.sh" <<'EOF'
if [ -n "${_MICROSEQ_OLD_QT_QPA_PLATFORM+x}" ]; then
  export QT_QPA_PLATFORM="$_MICROSEQ_OLD_QT_QPA_PLATFORM"
else
  unset QT_QPA_PLATFORM
fi
unset _MICROSEQ_OLD_QT_QPA_PLATFORM
EOF
```

Then run `conda deactivate && conda activate MicroSeq` once to apply.

## Logging vs Run Outputs

MicroSeq writes **runtime logs** and **pipeline artifacts** to different places:

* Logs: `logs/microseq_<session>.log` (and `logs/microseq_latest.log` symlink)
* Run outputs: your selected output/work directory (for example `*_microseq/`) containing artifacts such as `qc/`, `asm/`, `reads.fasta`, `hits.tsv`, and `hits_tax.tsv`
* Canonical output reference: [Output Artifacts Reference](output_artifacts.md)

Examples:

* CAP3 subprocess progress emitted through Python logging appears in `logs/microseq_<session>.log`.
* CAP3 per-sample captured streams remain in run artifacts like `asm/<sample>/cap3.stdout.txt` and `asm/<sample>/cap3.stderr.txt` when produced by compare/selection paths.

### Abbreviated paired single-run output tree

For a complete “what each file means” table, use [Output Artifacts Reference](output_artifacts.md). The shortened tree below is a quick GUI-oriented landmark map:

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

## Select Input 
* First you click Browse...and pick in this case a folder with ab1 files. In order to access folder selection it's impertative you click cancel first then choose folder will pop up for it. Then you can go down into the nested folders until you find the desired folder which directly in that file houses the .ab1 files you want to blast against. You can also type in the path in the search bar as well. 
* During execution the status bar shown the current stage and ETA; log pane streams MicroSeq's standard logging I have embedded in it.


## Configuring The Run
* In a future update I will have all the dials adjustable for now they are the default settings so keep that in mind.
* Use the **Assembly** dropdown to choose **Single** (default) or **Paired** runs. Paired mode enables forward/reverse regex fields plus well enforcement, duplicate policy controls, and assembler selection controls.
* If you want the BIOM file you will need to tick the BIOM box when clicking full pipeline run note you will need to supply a metadata table (GUI will warn if missing).
* Choose Fast(megablast) for routine Sanger reads; Comprehensive (blastn) for divergent amplicon reads. Note I do have it set up where if the coverage or percent ID doesn't reach 90% blastn will rerun instead given the nature of megablast algorithm.(Fast = NCBI Megablast large 28-bp word, tuned for ≥95 % identity hits between nearly identical Sanger traces; Comprehensive = standard blastn 11-bp word, slower but ∼5× more sensitive to <90 % identity or divergent/gapped matches; MicroSeq auto-falls back to blastn whenever the initial Megablast hit covers <90 % of the read or shows <90 % identity).
* Primer token detection is available via the primer-set dropdown. Toggle **Advanced regex** when you want to enter explicit forward/reverse regex patterns instead of presets.
* Primer **pairing labels** and primer **sequence trimming** are separate controls. Pairing labels (27F/1492R tokens) only determine forward/reverse matching; they do not edit read bases.
* Primer trimming now supports three modes: **Off**, **Detect** (report hits, no clipping), and **Clip** (trim matched primer sequence). Stage can be set to **Pre-quality** or **Post-quality**.
* Recommended default for standard 16S colony-PCR + Sanger workflows (PCR with 27F/1492R and sequencing with the same primers): keep primer trim mode **Off** and rely on **quality trimming** as the primary cleanup step.
* Use **Detect** to monitor potential primer-like hits without changing sequences; use **Clip** only when a known synthetic flanking sequence is expected to appear among the base-called bases (for example, vector backbone, PCR 5' overhangs/adapters/barcodes, or the opposite primer site in short amplicons), or when troubleshooting confirms clipping is beneficial.
* You can now paste **custom forward/reverse primer sequences** directly in the GUI (one per line). This is independent from pairing token presets and overrides config-based primer lists for the run.

* **Primer preview** now uses the same GUI primer stage/**Known synthetic flank preset**/custom-sequence settings as pipeline execution, but forces detect-only mode so you can validate hits before clipping.
* **Assembler selection** (paired mode) supports: **CAP3 default (legacy paired pipeline)**, **All assemblers (compare + pick best contig)**, or any single registered assembler from the dropdown.
* Compare assemblers produces `asm/compare_assemblers.tsv`, now shown in the GUI **Compare Assemblers** output tab for side-by-side engine review. When running full pipeline with **All assemblers**, MicroSeq reuses this comparison and ranks candidates by status (**assembled** > **merged** > others), then length, then deterministic tiebreak for downstream BLAST payloads.
* In compare-driven **Selected/All assemblers** full-pipeline mode, MicroSeq now writes `qc/pairing_report.tsv`; when a compared backend is CAP3, the **Use per-base quality scores** toggle is also applied to CAP3 input generation.
* **Mode semantics (authoritative):**
  * **CAP3 default (legacy paired pipeline):** runs `merge_two_reads` first (configured overlap engine/strategy), then falls back to CAP3 when merge is not accepted (subject to quality-policy exceptions).
  * **Single selected assembler:** runs compare flow constrained to the selected `assembler_id` only; no automatic global CAP3 rescue is injected just because that backend produced no payload.
  * **All assemblers:** runs every registered backend, writes `asm/compare_assemblers.tsv`, and selects a per-sample winner deterministically for BLAST payload generation.
* **`merge_two_reads` overlap engines (what they are):**
  * `merge_two_reads:ungapped` = placement-driven, end-anchored **ungapped sliding** overlap (fast; does not introduce indels in the overlap candidate itself).
  * `merge_two_reads:biopython` = Biopython pairwise-alignment backend that can model **gapped** overlaps (indel-aware).
  * `merge_two_reads:edlib` = edlib edit-distance/path backend that can model **gapped** overlaps (indel-aware).
  These are separate engines under the same `merge_two_reads` method (not aliases of each other). See also: `docs/workflow_resolution.md` section **“Working examples: ungapped sliding vs gapped alignment backends.”**
* Primer policy in paired mode is explicit in reports: mode (`off|detect|clip`), stage (`pre_quality|post_quality`), and source (`off|preset|custom|mixed`).
* In trim controls, **Known synthetic flank preset** applies only to Detect/Clip sequence scanning; pairing-token presets such as **27F/1492R** remain separate and only affect forward/reverse pairing.
* **Duplicate policy** lets you decide whether duplicate orientations error, keep first/last, merge, or keep separate (runs CAP3 per duplicate). The last selection is persisted between sessions.
* **Enforce well codes** requires A1-H12 plate tokens in filenames; files with missing or mismatched wells are skipped but reported.

### Primer trimming: what it means for AB1/Sanger reads
In dye-terminator Sanger, the sequencing primer itself is not represented as base-called sequence; called bases begin downstream of the primer binding site. In most Sanger runs, cleanup is therefore primarily **quality trimming** (low-Q leading/trailing bases), not clipping a sequencing primer sequence.

MicroSeq primer trimming is intended for cases where a known synthetic flanking sequence is actually present in called bases (for example, vector backbone before an insert, PCR 5' overhangs/adapters/barcodes, or the opposite primer site when reads span short amplicons).

Decision guide:
* **Off**: Default for standard 16S colony-PCR + Sanger using 27F/1492R as sequencing primers.
* **Detect**: Diagnostic mode to report primer/vector-like hits without modifying reads.
* **Clip**: Use only when a specific synthetic flank is expected in called bases and removal is confirmed to improve downstream results.

Stage guide:
* Prefer **Post-quality** to reduce false matches from noisy leading bases.
* Use **Pre-quality** only when you expect a high-confidence synthetic prefix that should be removed before quality heuristics.

Examples:
* **Example A (standard full-length 16S isolate)**: PCR with 27F/1492R, sequence with 27F for forward and 1492R for reverse. Use primer trim **Off** and keep quality trim enabled.
* **Example B (short amplicon reaches opposite end)**: A ~250-600 bp amplicon where reads can extend into opposite-end primer/overhang sequence. Use primer trim **Clip** with custom sequences, usually at **Post-quality**.
* **Example C (plasmid/clone sequencing)**: Sequencing with a vector primer where called bases begin in vector backbone before entering insert. Use **Detect** first, then **Clip** for the confirmed vector flank.

## Execution What do the buttons do??

* Button Stages Performed 
* **Run QC** Ab1 -> Fastq -> trim -> pass/fail 
* **Run Blast** BLAST only existing fasta/contigs 
* **Full Pipeline** (BIOM unchecked) QC + BLAST + Taxonomy (Note:Assembly for full length will be added in an update so I may just include a checkbox for partial 16s thenclicking FullPipeline as a workaround.)
* **Assembly** Run CAP3 in single or paired mode with the selected duplicate policy if in paired mode, primer patterns, and optional well enforcement. 
* **Preview-Pairs** (Paired Mode) shows which files mapped to Forward/Reverse, highlights missing mates, and lists suggested regex patterns to try before running CAP3. 
* **Post-BLAST** Collapse hits to best one + merge metadata to BIOM file that is outputted 

During execution the status bar shown the current stage and ETA; log pane streams MicroSeq's standard logging I have embedded in it. 

## Primer Pairing Control Panel real-world GUI workflow examples to consider 

Think of the paired-assembly widgets as a control panel that flips a handful of switches and routes that state into how tokens are shown, how regexes are built, how files are scanned, and how the pairing preview behaves. The table below summarises the switches before walking through concrete scenarios you can replicate with real folders of AB1/FASTQ files.

| Switch | Where to set it | Behaviour |
| --- | --- | --- |
| Assembly mode | **Assembly** dropdown (**Single** vs **Paired**) | Paired mode reveals primer/pairing controls and enables the **Preview pairs** flow; single mode hides paired-only controls. |
| Assembler mode | **Assembler selection** dropdown (paired mode) | Choose one of three routes: **CAP3 default** (merge_two_reads first, then CAP3 fallback when applicable), **Single selected assembler** (run only that backend; no automatic CAP3 rescue), or **All assemblers** (run all and deterministically select winner for BLAST inputs). |
| Primer set preset | **Primer set** dropdown (e.g., `16S (27F/1492R)`,`Custom`) | Prefills forward/reverse token fields and persists the choice via QSettings. Changing tokens by hand flips the preset to **Custom**. |
| Token vs regex override | **Advanced regex** checkbox | Unchecked uses token fields (`27F, 515F` -> `27F|515F`); checked exposes regex fields that bypass tokens entirely so you can enter expressions like `(?i)(27f|8f)` . |
| Enforce same well | **Enforce well codes** checkbox | When enabled, pairing requires both directions to share the same plate position (A1-H12 by default). |
| Auto-detected vs manual tokens | **Detect tokens** button vs typing vs presets | Detection scans filenames for `F`/`R` tokens (e.g., `10__27F.ab1`) and fills the token fields; manual edits and presets drive pairing if you do not run detection. |

### Scenario 1 | preset primers, paired mode, quick preview
**Goal:** You have a directory such as `plates/run42/` that already uses standard 16S primer names in filenames (`10__27F.ab1`, `10__1492R.ab1`, …) and you want to check pairing before launching CAP3.

1. Switch **Assembly** to **Paired**. The advanced regex checkbox appears but remains off.
2. Keep the default primer preset **“16S (27F/1492R)”**. Forward tokens show something like `27F, 8F, 515F`; reverse tokens show `1492R, 806R`.
3. Click **Preview pairs**. The GUI builds regexes from the tokens (e.g., `27F|8F|515F`) and scans your chosen folder. If **Enforce well codes** is ticked, only matches with the same well (e.g., `B10__27F.ab1` + `B10__1492R.ab1`) are paired. The dialog lists found pairs, missing mates, and suggested regex tweaks if anything is skipped.

### Scenario 2 | manual token tweaks that become “Custom”
**Goal:** You want minor tweaks (e.g., add `27Fmod`) without switching to regex overrides.

1. Stay in **Paired** mode and edit the forward token field to `27F, 27Fmod`.
2. The preset automatically flips to **Custom** so QSettings saves your edits.
3. Run **Preview pairs**. Pairing now uses `27F|27Fmod` for the forward side while keeping the existing reverse tokens. This is a quick way to reflect lab-specific primer nicknames without rewriting regexes.

### Scenario 3 | advanced regex for mixed naming conventions
**Goal:** Your filenames mix multiple primer conventions (e.g., `V1F`, `V2F`, `1492R`, `926R`) and tokens are too limiting.

1. Tick **Advanced regex** in paired mode. Regex fields become visible and editable.
2. Enter explicit patterns, for example `(?i)(\d{2}F|V[12]F)` for forward and `(?i)(1492R|926R|806R)` for reverse.
3. **Preview pairs** now routes those raw regexes into the detector; token fields are ignored while the checkbox is on. Untick the box anytime to return to token-based patterns.

### Scenario 4 | detect tokens from real filenames, then preview
**Goal:** You dropped 96 AB1 files into `runs/plateA/` and want MicroSeq to infer likely tokens instead of typing them.

1. Browse to the folder (or any file inside it) so `_infile` is set, then click **Detect tokens**.
2. The GUI scans `.fasta`, `.fastq`, and `.ab1` names for chunks ending in `F` or `R` (e.g., `27F`, `341F`, `1492R`). Top candidates are filled into the token fields and the preset flips to **Custom**.
3. Hit **Preview pairs** to see how those inferred tokens pair the files. This is ideal when collaborators send you plate exports with unfamiliar primer labels.

### Scenario 5 | compare well enforcement vs cross-well pairing
**Goal:** Your filenames contain wells (`B10__27F.ab1`, `H03__1492R.ab1`) and you want to see how strict well matching changes the preview.

1. In paired mode, set your tokens or regex (preset/detected/manual all work).
2. Leave **Enforce well codes** unchecked and click **Preview pairs** to view loose pairing (any sample ID match counts, regardless of well).
3. Tick **Enforce well codes** and preview again. Only pairs that share the same well survive; the dialog calls out cross-well mismatches so you can fix naming before running CAP3.


## Example of a Pilot Run Using Data From Quintara in single mode located in `test/reverse_orientation_only_ab1_demo_run`  
I click the browse option on the left side of MicroSeq and I go inside of MicroSeq to the test directory and search for the reverse_orientation_only_ab1_demo_run directory. 
This folder specifically that has the ab1 files I need directly inside it and no other folders. I exectued in clicking Full Pipeline it will go through the stages described above.
The run will go through a percentage bar that once it gets to 100% you know it will be finished and it also shows the eta and % of the run on the bottom left corner of the GUI. The output of the run will be located in the same parent directory or folder it was in with the same name with 'microseq' appended to it so you know MicroSeq was being used for the run and the results are in that specific directory. 

I will go through how the tree structure looks and the run...

```
reverse_orientation_only_ab1_demo_run_microseq/
├── raw_ab1/              # original chromatograms
├── raw_fastq/            # AB1->FASTQ before trimming
├── passed_qc_fastq/
├── failed_qc_fastq/
├── qc/
│   ├── *_avg_qual.txt    # Individual output metric avg Q for each sequence 
│   ├── trim_summary.tsv  # Summary output metrics for avg Q, length, and primer-trim telemetry for each file
│   └── trimmed.fasta     # Placed here for providence - same file as reads.fasta
├── reads.fasta           # merged QC-passed reads
├── hits.tsv              # BLAST results
└── hits_tax.tsv          # taxonomy-joined results
```
> Note: duplicate FASTA files > `qc/trimmed.fasta` and `reads.fasta` contain the exact same set of QC-passed reads. The first is emitted by the Trim step and kept inside qc/ for provenance; the second is copied to the top-level so downstream tools (assembly, BLAST, etc.) can find a single “canonical” FASTA without digging into sub-folders.

* raw_ab1/ inputs retained for provenance (copies of original ab1 files) will include option to symlink later. 
* raw_fastq/ FASTQ files produced from the ab1 files before trimming 
* qc/ per-read stats and trimmed FASTA quality‑control summary files:
  - *_avg_qual.txt is the per‑read length and Phred average value.
  - trim_summary.tsv is the combined metrics file for all reads, with core quality/length metrics and primer-trim telemetry columns.
  - trimmed.fasta is the final FASTA used for assembly or BLAST that you can manually hand off to NCBI website for example.
  - Merge confidence is primarily determined by overlap-local metrics (length/identity/quality policy). 
  - CLI users can enforce retention with `--min-reads-kept` to fail files that lose too many reads during QC. 
* passed_qc_fastq/ quality filter outcome with default at (≥ Q20); these are individual FASTQ files that passed that metric. 
* failed_qc_fastq/ quality filter outcome with default at (≥ Q20); these are individual files that failed that metric whose average quality was below the threshold. 
* reads.fasta the merged FASTA of all QC‑passed reads (made from passed_qc_fastq/). Created during the pipeline when FASTQs are converted to a single FASTA that is used for blasting by MicroSeq. 
* hits.tsv `microseq blast` output 
* hits_tax.tsv This file has the taxonomy needed to interpret the blast results by adding the taxonomy column using `microseq add_taxonomy`

# QC telemetry vs merge decisions in MicroSeq

MicroSeq trims reads first (sliding‑window or Trimmomatic). `trim_summary.tsv` records quality/length metrics and primer-trim telemetry, while merge decisions are based on overlap-local evidence.

### What MicroSeq reports

`trim_summary.tsv` includes per-file read count, average length, average Phred quality, QC status, and primer-trim telemetry.

Primer-trim telemetry columns (when primer trim is enabled):

* `primer_bases_trimmed` (legacy, avg trim length per trimmed read).
* `primer_trim_len_avg` (explicit avg trim length).
* `primer_bases_trimmed_total` (total bases removed by primer trimming).
* `primer_mismatch_avg`, `primer_hit`, `primer_offset_avg`.
* `primer_orientation`, `primer_orientation_source`, `primer_iupac_mode`.

Primer detect mode writes `qc/primer_detect_report.tsv` with per-file hit-rate telemetry (`reads_scanned`, `reads_with_primer_hit`, `hit_rate`, mismatch/offset averages, and matched primers).

### Primer-trim config migration

Legacy config used `primer_trim.enabled: true|false`. MicroSeq now normalises this to `primer_trim.mode` automatically:

* `enabled: false` -> `mode: off`
* `enabled: true` -> `mode: clip`

Supported values are `mode: off|detect|clip` and `stage: pre_quality|post_quality`.

If an older config references `primer_trim.preset: 16S_27F_1492R`, MicroSeq now migrates safely at runtime: the removed preset is cleared, and users should provide custom synthetic flank sequences (or choose one of the example synthetic-flank presets) for Detect/Clip.

### Primer trim usage examples

* **Example A (standard 16S amplicon):** set primer trim mode to **Off**. Pairing tokens like `27F/1492R` still work for paired assembly matching.
* **Example B (known synthetic construct):** set mode to **Detect** with either custom synthetic flank sequences or a **Known synthetic flank preset**.
* **Example C (validated synthetic flank clipping):** set mode to **Clip** and provide validated custom/synthetic flank sequences before production use.

### Why this policy exists

For long Sanger reads, MicroSeq prioritizes overlap-local evidence for merge decisions.

Because assembly trust is an overlap-local decision, MicroSeq prioritizes:

* end-anchored overlap geometry,
* overlap length + identity thresholds, and
* overlap-quality policy (`warning` vs `blocking`). In `blocking` mode, overlap quality must be present and meet threshold; missing QUAL evidence is treated as non-feasible for merge gating.

Overlap-local thresholds and policy determine merge gating.



## Example of a Paired AB1 Run (full pipeline with CAP3 assembly)

Below is a real paired-mode walkthrough using the demo data in `tests/paired_ab1_demo_run/`. It assumes filenames already encode primer tokens (e.g., `27F`/`1492R`) and plate wells, so the pairing preview finds mates without extra regex tweaks.

1. Click **Browse...**, select the `tests/paired_ab1_demo_run/10292025_1080497/` folder (or any file inside it), and switch **Assembly** to **Paired**.
2. Keep the preset **16S (27F/1492R)** tokens, leave **Enforce well codes** off (sample ID matching is enough here), and choose the **GG2** database.
3. Click **Full pipeline**. MicroSeq will stage the AB1s, run QC/trim, pair and assemble with CAP3, then BLAST the contigs against GG2.

Expected outputs are written to `tests/paired_ab1_demo_run/10292025_1080497_microseq/`:

```
10292025_1080497_microseq/
├── raw_ab1/                 # copies of the original chromatograms
├── raw_fastq/               # AB1->FASTQ before trimming
├── qc/
│   ├── *_avg_qual.txt       # per-read length and avg Phred
│   ├── trim_summary.tsv     # combined QC + telemetry metrics across reads
│   ├── pairing_report.tsv   # which forward/reverse files paired (includes well codes if enforced)
│   └── paired_fasta/        # per-read trimmed FASTA files before CAP3
├── passed_qc_fastq/         # individual FASTQs that passed QC (per primer direction)
├── failed_qc_fastq/         # any reads failing QC
├── asm/
│   ├── <sample_id>/
│   │   ├── *_paired.fasta           # CAP3 contigs for that sample
│   │   ├── *_paired.fasta.cap.*     # CAP3 provenance (contigs, qual, links, ace)
│   │   └── *_paired.fasta.cap.singlets # reads that did not assemble into a contig
│   ├── contig_map.tsv        # input->contig traceability (which FASTQ/FASTA formed each contig)
│   └── paired_contigs.fasta  # all per-sample contigs concatenated
├── hits.tsv                  # BLAST results on assembled contigs
├── hits_tax.tsv              # taxonomy-joined results
└── reads.fasta               # canonical FASTA of QC-passed reads (pre-assembly)
```

Key files to highlight in the GUI doc:

* `qc/pairing_report.tsv` shows which forward/reverse files paired, along with any missing mates or well conflicts.
* `asm/*_paired.fasta.cap.singlets` is expected from CAP3; it lists input reads that could not be merged into a contig at the overlap/identity thresholds.
* `asm/contig_map.tsv` links each contig back to the input reads used, which is helpful when auditing assemblies or troubleshooting mismatched pairs.

## Troublehshooting what the different Status sign means in Assembly Summary Tab
For `CAP3 default (legacy paired pipeline)`, MicroSeq runs `merge_two_reads` first using the configured overlap engine/strategy, and if that merge path is not accepted it falls back to CAP3 using the selected CAP3 profile (strict, relaxed, diagnostic), except where quality-policy routing keeps singlets/no-payload by design.

If you get `ambiguous_overlap`, MicroSeq found near-tied top feasible overlap candidates and intentionally refused to force a unique merge. For canonical decision rules (including the **“Candidate generation in 3 stages”** walkthrough, feasibility gates, top-1 vs top-2 tie checks, and policy outcomes like `topk`/`best_guess`), see **Workflow Resolution Funnel -> “Assemble -> Validate handoff: how `ambiguous_overlap` is decided”** in [`docs/workflow_resolution.md`](workflow_resolution.md). 

Status glossary cross-link: for the full trigger-to-status routing matrix (`merged`, `ambiguous_overlap`, `high_conflict`, `quality_low`, `cap3_unverified`) and the exact artifact fields to audit, see **Workflow Resolution Funnel -> “Assemble/validate status routing matrix”** in [`docs/workflow_resolution.md`](workflow_resolution.md).

## Explaning what each of the columns means in the Tabs Assembly Summary, Blast Inputs, Diagnostics, Compare Assemblers 

### Assembly Summary Tab 
`sample_id`: The sample key used throughout paired assembly outputs. 
`status`: Finaly summary status for that sample. In legacy Cap3 mode this can be cap3-derived (`assembled`, `singlets_only`, `cap3_no_output`) and may be overridden by overlap audit `overlap_*` or for example as I suggested earlier `ambiguous_overlap`.
`assembler`: Which assembler path was selected for summary (ie `cap3:relaxed, merge_two_reads:*`) etc, in the selected/all mode, or `cap3/merge_two_reads` in legacy report parsing. 
`contig_len`: Length of selected/produced contig payload (max length here when multple records) 
`blast_payload`: This is what is being passed to BLAST (contig, singlet, no_payload, pair_missing depending on mode/results). 
`selected_engine`: Engine that produced the chosen result for example overlap engine for `merge_two_reads` rows, `cap3` for cap3 results etc. 
`configured_engine`: Configured overlap engine (legacy parse path populates; selected/all summary currently leaves blank will adjust this later on). 
`merge_status`: Merge prepass status (`merged`, `identity_low`, `overlap_too_short`, `quality_low`, `ambiguous_overlap`, `not_end_anchored`, `high_conflict`)
`merge_overlap_len`: Overlap length measured by merge pre-pass (legacy parse path).
`merge_identity`: Overlap identity from merge pre-pass. 
`overlaps_saved / overlaps_removed`: Parsed from CAP3 `.cap.info` overlap counters. 
`primer_mode / primer_stage`: Run time primer trimming policy (off/detect/clip and pre/post quality stage). 

### BLAST Inputs Tab 
`sample_id`: Sample key for payload row. 
`blast_payload`: Type of sequence payload sent to BLAST(`contig`, `singlet`, `no_payload`, `pair_missing`)
`reason`: Why that payload choice happened. Legacy examples `contigs_present`, `singlets_only`, `cap3_no_output`, `pair_missing`. Selected/all compare examples: `selected_payload`, `winner_no_payload`, `selected_backend_no_payload`, `pair_missing`. 
`payload_ids`: Mapping of rewritten FASTA ids to original ids, eg `sample|contig|cap3_relaxed1=contig1`

### Diagnostic Tab 
`sample_id`: Sample key audited. 
`overlap_len`: Chosen overlap length. 
`overlap_identity`: Identity of selected overlap candidate
`overlap_quality`: Mean overlap quality 
`orientation`: Selected overlap orientation (`revcomp` / `forward`)
`status`: Overlap audit status (`ok`, `overlap_too_short`, `overlap_identity_low`, `overlap_quality_low`, `ambiguous_overlap`, `not_end_anchored`)
`best_identity`: Best identity value seen among condidates used in diagnostic context 
`best_identity_orientation`: Orientation associated with best identity candidate. 
`anchoring_feasible`: Whether at least one candidate can satify end anchor + thresholds. 
`end_anchored_possible`: Whether any candidate is end anchored at all. 
`fwd_best_identity` / `revcomp_best_identity`: Best identity by orientation within feasibility-space (`overlap_len >= min_overlap`).
`fwd_best_overlap_len` / `revcomp_best_overlap_len`: Overlap length paired to the feasible orientation best-identity values.
`fwd_anchor_feasible` / `revcomp_anchor_feasible`: Whether each orientation has any end-anchored candidate that passes thresholds.
`identity_delta_revcomp_minus_fwd`: `revcomp_best_identity - fwd_best_identity` (feasible-space).
`selected_vs_best_identity_delta`: Best-identity minus selected overlap identity for the chosen engine.
`top_candidate_count`: Number of near-tied top feasible candidates using ambiguity epsilons.
`top2_identity_delta` / `top2_overlap_len_delta` / `top2_quality_delta`: Absolute deltas between top-1 and top-2 feasible candidates.
`tie_reason_code`: Why a near tie was flagged (`len_equal`, `mismatch_equal`, `identity_eps`, `quality_eps`).
`fwd_best_identity_any` / `revcomp_best_identity_any`: Best identity by orientation without `min_overlap` filtering (raw evidence-space).
`fwd_best_overlap_len_any` / `revcomp_best_overlap_len_any`: Overlap lengths paired to the raw evidence-space best-identity values.
`pretrim_best_identity` / `pretrim_best_overlap_len` / `pretrim_status`: Best-overlap diagnostics recomputed from pre-primer-trim paired FASTA.
`posttrim_best_identity` / `posttrim_selected_overlap_len` / `posttrim_status`: Post-trim best identity plus the selected overlap length and status in the final run.
`ambiguity_identity_delta_used` / `ambiguity_quality_epsilon_used`: Effective thresholds used to decide ambiguous/near-tie behavior.
`primer_trim_bases_fwd` / `primer_trim_bases_rev`: Total primer-trimmed bases per orientation from `qc/primer_trim_report.tsv`.
`selected_engine`: Overlap engine selected for audit row.
`fallback_used`: yes/no if cascade/all strategy had to move off first engine.
`overlap_engine` (UI label; file writes configured_engine): Configured default overlap engine from config (auto, ungapped, etc.).

### Compare Assemblers Tab 
- Columns shown:
- sample_id, assembler_id, assembler_name, status, selected_engine, contig_len, warnings, diag_code_for_machine, diag_detail_for_human, cap3_contigs_n, cap3_singlets_n, warnings.
- Additional compare TSV-only fields (not shown as table columns) and used by compare-row details/open actions:
  - cap3_info_path
  - cap3_stdout_path (legacy/backward-compatible)
  - cap3_stderr_path (legacy/backward-compatible)
  - tool_name
  - tool_stdout_path
  - tool_stderr_path
  - selection_trace_path
  - winner_reason
  - payload_fasta
- Underlying compare TSV also includes dup_policy and payload_fasta (not displayed in current GUI table).
- Meanings
    - sample_id
        - Sample being evaluated by each backend.
    - assembler_id
        - Internal id like merge_two_reads:biopython, cap3:strict, etc.
    - assembler_name
        - Human-readable label shown in dropdown/table.
    - status
        - Per-backend run outcome:
            - merge backend: merged, identity_low, overlap_too_short, quality_low, ambiguous_overlap, not_end_anchored, high_conflict, etc.
            - CAP3 backend: assembled or cap3_no_output; error on exception.
    - selected_engine
        - Merge overlap engine used (merge rows) or cap3 (CAP3 rows).
    - contig_len
        - Max contig length produced by that backend row (blank if no payload).
    - warnings
        - Merge warning text or caught exception text for failed run.
     - dup_policy
        - Duplicate handling policy used by that backend for this sample (for example `keep_first`, `keep_longest`, etc.).
    - diag_code_for_machine
        - Compact machine-readable diagnosis code for the row result. Examples include `pair_missing`, `merge_identity_low`, `cap3_contigs_present`, `cap3_singlets_only`, `cap3_nonzero_exit_with_output`, `cap3_nonzero_exit_no_output`, and `exception`.
    - diag_detail_for_human
        - Human-readable diagnostic message matching the machine code above (often includes return code/profile/input filename and contig/singlet counts for CAP3 rows).
    - cap3_contigs_n
        - Number of CAP3 contig records detected for that compare row.
    - cap3_singlets_n
        - Number of CAP3 singlet records detected for that compare row.
    - cap3_info_path
        - Path to CAP3 `.cap.info` artifact for that compare row (when present).
    - cap3_stdout_path
        - Legacy CAP3-specific stdout path retained for compatibility with historical compare TSVs.
    - cap3_stderr_path
        - Legacy CAP3-specific stderr path retained for compatibility with historical compare TSVs.
    - tool_name
        - Backend/tool label used to generate row-scoped process diagnostics for this compare row.
    - tool_stdout_path
        - Generic stdout/diagnostic artifact path for the backend selected by this row (CAP3 and non-CAP3 backends).
    - tool_stderr_path
        - Generic stderr artifact path for the backend selected by this row (empty/diagnostic placeholder for pure-Python backends).
    - selection_trace_path
        - Path to per-sample winner-selection trace artifact written under `asm/selection_trace/<sample_id>.selection_trace.tsv`.
        - Contains every candidate row considered for that sample, derived ranking keys (`success_rank`, `payload_rank`, ambiguity penalty, normalized length, tie-break key), final sorted order, winner flag, and winner reason.
    - winner_reason
        - Concise human-readable explanation of why the winner row was selected (populated on the selected/winner row; blank for non-winner rows).
    - payload_fasta
        - Path to the backend payload FASTA produced by that row; populated when a payload exists and used by compare-row Details/open-file actions.

- Compare tab Details panel is now row-scoped by `(sample_id, assembler_id)` instead of only `sample_id`, so selecting two rows for the same sample but different backends shows backend-specific reasons/artifacts instead of collapsing them together.

- In all/selected modes, winners are chosen per sample with deterministic ranking keys (in order):
    - success tier (`assembled` and `merged` are equivalent success outcomes and appear as green/success states in compare status semantics),
    - payload rank (`contig` > `contig_alt` > `singlet` > `none`),
    - ambiguity penalty (`0` before `1`),
    - normalized length (descending),
    - assembler id tie-break (ascending).
- Example tie-break: if two rows for one sample both have payload_kind=`contig`, ambiguity_flag=`0`, and payload_max_len=`1200`, the row with lexicographically smaller `assembler_id` wins.
- That winner selection affects what appears in Assembly Summary/BLAST Inputs for that sample.


### What does anchoring feasability even mean? 
Lets cover some core terms first. 

Overlap: This means I have a forward (F) and reverse (R) reads from the same sample. The assembler tries to find where they cover the same DNA region. That shared region is the overlap. 

Candidate: A candidate is one possible overlap alignment the algorithm found. There can be mutlple alignments because:
* alignment can shift slightly
* orentation may differ (forward versus reversus complimentation {revcomp})
* different engines may score differently. 

The pipeline will then score/filter these candidates and picks one that is best or marks is ambiguous because the overlap would be good from a top 2. 

End anchored: This means the overlap is not just somewhere in the middle, so it behaves like a read merge should where one read's end meets the other read. 
If an overlap exists but doesn't satisfy this boundary behavior then the status will be marked as `not_end_anchored` to denote that merge quality constraint. 

Geometrically 
This refers to the shape/placement of the alignment so think here the position of the overall alignment not just percentage identity. So geometric feasability refers to the alignment being biologically plausible end overlap not just an arbitrary internal match. 

### What metrics are being used to assess alignment?
From config defaults (config/config.yaml) and overlap logic:

Minimum overlap length: min_overlap = 100

Minimum identity: min_identity = 0.8 (80%)

Minimum overlap quality: min_quality = 20.0 (mean Phred in overlap, when quality is available which is the default)

Quality policy: warning or blocking (blocking enforces quality threshold if quality is expected)

So a candidate is considered “good/feasible” when it clears the required gates (length + identity + geometry, plus quality).

Putting it together 
end_anchored_possible = yes
then at least one candidate had valid end anchor geometry alignment.

anchoring_feasible = yes
then at least one candidate was end anchored and met threshold checks (length/identity/(quality policy)).

ambiguous_overlap
then multiple top candidates are effectively tied, so selector refuses to choose one confidently.

### Going into more detail on the algorithms MicroSeq uses for assemblies 

MicroSeq’s Assembly step has two distinct problem settings:

Single-read assembly: assemble/clean one-direction reads (or read sets) into contigs/singlets.

Paired-read assembly: assemble a forward (F) and reverse (R) read from the same sample into a single contig (or report why that failed).

In the GUI, paired mode additionally exposes an Assembler selection that routes the paired workflow into one of multiple backends (legacy CAP3 path, a user-selected backend, or “run-all + deterministically pick best”).

CAP3 assembler (Overlap–Layout–Consensus, quality-aware)

What it is
CAP3 is a classical OLC (Overlap–Layout–Consensus) assembler designed for Sanger-style reads. Its core workflow is:

Preprocessing / clipping: CAP3 can clip low-quality regions at the ends of reads.

Overlap detection: it computes overlaps between reads (optionally using base quality values in scoring).

Layout + alignment: it constructs read layouts and multiple sequence alignments for contig construction.

Consensus calling: it generates consensus sequences and uses quality values to reduce consensus errors.

Forward–reverse constraints: CAP3 can use direction/orientation constraints to correct errors and link contigs.

Why it matters in MicroSeq

CAP3 is the “legacy paired pipeline” default and the canonical contig builder for paired reads when you want OLC behavior.

When Use per-base quality scores is enabled, MicroSeq can feed quality information into CAP3’s overlap/consensus steps (CAP3 explicitly supports using base quality values for overlaps, alignments, and consensus generation).

Two-read overlap merge (deterministic F/R merge, “merge_two_reads” backend)

What it is
This is a specialized assembler for the paired-read case: rather than assembling a graph of many reads, it tries to merge exactly two reads (F and reverse-complemented R) into one contig by finding and validating an overlap.

Algorithmically, it’s the same class of method as paired-end read mergers used in short-read pipelines (find an overlap, score candidates, pick the best merge, and emit a merged sequence or a failure mode), even though MicroSeq’s inputs are typically Sanger-derived rather than Illumina.

Core mechanics
Given forward read 𝐹
F and reverse read 𝑅
R:Reverse-complement 𝑅
R so both are in the same orientation.

Generate overlap candidates (multiple possible offsets/orientations/alignments depending on engine).

Apply feasibility gates (the “geometry” and threshold logic described):

End-anchored geometry: overlap must behave like a plausible end-to-end merge (read ends meet as expected for forward/reverse Sanger geometry), not an arbitrary internal match.

Minimum overlap length and minimum identity gates.

Quality gate (when quality is available/expected): compute mean Phred in the overlap and require it under a configured policy (warning vs blocking).

Select a best candidate deterministically:

rank by status/feasibility, then by identity/length/quality as configured,

detect near-ties and emit ambiguous_overlap if multiple candidates are effectively equivalent.

Consensus construction in the overlap:

if bases agree then emit base,

if bases disagree then prefer higher-quality base (or mark/handle conflicts per policy),

track high-conflict cases.

Why it matters in MicroSeq

It provides a fast, explicit, interpretable paired-read merge path (with clear failure categories like overlap too short, identity too low, quality too low, ambiguous overlap).

In “compare-driven” operation (Selected / All assemblers), MicroSeq can run this backend side-by-side with CAP3 and pick the best payload deterministically.

*Scientific references for the method class*
Overlap-based merging of paired reads is a standard approach in sequencing pipelines; FLASH and PEAR are widely cited examples of the overlap-then-merge paradigm (candidate overlap scoring + merge selection).

Alignment / overlap engines used to generate candidates

MicroSeq’s paired merging logic relies on an overlap engine to propose/score candidate alignments. Different engines trade off speed, sensitivity to indels, and determinism. The key families are:

Dynamic-programming gapped alignment (Needleman–Wunsch / Smith–Waterman style)

When you allow mismatches and indels, classical dynamic programming (DP) provides an optimal alignment under a scoring model:

Needleman–Wunsch is the canonical global alignment DP formulation.

Smith–Waterman is the canonical local alignment DP formulation.

In practice, MicroSeq can use Biopython’s alignment machinery as an implementation substrate for DP-based pairwise alignment and scoring.

Exact edit-distance alignment (Edlib)

Edlib is a high-performance library for exact sequence alignment using edit distance, optimized for speed while still returning exact results under the edit-distance model. It’s useful when you want deterministic exactness under an indel-aware distance metric (particularly for near-identical overlaps).

Quality scores as a first-class signal (Phred)

When per-base quality is available, MicroSeq’s overlap evaluation and conflict resolution can incorporate Phred quality scores, which represent log-scaled error probabilities:

Q=−10log10 P(error)

Phred base calling and the error-probability interpretation are described in the original Genome Research papers.

This matters in two places:

Overlap feasibility: mean Phred in the overlap can be used as a gate or warning signal.

Consensus resolution: when F and R disagree, quality-guided selection reduces error propagation into the merged contig (and CAP3 itself is explicitly quality-aware for overlap scoring and consensus).

References used here in my thinking. 

Huang X, Madan A. CAP3: A DNA sequence assembly program. Genome Research (1999).

Magoč T, Salzberg SL. FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics (2011).

Zhang J, Kobert K, Flouri T, Stamatakis A. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics (2014).

Šošić M, Šikić M. Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. Bioinformatics (2017).

Needleman SB, Wunsch CD. A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol (1970).

Smith TF, Waterman MS. Identification of common molecular subsequences. J Mol Biol (1981).

Cock PJA et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics (2009).

Ewing B, Hillier L, Wendl MC, Green P. Base-calling of automated sequencer traces using phred. I. Accuracy assessment. Genome Research (1998).

Ewing B, Green P. Base-calling of automated sequencer traces using phred. II. Error probabilities. Genome Research (1998)

<<<<<<< ours


- Process log contract: compare rows now write deterministic per-row process logs named with both sample + assembler labels (collision-safe suffixing if needed). In the GUI Details panel for Compare rows, use **Open tool stdout/stderr** to open backend process artifacts and **Open selection trace** to inspect winner-ranking decisions; for older runs the GUI falls back to `cap3_stdout_path` / `cap3_stderr_path`.

=======
>>>>>>> theirs

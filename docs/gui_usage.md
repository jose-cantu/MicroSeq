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



## Select Input 
* First you click Browse...and pick in this case a folder with ab1 files. In order to access folder selection it's impertative you click cancel first then choose folder will pop up for it. Then you can go down into the nested folders until you find the desired folder which directly in that file houses the .ab1 files you want to blast against. You can also type in the path in the search bar as well. 
* During execution the status bar shown the current stage and ETA; log pane streams MicroSeq's standard logging I have embedded in it.


## Configuring The Run
* In a future update I will have all the dials adjustable for now they are the default settings so keep that in mind.
* Use the **Assembly** dropdown to choose **Single** (default) or **Paired** CAP3 runs. Paired mode enables forward/reverse regex fields plus well enforcement and duplicate policy controls.
* If you want the BIOM file you will need to tick the BIOM box when clicking full pipeline run note you will need to supply a metadata table (GUI will warn if missing).
* Choose Fast(megablast) for routine Sanger reads; Comprehensive (blastn) for divergent amplicon reads. Note I do have it set up where if the coverage or percent ID doesn't reach 90% blastn will rerun instead given the nature of megablast algorithm.(Fast = NCBI Megablast—large 28-bp word, tuned for ≥95 % identity hits between nearly identical Sanger traces; Comprehensive = standard blastn—11-bp word, slower but ∼5× more sensitive to <90 % identity or divergent/gapped matches; MicroSeq auto-falls back to blastn whenever the initial Megablast hit covers <90 % of the read or shows <90 % identity).
* Primer token detection is available via the primer-set dropdown. Toggle **Advanced regex** when you want to enter explicit forward/reverse regex patterns instead of presets.
* **Duplicate policy** lets you decide whether duplicate orientations error, keep first/last, merge, or keep separate (runs CAP3 per duplicate). The last selection is persisted between sessions.
* **Enforce well codes** requires A1–H12 plate tokens in filenames; files with missing or mismatched wells are skipped but reported.

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
| Assembly mode | **Assembly** dropdown (**Single** vs **Paired**) | Paired mode reveals primer/pairing controls and enables the **Preview pairs** flow; single mode hides them and runs CAP3 without pairing. |
| Primer set preset | **Primer set** dropdown (e.g., `16S (27F/1492R)`, `ITS (ITS1/ITS4)`, `Custom`) | Prefills forward/reverse token fields and persists the choice via QSettings. Changing tokens by hand flips the preset to **Custom**. |
| Token vs regex override | **Advanced regex** checkbox | Unchecked uses token fields (`27F, 515F` → `27F|515F`); checked exposes regex fields that bypass tokens entirely so you can enter expressions like `(?i)(27f|8f)` . |
| Enforce same well | **Enforce well codes** checkbox | When enabled, pairing requires both directions to share the same plate position (A1–H12 by default). |
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
├── raw_fastq/            # AB1→FASTQ before trimming
├── passed_qc_fastq/
├── failed_qc_fastq/
├── qc/
│   ├── *_avg_qual.txt    # Individual output metric avg Q for each sequence 
│   ├── trim_summary.tsv  # Summary output metrics for avg Q for each file combined
│   └── trimmed.fasta     # Placed here for providence - same file as reads.fasta
├── reads.fasta           # merged QC-passed reads
├── hits.tsv              # BLAST results
└── hits_tax.tsv          # taxonomy-joined results
```
> Note: duplicate FASTA files > `qc/trimmed.fasta` and `reads.fasta` contain the exact same set of QC-passed reads. The first is emitted by the Trim step and kept inside qc/ for provenance; the second is copied to the top-level so downstream tools (assembly, BLAST, etc.) can find a single “canonical” FASTA without digging into sub-folders.

* raw_ab1/ inputs retained for provenance (copies of original ab1 files) will include option to symlink later. 
* raw_fastq/ FASTQ files produced from the ab1 files before trimming 
* qc/ per-read stats and trimmed FASTA quality‑control summary files:
  - *_avg_qual.txt is the per‑read length and Phred averages.
  - trim_summary.tsv is the combined metrics for all files.
  - trimmed.fasta is the final FASTA used for assembly or BLAST that you can manually hand off to NCBI website for example.
* passed_qc_fastq/ quality filter outcome with default at (≥ Q20); these are individual FASTQ files that passed that metric. 
* failed_qc_fastq/ quality filter outcome with default at (≥ Q20); these are individual files that failed that metric whose average quality was below the threshold. 
* reads.fasta the merged FASTA of all QC‑passed reads (made from passed_qc_fastq/). Created during the pipeline when FASTQs are converted to a single FASTA that is used for blasting by MicroSeq. 
* hits.tsv `microseq blast` output 
* hits_tax.tsv This file has the taxonomy needed to interpret the blast results by adding the taxonomy column using `microseq add_taxonomy` 




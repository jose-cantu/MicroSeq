---
layout: page
title: Gui Walkthrough
permalink: /gui/
---

## What this page is for

Use this page when you want to run MicroSeq through the graphical interface on a local workstation.

This page is operational. It focuses on how to launch the GUI, choose inputs, configure a run, interpret the main on-screen controls, and know where to look during and after a run.

For deeper references, use these pages:

* [Output Artifacts Reference](output_artifacts.md) for the canonical file-by-file output table
* [Paired CAP3 Assembly](paired_assembly.md) for paired-mode semantics, BLAST handoff meaning, and CAP3 interpretation
* [Workflow Resolution Funnel](workflow_resolution.md) for canonical status semantics such as `ambiguous_overlap`
* [GUI Table Reference](gui_table_reference.md) for detailed column-by-column meanings in GUI tables
* [MicroSeq design: algorithm and architecture choices for Sanger assembly](design_algorithm_choices.md) for algorithm tradeoffs
* [GUI Color Legend](gui_color_legend.md) for color meanings used in table views

## Launch

Open a terminal, activate the MicroSeq conda environment, and launch the GUI:

```bash
conda activate MicroSeq
microseq-gui
```

When the window opens, the GUI shows the current run controls, including database choice, BLAST controls, assembly mode, primer pairing labels, primer trimming controls, duplicate policy, and well enforcement.

## Platform-specific GUI stability notes

### Linux / Wayland stability mode

If you are on Linux Wayland and see maximize crashes such as:

```text
xdg_surface buffer ... does not match configured state
```

MicroSeq defaults to `QT_QPA_PLATFORM=xcb` for stability.

Override options:

* `MICROSEQ_QT_BACKEND=wayland microseq-gui` to opt into native Wayland
* `MICROSEQ_QT_BACKEND=xcb microseq-gui` to force XWayland
* `MICROSEQ_QT_BACKEND=offscreen microseq-gui` for headless diagnostics

The GUI also persists normal geometry and restores maximized state only after first show to reduce Wayland configure/resize timing hazards.

To make this automatic inside the MicroSeq conda environment, add activate/deactivate hooks:

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

Then run:

```bash
conda deactivate
conda activate MicroSeq
```

## Choose input

Click **Browse...** and select cancel so you can choose the folder you want or files you want to process.

Practical notes:

* For AB1 runs, point the GUI at the folder that directly contains the `.ab1` files you want to use.
* You can browse graphically or type a path directly.
* In paired mode, the GUI stages the selected files into a temporary working area so the paired assembly path sees a consistent input set.
* During execution, the **status bar** shows the current stage and progress information, and the **log pane** streams MicroSeq logging in real time.

## Configure the run

## Main run controls

| Control | What it does |
| --- | --- |
| **DB** | Chooses the reference database used for BLAST and taxonomy, such as `gg2`, `silva`, or `ncbi`. |
| **BLAST mode** | Chooses **Fast (megablast)** or **Comprehensive (blastn)**. |
| **Assembly** | Chooses **Single** or **Paired** mode. Paired mode exposes forward/reverse pairing controls. |
| **Assembler selection** | In paired mode, selects the paired backend route: legacy CAP3 default, a single selected backend, or all assemblers for compare-and-pick behavior. |
| **Primer pairing labels** | Controls how forward and reverse files are recognized from filenames. These labels affect pairing only, not base clipping. |
| **Primer trim mode** | Chooses whether primer-like sequence search is **Off**, **Detect**, or **Clip**. |
| **Primer trim stage** | Chooses whether primer-like sequence search happens **Pre-quality** or **Post-quality**. |
| **Duplicate policy** | Controls what happens when multiple forward or reverse files map to the same sample. |
| **Enforce well codes** | Requires matching plate well labels such as `A1` to `H12` in paired mode. |
| **BIOM** | Includes BIOM output during post-BLAST processing when that output is desired. |

### BLAST controls: adjustable ranges and defaults

The BLAST dials are adjustable in the GUI.

| Dial | Allowed range | Default |
| --- | --- | --- |
| Identity (%) | 50 to 100 | 97 |
| Q-cov (%) | 10 to 100 | 80 |
| Max hits | 1 to 500 | 5 |
| Threads | 1 to 32 | 4 |

### Fast vs Comprehensive

MicroSeq exposes two BLAST algorithms in the GUI:

* **Fast (megablast)** is the normal and default path for routine Sanger work, especially when you expect close matches.
* **Comprehensive (blastn)** is the more sensitive path for divergent amplicons, harder cases, novelty screening, or situations where sensitivity matters more than speed.

MicroSeq can also fall back to `blastn` when the initial megablast result is weak, such as cases where the top hit is under 90% identity or 90% coverage.

### Recommended starting settings

The tables below are **recommended starting settings**, not hard-coded GUI defaults.

#### Full-length 16S Sanger

Use this when you expect near full-length isolate reads and want routine classification.

| Setting | Recommended start |
| --- | --- |
| Algorithm | Fast (megablast) |
| Identity (%) | 97 |
| Q-cov (%) | 85 to 90 |
| Max hits | 3 to 5 |
| Threads | 4 for a typical local run, or match your HPC allocation |

Why:

* Full-length reads can support stricter coverage expectations.
* `97` matches the GUI default and common species-level convention already used in the docs.
* Keeping a few hits preserves auditability without cluttering the output.

#### Partial-length 16S Sanger

Use this when only part of the 16S region is present and total query coverage is expected to be lower.

| Setting | Recommended start |
| --- | --- |
| Algorithm | Fast (megablast) first |
| Identity (%) | 97 |
| Q-cov (%) | 70 to 80 |
| Max hits | 5 to 10 |
| Threads | Same logic as full-length; hardware-dependent |

Why:

* Shorter amplicons usually need relaxed coverage before relaxed identity.
* Shorter regions produce more tied or near-tied hits, so keeping more than one hit is more informative.
* This preserves specificity while acknowledging that full-length query coverage is not expected.

#### Harder or more divergent cases

Use **Comprehensive (blastn)** up front when:

* you expect a divergent amplicon
* you suspect novelty
* the read is short enough that megablast may be too brittle
* you care more about sensitivity than speed

For these cases, use starting settings like:

**Full-length but harder / more divergent case**

| Setting | Recommended start |
| --- | --- |
| Algorithm | Comprehensive (blastn) |
| Identity (%) | 95 to 97 |
| Q-cov (%) | 80 to 90 |
| Max hits | 5 |
| Threads | Hardware-dependent |

**Partial-length and harder / more divergent case**

| Setting | Recommended start |
| --- | --- |
| Algorithm | Comprehensive (blastn) |
| Identity (%) | 95 to 97 |
| Q-cov (%) | 60 to 75 |
| Max hits | 10 |
| Threads | Hardware-dependent |

### Practical guidance

* Lower **Q-cov** before lowering **identity** for partial-length reads.
* Raise identity toward **99** only when you want stricter species-level filtering.
* Use **max hits = 1** only for a strict top-hit production view; it is not the best exploratory default.
* Treat **threads** as a compute dial, not a biology dial.

## Primer trimming in practice

Primer pairing labels and primer trimming are separate controls.

* **Primer pairing labels** tell MicroSeq how to recognize forward and reverse files from filenames.
* **Primer trimming** controls whether MicroSeq searches for and optionally clips known synthetic flank sequence within the called bases.

### What primer trimming means for Sanger reads

In dye-terminator Sanger sequencing, the sequencing primer itself is usually not represented in the called sequence. In most standard Sanger runs, cleanup is therefore mainly **quality trimming**, not clipping of the sequencing primer.

MicroSeq primer trimming is intended for cases where a known synthetic flank sequence is actually present in the called bases, such as:

* vector backbone before an insert
* PCR overhangs, adapters, or barcodes
* opposite-end primer sequence in short amplicons

### Decision guide

| Mode | When to use it |
| --- | --- |
| **Off** | Default for standard 16S colony-PCR + Sanger workflows using primers such as `27F` and `1492R` as sequencing primers. |
| **Detect** | Report primer-like or synthetic-flank hits without changing the reads. |
| **Clip** | Remove a validated synthetic flank sequence when you have evidence that clipping improves downstream interpretation. |

### Stage guide

| Stage | When to use it |
| --- | --- |
| **Post-quality** | Preferred in most cases because it reduces false matches from noisy leading bases. |
| **Pre-quality** | Use only when you expect a strong synthetic prefix that should be removed before quality heuristics are applied. |

### Practical examples

* **Example A: standard full-length 16S isolate**
  * PCR with `27F/1492R`
  * sequence with the same primers
  * recommended setting: primer trim **Off**

* **Example B: short amplicon reaching the opposite end**
  * the read may extend into opposite-end primer or synthetic sequence
  * recommended setting: primer trim **Detect** first, then **Clip** only if confirmed

* **Example C: plasmid or clone sequencing**
  * the called bases begin in vector backbone before entering the insert
  * recommended setting: **Detect** first, then **Clip** the confirmed vector flank

## What each button does

| Button | What it does |
| --- | --- |
| **Run QC** | AB1 -> FASTQ conversion -> quality trimming/QC -> pass/fail outputs |
| **Run Blast** | BLAST only on an existing FASTA or assembly output |
| **Full Pipeline** | QC/trim -> assembly if applicable -> BLAST -> taxonomy join -> optional BIOM/post-BLAST outputs |
| **Assembly** | Runs assembly only, using the selected single or paired mode and the current assembly-related settings |
| **Preview pairs** | In paired mode, shows how files map to forward and reverse, highlights missing mates or well conflicts, and helps diagnose filename-based pairing |
| **Post-BLAST** | Collapses BLAST/taxonomy output into downstream reporting artifacts such as BIOM and related tables |

## Two common GUI recipes

### Recipe 1: routine full-length 16S paired isolate run

Use this for standard paired AB1 isolate sequencing where filenames already encode forward/reverse primer labels.

1. Launch the GUI.
2. Click **Browse...** and select the folder containing the paired `.ab1` files.
3. Set **Assembly** to **Paired**.
4. Keep the common primer pairing labels, such as `27F` and `1492R`, or choose the preset that matches your filenames.
5. Keep primer trim **Off** unless you are deliberately checking for synthetic flank sequence.
6. Start with these BLAST settings:
   * **Algorithm**: Fast (megablast)
   * **Identity**: 97
   * **Q-cov**: 85 to 90
   * **Max hits**: 3 to 5
7. Use **Duplicate policy = error** unless you have a deliberate reason to keep or merge duplicate directions.
8. Turn on **Enforce well codes** only when you are working from plate exports and cross-well swaps are a real risk.
9. Click **Preview pairs** if you want to verify pairing before the run.
10. Click **Full Pipeline**.

What to check first afterward:

* `hits_tax.tsv` for the main taxonomic result
* `qc/pairing_report.tsv` to confirm forward/reverse pairing
* `asm/assembly_summary.tsv` to see whether a contig or singlet was used
* `asm/blast_inputs.tsv` to see what sequence actually went to BLAST

### Recipe 2: partial-length or harder / more divergent case

Use this when the 16S region is shorter, the organism may be divergent, or you want a more sensitive first pass.

1. Launch the GUI and choose the input folder.
2. Use **Single** or **Paired** assembly according to the data you actually have.
3. Keep primer trimming **Off** unless you have a concrete synthetic-flank reason to use **Detect** or **Clip**.
4. Start with these BLAST settings for a shorter routine partial 16S read:
   * **Algorithm**: Fast (megablast) first
   * **Identity**: 97
   * **Q-cov**: 70 to 80
   * **Max hits**: 5 to 10
5. Switch to **Comprehensive (blastn)** up front when:
   * the region is short enough that megablast may be brittle
   * you expect divergence
   * you suspect novelty
   * the first pass is weak and you want higher sensitivity
6. For harder cases, start with:
   * **Algorithm**: Comprehensive (blastn)
   * **Identity**: 95 to 97
   * **Q-cov**: 60 to 75 for shorter regions, or 80–90 for harder full-length cases
   * **Max hits**: 5 to 10
7. Run **Full Pipeline** or the specific stage you need.

Interpretation rule of thumb:

* Lower coverage before lowering identity for shorter reads.
* Keep more than one hit when you expect tied or near-tied partial-length matches.

## Where to look during and after a run

### During the run

Use the GUI itself as the first monitor:

* **Status bar**: current stage and progress
* **Log pane**: streaming runtime messages
* **Preview pairs**: pairing audit before paired assembly
* **Compare Assemblers** tab when relevant: side-by-side backend results for compare-driven paired runs

### After the run

MicroSeq separates **runtime logs** from **run artifacts**:

* **Logs**: `logs/microseq_<session>.log` and `logs/microseq_latest.log`
* **Run outputs**: the selected output/work directory, often ending in `_microseq/`

The first files most users check are:

| File | Why check it first |
| --- | --- |
| `hits_tax.tsv` | Main taxonomy-interpretable result |
| `qc/trim_summary.tsv` | Read quality, length, and primer-trim telemetry |
| `qc/pairing_report.tsv` | Forward/reverse pairing and well behavior in paired mode |
| `asm/assembly_summary.tsv` | Sample-level assembly outcome |
| `asm/blast_inputs.tsv` | What sequence was actually sent to BLAST, and why |

The first GUI table columns most users need are:

| Table | Columns to look at first |
| --- | --- |
| **Assembly Summary** | `sample_id`, `status`, `contig_len`, `blast_payload` |
| **BLAST Inputs** | `sample_id`, `blast_payload`, `reason`, `payload_ids` |
| **Diagnostics** | `sample_id`, `status`, `overlap_len`, `overlap_identity`, `orientation` |
| **Compare Assemblers** | `sample_id`, `assembler_id`, `status`, `contig_len`, `warnings` |

For the full column reference, use [GUI Table Reference](gui_table_reference.md).

For the full output tree and exact meaning of each file, use [Output Artifacts Reference](output_artifacts.md).

## Short troubleshooting and links to deeper docs

| Situation | First thing to check | Read next |
| --- | --- | --- |
| GUI does not launch cleanly on Linux / Wayland | Use the stability notes and `MICROSEQ_QT_BACKEND` settings | This page: **Platform-specific GUI stability notes** |
| Paired run finds no valid pairs | Check primer pairing labels, Preview pairs, and well enforcement | [Paired CAP3 Assembly](paired_assembly.md) |
| A paired sample ends up as `singlets_only` or `cap3_no_output` | Check `asm/assembly_summary.tsv`, `qc/overlap_audit.tsv`, and the sample `.cap.info` file | [Paired CAP3 Assembly](paired_assembly.md) |
| You see `ambiguous_overlap` and want the exact meaning | Do not interpret it from the GUI page alone | [Workflow Resolution Funnel](workflow_resolution.md) |
| You are unsure what a file in the output directory means | Use the canonical output table | [Output Artifacts Reference](output_artifacts.md) |
| You want the exact meaning of GUI table columns | Use the dedicated table reference page | [GUI Table Reference](gui_table_reference.md) |
| You want the color meaning of a row or status cell | Use the GUI color legend | [GUI Color Legend](gui_color_legend.md) |
| You want the algorithm rationale behind CAP3 vs staged overlap behavior | Use the design page | [MicroSeq design: algorithm and architecture choices for Sanger assembly](design_algorithm_choices.md) |

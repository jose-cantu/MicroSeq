---
layout: page
title: Paired CAP3 assembly
permalink: /paired-assembly/
---

> Use this page to understand what paired mode means biologically and how MicroSeq decides what sequence, if any, should move forward to BLAST.
>
> For command examples and wrapper usage, see [CLI Workflows](cli_workflow_tutorial.md).
> For button-by-button GUI use, see [GUI Walkthrough](gui_usage.md).
> For the canonical artifact table, see [Output Artifacts Reference](output_artifacts.md).
> For ambiguity routing and sample-level resolution states, see [Workflow Resolution Funnel](workflow_resolution.md).
> For algorithm and design tradeoffs, see [MicroSeq design: algorithm and architecture choices for Sanger assembly](design_algorithm_choices.md).

## Scope and cross-links

This page is the semantic reference for paired forward/reverse assembly in MicroSeq.
It focuses on the biological meaning of pairing, contigs, singlets, overlap failure,
and BLAST handoff. It does **not** try to teach the wrapper flow or the GUI.

In microbiology terms, paired mode starts from one practical question:
do these forward and reverse reads from the same sample support one defensible
consensus sequence, or do they only support partial evidence?

## Paired-mode semantic model

In paired mode, MicroSeq normalizes the input reads, groups files into forward/reverse
pairs, runs paired assembly logic, and then chooses the sequence output that will
represent that sample downstream.

At the sample level, the main outcomes are:

| Outcome | What it means biologically | What MicroSeq does downstream |
| --- | --- | --- |
| `contig` | Forward and reverse reads support one merged consensus sequence. | Sends the contig to BLAST. |
| `singlet` | A valid pair existed, but assembly did not yield a defended merged contig; at least one standalone sequence remains usable. | Sends singlet sequence(s) to BLAST according to handoff policy. |
| `no_payload` | The sample reached assembly, but no usable sequence output survived for downstream search. | Records the sample in the manifest; sends no sequence to BLAST. |
| `pair_missing` | Only one orientation survived pairing detection or well enforcement, so no true paired assembly was possible. | Records the sample in the manifest; sends no sequence to BLAST. |

For microbiology workflows, the key practical idea is this:
paired mode is not just “run CAP3.” It is a structural decision layer that asks
whether two directional reads support a single sample-level sequence interpretation.

## Pairing rules: tokens/primer label, wells, duplicate policy

MicroSeq pairs files deterministically from filenames.
It strips known forward/reverse tokens/primer labels from the name, resolves a sample ID,
and optionally checks plate well codes. This favors reproducibility over guesswork.

### Core pairing rules

| Concept | Rule | Why it matters |
| --- | --- | --- |
| Sample ID | Forward and reverse files pair only if the remaining sample ID matches after token/primer-label stripping. | Prevents accidental cross-sample pairing. |
| Token/primer-label position | Primer labels can appear at the front, middle, or end of the filename. | Lets labs keep existing naming conventions. |
| Well codes | Wells matter only when well enforcement is enabled. | Useful for plate exports where cross-well swaps are risky. |
| Duplicate policy | Controls what happens when multiple forwards or reverses map to the same sample. | Makes duplicate handling explicit instead of silent. |

### Filename examples

These pair successfully because the sample ID resolves to `KD001` in both files:

* `KD001_27F.ab1` ↔ `KD001_1492R.ab1`
* `27F_KD001_A01.ab1` ↔ `1492R_KD001_A01.ab1`
* `KD001_A01_27F.ab1` ↔ `KD001_A01_1492R.ab1`

These do **not** pair because the resolved sample IDs differ:

* `KD001_27F.ab1` ↔ `KD004_1492R.ab1`

These only pair when well enforcement is **off**:

* `KD001_A01_27F.ab1` ↔ `KD001_B01_1492R.ab1`

### What duplicate policy means

| Policy | Meaning |
| --- | --- |
| `error` | Stop and force the duplicate to be resolved explicitly. |
| `keep-first` | Keep the first matching file per orientation. |
| `keep-last` | Keep the last matching file per orientation. |
| `merge` | Merge duplicate same-orientation inputs before downstream paired assembly. |
| `keep-separate` | Preserve duplicate branches as separate paired hypotheses. |

For most routine isolate work, `error` is the safest default because it prevents
silent mixing of repeated reactions or mislabeled files.

## Assembly -> BLAST handoff semantics

This is the most important downstream contract in paired mode:
`asm/blast_inputs.fasta` is what BLAST actually saw, and `asm/blast_inputs.tsv`
is the audit trail explaining why each sequence was selected.

### Main handoff artifacts

| Artifact | Meaning |
| --- | --- |
| `asm/blast_inputs.fasta` | Final sequence records sent to BLAST. |
| `asm/blast_inputs.tsv` | Per-sample manifest linking sequence choice to assembly outcome. |

### Contig -> singlet fallback order

For each paired sample, MicroSeq chooses the first available sequence output in this order:

1. CAP3 contigs (`*.cap.contigs`)
2. CAP3 singlets (`*.cap.singlets`)
3. No usable sequence output

### FASTA header rewrite

BLAST query IDs are rewritten so the sample and payload type remain visible:

```text
sampleA|contig|cap3_c1
sampleA|singlet|cap3_s1
```

### `blast_inputs.tsv` columns

| Column | Meaning |
| --- | --- |
| `sample_id` | Paired sample key. In `keep-separate` mode this can include branch suffixes such as `_1`, `_2`, and so on. |
| `blast_payload` | One of `contig`, `singlet`, `no_payload`, or `pair_missing`. |
| `payload_ids` | Mapping from rewritten BLAST IDs back to original CAP3 or read IDs. |
| `reason` | Why that output was chosen. |

### `blast_payload` and `reason` taxonomy

| `blast_payload` | `reason` | Interpretation |
| --- | --- | --- |
| `contig` | `contigs_present` | A merged consensus sequence was available and used. |
| `singlet` | `singlets_only` | No contig survived, but at least one standalone sequence remained usable. |
| `no_payload` | `cap3_no_output` | Assembly produced no usable sequence output for BLAST. |
| `pair_missing` | `pair_missing` | No valid forward/reverse pair existed after pairing rules were applied. |

### How `pair_missing` is assigned

MicroSeq uses `qc/pairing_report.tsv` together with the staged paired inputs.
If a sample has only one surviving orientation after token/primer-label detection
and optional well enforcement, it is marked as `pair_missing`. The sample remains visible in
`asm/blast_inputs.tsv`, but no sequence is sent to BLAST.

For microbiology interpretation, this matters because `pair_missing` is not a
taxonomic failure. It is a structural or input-state failure: the sample never had
a valid paired assembly opportunity.

## CAP3 diagnostics and overlap-failure interpretation

When a sample does not produce a contig, the main question is whether the problem
came from the reads themselves, from overlap geometry, or from an assembly gate
that was too strict. The three most useful places to look are:

* `asm/assembly_summary.tsv`
* `qc/overlap_audit.tsv`
* `asm/<sample>/*.cap.info`

### Common paired outcomes

| Status or outcome | Practical meaning | First place to look |
| --- | --- | --- |
| `assembled` | A contig was produced and selected. | `asm/assembly_summary.tsv` |
| `singlets_only` | The sample had usable sequence evidence, but not a defended merged contig. | `*.cap.singlets`, `qc/overlap_audit.tsv` |
| `pair_missing` | One direction was absent after pairing logic. | `qc/pairing_report.tsv` |
| `cap3_no_output` | CAP3 did not produce usable contigs or singlets for handoff. | `*.cap.info` |
| `ambiguous_overlap` | More than one plausible overlap candidate was near-tied, so MicroSeq refused to force a single merge. | `qc/overlap_audit.tsv`, `workflow_resolution.md` |
| `quality_low` | Overlap existed, but quality evidence did not support a safe merge. | `qc/overlap_audit.tsv`, QUAL files |

### CAP3 overlap-removal gates

CAP3 overlap acceptance is a conjunction of multiple filters, so it helps to know
which gate actually failed before relaxing parameters.

* **Clipping range failures (`-y`)**: CAP3 can miss a real overlap if clipping removed too much useful end sequence.
* **Difference score (`-b/-d`)**: Too many quality-weighted mismatches can cause CAP3 to reject the overlap.
* **Difference cap (`-e`)**: CAP3 rejects overlaps with more differences than the allowed mismatch tolerance.
* **Similarity score (`-s`)**: Low quality-weighted similarity can reject an overlap even when some sequence similarity exists.
* **Hard gates (`-o`, `-p`)**: Minimum overlap length and minimum percent identity still must be met.

### Why QUAL propagation matters

QUAL is required for correct CAP3 scoring.
Without `.qual` files, CAP3 treats every base as Q=10, which distorts clipping,
mismatch penalties, and similarity scoring. For Sanger data, that can change
whether two reads are judged mergeable at all.

### Audit-driven relaxation strategy

Use the audit to target the gate that failed instead of relaxing everything at once:

* If overlap length is just below threshold, adjust overlap length first.
* If the overlap looks real but clipping removed too much end sequence, revisit clipping before lowering identity.
* If overlap length and identity look acceptable but CAP3 still rejects the merge, inspect quality-weighted gates before broad relaxation.

### Validation criteria for rescued contigs

If you rescue a borderline case by relaxing a gate, the strongest validation pattern is:

* forward and reverse singlet top hits still agree taxonomically,
* the rescued contig BLAST result agrees with both singlets,
* and the overlap region does not show implausible clusters of disagreement.

That keeps the rescue biologically interpretable instead of merely computationally permissive.

## Short FAQ on pairing rules

### How can I tell whether files were auto-paired correctly?

Use `qc/pairing_report.tsv`.
MicroSeq pairs files by resolved sample ID after stripping forward/reverse tokens/primer labels,
with wells participating only when well enforcement is enabled.

### What does `pair_missing` mean in practice?

It means one orientation was missing after pairing logic.
It does **not** mean BLAST failed. It means the sample never reached true paired assembly.

### What is a singlet?

A singlet is a standalone sequence record that survived assembly processing when
a defended merged contig was not produced. In practice, it means the sample still
has usable sequence evidence, but the forward/reverse pair did not support one
clean consensus contig.

### When should I care about well enforcement?

Use it when filenames come from plate-based exports and cross-well mismatches are a real risk.
Leave it off when sample IDs are already unique and wells are not biologically meaningful.

## Related docs

* [CLI Workflows](cli_workflow_tutorial.md) for wrapper commands and headless runs
* [GUI Walkthrough](gui_usage.md) for button-by-button GUI usage
* [Output Artifacts Reference](output_artifacts.md) for the canonical file table
* [Workflow Resolution Funnel](workflow_resolution.md) for ambiguity routing and sample-level resolution states
* [MicroSeq design: algorithm and architecture choices for Sanger assembly](design_algorithm_choices.md) for algorithm tradeoffs

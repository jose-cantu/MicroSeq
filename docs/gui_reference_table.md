---
layout: page
title: GUI Table Reference
permalink: /gui-table-reference/
---

## What this page is for

This page is the detailed column-by-column reference for the main GUI tables.

Use this page when you want the exact meaning of the values shown in:

* **Assembly Summary**
* **BLAST Inputs**
* **Diagnostics**
* **Compare Assemblers**

This page is a table reference, not the canonical source for output-file meanings or the canonical source for semantic routing rules.

Use these pages alongside it:

* [Gui Walkthrough](gui_usage.md) for operating the interface
* [Output Artifacts Reference](output_artifacts.md) for the meaning of files on disk
* [Paired CAP3 Assembly](paired_assembly.md) for paired-mode semantics and BLAST handoff meaning
* [Workflow Resolution Funnel](workflow_resolution.md) for canonical status semantics such as `ambiguous_overlap`
* [MicroSeq design: algorithm and architecture choices for Sanger assembly](design_algorithm_choices.md) for algorithm tradeoffs

## Assembly Summary

The Assembly Summary table is the fastest sample-level overview of what happened in the assembly stage and what kind of sequence output is being handed downstream.

### Key columns

| Column | Meaning |
| --- | --- |
| `sample_id` | The sample key used throughout paired assembly outputs. |
| `status` | Final sample-level assembly status for that row. In legacy CAP3 mode this may be CAP3-derived (`assembled`, `singlets_only`, `cap3_no_output`) and may also reflect overlap-driven outcomes in workflows that synchronize merge decisions. |
| `assembler` | Which assembler path or backend produced the summarized row, such as `cap3:strict` or a `merge_two_reads:*` backend. |
| `contig_len` | Length of the selected or produced contig sequence output. When multiple records exist, this is usually the maximum length for that row. |
| `blast_payload` | What kind of sequence output is being passed to BLAST: `contig`, `singlet`, `no_payload`, or `pair_missing`. |

### Additional columns

| Column | Meaning |
| --- | --- |
| `selected_engine` | The engine that produced the chosen result, such as a specific overlap backend or `cap3`. |
| `configured_engine` | The configured overlap engine for the run. Some summaries populate this directly; some compare-driven summaries may leave it blank. |
| `merge_status` | Merge prepass outcome, such as `merged`, `identity_low`, `overlap_too_short`, `quality_low`, `ambiguous_overlap`, `not_end_anchored`, or `high_conflict`. |
| `merge_overlap_len` | Overlap length reported by the merge prepass. |
| `merge_identity` | Overlap identity reported by the merge prepass. |
| `overlaps_saved` | Number of overlaps CAP3 kept according to parsed `.cap.info` counters. |
| `overlaps_removed` | Number of overlaps CAP3 removed according to parsed `.cap.info` counters. |
| `primer_mode` | Primer trimming policy used for the run, such as `off`, `detect`, or `clip`. |
| `primer_stage` | Whether primer trimming or detection was applied `pre_quality` or `post_quality`. |

### Reading `status` correctly

The `status` column is useful in the GUI, but its canonical meaning lives elsewhere:

* use [Workflow Resolution Funnel](workflow_resolution.md) for exact routing semantics such as `ambiguous_overlap`
* use [Paired CAP3 Assembly](paired_assembly.md) for paired-mode interpretation such as `contig`, `singlet`, `pair_missing`, and CAP3 failure interpretation

## BLAST Inputs

The BLAST Inputs table explains what sequence was actually sent to BLAST and why that output was chosen for each sample.

### Main displayed columns

| Column | Meaning |
| --- | --- |
| `sample_id` | Sample key for that BLAST-input row. |
| `blast_payload` | Type of sequence sent to BLAST: `contig`, `singlet`, `no_payload`, or `pair_missing`. |
| `reason` | Why that payload type was selected. Examples include `contigs_present`, `singlets_only`, `cap3_no_output`, `pair_missing`, `selected_payload`, or `winner_no_payload`. |
| `payload_ids` | Mapping of rewritten BLAST FASTA IDs back to the original CAP3 or read IDs. |

### Common `blast_payload` interpretations

| `blast_payload` | Practical meaning |
| --- | --- |
| `contig` | A merged consensus sequence was available and used. |
| `singlet` | No merged contig survived, but at least one standalone sequence remained usable. |
| `no_payload` | Assembly reached a point where no usable sequence output remained for BLAST. |
| `pair_missing` | No valid forward/reverse pair existed after pairing logic. |

### Handoff logic to remember

In paired mode, MicroSeq chooses the first available output in this order:

1. CAP3 contigs
2. CAP3 singlets
3. no usable sequence output

### Additional underlying fields you may see in the TSV

Depending on the run path and synchronization contract, the underlying `asm/blast_inputs.tsv` can also include fields such as:

| Column | Meaning |
| --- | --- |
| `payload_entity_n` | Number of concrete sequence entities emitted for that row. |
| `hypothesis_map` | Mapping from `qseqid` to structural hypothesis ID. |
| `source_id_map` | Mapping from `qseqid` back to the original source sequence ID. |

For the canonical schema contract, use [Workflow Resolution Funnel](workflow_resolution.md).

## Diagnostics

The Diagnostics table is where you look when the question is not just “what happened,” but “why did the overlap decision go this way?”

### Core diagnostic columns

| Column | Meaning |
| --- | --- |
| `sample_id` | Sample key being audited. |
| `status` | Diagnostic overlap outcome, such as `ok`, `overlap_too_short`, `overlap_identity_low`, `overlap_quality_low`, `ambiguous_overlap`, or `not_end_anchored`. |
| `overlap_len` | Selected overlap length. |
| `overlap_identity` | Identity of the selected overlap candidate. |
| `overlap_quality` | Mean overlap quality when quality information is available. |
| `orientation` | Selected overlap orientation, usually `forward` or `revcomp`. |

### Feasibility and orientation columns

| Column | Meaning |
| --- | --- |
| `best_identity` | Best identity value seen among relevant candidates in the diagnostic context. |
| `best_identity_orientation` | Orientation associated with that best-identity candidate. |
| `anchoring_feasible` | Whether at least one candidate was end-anchored **and** passed the active feasibility gates. |
| `end_anchored_possible` | Whether any candidate was end-anchored at all, even before full feasibility filtering. |
| `fwd_best_identity` | Best identity among feasible forward-orientation candidates. |
| `revcomp_best_identity` | Best identity among feasible reverse-complement candidates. |
| `fwd_best_overlap_len` | Overlap length paired to the best feasible forward identity. |
| `revcomp_best_overlap_len` | Overlap length paired to the best feasible reverse-complement identity. |
| `fwd_anchor_feasible` | Whether the forward orientation had any end-anchored feasible candidate. |
| `revcomp_anchor_feasible` | Whether the reverse-complement orientation had any end-anchored feasible candidate. |
| `identity_delta_revcomp_minus_fwd` | Difference between reverse-complement and forward best-identity values in feasible space. |
| `selected_vs_best_identity_delta` | Difference between the best-identity candidate and the candidate actually selected. |

### Tie and ambiguity columns

| Column | Meaning |
| --- | --- |
| `top_candidate_count` | Number of near-tied top feasible candidates under the ambiguity thresholds. |
| `top2_identity_delta` | Absolute identity difference between top-1 and top-2 feasible candidates. |
| `top2_overlap_len_delta` | Absolute overlap-length difference between top-1 and top-2 feasible candidates. |
| `top2_quality_delta` | Absolute overlap-quality difference between top-1 and top-2 feasible candidates. |
| `tie_reason_code` | Why a near tie was flagged, such as `len_equal`, `mismatch_equal`, `identity_eps`, or `quality_eps`. |
| `ambiguity_identity_delta_used` | Effective identity threshold used to decide ambiguity under missing-quality conditions. |
| `ambiguity_quality_epsilon_used` | Effective quality epsilon used to decide near-tie ambiguity when quality exists. |

### Any-space and pre/post-trim comparison columns

| Column | Meaning |
| --- | --- |
| `fwd_best_identity_any` | Best forward identity without the `min_overlap` feasibility filter. |
| `revcomp_best_identity_any` | Best reverse-complement identity without the `min_overlap` feasibility filter. |
| `fwd_best_overlap_len_any` | Overlap length associated with the forward best identity in raw evidence space. |
| `revcomp_best_overlap_len_any` | Overlap length associated with the reverse-complement best identity in raw evidence space. |
| `pretrim_best_identity` | Best overlap identity recomputed from pre-primer-trim paired inputs. |
| `pretrim_best_overlap_len` | Best overlap length recomputed from pre-primer-trim paired inputs. |
| `pretrim_status` | Diagnostic status before primer-trim effects are applied. |
| `posttrim_best_identity` | Best overlap identity after trimming. |
| `posttrim_selected_overlap_len` | Selected overlap length after trimming. |
| `posttrim_status` | Final diagnostic status after trimming. |

### Engine and fallback columns

| Column | Meaning |
| --- | --- |
| `selected_engine` | Overlap engine used for the diagnostic row. |
| `fallback_used` | Whether a cascade or all-engine strategy had to move away from the first engine choice. |
| `overlap_engine` | UI label for the configured default overlap engine. In file output this may appear as `configured_engine`. |

### Primer-trim-linked columns

| Column | Meaning |
| --- | --- |
| `primer_trim_bases_fwd` | Total primer-trimmed bases for the forward read, derived from primer-trim telemetry. |
| `primer_trim_bases_rev` | Total primer-trimmed bases for the reverse read, derived from primer-trim telemetry. |

### What `anchoring_feasible` means

This is one of the most important diagnostics columns to understand.

A candidate overlap can look sequence-similar but still be a poor merge candidate if it does not behave like a biologically plausible end-to-end stitching overlap.

Use these distinctions:

| Field | Meaning |
| --- | --- |
| `end_anchored_possible` | At least one candidate could be placed as an end-anchored overlap. |
| `anchoring_feasible` | At least one end-anchored candidate also passed the active length, identity, and quality gates. |

So:

* `end_anchored_possible = yes` means the geometry exists
* `anchoring_feasible = yes` means the geometry exists **and** the candidate passed the feasibility thresholds

For the canonical ambiguity and overlap-routing rules, use [Workflow Resolution Funnel](workflow_resolution.md).

## Compare Assemblers

The Compare Assemblers table is the side-by-side backend view used when MicroSeq runs more than one paired backend or profile and then chooses a winner deterministically.

### Main columns

| Column | Meaning |
| --- | --- |
| `sample_id` | Sample being evaluated by each backend row. |
| `assembler_id` | Internal backend identifier, such as `merge_two_reads:biopython` or `cap3:strict`. |
| `assembler_name` | Human-readable label shown in the GUI. |
| `status` | Per-backend outcome for that row. |
| `selected_engine` | Engine actually used for that backend row. |
| `contig_len` | Maximum usable contig length produced by that backend row, if any. |
| `warnings` | Warning text, merge warning text, or exception detail associated with that row. |

### Additional compare columns

| Column | Meaning |
| --- | --- |
| `dup_policy` | Duplicate-handling policy used for that row. |
| `diag_code_for_machine` | Compact machine-readable diagnostic code for the row result. |
| `diag_detail_for_human` | Human-readable explanation of the row result. |
| `cap3_contigs_n` | Number of CAP3 contig records detected for that row. |
| `cap3_singlets_n` | Number of CAP3 singlet records detected for that row. |
| `cap3_info_path` | Path to the CAP3 `.cap.info` artifact when present. |
| `cap3_stdout_path` | Legacy CAP3-specific stdout path retained for compatibility. |
| `cap3_stderr_path` | Legacy CAP3-specific stderr path retained for compatibility. |
| `tool_name` | Backend or tool label used for process diagnostics. |
| `tool_stdout_path` | Generic stdout or diagnostic artifact path for that row. |
| `tool_stderr_path` | Generic stderr artifact path for that row. |
| `selection_trace_path` | Path to the per-sample winner-selection trace written under `asm/selection_trace/`. |
| `winner_reason` | Short explanation of why that row was selected as the winner. |
| `payload_fasta` | Path to the backend sequence-output FASTA produced by that row. |

### How to read `status` in Compare Assemblers

The `status` column is backend-specific.

Examples:

* merge backend rows may show `merged`, `identity_low`, `overlap_too_short`, `quality_low`, `ambiguous_overlap`, `not_end_anchored`, or `high_conflict`
* CAP3 backend rows may show `assembled` or `cap3_no_output`
* exception paths may show error or warning text in the diagnostics fields

### How winners are chosen

In all-assembler or compare-driven modes, winners are chosen per sample with deterministic ranking keys. The current logic favors:

1. successful outcomes
2. stronger payload kinds such as `contig` over `singlet`
3. lower ambiguity burden
4. longer normalized length
5. deterministic backend tie-break rules

That winner then affects what appears downstream in Assembly Summary and BLAST Inputs.

For the canonical structural and resolution logic, use:

* [Paired CAP3 Assembly](paired_assembly.md)
* [Workflow Resolution Funnel](workflow_resolution.md)

### Backend glossary

Use this quick glossary when reading backend rows:

| Backend | Meaning |
| --- | --- |
| `cap3:*` | CAP3-based paired assembly profile |
| `merge_two_reads:ungapped` | End-anchored ungapped overlap merge |
| `merge_two_reads:biopython` | Gapped pairwise-alignment backend using Biopython machinery |
| `merge_two_reads:edlib` | Gapped edit-distance backend using edlib |

Use [MicroSeq design: algorithm and architecture choices for Sanger assembly](design_algorithm_choices.md) for the design rationale behind staged assembly and CAP3 fallback.

---
layout: page
title: Workflow Resolution Funnel
---

# Updated MicroSeq workflow mental model

MicroSeq now follows a **4-stage funnel** so the default user view stays simple while preserving audit-grade details.

## 1) Ingest (Trace QC)

MicroSeq computes AB1 trace-QC metrics (signal/noise, mixed-signal proxies, vendor tags when present).
Threshold flags are optional (`--trace-qc-flags` or `trace_qc.enable_flags`).

Mixture escalation is config-driven:

- `trace_qc.enable_mixture_inference`
- `trace_qc.mixture_suspect_threshold`

## 2) Assemble (Structural hypotheses)

MicroSeq generates multiple hypotheses **only when there is a structural tie** (for example, ambiguous overlap outcomes).
If there is only one structural path, the sample is structurally unambiguous.


## Assemble → Validate handoff: how `ambiguous_overlap` is decided

### Candidate generation in 3 stages

#### Stage 1  Generate candidate overlaps (the “sliding placements” step)

Stage 1  Generate candidates: MicroSeq tries multiple orientations (`R` as-is vs `revcomp(R)`) and multiple relative placements (sliding the reverse read along the forward read). Each placement defines an overlap region. For that overlap, MicroSeq computes overlap length, mismatches, identity, and (when QUAL is available) an overlap-quality score. Candidates that do not connect the reads at their ends are rejected by the end-anchoring rule (`anchor_tolerance_bases` allows small end drift due to trimming).

Definition: a **candidate** is one proposed end-anchored overlap between `F` and (`R` or `revcomp(R)`), characterized by relative placement (and, for gapped backends, an alignment path), from which `overlap_len`, `mismatches`, `identity`, and `overlap_quality` are computed.

- It evaluates both orientations for the reverse read (`forward` and `revcomp`).
- Ungapped intuition: slide one read across the other (`... left, center, right ...`), where each
  placement implies a different overlap slice.
- Gapped backends may represent differences with an indel path, then anchoring is checked using
  end-anchored placement rules.

> Important: this is not “use only a 30 bp overlap.”
> `anchor_tolerance_bases` is a forgiveness margin for **end placement** (how close overlap must be
> to read ends after trimming), not the overlap length itself. Overlaps can still be much longer
> (for example 100+ bp).

#### Stage 2  Feasibility filtering (hard gates)

From all generated candidates, MicroSeq keeps only feasible candidates:

- `overlap_len >= min_overlap`
- `identity >= min_identity`
- quality gate only when `quality_mode=blocking`:
  - `overlap_quality` must exist and be `>= min_quality`
- end-anchored requirement (`end_anchored=True`)

Under current repo defaults (`min_overlap=100`, `min_identity=0.8`, `min_quality=20.0`,
`quality_mode=warning`), quality is advisory by default, while overlap length/identity and
end-anchoring remain hard feasibility checks.

#### Stage 3  Final decision (unique best vs `ambiguous_overlap`)

Feasible candidates are ranked deterministically (longer overlap first, then fewer mismatches, then
overlap quality, then identity), and ambiguity is evaluated only on the top-2 candidates.

Tie rule for `ambiguous_overlap`:

1. top-1 and top-2 have the same `overlap_len`
2. top-1 and top-2 have the same `mismatches`
3. and either:
   - both have quality and `|q1 - q2| <= ambiguity_quality_epsilon`, or
   - quality is unavailable and `|id1 - id2| <= ambiguity_identity_delta`

If all hold, selector returns `status=ambiguous_overlap`; otherwise top-1 is the unique winner.

### Working examples: ungapped sliding vs gapped alignment backends

Legend: `|` = match, `x` = mismatch, `-` = gap (indel introduced by aligner).
These diagrams illustrate candidate *types*; exact mismatch/score values can vary by backend and scoring settings.

#### A) Ungapped, end-anchored sliding candidates (placement-driven)

For ungapped candidates, the generator effectively slides one sequence relative to the other and
scores each valid end-anchored placement. For a given orientation and placement, you get one ungapped
candidate.

Toy example (repetitive tail so nearby placements can look similarly good):

- `F   = ATATATATATGGGG`
- `Rrc = ATATATATCCCCCC`

Candidate A (placement 1: `Rrc` starts aligned at the beginning of `F`):

- `F_overlap = ATATATATAT`
- `R_overlap = ATATATATCC`

```text
F:    ATATATATATGGGG
Rrc:  ATATATATCCCCCC
      ||||||||xx
```

Candidate B (placement 2: `Rrc` starts one base later in `F`):

- `F_overlap = TATATATATG`
- `R_overlap = ATATATATCC`

```text
F_overlap:  TATATATATG
R_overlap:  ATATATATCC
            x|||||||xx
```

In repetitive sequence, A and B can be near-tied in identity/mismatches/quality, which is exactly
how ambiguous overlap candidates arise in Stage 1 before tie-resolution in Stage 3.

#### B) Gapped backend candidates (Biopython/edlib): indels are allowed

Gapped backends can explain insertion/deletion differences with gaps, instead of forcing a mismatch
cascade.

Example:

- `F   = ATGCCCTTAG`
- `Rrc = ATGCCTTAG` (one base shorter around `C` run)

Ungapped intuition (no indels allowed):

In ungapped mode, length differences manifest as forced mismatches/edge penalties (the trailing `-`
below is only to keep the visualization aligned):

```text
F:    ATGCCCTTAG
Rrc:  ATGCCTTAG-
      ||||x||||-
```

Gapped alignment can model this as one indel event:

```text
F:    ATGCCCTTAG
Rrc:  ATG-CCTTAG
      ||| ||||||
```

This can produce a stronger candidate than any ungapped placement for the same read pair.

#### C) Repeats/homopolymers can create near-equivalent gapped interpretations

In low-complexity sequence, an indel can often be represented in nearby positions with similar
scores. Conceptually:

```text
Path 1: A AAAAGGGTT
        - AAAAGGGTT

Path 2: AAAA AGGGTT
        AAA- AGGGTT
```

Both describe the same underlying deletion in a homopolymer context.

Implementation note: in practice MicroSeq typically emits one best alignment per backend per orientation
(and sometimes a small bounded set), rather than enumerating a combinatorial set of all optimal
paths. So many practical near-ties come from relative-placement/orientation near-ties, while gapped
backends still explain why mismatch/indel trade-offs can differ from ungapped sliding.

#### D) Why end-anchoring still matters

End-anchoring constrains relative placement (overlap must reach the expected stitching ends within tolerance so reads can connect into one contiguous sequence). Within that
constraint:

- ungapped mode explores placement shifts,
- gapped mode can change mismatch/indel trade-offs and effective overlap metrics.

So candidate sets can differ by backend even for the same sample and orientation.


Plain-language end-anchoring rule:

- **Accept** when overlap reaches the expected stitching ends (within `anchor_tolerance_bases`) so the reads connect into one contiguous sequence.
- **Reject** when a high-identity block is internal-only and does not connect those stitching ends.

```text
ACCEPT (end-anchored):
F:      AAAAACCCCCGGGGG
                 |||||||
Rrc:          CCCCCGGTTTTT
             overlap reaches the expected stitching end (or is within tolerance)

REJECT (internal-only match):
F:      AAAAACCCCCGGGGGTTTTT
                 |||||||
Rrc:      XXYYCCCCCGGZZWW
         good local match, but does not connect stitching ends
```

### Worked example (ambiguous)

Use the same feasible-candidate context as above (passes Stage 2 gates):

- Candidate A: `overlap_len=120`, `mismatches=2`, `overlap_quality=34.20`, `identity=0.9833`
- Candidate B: `overlap_len=120`, `mismatches=2`, `overlap_quality=34.16`, `identity=0.9833`

Configured ambiguity thresholds:

- `ambiguity_quality_epsilon = 0.10`
- `ambiguity_identity_delta = 0.0025`

Checks:

- same length? yes (`120 == 120`)
- same mismatches? yes (`2 == 2`)
- quality near-tie? yes (`|34.20 - 34.16| = 0.04 <= 0.10`)

Result: `ambiguous_overlap`.

### Counterexample (not ambiguous)

- A: `overlap_len=120`, `mismatches=2`
- B: `overlap_len=119`, `mismatches=2`

Different overlap length breaks the ambiguity rule immediately, so A is selected as the unique best
candidate.

### What ambiguous policies do (`merge_two_reads`)

After `ambiguous_overlap`, policy controls output behavior:

- `strict`: keep ambiguous outcome; no forced merge payload.
- `singlets`: emit singlets (`ambiguous_overlap_singlets`).
- `best_guess`: force top-1 candidate and emit one consensus (`merged_best_guess`).
- `topk`: emit `alt1..altK` consensus sequences (`ambiguous_topk`), where `K=ambiguous_top_k`
  (bounded by available feasible candidates).

When `topk` is used, each emitted alternative is tracked as a structural branch through
`hypothesis_map` (`qseqid -> structural_hypothesis_id`). Those branches are then independently
validated by BLAST/taxonomy and collapsed into sample-level resolution state.

## 3) Validate (BLAST + taxonomy)

MicroSeq runs BLAST/taxonomy against payloads and ranks best hit per hypothesis deterministically:

1. `bitscore` (desc)
2. `pident` (desc)
3. `qcovhsp` (desc)
4. `evalue` (asc)
5. `qseqid` (asc tie-break)

Taxonomy agreement can be evaluated at a configured rank (`species` by default) using parsed lineage tokens (`k__`, `p__`, `c__`, `o__`, `f__`, `g__`, `s__`).

## 4) Resolve (sample-level state)

At sample level, MicroSeq assigns:

- `unambiguous`
- `resolved_by_evidence`
- `needs_review`

And stores:

- `resolution_state`
- `resolved_hypothesis`
- `resolution_reason`

These fields are emitted to review/summary TSV outputs so the UI can stay quiet by default.

## What changed for users (practical impact)

The latest contract updates are mainly about **predictability** and **traceability**:

1. **`hypothesis_map` now has one meaning everywhere**
   - It maps `qseqid -> structural_hypothesis_id` only.
   - This prevents provenance IDs from being mixed into structural decision logic.

2. **`source_id_map` was added for provenance**
   - It maps `qseqid -> original source sequence id`.
   - You can now audit where each BLAST input came from without overloading hypothesis logic.

3. **Payload size and structural ambiguity are now separated**
   - `payload_entity_n` tracks how many concrete sequence entities were emitted.
   - `structural_hypothesis_n` tracks decision branches.
   - Result: a sample can have multiple payload entities and still be structurally unambiguous.

4. **Multi-entity payloads are advisory by default**
   - Non-`contig_alt` multi-entity payloads add `multi_payload` to `warning_flags`.
   - This is non-blocking unless other safety/review conditions escalate.

5. **Rank-aware missing taxonomy is explicit**
   - If hits exist but no usable label can be extracted at the configured rank, resolution is `rank_missing`.
   - This is more informative than folding these cases into generic ambiguity.

6. **Review vs advisory is explicit**
   - `review_action` + `review_reason` govern queue inclusion.
   - `advisory_reason` summarizes non-blocking signals for UI/triage.

## Decision table

| Condition | resolution_state | resolution_reason |
| --- | --- | --- |
| Single structural hypothesis and all required evidence present | `unambiguous` | `single_hypothesis` |
| Multiple structural hypotheses, all validated, same taxonomy label at configured rank, no blocking safety flags | `resolved_by_evidence` | `hypotheses_agree_<rank>` |
| Taxonomy disagreement across validated hypotheses | `needs_review` | `ambiguous_taxonomy` |
| Structural hypotheses > hypotheses with hits | `needs_review` | `partial_hits` |
| Trace QC `FAIL` (sticky) | `needs_review` | `trace_fail` |
| Safety escalation (e.g. high conflict) | `needs_review` | safety flag value |

## No-hit and missing-payload policy

MicroSeq drives review queue population from blast-input contract rows, then joins taxonomy evidence when present.

| Contract/evidence condition | review behavior |
| --- | --- |
| `blast_payload=pair_missing` | `needs_review`, reason `pair_missing` |
| `blast_payload=no_payload` | `needs_review`, reason `no_payload` |
| structural hypotheses exist but zero taxonomy hits | `needs_review`, reason `no_hits` |
| taxonomy file unavailable / not parseable | `needs_review`, reason `taxonomy_missing` |

## Trace escalation rules (paired samples)

Sample-level trace status is computed from F/R statuses using:

`FAIL > WARN > PASS > NA`

- `trace_fail` => sticky safety escalation + `needs_review`
- `trace_warn` => non-blocking advisory in `warning_flags`/`advisory_reason`
- mixture inference (when enabled) can set `review_reason=mixture_suspected`

## UI behavior

- **Default view:** show `unambiguous` + `resolved_by_evidence`
- **Expert view:** include `needs_review` and expandable hypotheses

This keeps output Geneious-like for routine samples while preserving full provenance.

## Assemble/validate status routing matrix

| Trigger condition | Status label emitted (`merged`, `ambiguous_overlap`, `high_conflict`, `quality_low`, `cap3_unverified`) | IUPAC involvement (`yes/no`, only after a single structural path is selected) | Routing/next action | Primary artifact field(s) to inspect (`merge_status`, `status`, `hypothesis_map`, `review_reason`, `warning_flags`) |
| --- | --- | --- | --- | --- |
| One unique top feasible overlap candidate (passes overlap + anchoring gates) | `merged` | no (only after a single structural path is selected) | Emit merged sequence payload; continue to validate/rank hits on that single path | `merge_status`; sample `status` |
| Top-1 and top-2 feasible candidates are tie-equivalent under ambiguity thresholds | `ambiguous_overlap` | no (only after a single structural path is selected) | Branch hypotheses (`best_guess`/`topk`) or hold as strict ambiguity for review routing | `merge_status`; `hypothesis_map`; `review_reason` |
| Overlap exists but disagreement burden exceeds configured high-confidence conflict guardrail | `high_conflict` | no (only after a single structural path is selected) | Route by configured conflict action (typically CAP3 fallback or explicit review escalation) | `merge_status`; sample `status`; `warning_flags`; `review_reason` |
| Candidate is structurally feasible but blocked by quality policy (`quality_mode=blocking`) | `quality_low` | no (only after a single structural path is selected) | Do not accept fast merge; route to CAP3 fallback or review path per policy | `merge_status`; sample `status`; `warning_flags` |
| CAP3 fallback output fails verification contract (for example, missing source-read representation) | `cap3_unverified` | no (only after a single structural path is selected) | Keep singlet/no-payload fallback and escalate to review if required | sample `status`; `review_reason`; `warning_flags` |

## Why this method instead of manual IUPAC consensus curation

Older Sanger workflows often resolve ambiguous overlap or mixed-base positions by creating a single IUPAC consensus, then manually curating traces in a GUI tool.
That approach is useful for expert review, but it has trade-offs for reproducible, batch-scale pipelines.

MicroSeq intentionally uses a **hypothesis + evidence resolution** method first, with manual review as a targeted fallback:

1. **Preserves uncertainty explicitly**
   - Instead of collapsing early into one IUPAC sequence, MicroSeq keeps structural alternatives as hypotheses.
   - This avoids hiding meaningful disagreement before validation.

2. **Deterministic at scale**
   - Resolution uses a fixed ranking tuple and contract fields, so repeated runs are explainable and reproducible.
   - Manual consensus editing can vary by operator and session.

3. **Separates routine from exceptional cases**
   - When hypotheses agree taxonomically and safety checks pass, MicroSeq resolves automatically (`resolved_by_evidence`).
   - Only the small subset with conflicts/no-hit/trace-fail is escalated to `needs_review`.

4. **Better auditability for regulated/shared labs**
   - Review queue, reasons, warnings, and structural-vs-hit counts are emitted as machine-readable artifacts.
   - Manual curation can still be done, but now it is focused and documented per sample.

5. **Still compatible with manual trace review**
   - This design does **not** prohibit classic chromatogram curation.
   - It reduces manual burden by routing analysts only to samples where evidence actually disagrees or quality/safety thresholds fail.

In short: IUPAC/manual consensus remains a valuable expert tool, but MicroSeq’s default model is optimized for reproducibility, throughput, and explicit decision provenance.

## Schema contracts

### `asm/blast_inputs.tsv`

Required columns:

- `sample_id`, `blast_payload`, `payload_ids`, `reason`
- `payload_kind`, `payload_n`, `payload_entity_n`, `payload_max_len`
- `ambiguity_flag`, `safety_flag`, `decision_source`
- `review_reason`, `warning_flags`
- `structural_hypothesis_n`, `hypotheses_with_hits_n`, `missing_hits_n`
- `hypothesis_map` (`qseqid=structural_hypothesis` pairs)
- `source_id_map` (`qseqid=original_source_id` pairs)

Interpretation notes:

- `payload_n` / `payload_entity_n` are payload counts, **not** structural branch counts.
- `structural_hypothesis_n` is the number of structural alternatives considered by resolution logic.
- `warning_flags` contains all non-blocking warnings (semicolon-separated, deduplicated, sorted).

### `asm/assembly_summary.tsv`

In addition to assembly/reporting columns, resolution fields are synchronized post-taxonomy:

- `resolution_state`, `resolved_hypothesis`, `resolution_reason`
- `review_action`, `review_reason`, `advisory_reason`, `warning_flags`
- `structural_hypothesis_n`, `hypotheses_with_hits_n`, `missing_hits_n`
- `trace_status`, `trace_status_f`, `trace_status_r`, `trace_flags`

### `qc/review_queue.tsv`

Primary output contract:

- `sample_id`
- `review_action` *(renamed from status to avoid cross-table ambiguity)*
- `review_reason`, `advisory_reason`
- `warning_flags`
- `structural_hypothesis_n`, `hypotheses_with_hits_n`, `missing_hits_n`
- `top_labels`
- `resolution_state`, `resolved_hypothesis`, `resolution_reason`
- `trace_status`, `trace_flags`

Behavior contract:

- `review_action=queue` only for actionable review cases.
- `review_reason` is populated only when `review_action=queue`.
- `advisory_reason` summarizes one prioritized non-blocking signal.

## Worked example (happy path)

If sample hypotheses (`hyp1`, `hyp2`, `hyp3`) all hit the same species with near-identical stats:

- Structural ambiguity is retained internally,
- taxonomy agreement collapses the sample to one resolved result,
- `resolution_state=resolved_by_evidence`,
- no manual review unless safety/trace flags escalate.

This is the intended “hidden complexity, simple default output” behavior.

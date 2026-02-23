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

# MicroSeq design: algorithm and architecture choices for Sanger assembly

This document explains why MicroSeq uses a staged assembly strategy for paired Sanger reads and how that compares with common alternatives.

## Current strategy (implemented)

1. **Fast overlap merge first** (`merge_two_reads`):
   - End-anchored, ungapped overlap search across forward/revcomp orientations.
   - Candidate ranking by overlap length, mismatch count, overlap-quality tie-breakers.
   - Explicit **ambiguous overlap** detection when top candidates are not uniquely best.
2. **Policy-aware merge decision**:
   - Terminal hard policy in strict mode (`quality_low` in `blocking`).
   - **High-confidence conflict guardrail** (`high_conflict_q_threshold`, `high_conflict_action`).
   - IUPAC ambiguity for low-confidence quality ties.
3. **CAP3 fallback** for non-merged fast-path outcomes (including ambiguous/high_conflict when routed).
4. **Post-CAP3 validation**:
   - Accept CAP3 contig only when both source reads are represented.
   - Mark failures as `cap3_unverified` and keep singlet fallback.

## Why this architecture

### Strengths
- **Speed on clean data**: most pairs resolve in the fast path.
- **Robust rescue on messy tails/indels**: CAP3 handles difficult Sanger edge cases.
- **Auditability**: explicit statuses (`merged`, `ambiguous_overlap`, `high_conflict`, `quality_low`, `cap3_unverified`) are easier to reason about than opaque sensitivity knobs.
- **Information retention**: IUPAC retains more evidence than blanket `N` at low-confidence ties.

### Tradeoffs
- More moving parts than a single monolithic assembler.
- Requires clear telemetry and tests to avoid state drift between fast and fallback paths.

## Comparison with common approaches

### 1) Single-pass greedy overlap assembler
**Pros**
- Simple operational model (one engine).
- Mature tools exist.

**Cons**
- Pays full computational cost even when easy fast merges dominate.
- Harder to expose explicit failure classes for automated QC routing.

### 2) Strict expected-errors gating before assembly
**Pros**
- Easy global quality filter.

**Cons**
- Length-coupled rejection can be too punitive for long Sanger reads.
- May discard salvageable read pairs where overlap-local evidence is strong.

### 3) Reference-guided assembly
**Pros**
- Can resolve ambiguous de novo overlaps.
- Useful for targeted loci with trusted references.

**Cons**
- Introduces reference bias.
- Requires curation/versioning of reference sets and may hide novel variation.

## When to use reference-guided mode later
Add as an optional branch when:
- overlap ambiguity is frequent,
- a validated locus reference panel exists,
- or downstream analysis tolerates reference bias.

Keep de novo staged assembly as the default for broad, transparent Sanger workflows.

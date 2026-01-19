---
layout: page
title: MicroSeq's Design: Questions & Answers 
---

# MicroSeq's Design: Questions and Answers 

## What are the advantages of using MicroSeq for blasting vs just grabbing my Sanger sequence and blasting it on the NCBI website?

**Short answer:** MicroSeq is a reproducible end‑to‑end pipeline, while a website BLAST is a single one‑off query.

**Why MicroSeq can be better in practice:**

* **Starts from raw AB1 or FASTQ and standardizes QC/trim.** MicroSeq turns raw traces into QC‑passed FASTA, making the BLAST step consistent across samples and runs.
* **Batch‑friendly and reproducible.** Same config, same settings, same outputs across a full plate or multiple runs.
* **Multiple curated databases with one switch.** Choose Greengenes2, SILVA, or NCBI 16S.
* **Paired assembly before BLAST.** If you have forward/reverse Sanger reads, MicroSeq can assemble contigs before searching.
* **Downstream‑ready outputs.** BLAST hits can be immediately joined to taxonomy and exported to BIOM/CSV for QIIME‑2/phyloseq workflows.

If you want a one‑off quick check, NCBI is fine. If you want a repeatable, batchable workflow with QC, assembly, taxonomy joins, and BIOM output, MicroSeq saves a lot of manual steps.

---

## Why are there dials (identity, coverage, mode) for BLAST? Doesn’t that make it more confusing?

**Short answer:** the dials exist because your data and goals vary, and identity/coverage thresholds are not one‑size‑fits‑all.

MicroSeq exposes those parameters explicitly so you can:

* Match the stringency to your biology (e.g., species‑level vs genus‑level resolution).
* Handle partial reads or low‑quality Sanger traces without silently losing all hits.
* Switch between **fast** (megablast) and **comprehensive** (blastn) search strategies.

Defaults (97% identity / 80% coverage / megablast) are meant to be a sane starting point, not a rigid rule. The dials make the logic transparent instead of hiding it.

---

## What if you get no hits at 97% identity and 80% coverage in MicroSeq?

MicroSeq’s design anticipates this:

* **It automatically falls back to a more sensitive search.** If megablast finds no close hit or the top hit is below 90% identity/coverage, MicroSeq reruns with blastn.
* **You can run a relaxed search and filter later.** This gives you a broader hit list first, then lets you tighten cutoffs based on the results.
* **QC outputs are preserved.** You can inspect the trimmed FASTA and compare to NCBI BLAST manually if needed.

In short: “no hits” is a signal to either (1) loosen thresholds, (2) use the more sensitive mode, or (3) inspect QC and trimming outcomes.

---

## Why is pairing automatic and tied to naming? Shouldn’t there be a simpler way?

**Short answer:** pairing must be deterministic, and filename tokens are the only reliable metadata available in typical Sanger workflows.

AB1/FASTQ files do **not** contain a universal field that says “this is the reverse mate for sample X.” So MicroSeq pairs reads using the only stable identifier it can trust at scale: filenames (sample IDs + primer tokens, optionally wells).

**What MicroSeq does to make this easier:**

* Auto‑detects primer tokens and suggests regex patterns when pairing fails.
* Provides a pairing preview so you can see exactly what matched (and what didn’t).
* Supports both token‑based and regex‑based pairing so you can adapt to different lab naming conventions.

If filenames don’t encode sample identity and direction, any pairing attempt would be guesswork. MicroSeq’s design choice favors reproducibility and safety over silent, error‑prone inference.

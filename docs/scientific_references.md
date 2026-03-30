---
layout: page
title: Scientific References & Further Reading
permalink: /scientific-references/
---

## What this page is for

This page collects the research literature that directly supports MicroSeq’s methods, terminology, and design decisions.

It is not an exhaustive bibliography of everything related to Sanger sequencing or 16S analysis.
Instead, it is a curated reading list showing which papers matter most for understanding how MicroSeq works and why certain choices were made.

## Papers cited directly in the docs

| Topic | Citation | Why it matters to MicroSeq | Used in docs |
| --- | --- | --- | --- |
| Sanger base-calling / Phred | Ewing B, Hillier L, Wendl MC, Green P. Base-calling of automated sequencer traces using phred. I. Accuracy assessment. *Genome Research* (1998). | Grounds Phred-style quality interpretation and why per-base quality matters in Sanger workflows. | `gui_usage.md`, `paired_assembly.md` |
| Sanger base-calling / error probabilities | Ewing B, Green P. Base-calling of automated sequencer traces using phred. II. Error probabilities. *Genome Research* (1998). | Supports the interpretation of Phred as an error-probability scale. | `paired_assembly.md` |
| CAP3 | Huang X, Madan A. CAP3: A DNA sequence assembly program. *Genome Research* (1999). | Primary literature for the CAP3 assembler used in MicroSeq’s paired-read workflow. | `paired_assembly.md`, `design_algorithm_choices.md` |
| Paired-read overlap merging | Magoč T, Salzberg SL. FLASH: fast length adjustment of short reads to improve genome assemblies. *Bioinformatics* (2011). | Helpful for understanding the general overlap-then-merge paradigm. | `design_algorithm_choices.md` |
| Paired-read overlap merging | Zhang J, Kobert K, Flouri T, Stamatakis A. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. *Bioinformatics* (2014). | Useful comparison point for overlap-based merging logic. | `design_algorithm_choices.md` |
| Edit-distance alignment | Šošić M, Šikić M. Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. *Bioinformatics* (2017). | Relevant to indel-aware overlap backends and exact edit-distance alignment behavior. | `design_algorithm_choices.md`, `gui_table_reference.md` |
| Global alignment | Needleman SB, Wunsch CD. A general method applicable to the search for similarities in the amino acid sequence of two proteins. *J Mol Biol* (1970). | Foundational global alignment framework behind pairwise alignment thinking. | `design_algorithm_choices.md` |
| Local alignment | Smith TF, Waterman MS. Identification of common molecular subsequences. *J Mol Biol* (1981). | Foundational local alignment framework behind pairwise alignment reasoning. | `design_algorithm_choices.md` |
| Biopython | Cock PJA et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* (2009). | Relevant to Biopython-backed alignment and sequence-processing components. | `design_algorithm_choices.md` |



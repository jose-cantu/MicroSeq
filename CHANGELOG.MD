# Changelog 

All notable changes in MicroSeq I will keep track and documented here for the forseable future.

Current MicroSeq Setup: Version 1.0 

Any future version with changes/features will be below. 

Note due to feedback from colleagues in the lab I'm in I will make some revision and will call the features/adjustments in 1.1 - so stay tuned. =) 

----------------------------------------------------------

## [Unreleased] 



### Added Feature 

- Compare-all workflow now writes and auto-loads `asm/compare_assemblers.tsv` into the GUI results tables.
- Assembly reporting now surfaces primer source classification (`off|preset|custom|mixed`) alongside primer mode/stage.
- GUI output tabs now include **Compare Assemblers** to display `asm/compare_assemblers.tsv` after comparison runs.
- Primer-trim configuration now supports `mode: off|detect|clip` and `stage: pre_quality|post_quality`, with migration from legacy `primer_trim.enabled`.
- Added primer detect-only reporting (`qc/primer_detect_report.tsv`) and GUI **Primer preview** action for dry-run primer hit telemetry.
- Paired assembly summary now includes canonical cross-path fields (`assembler`, selected/configured engine, overlap metadata, `contig_len`, payload, primer mode/stage).
- Added reporting toggles in config (`reporting.emit_engine_audit`, `reporting.emit_per_sample_merge_report`) to keep debug artifacts optional.
- Primer preview now consumes the same GUI primer policy overrides (stage/preset/custom sequences) and forces detect mode for parity with run behavior.
- Detect mode now avoids writing primer-trim FASTQ copies, emitting only detect telemetry reports by default.
- Added GUI custom primer sequence entry (forward/reverse, one per line) with validation for IUPAC DNA alphabet and minimum length.
- Canonical `assembly_summary.tsv` now derives selected/configured engine fields from merge/config truth with audit fields separated as `audit_*`.

- Session-based logging; export `MICROSEQ_SESSION_ID` to aggregate logs.
- GUI toggle for speed vs. sensitivity BLAST algorithms.
- Fallback to classic `blastn` when identity or coverage falls below 90%.
- Added pytests covering the new logging behaviour.
- New `recommend-threads` command suggests a safe default for `--threads`.
- CLI `--session-id` overrides `MICROSEQ_SESSION_ID`.
- `microseq-setup` templates for zsh and fish shells.
- Cross-platform CI verifies Linux, macOS and Windows.
- Paired CAP3 assembly end-to-end: CLI/GUI support for `--mode paired`, primer presets or custom `--fwd-pattern/--rev-pattern`, duplicate policies (error/keep-first/keep-last/merge/keep-separate), plate well enforcement, pairing preview diagnostics, and per-sample contig outputs with `contig_map.tsv` manifest 
- GUI primer set dropdown plus advanced regex override for pairing. primer token detection and suggestions when no pairs are found. 
- Paired assembly now normalizes FASTA/FASTQ/AB1 inputs to per-read FASTA, re-writes contig headers with sample-aware IDs, and records detector/well metadata in pairing reports. 
- Documented paired vs single-mode CAP3 usage with real-world scenarios used under tests for data samples (plate runs, ad-hoc single mode) in README and reference in GUI/doc tutorials.
- Trim summaries now report expected error telemetry (EE/kb) and qeff based lebels. 
- Paired assembly now supports CAP3 profiling presets (strict/diagnostic/relaxed), QUAL usuage, overlap auditing, and blast input manifests (`asm/blast_inputs.fasta` + `asm/blast_inputs.tsv`) 
- GUI paired assembly panels now include CAP3 profile selection, BLAST input toggles, overlap audit display, and summary/detail views for paired outputs. 

### Changes Made 
- Simplified `microseq-setup` wizard, removing obsolete logging mode.
- Updated README with new logging instructions.
- Default log `backup_count` set to 0 and new `session_env` config key.
- Algorithms list now shows correctly via wrapper on second run.
- Logging files now open in append mode for multi-process runs.
- Buffering tool detection prefers `unbuffer` or `stdbuf` depending on OS.
- Added note about log files on NFS in the README.
- Relaxed plate well regex detection in handle wells adjacent to underscores/dashes and added validation coverage. 
- Added default contig map manifest creation during paired contig merges and rewrote cap3 with sample-aware identifiers to reference back to files being used in assembly 
- Paired assembly task scheduling now skips empty task lists and respects duplicate policies when creating cap3 jobs. 
- Removed outdated GUI warnings about blastn instability and pipeline workarounds that no longer apply
- QC trim summaries include EE/kb telemetry, qeff labels, and align per-read stats with expected-error calculations.
- Included vsearch replicate collapse per sample and reference chimera filtering in the full pipeline that you can use optionally (like a toggle in the gui to include it), with config backed references per DB (SILVA & GG2) w/ recommendations to use SILVA as the default. 
- vsearch docs and tests added for derep / uchime plus replicate weights tsv applied in BIOM. 
- Added CLI workflow scripts that will make it easier to use the CLI and provided examples in the /docs how to use them in the test sample datasets. 
- Paired contigs outputs now support `paired_contigs.fasta` which is contigs only or `blast_inputs.fasta` (which is the contigs + singlets fallback) for BLAST, with payload tracebility in `asm/blast_inputs.tsv`. 

### Fixed Bugs 
- Percent signs escape properly for `argparse`.
- Second BLAST run now wraps the command correctly.
- Patched blastdbcmd so that it can see the fullpath and avoid MemoryMap Error on wizard command microseq-setup
- Streamed BLAST output correctly on all platforms using `stdbuf` or `unbuffer`.
- Prevented paired assembly helpers from importing with indentation errors; corrected fwd/rev parameter annotations. 
- Fixed duplicate policy UI wiring (advanced regex checkbox, reverse token edits) and missing tokens dialog so GUI control renders and presists correctly. 
- Corrected pairing metadata handling so wells and computed pair keys are retained. Ensured missing-well tracking initiates before file processing. Duh. 
- Fixed paired preview helper signature and added well-pattern parameters; primer token scanning uses contiguous regex to avoid missing matches. 
- Resolved GUI dry run toggle glitches and missing message text in dialog box. 
- Perserved original `qseqid` while adding `sample_id` during taxonomy join so downstream merges do not drop identifiers. 


---

Example Format:
## [1.0] - 2025-05-25
### Fixed 
- Fixed GUI settings file with wrong JSON key 

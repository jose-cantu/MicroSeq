# Changelog 

All notable changes in MicroSeq I will keep track and documented here for the forseable future.

Current MicroSeq Setup: Version 1.0 

Any future version with changes/features will be below. 

Note due to feedback from colleagues in the lab I'm in I will make some revision and will call the features/adjustments in 1.1 - so stay tuned. =) 

----------------------------------------------------------

## [Unreleased] 



### Added Feature 

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



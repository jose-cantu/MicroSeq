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

### Changes Made 
- Simplified `microseq-setup` wizard, removing obsolete logging mode.
- Updated README with new logging instructions.
- Default log `backup_count` set to 0 and new `session_env` config key.
- Algorithms list now shows correctly via wrapper on second run.


### Fixed Bugs 
- Percent signs escape properly for `argparse`.
- Second BLAST run now wraps the command correctly.
- Patched blastdbcmd so that it can see the fullpath and avoid MemoryMap Error on wizard command microseq-setup 


---

Example Format:
## [1.0] - 2025-05-25
### Fixed 
- Fixed GUI settings file with wrong JSON key 



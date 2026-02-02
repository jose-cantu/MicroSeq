# MicroSeq/src/microseq_tests/assembly/cap3_profiles.py 

from __future__ import annotations 

"""CAP3 profile defaults used to build command line arguments."""

from typing import Sequence

CAP3_PROFILES: dict[str, list[str]] = {
    # HIGH CONFIDENCE: No clipping (-k 0), rely on MicroSeq's own trimming.
    "strict":     ["-o", "100", "-p", "95", "-s", "900", "-k", "0", "-y", "100"],

    # DIAGNOSTIC: Same as strict but enables CAP3's internal clipping (-k 1).
    # This helps identify if CAP3's clipping engine solves an assembly failure.
    "diagnostic": ["-o", "100", "-p", "95", "-s", "900", "-k", "1", "-y", "100"],

    # RESCUE: Modestly relaxes overlap length/identity based on UGENE defaults.
    "relaxed":    ["-o", "40",  "-p", "90", "-s", "900", "-k", "0", "-y", "100"],
}

def resolve_cap3_profile(profile: str, extra_args: Sequence[str] | None = None) -> list[str]:
    """Return CAP3 args for the named profile, optionally appending extra args."""
    if profile not in CAP3_PROFILES:
        raise ValueError(f"Unknown CAP3 profile: {profile}")
    args = list(CAP3_PROFILES[profile])
    if extra_args:
        args.extend(extra_args)
    return args

# src/microseq_tests/assembly/pairing.py 

"""
Main purpose is to scan a folder full of DNA sequence files which often come in pairs and identify only complete pairs. A complete pair meaning a forward and a reverse file for the same sample. It discards any samples that are missing either the F or the R file. Its meant to be a full registery system to find the correct file regardless of naming convention if forwhatever reason it doesn't then you can add your custom regex to the resitery below. 
"""

from __future__ import annotations # Postpones evaluation of type annotations (PEP 563) so they are no longer evaluated at function definition time - treated as string instead first 
import re # For regular expressions (text pattern matching) 
from pathlib import Path # Handling file systems 
import logging # Print warning messages 
from collections import defaultdict # Creating specialized dictionaries dealing with None Values 
from typing import Callable, Sequence, Union # FOr creating speicific type hints 
from enum import Enum # Creating enumerable, constant values 

# Inheriting from 'str' and "Enum" allows members to be compared directly to strings. For example, DupPolicy.Error = "error" will be True. 
class DupPolicy(str, Enum):
    """Defines the policy for handling duplicate files for the same sample/orientation.""" 
    ERROR = "error"
    KEEP_FIRST = "keep-first"
    KEEP_LAST = "keep-last"
    MERGE = "merge" 
    KEEP_SEPARATE = "keep-separate" 

"""
Here creating a type alias named Detector It defines the blueprint for any function that takes a string filename and returns a tuple of string, string or none. Making the code easier to read for future me. =) 
"""

Detector = Callable[[str], tuple[str, str | None]] # Calling a function ["string filename" and returning a tuple that has "string sampleID and F/R orientation OR None"  

# ----------- Detectors ---------------------------
# Pre-compiled regular expressions that will be used to find primer tokens within filenames automatically note I added re.I for case-insensitive. 

#_MID_RX: Finds a word like "_27F_" or "-1492R-" in the middle of a filename. 
_MID_RX = re.compile(r"[_\-]([A-Za-z0-9]+[FR])[_\-]", re.I) 
# _PREFIX_RX Finds a word like "27F_" or "1492R-" at the very start of filename.
_PREFIX_RX = re.compile(r"^([A-Za-z0-9]+[FR])[_\-]", re.I) 
# _SUFFIX_RX Finds a word like "_27F.fasta" or "-1492R.fasta" at the end of a filename 
_SUFFIX_RX = re.compile(r"[_\-]([A-Za-z0-9]+[FR])\.[^.]+$", re.I)  

def mid_token_detector(name: str):
    """Detects primer tokens in the middle of a filename."""
    # Search the filename for the middle-token pattern. 
    if m := _MID_RX.search(name):
        # If a match is found, extract the token/word (e.g. "27F") and orientation ("F") plus convert always to uppercase to standardize. 
        tok, orient = m[1], m[1][-1].upper() 
        # The sample ID is assumed to be everything befoer the token.
        sid = name.split(f"{tok}", 1)[0].rstrip("_-")
        return sid, orient 
    # If no match found, return the filename stem and None for orientation. 
    return Path(name).stem, None 


def prefix_detector(name: str):
    """Detects primer tokens at the start of a filename"""
    # Match the filename for the prefix-token pattern I've set up 
    if m := _PREFIX_RX.match(name):
        # Extract the token (Primer name) and orientation 
        tok, orient = m[1], m[1][-1].upper() 
        # The sample ID is everything after the token (hence prefix = 1st) 
        sid = name.split(f"{tok}", 1)[1].lstrip("_-").split(".",1)[0] 
        return sid, orient 
    # If no match, return the default 
    return Path(name).stem, None 

def suffix_detector(name: str):
    """Detects primer tokens at the end of a filename"""
    # Search the filename for the suffix pattern 
    if m := _SUFFIX_RX.search(name):
        # Extract the token (the primer name) and orientation
        tok = m[1] # example 27F or 1492R 
        orient = tok[-1].upper() # F or R    
        # The sample ID is everything before the token. Handling both dash and underscore scenarios
        sid = name[:m.start()].rstrip("_-") 
        return sid, orient 
    # If no match then return default (Just sample ID with no orient/token) 
    return Path(name).stem, None 


# This list acts as a registry for all available detector functions,
# The order is important, as they will be tried sequentially. 
DETECTORS: list[Detector] = [
    mid_token_detector,
    prefix_detector,
    suffix_detector,
]

def _strip_token(name: str, token: str) -> str:
    """Remove token from name and tidy separators/extension."""
    base = name
    token_idx = base.upper().find(token.upper())
    if token_idx >= 0:
        before = base[:token_idx].rstrip("_-")
        after = base[token_idx + len(token):].lstrip("_-.")
        cleaned = f"{before}_{after}" if before and after else before or after
    else:
        cleaned = base

    return Path(cleaned).stem


def make_pattern_detector(fwd_pattern: str, rev_pattern: str) -> Detector:
    """Build a detector that searches for custom forward/reverse regex tokens."""

    fwd_rx = re.compile(f"({fwd_pattern})", re.I)
    rev_rx = re.compile(f"({rev_pattern})", re.I)

    def detector(name: str) -> tuple[str, str | None]:
        if m := fwd_rx.search(name):
            token = m[1]
            return _strip_token(name, token), "F"
        if m := rev_rx.search(name):
            token = m[1]
            return _strip_token(name, token), "R"
        return Path(name).stem, None

    return detector


# --------------- Public Helpers ---------------------------
def extract_sid_orientation(name: str, *, detectors: Sequence[Detector] | None = None) -> tuple[str, str | None]: 
    """Try each detector until one recognises a primer token."""
    for det in DETECTORS:
        sid, orient = det(name)
        if orient in ("F", "R"):
            return sid, orient
    # Loop exhausted -> no primer recognized 
    return Path(name).stem, None 

def group_pairs(
    folder: Union[str, Path], 
    dup_policy: DupPolicy = DupPolicy.ERROR,
    *, fwd_pattern: str | None = None, rev_pattern: str | None = None, detectors: Sequence[Detector] | None = None
    ) -> dict[str, dict[str, Union[Path, list[Path]]]]:
    """
    Groups Forward/Reverse reads from a folder into pairs 
    """
    # Special dictionary will automatically create a nested dictionary 
    # first time a new sample ID (sid) is encountered.
    """
    # how it looks like:
    # {
  'SampleA': {  # Outer key (str)
    'F': Path('/path/to/file_F.fasta'),  # Inner key (str) and value (Path)
    'R': Path('/path/to/file_R.fasta')
  },
   """ 
    pairs: dict[str, dict[str, Union[Path, list[Path]]]] = defaultdict(dict) 
     
    active_detectors: Sequence[Detector]
    if detectors is not None:
        active_detectors = detectors
    elif fwd_pattern and rev_pattern:
        active_detectors = [make_pattern_detector(fwd_pattern, rev_pattern), *DETECTORS]
    else:
        active_detectors = DETECTORS


    # Use iterdir() to look at every item in the folder. 
    for p in Path(folder).iterdir():
        # Check the file extension in a case-insensitive way. 
        # If it's not a recognized FASTA extension, skip to the next file. 
        if p.suffix.lower() not in {".fasta", ".fa", ".fna"}:
            continue

        # Get the sample ID and orientation ('F', 'R', or None)
        sid, orient = extract_sid_orientation(p.name, detectors=active_detectors)

        # If no valid orientation was found, skip to the next file. 
        if orient not in ("F", "R"):
            continue 

        # This here is a get or set operation. It looks for 'pairs[sid][orient]
        # If slot exists, it returns its value (bucket) 
        # If slot not exists, sets to None and returns None. 
        # Prevents KeyError and initialized the slot in one step. 
        bucket = pairs[sid].get(orient)  

        # If buck is 'None', it means this is the first time 
        # For this combination sampleID and orientation. 
        if bucket is None: 
            # store the file path in its designated slot. 
            pairs[sid][orient] = p 
        else:
            # A file for this sample/orientation already exists, so we have a duplicate in place. Now is the time for user-selected policy. 

            # Policy 1: Raise an error and stop the program MicroSeq. Safe 
            # Default 
            if dup_policy == DupPolicy.ERROR:
                raise ValueError(f"Duplicate {orient} read for sample {sid}: found {p} but {bucket} already exists.") 


            # Policy 2: Log a warning but make no changes, keep first only 
            elif dup_policy == DupPolicy.KEEP_FIRST:
                logging.warning("Duplicate %s/%s ignored due to 'keep-first' policy: %s", sid, orient, p) 

            # Policy 3: Log a warning and overwrite the old path with the new one 
            elif dup_policy == DupPolicy.KEEP_LAST:
                logging.warning("Duplicate %s/%s overwriting previous entry %s due to 'keep-last' policy.", sid, orient, bucket) 
                pairs[sid][orient] = p 

            # Policy 4: Collect all duplicate paths into a list 
            elif dup_policy in (DupPolicy.MERGE, DupPolicy.KEEP_SEPARATE):
                # Step 1: Prepare the list. If 'bucket` is already a list, use it 
                # If 'bucket' (the old stored path) is a single path object, create a new list containing it. 
                lst = bucket if isinstance(bucket, list) else [bucket] 
                # Append the path of the new duplicate file to the list 
                lst.append(p)
                # Put the updated list back into the dictionary 
                pairs[sid][orient] = lst 
    # Finally filtering pairs keeping sample IDs ('s') have "F' and "R' keys 
    return {s: d for s, d in pairs.items() if {"F", "R"} <= d.keys()}

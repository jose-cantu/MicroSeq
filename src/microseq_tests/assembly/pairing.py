# src/microseq_tests/assembly/pairing.py 

"""
Main purpose is to scan a folder full of DNA sequence files which often come in pairs and identify only complete pairs. A complete pair meaning a forward and a reverse file for the same sample. It discards any samples that are missing either the F or the R file. Its meant to be a full registery system to find the correct file regardless of naming convention if forwhatever reason it doesn't then you can add your custom regex to the resitery below. 
"""

from __future__ import annotations
import re # For regular text expressions (text pattern matching) 
from pathlib import Path # Cool way of handling file systems 
import logging # For printing warning messages 
from collections import defaultdict # Creating special dictionaries dealing with None Values 
from typing import Callable, Union # For creating specific type hints here...
from enum import Enum # For creating enumerable, constant values 

# Inheriting from 'str' and 'Enum' allows members to be compared directly to strings. For exmaple, DupPolicy.ERROR == "error" will be True. 
class DupPolicy(str, Enum):
    """Defines the policy for handling duplicate files for the same sample/orientation.""" 
    ERROR = "error" 
    KEEP_FIRST = "keep-first"
    KEEP_LAST = "keep-last"
    MERGE = "merge" 

"""
Here creating a type alias named Detector It defines the blueprint for any function that takes a string filename and returns a tuple of string, string or none. Making the code easier to read for future me. =) 
"""

Detector = Callable[[str], tuple[str, str | None]] 


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
        orient = m[1], m[1][-1].upper()   
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

# --------------- Public Helpers ---------------------------
def extract_sid_orientation(name: str) -> tuple[str, str | None]: 
    """Try each detector until one recognises a primer token."""
    for det in DETECTORS:
        sid, orient = det(name)
        if orient in ("F", "R"):
            return sid, orient
    # Loop exhausted -> no primer recognized 
    return Path(name).stem, None 

def group_pairs(
    folder: Union[str, Path], 
    dup_policy: DupPolicy = DupPolicy.ERROR 
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
     
    # Use iterdir() to look at every item in the folder. 
    for p in Path(folder).iterdir():
        # Check the file extension in a case-insensitive way. 
        # If it's not a recognized FASTA extension, skip to the next file. 
        if p.suffix.lower() not in {".fasta", ".fa", ".fna"}:
            continue

        # Get the sample ID and orientation ('F', 'R', or None)
        sid, orient = extract_sid_orientation(p.name)

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
            elif dup_policy == DupPolicy.MERGE:
                # Step 1: Prepare the list. If 'bucket` is already a list, use it 
                # If 'bucket' (the old stored path) is a single path object, create a new list containing it. 
                lst = bucket if isinstance(bucket, list) else [bucket] 
                # Append the path of the new duplicate file to the list 
                lst.append(p)
                # Put the updated list back into the dictionary 
                pairs[sid][orient] = lst 
    # Finally filtering pairs keeping sample IDs ('s') have "F' and "R' keys 
    return {s: d for s, d in pairs.items() if {"F", "R"} <= d.keys()}

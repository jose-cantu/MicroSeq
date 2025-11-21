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
from Bio import SeqIO 

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
# Match plate positions (A01-H12) that may sit next to dashes/underscores 
# rather than relying on word boundaries 
_WELL_RX = re.compile(r"(?i)(?<![A-Z0-9])([A-H](?:0?[1-9]|1[0-2]))(?![A-Z0-9])") 

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
            sid = name.split(token, 1)[0].rstrip("_-")
            return sid or Path(name).stem, "F"
        if m := rev_rx.search(name):
            token = m[1]
            sid = name.split(token, 1)[0].rstrip("_-")
            return sid or Path(name).stem, "R" 
        return Path(name).stem, None

    return detector

def _extract_well(name: str, *, pattern: str | re.Pattern[str] | None = None) -> str | None: 
    """ Return a normalized plate well code (example "A07") if present."""

    if pattern is None:
        rx = _WELL_RX
    elif isinstance(pattern, str):
        rx = re.compile(pattern, re.I)
    else:
        rx = pattern 

    if m := rx.search(name):
        raw = m[1]
        letter, num = raw[0].upper(), int(raw[1:])
        return f"{letter}{num:02d}"

    return None 

def _pair_key(sid: str, well: str | None, enforce_same_well: bool) -> str: 
    """Derive the key used to bucket forward/reverse reads."""

    if enforce_same_well and well:
        return well 
    return sid 

# --------------- Public Helpers ---------------------------
def extract_sid_orientation(name: str, *, detectors: Sequence[Detector] | None = None) -> tuple[str, str | None]: 
    """Try each detector until one recognises a primer token."""
    active_detectors = detectors or DETECTORS

    for det in active_detectors: 
        sid, orient = det(name)
        if orient in ("F", "R"):
            return sid, orient
    # Loop exhausted -> no primer recognized 
    return Path(name).stem, None 

def _detect_sid_orientation(
    name: str, detectors: Sequence[Detector]
) -> tuple[str, str | None, str | None]:
    """Return (sid, orientation, detector_name) for the first detector that matches."""

    for det in detectors:
        sid, orient = det(name)
        if orient in ("F", "R"):
            det_name = getattr(det, "__name__", det.__class__.__name__)
            return sid, orient, det_name

    return Path(name).stem, None, None

def group_pairs(
    folder: Union[str, Path], 
    dup_policy: DupPolicy = DupPolicy.ERROR,
    *, fwd_pattern: str | None = None, rev_pattern: str | None = None, detectors: Sequence[Detector] | None = None, return_metadata: bool = False, enforce_same_well: bool = False, well_pattern: str | re.Pattern[str] | None = None
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
    meta: dict[str, dict[str, Union[str, list[str]]]] = defaultdict(dict)
     
    active_detectors: Sequence[Detector]
    if detectors is not None:
        active_detectors = detectors
    elif fwd_pattern and rev_pattern:
        active_detectors = [make_pattern_detector(fwd_pattern, rev_pattern), *DETECTORS]
    else:
        active_detectors = DETECTORS

    well_rx = _WELL_RX if well_pattern is None else well_pattern 

    def _store_entry(key: str, sid: str, orient: str, path: Path, detector_name: str | None, well: str | None) -> None:
        bucket = pairs[key].get(orient)
        meta_bucket = meta[key].get(orient)

        if bucket is None:
            pairs[key][orient] = path
            meta[key][orient] = detector_name or "unknown"
            if well:
                meta[key][f"well_{orient}"] = well 
            return

        if dup_policy == DupPolicy.ERROR:
            raise ValueError(f"Duplicate {orient} read for sample {sid}: found {path} but {bucket} already exists.")

        if dup_policy == DupPolicy.KEEP_FIRST:
            logging.warning("Duplicate %s/%s ignored due to 'keep-first' policy: %s", sid, orient, path)
            return

        if dup_policy == DupPolicy.KEEP_LAST:
            logging.warning("Duplicate %s/%s overwriting previous entry %s due to 'keep-last' policy.", sid, orient, bucket)
            pairs[key][orient] = path
            meta[key][orient] = detector_name or "unknown"
            if well:
                meta[key][f"well_{orient}"] = well 
            return

        if dup_policy in (DupPolicy.MERGE, DupPolicy.KEEP_SEPARATE):
            lst = bucket if isinstance(bucket, list) else [bucket]
            lst.append(path)
            pairs[key][orient] = lst
            meta_lst = meta_bucket if isinstance(meta_bucket, list) else [meta_bucket]
            meta_lst.append(detector_name or "unknown")
            meta[key][orient] = meta_lst
            if well:
                wells = meta[key].get(f"well_{orient}")
                if wells is None:
                    meta[key][f"well_{orient}"] = [well] 
                elif isinstance(wells, list):
                    wells.append(well)
                else:
                    meta[key][f"well_{orient}"] = [wells, well] 


    path = Path(folder)
    fasta_exts = {".fasta", ".fa", ".fna"}

    missing_well: list[str] = []

    if path.is_file():
        if path.suffix.lower() not in fasta_exts:
            raise ValueError(f"Unsupported FASTA input: {path}")
  
        record_counts: dict[tuple[str, str], int] = defaultdict(int)
        tmp_dir = path.parent / f".{path.stem}_paired_records"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        for record in SeqIO.parse(path, "fasta"):
            sid, orient, det_name = _detect_sid_orientation(record.id, active_detectors) 
            if orient not in ("F", "R"):
                continue
            # if well is missing doc it 
            well = _extract_well(record.id, pattern=well_rx) if enforce_same_well else None 
            if enforce_same_well and not well:
                missing_well.append(record.id)
                continue 

            record_counts[(sid, orient)] += 1
            rec_idx = record_counts[(sid, orient)]
            safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", record.id)
            rec_path = tmp_dir / f"{safe_id}_{rec_idx}.fasta"
            SeqIO.write([record], rec_path, "fasta")
            key = _pair_key(sid, well, enforce_same_well) 
            _store_entry(key,sid, orient, rec_path, det_name, well)

    if path.is_dir(): 
        for p in path.iterdir():
            # Check the file extension in a case-insensitive way.
            # If it's not a recognized FASTA extension, skip to the next file.
            if p.suffix.lower() not in fasta_exts:
                continue

            # Get the sample ID and orientation ('F', 'R', or None)
            sid, orient, det_name = _detect_sid_orientation(p.name, detectors=active_detectors)

            # If no valid orientation was found, skip to the next file.
            if orient not in ("F", "R"):
                continue

            well = _extract_well(p.name, pattern=well_rx) if enforce_same_well else None 
            if enforce_same_well and not well: 
                missing_well.append(p.name)
                continue 

            key = _pair_key(sid, well, enforce_same_well)
            _store_entry(key, sid, orient, p, det_name, well) 
    else:
        raise ValueError(f"Input path must be a FASTA file or directory: {path}")

        # Finally filtering pairs keeping sample IDs ('s') have "F' and "R' keys
    paired_only = {s: d for s, d in pairs.items() if {"F", "R"} <= d.keys()}
    meta_only = {s: meta[s] for s in paired_only}
    if enforce_same_well and missing_well:
        meta_only["_missing_well"] = {"files": missing_well}

    if return_metadata:
        return paired_only, meta_only

    return paired_only


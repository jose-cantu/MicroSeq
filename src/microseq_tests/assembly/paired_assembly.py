# src/microseq_tests/assembly/paired_assembly.py 

"""
Utilities for running CAP3 assembly on forward/reverse pairs.
""" 
from __future__ import annotations # Postpones evaluation of type annotations (PEP 563) so they are no longer evaluated at function definition time - treated as string instead first 
import logging # print warning messages  
from os import PathLike
import subprocess 
from pathlib import Path 
from typing import Iterable, Sequence

from microseq_tests.utility.utils import load_config 

from .pairing import DupPolicy, group_pairs 

L = logging.getLogger(__name__) 

def _iter_paths(value: Path | list[Path]) -> Iterable[Path]:
    """Normalize a stored path entry such as a list into an iterator."""
    if isinstance(value, list):
        for item in value: 
            yield Path(item)
    else:
            yield Path(value)  

def _as_path_list(value: Path | list[Path]) -> list[Path]:
    """Return a list of `path` objects perserving insertion order."""
    return [Path(p) for p in _iter_paths(value)]

def _build_keep_separate_pairs(forward: list[Path], reverse: list[Path]) -> list[tuple[Path, Path]]:
    """Produce forward/reverse pairings for the keep-separate policy."""
    if not forward or not reverse: 
        raise ValueError("Forward and reverse inputs must not be empty when paring reads.") 

    f_count = len(forward)
    r_count = len(reverse)

    # If both orientations contains multiple primers, require the counts to align
    # so that we can preserve primer ordering w/o inventing new orientations 
    if f_count >1 and r_count > 1 and f_count != r_count:
        raise ValueError("Mismatched duplicate counts detected while using 'keep-separate' policy:" 
        f"{f_count} forward vs {r_count} reverse reads." 
        ) 

    pairs: list[tuple[Path, Path]] = []
    if f_count == r_count:
        pairs = list(zip(forward, reverse))
    elif f_count > r_count:
        if r_count != 1: 
            raise ValueError("Unable to align forward duplicates with revcerse reads under 'keep-separate' policy.")
        pairs = [(f_path, reverse[0]) for f_path in forward]
    else: # r_count > f_count 
        if f_count != 1:
            raise ValueError(
                "Unable to align reverse duplicates with forward reads under 'keep-separate' policy.")
        pairs = [(forward[0], r_path) for r_path in reverse] 

    return pairs 



def _write_combined_fasta(sources: Iterable[Path], destination: Path) -> None: 
    """ Here I will combine multiple FASTA files into 'destination'."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as fout:
        for src in sources:
            data = Path(src).read_text(encoding="utf-8")
            # check if the data file is empty or not 
            if not data:
                # if its empty skip and continue 
                continue 
            # write data to file 
            fout.write(data)
            # make sure it adheres to FASTA format and check if newline if not newline will make one for next sample 
            if not data.endswith("\n"):
                fout.write("\n") 

def assemble_pairs(input_dir: PathLike, output_dir: PathLike, *, dup_policy: DupPolicy = DupPolicy.ERROR, cap3_options: Sequence[str] | None = None, fwd_pattern: str | None = None, rev_pattern: str | None = None  
                   ) -> list[Path]:
    """
    Run CAP3 assemblies for each forward and reverse pair discovered in ``input_dir``. 
    
    Parameters:
    ----------
    input_dir:
        Directory containiing per-sample FASTA files. Files are paired via 
        via :func:`microseq_tests.assembly.pairing.groups_pairs`.
    output_dir:
        Destination directory for per-sample CAP3 runs. Each sample receives a subdirectory containing 
        the concatenated FASTA and CAP3 outputs.
    dup_policy:
        Policy for handling duplicate forward/reverse files when pairing. 
        See :class: `microseq_tests.assembly.pairing.DupPolicy`.
    cap3_options:
        Optional additional command-line arguments appended to CAP3 call. 

    Return:
    ------
    
    list[Pathlib.Path]

    A list of paths to the generated ``*.cap.contigs`` files, ordered by
        sample identifier. 
        list[Pathlib.Path]

    """
    cfg = load_config() 
    cap3_exe = cfg["tools"]["cap3"]

    in_dir = Path(input_dir).resolve() 
    out_dir = Path(output_dir).resolve() 
    out_dir.mkdir(parents=True, exist_ok=True) 

    pairs = group_pairs(in_dir, dup_policy=dup_policy, fwd_pattern=fwd_pattern, rev_pattern=rev_pattern)
    if not pairs:
        L.info("No paired reads detected in input path: %s", in_dir)
        # Exit out since nothing here to do 
        return [] 
    
    # goign to collect the samples started with an empty list 
    contig_paths: list[Path] = [] 
    # Will process the ordered sampled paired list one by one 
    for sid in sorted(pairs):
        # Find the location of the orientation based on the file name 
        entries = pairs[sid]

        tasks: list[tuple[str, list[Path]]] = [] 
        if dup_policy == DupPolicy.KEEP_SEPARATE:
            f_sources = _as_path_list(entries["F"])
            r_sources = _as_path_list(entries["R"])
            primer_pairs = _build_keep_separate_pairs(f_sources, r_sources)
            for idx, (fwd, rev) in enumerate(primer_pairs, start=1):
                sample_key = sid if len(primer_pairs) == 1 else f"{sid}_{idx}"
                tasks.append((sample_key, [fwd, rev]))
            else:
                sources = list(_iter_paths(entries["F"])) + list(_iter_paths(entries["R"]))
                tasks.append((sid, sources))

            for sample_key, sources in tasks: 
                sample_dir = out_dir / sample_key 
                sample_dir.mkdir(parents=True, exist_ok=True) 
                sample_fasta = sample_dir / f"{sample_key}_paired.fasta" 

                _write_combined_fasta(sources, sample_fasta) 

                cmd = [cap3_exe, sample_fasta.name]
                if cap3_options:
                    cmd.extend(cap3_options)
                L.info(" Run CAP3 (apired) %s: %s", sample_key, "".join(cmd)) 

                try:
                    subprocess.run(
                        cmd, 
                        check=True,
                        cwd=sample_dir,
                        stderr=subprocess.PIPE,
                        text=True,
                    )
                except subprocess.CalledProcessError as exc:
                    L.error( "CAP3 failed for sample %s (exit %s):\n%s", sample_key, exc.returncode, exc.stderr
                    )
                    raise # stop the run here on raise 
                contig_path = sample_dir / f"{sample_key}_paired.fasta.cap.contigs"
                if not contig_path.exists():
                    raise FileNotFoundError(contig_path)

                contig_paths.append(contig_path)
                L.info("Cap3 paired assembly finished for %s and for contigs: %s", sample_key, contig_path) 
    
    return contig_paths 

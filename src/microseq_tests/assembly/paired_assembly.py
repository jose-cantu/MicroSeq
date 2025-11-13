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

def assemble_pairs(input_dir: PathLike, output_dir: PathLike, *, dup_policy: DupPolicy = DupPolicy.ERROR, cap3_options: Sequence[str] | None = None 
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

    pairs = group_pairs(in_dir, dup_policy=dup_policy)
    if not pairs:
        L.info("No paired reads detected in input directory: %s", in_dir)
        # Exit out since nothing here to do 
        return [] 
    
    # goign to collect the samples started with an empty list 
    contig_paths: list[Path] = [] 
    # Will process the ordered sampled paired list one by one 
    for sid in sorted(pairs):
        # Find the location of the orientation based on the file name 
        entries = pairs[sid] 
        # create sub directory for each assembly 
        sample_dir = out_dir / sid 
        # safty check 
        sample_dir.mkdir(parents=True, exist_ok=True) 
        sample_fasta = sample_dir / f"{sid}_paired.fasta" 
        # grab all F/R 
        sources = list(_iter_paths(entries["F"])) + list(_iter_paths(entries["R"]))
        # combine each F/R pair 
        _write_combined_fasta(sources, sample_fasta)
        
           # writing out manual commmand instructions and run it and raise error if failure
           # here only the file name is needed since I change workding directory to the sample directory its in! 
           cmd = [cap3_exe, sample_fasta.name] 
           # extend cap3 options commands
           if cap3_options:
               cmd.extend(cap3_options)
           L.info(" Run CAP3 (paired) %s: %s", sid, " ".join(cmd)) 
           # now try to run or error and cancel run 
           try: 
               subprocess.run(
                    cmd, # run the command 
                    check=True, # check if ran complete 
                    cwd=sample_dir, # change working direcotry 
                    stderr=subprocess.PIPE, # pipe the error captured in except process 
                    text=True # return logging in text rather than byte form making it readable 
                ) 
           except subprocess.CalledProcessError as exc:
                L.error(
                    "CAP3 failed for sample %s (exit %s):\n%s", sid, exc.returncode, exc.stderr
                )
                raise # raise this issue and exit 

           contig_path = sample_dir / f"{sample_fasta}.cap.contigs" 
           if not contig_path.exists():
                raise FileNotFoundError(contig_path)

           contig_paths.append(contig_path)
           L.info("Cap3 paired assembly finished for %s and for contigs: %s", sid, contig_path) 

    return contig_paths 


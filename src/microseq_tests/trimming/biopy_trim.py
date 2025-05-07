# src/microseq_tests/trimming/biopy_trim.py 
from __future__ import annotations
from pathlib import Path
from typing import List, Optional, Iterable 
from Bio import SeqIO                         
import logging 
from importlib import import_module 

L = logging.getLogger(__name__)


# ------- helper functions ----------------------
def _tick_safe(bar, inc: int = 1) -> None: 
    """ Calls bar.update(inc) only if bar exposes that attribute."""
    if hasattr(bar, "update"):
        bar.update(inc)  # safety guard 


# ----------------------------------------------------------------------
# Single‑read sliding‑window trim
# ----------------------------------------------------------------------
def dynamic_trim(record, win: int = 5, q: int = 20):
    """
    Return trimmed SeqRecord or None if the read never reaches quality q. 
    """
    quals: List[int] = record.letter_annotations["phred_quality"]
    n = len(quals)
    if n < win:
        return None                            # too short

    # 5′ → 3′
    left = next(
        (i for i in range(n - win + 1)
         if sum(quals[i:i + win]) / win >= q),
        None,
    )
    if left is None:
        return None

    # 3′ → 5′
    right = next(
        (j for j in range(n, win - 1, -1)
         if sum(quals[j - win:j]) / win >= q),
        None,
    )
    if right is None or right <= left:
        return None

    return record[left:right]                 # trimmed SeqRecord


# ────────────────────────────────────────────────────────────────────────
# helper that does the real trimming work (keeps outer function compact)
# ────────────────────────────────────────────────────────────────────────
def _trim_all(
    fastqs: Iterable[Path],
    *,
    window_size: int,
    per_base_q: int,
    file_q_threshold: float,
    passed_dir: Path,
    failed_dir: Path,
    stats_root: Path,
    comb: Optional[open],
    bar,
) -> None:
    """
      Core loop: trims each FASTQ, writes per-read stats, dispatches files
    to passed_dir / failed_dir, and updates the progress bar once
    per input file.
    """ 
    for fq in fastqs:
        base         = fq.stem
        stats_path   = stats_root / f"{base}_avg_qual.txt"
        trimmed_path = passed_dir / f"{base}_trimmed.fastq"
        
        reads = bases = qsum = 0
        trimmed_recs = []
        
        for rec in SeqIO.parse(fq, "fastq"):
            trimmed = dynamic_trim(rec, window_size, per_base_q)
            if not trimmed:
                continue
            ph = trimmed.letter_annotations["phred_quality"]
            reads += 1
            bases += len(trimmed)
            qsum  += sum(ph)
            trimmed_recs.append(trimmed)
            
        avg_q   = qsum  / bases if bases else 0
        avg_len = bases / reads if reads else 0
        
        # per-read stats
        with open(stats_path, "w") as fh:
            for r in trimmed_recs:
                ph = r.letter_annotations["phred_quality"]
                fh.write(f"{r.id}\t{len(r)}\t{sum(ph)/len(ph):.2f}\n")
                
        # pass / fail per file
        if avg_q < file_q_threshold:
            (failed_dir / fq.name).write_bytes(fq.read_bytes())
            (failed_dir / stats_path.name).write_bytes(stats_path.read_bytes())
            L.info("[FAIL] %s  (avgQ %.2f)", fq.name, avg_q)
        else:
            SeqIO.write(trimmed_recs, trimmed_path, "fastq")
            L.info("[PASS] %s → %s (avgQ %.2f)", fq.name, trimmed_path, avg_q)
            if comb:
                comb.write(f"{fq.name}\t{reads}\t{avg_len:.1f}\t{avg_q:.2f}\n")
                
        _tick_safe(bar) # exactly one tick per FASTQ  
        
# ────────────────────────────────────────────────────────────────────────
# public folder-level driver
# ────────────────────────────────────────────────────────────────────────
def trim_folder(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    window_size: int = 5,
    per_base_q: int = 20,
    file_q_threshold: float = 20.0,
    combined_tsv: str | Path | None = None,
    threads: int = 1,                          # kept for API parity; unused
    **kwargs,
) -> None:
    """
    Quality-trim every *.fastq in input_dir and write results under
    output_dir:

      ├── passed_qc_fastq/
      ├── failed_qc_fastq/
      └── per-read stat files + combined TSV

    A single outer progress bar ticks once per FASTQ.  Works both with
    a real tqdm bar and with the dummy bar used by the test-suite’s
    monkey-patch.
    """
    input_dir  = Path(input_dir)
    output_dir = Path(output_dir)
    passed_dir = output_dir.parent / "passed_qc_fastq"
    failed_dir = output_dir.parent / "failed_qc_fastq"
    for p in (output_dir, passed_dir, failed_dir):
        p.mkdir(parents=True, exist_ok=True)
     
     # combined summary 
    comb: Optional[open] = None
    if combined_tsv:
        comb = open(combined_tsv, "a")
        if comb.tell() == 0:                   # new file → header
            comb.write("file\treads\tavg_len\tavg_q\n")
            
    fastqs = sorted(input_dir.glob("*.fastq"))
    
    # late-bind progress helper so pytest can monkey-patch it
    prog = import_module("microseq_tests.utility.progress")
    cm_or_gen   = prog.stage_bar(len(fastqs), desc="trim", unit="file")
    
    def _run(bar):
        _trim_all(
                fastqs,
                window_size=window_size,
                per_base_q=per_base_q,
                file_q_threshold=file_q_threshold,
                passed_dir=passed_dir,
                failed_dir=failed_dir,
                stats_root=output_dir,
                comb=comb,
                bar=bar,   # inside _trim_all still call bar.update 
            )
    if hasattr(cm_or_gen, "__enter__"):  # proper-context manager 
        with cm_or_gen as bar:
            _run(bar) 

    else:                        # generator or bare bar 
        try:
            bar = next(cm_or_gen) # unwrap generator 
        except StopIteration:
            bar = cm_or_gen 
        _run(bar) 

        if hasattr(cm_or_gen, "close"):
            cm_or_gen.close() 
        if hasattr(bar, "close") and bar is not cm_or_gen:
            bar.close()
            
    if comb:
        comb.close() 

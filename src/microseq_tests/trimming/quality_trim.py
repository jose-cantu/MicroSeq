# microseq_tests/src/microseq_tests/trimming/quality_trim.py 

from __future__ import annotations 
import logging
import argparse
import subprocess 
from pathlib import Path
from typing import Iterable
from Bio import SeqIO 
import shutil 

from microseq_tests.utility.utils import load_config, setup_logging 
from microseq_tests.trimming.expected_errors import expected_errors, qeff_from_mee_per_kb 
from microseq_tests.trimming.biopy_trim import _mee_qc_label 

L = logging.getLogger(__name__)

PathLike = str | Path 

def quality_trim(input_file: PathLike, output_file: PathLike, *, threads: int=1, **kwargs, ) -> Path:
    """
    Run Trimmomatic single-end mode and return the absolute path to trimmed FASTAQ/FASTA so down stream steps (assembly, GUI, tests) can chain without guessing. 
    threads : int 
    Forwarded to Trimmomatics's threads option. 
    kwards: dict 
    Currnely unused for now will reserve it for more trimmomatic patterns.... `

        """
    cfg = load_config()
    trimm = cfg["tools"]["trimmomatic"]

    cmd = [ 
        trimm, "SE", "-threads", str(threads), "-phred33",
        str(input_file), str(output_file),    # <<< cast once here 
        "SLIDINGWINDOW:5:20",
        "MINLEN:200",
        ]

    L.info("RUN Trimmomatic: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        L.error("Trimmomatic failed (exit %s):\n%s", e.returncode, e.stderr)
        raise 

    out_path = Path(output_file).resolve()
    if not out_path.exists():
        raise FileNotFoundError(out_path)
    return out_path 

def _looks_like_fastq(fp: Path) -> bool:
    """
    Heuristic: treat file names with suffixws containing '.fastq' or '.fq' (case-insentitive) 
    as FASTQ canidates. Handles compressed names like *.fastq.gz b/c '.fastq' will be present amongst suffixes. 
    """
    if not fp.is_file():
        return False 
    suffixes = {s.lower() for s in fp.suffixes}
    return ".fastq" in suffixes or ".fq" in suffixes 

def _fastq_stats(path: Path) -> tuple[int, float, float, float, float, float, str]:
    """
    Stream a trimmed FASTQ and compute basic QC metrics:
        returns (reads, avg_len, avg_q, avg_mee), where avg_q is the base-wise mean Phred (Phred+33) 
    """
    reads = bases = qsum = 0 
    mee_sum = 0.0 
    for rec in SeqIO.parse(path, "fastq"):

        ph = rec.letter_annotations["phred_quality"]
        reads += 1 
        bases += len(rec) 
        qsum += sum(ph)
        mee = expected_errors(ph)
        mee_sum += mee 
    avg_q = qsum / bases if bases else 0 
    avg_len = bases / reads if reads else 0 
    avg_mee = mee_sum / reads if reads else 0 
    avg_mee_per_kb = (1000 * mee_sum / bases) if bases else 0 
    avg_qeff = qeff_from_mee_per_kb(avg_mee_per_kb) 
    mee_qc_label = _mee_qc_label(avg_mee_per_kb) 
    return reads, avg_len, avg_q, avg_mee, avg_mee_per_kb, avg_qeff, mee_qc_label 
     

def _iter_fastq_sources(path: Path) -> Iterable[Path]:
    """
    Yield source FASTQ files to trime. 
    - If 'path' is a directory then yeield files that look like FASTQ in lexicographic order. 
    - If 'path' is a file and exists then yeild it (no additional validation needed). 
    """
    if path.is_dir():
        for fp in sorted(path.iterdir()):
            if _looks_like_fastq(fp):
                yield fp 
    elif path.exists():
        yield path

def trim_fastq_inputs(input_path: PathLike, trim_dir: PathLike, *, summary_tsv: PathLike | None = None) -> Path: 
    """
    Trim FASTQ input that may be a single file or a directory of FASTQs. 

    For a directory, each FASTQ is trimmed individually and concatentated into ``trim_dir/trimmed.fastq`` in lexicographic order. 
    If ``summary_tsv`` is set, one per-file row (file, reads, avg_len, avg_q) is appended for each source. 
    Returns the combined trimmed FASTQ written to ``trim_dir/trimmed.fastq``. 
    """
    # converting inputs into Path objects for file handling 
    input_path = Path(input_path)
    trim_dir = Path(trim_dir) 
    
    # Ensure output directory exists if not create it 
    trim_dir.mkdir(parents=True, exist_ok=True) 

    # Define final output path where all trimmed data will be combined 
    out_fq = trim_dir / "trimmed.fastq" 

    # Find all source FASTQ files, either from a single file path or by scanning a direcotory 
    sources = list(_iter_fastq_sources(input_path))
    if not sources:
        raise FileNotFoundError(f"No FASTQ files found in {input_path}")
    
    # If a run is interupted via user and the combined output file already exists from the previous run delete to start over 
    if out_fq.exists():
        out_fq.unlink() 

    # Initialze a list to hold summary stats for each processed file 
    summary_rows : list[tuple[str, int, float, float, float, float, str]] = [] 

    # Iterate over FASTQ file found 
    for src in sources:
        # Using a temp file because quality_trim expects a single output destination 
        tmp_out = trim_dir / f"{src.name}.trimmed.fastq" 
        # Run acutaly trimming using Trimmomatic tool 
        trimmed = quality_trim(src, tmp_out) 

        # If user requested a summary report, calculate stats for trimmed data
        reads, avg_len, avg_q, avg_mee, avg_mee_per_kb, avg_qeff, mee_qc_label = _fastq_stats(trimmed)
        if summary_tsv:  
            summary_rows.append((src.name, reads, avg_len, avg_q, avg_mee, avg_mee_per_kb, avg_qeff, mee_qc_label))

        # Append contents of the temporary trimmed file to the final combined output file 
        # 'ab' mode opens the combined file for appending in binary mode 
        # 'rb' mode opens the temp file for reading in binary mode 
        with open(out_fq, "ab") as combined, open(trimmed, "rb") as fh:
            shutil.copyfileobj(fh, combined) 

        # Delete the temporary file now that its content has been safely moved to the combined file 
        # missing_ok=True handles cases where the script might be interrupted and restarted.... 
        trimmed.unlink(missing_ok=True) 

    # After processing all source files, write the summary stats if requested 
    if summary_tsv and summary_rows: 
        summary_fp = Path(summary_tsv) 
        # Ensure the directory for the summary files exists 
        summary_fp.parent.mkdir(parents=True, exist_ok=True)

        # Check if needed to write header (really only necessary if anppending to an exsiting file) 
        write_header = not summary_fp.exists()

        # For global avg sum handling
        total_reads = sum(reads for _, reads, _, _, _, _, _, _ in summary_rows) 
        total_bases = sum(reads * avg_len for _, reads, avg_len, _, _, _, _, _ in summary_rows) 
        total_qsum = sum(avg_q * reads * avg_len for _, reads, avg_len, avg_q, _, _, _, _ in summary_rows)
        total_mee_sum = sum(avg_mee * reads for _, reads, _, _, avg_mee, _, _, _ in summary_rows)
        combined_avg_len = total_bases / total_reads if total_reads else 0 
        combined_avg_q = total_qsum / total_bases if total_bases else 0 
        combined_avg_mee = total_mee_sum / total_reads if total_reads else 0 
        combined_avg_mee_per_kb = (1000 * total_mee_sum / total_bases) if total_bases else 0 
        combined_avg_qeff = qeff_from_mee_per_kb(combined_avg_mee_per_kb) 
        combined_mee_qc_label = _mee_qc_label(combined_avg_mee_per_kb)

        # Open summary file in append mode 
        with open(summary_fp, "a") as comb:
            if write_header:
                comb.write("file\treads\tavg_len\tavg_q\tavg_mee\tavg_mee_per_kb\tavg_qeff\tmee_qc_label\n")
            # Write each row of collected stats to the TSV file 
            for name, reads, avg_len, avg_q, avg_mee, avg_mee_per_kb, avg_qeff, mee_qc_label in summary_rows:
                comb.write(f"{name}\t{reads}\t{avg_len:.1f}\t{avg_q:.2f}\t{avg_mee:.3f}\t"
                           f"{avg_mee_per_kb:.3f}\t{avg_qeff:.2f}\t{mee_qc_label}\n"
                )
            comb.write(
                f"_combined\t{total_reads}\t{combined_avg_len:.1f}\t{combined_avg_q:.2f}\t"
                f"{combined_avg_mee:.3f}\t{combined_avg_mee_per_kb:.3f}\t{combined_avg_qeff:.2f}\t" 
                f"{combined_mee_qc_label}\n"
            )
    return out_fq 

def main():
    setup_logging() # Logging set up to capture all later messages aka Python's root logger (format, level, file handler) 
    parser = argparse.ArgumentParser(description="Quality trim using Trimmomatic")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ")
    parser.add_argument("-o", "--output", required=True, help="Output FASTQ")
    args = parser.parse_args() # validates flags being used here otherwise exits with proper flags to use 

    quality_trim(args.input, args.output)

if __name__ == "__main__":
    main() 

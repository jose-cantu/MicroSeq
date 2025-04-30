# src/microseq_tests/microseq.py
from __future__ import annotations
from tqdm.auto import tqdm 
import argparse, pathlib, logging, shutil   
import microseq_tests.trimming.biopy_trim as biopy_trim
# ── pipeline wrappers (return rc int, handle logging) ──────────────
from microseq_tests.pipeline import (
        run_ab1_to_fastq, run_fastq_to_fasta,) 
from microseq_tests.utility.utils import setup_logging, load_config 
from microseq_tests.trimming.quality_trim import quality_trim 
from microseq_tests.assembly.de_novo_assembly import de_novo_assembly 
from microseq_tests.blast.run_blast import run_blast
from microseq_tests.post_blast_analysis import run as postblast_run 
from microseq_tests.utility.add_taxonomy import run_taxonomy_join
# ── low-level helpers (do the actual work, return Path/Seq list) ──
from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq 
from microseq_tests.trimming.fastq_to_fasta import fastq_folder_to_fasta
from Bio import SeqIO 



def main() -> None:
    cfg = load_config() 
    ap = argparse.ArgumentParser(
    prog="microseq", description="MicroSeq QC-trim Fastq; optional CAP3 assembly; blastn search; taxonomy join; optional BIOM export")
    # adding global flags here ...... 
    ap.add_argument("--workdir", default=cfg.get("default_workdir","data"), help="Root folder for intermediate outputs (default: ./data) note: Yaml is placed as a 2ndary place for a shared repo project which can you modify and change without using workdir flag otherwise use --workdir to point where you want to set up your individual project")
    ap.add_argument("--threads", type=int, default=1, 
                    help="CPU threads for parallel stages")
    ap.add_argument("-v", "--verbose", action="count", default=0, help="-v: info, -vv: use for debugging") 
    sp = ap.add_subparsers(dest="cmd", required=True)

    # trimming sub command 
    p_trim = sp.add_parser("trim", help="Quality trimming via Trimmomatic")
    p_trim.add_argument("-i", "--input", required=True, help="FASTQ")
    p_trim.add_argument("-o", "--output", required=False, help="ignored when --workdir is used or you can specifiy your own output if you don't want the automated version workdir gives you")
    p_trim.add_argument("--sanger", dest="sanger", action="store_true", help="Use BioPython trim for abi files Input is the AB1 folder -> convert + trim; autodetected if omitted") 
    p_trim.add_argument("--no-sanger", dest="sanger", action="store_false", help="Force FASTQ mode") 
    p_trim.set_defaults(sanger=None) # default = auto-detect in this case 
    p_trim.add_argument("--link-raw", action="store_true", help="Symlink AB1 traces into workdir instead of copying") 
    
    # ── AB1 → FASTQ -------------------------------------------------------
    p_ab1 = sp.add_parser("ab1-to-fastq",
                          help="Convert ABI chromatograms to FASTQ")
    p_ab1.add_argument("-i", "--input_dir",  required=True, metavar="DIR",
                       help="Folder containing *.ab1 files")
    p_ab1.add_argument("-o", "--output_dir", required=True, metavar="DIR",
                       help="Folder to write *.fastq files")
    p_ab1.add_argument("--overwrite", action="store_true",
                       help="Re-create FASTQ even if it exists")

    # ── FASTQ → FASTA -----------------------------------------------------
    p_fq = sp.add_parser("fastq-to-fasta",
                         help="Merge all FASTQ in a folder into one FASTA")
    p_fq.add_argument("-i", "--input_dir",    required=True, metavar="DIR",
                      help="Folder with *.fastq / *.fq")
    p_fq.add_argument("-o", "--output_fasta", required=True, metavar="FASTA",
                      help="Output FASTA file path") 
    
    # assembly 
    p_asm = sp.add_parser("assembly", help="De novo assembly via CAP3")
    p_asm.add_argument("-i", "--input", required=True)
    p_asm.add_argument("-o", "--output", required=True) 

    # blast 
    db_choices = list(cfg["databases"].keys())    # e.g. here ['gg2', 'silva', 'ncbi16s']
    p_blast = sp.add_parser("blast", help="Blast search against 16S DBs")
    p_blast.add_argument("-i", "--input", required=True)
    p_blast.add_argument("-d", "--db", choices=db_choices, required=True)
    p_blast.add_argument("-o", "--output", required=True)
    p_blast.add_argument("--identity", type=float, default=97.0, help="percent-identity threshold (default: %(default)s) you can adjust value based on needs of project")
    p_blast.add_argument("--qcov", type=float, default=80.0, help="query coverage %% (default: %(default)s) again you can adjust value based on needs of project")
    p_blast.add_argument("--max_target_seqs", type=int, default=5, help="How many DB hits to retain per query (passed to BLAST")

                         
    # taxonomy join after postblast (GG2 only) 
    p_tax = sp.add_parser("add_taxonomy", help="Append Greengenes2 taxon names to a BLAST table")
    p_tax.add_argument("-i", "--hits", required=True, help="Blast merge table (needs sseqid & qseqid)") 
    p_tax.add_argument("-t", "--taxonomy", required=True, help="Greengenes2 taxonomy.tsv") 
    p_tax.add_argument("-o", "--output", required=True, help="Output CSV/TSV == Perferable do TSV please!")


    # postblast BIOM 
    p_BIOM = sp.add_parser("postblast", help="Join BLAST + metadata -> BIOM(+CSV)")
    p_BIOM.add_argument("-b", "--blast_file", required=True, help="BLAST hits TSV produced by MicroSeq blast =)")
    p_BIOM.add_argument("-m", "--metadata", required=True, help="Metadata TSV (must have the sample_id column)")
    p_BIOM.add_argument("-o", "--output_biom", required=True, help="Output .biom path; .csv written alongside")
    p_BIOM.add_argument("--sample-col", help="Column in metadata to treat as sample_id helps MicroSeq known which column to treat as such if not sample_id itself") 
  
    # parse out arguments 
    args = ap.parse_args()

    LEVEL = {0: logging.WARNING, 1: logging.INFO}.get(args.verbose, logging.DEBUG)
    setup_logging(level=LEVEL)  # reusing helper, but expose by level 

    # createing the directory for workdir
    workdir_arg = args.workdir or cfg.get("default_workdir", "data")
    workdir = pathlib.Path(workdir_arg).expanduser().resolve() 
    for sub in ("raw", "raw_fastq", "qc", "asm", "blast", "biom", "passed_qc_fastq", "failed_qc_fastq"):
        (workdir / sub).mkdir(parents=True, exist_ok=True) 
   
   # use workdir in every brnach 
    if args.cmd == "trim":
        inp = pathlib.Path(args.input)
        if args.sanger is None:
            args.sanger = (inp.suffix.lower() == ".ab1") or any(inp.glob("*.ab1"))
        
        if args.sanger:
            # -- preparing raw_ab1 folder copy or symlink --- prefernece for symlink here for me --- 
            dst = workdir / "raw_ab1"

            # Here you clean up leftovers from a previous run so the command is idempotent assuming you use the same directory... (rare scenario) 

            if dst.exists():
                if dst.is_symlink() or dst.is_file():
                    dst.unlink()     # removes stale symlink or stray file 
                else:
                    shutil.rmtree(dst)  # removes old directory tree 

            if args.link_raw:
                # symlink handles file vs directory 
                dst.symlink_to(inp.resolve(), target_is_directory=inp.is_dir())
            else: 
                # physical copy now 
                if inp.is_dir():
                    # copy folder into raw_ab1/ 
                    shutil.copytree(inp, dst, dirs_exist_ok=True)
                else: 
                    # single AB1 file 
                    dst.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(inp, dst / inp.name)
            
            # --- convert AB1 to FASTQ --------------- 
            raw_fastq_dir = workdir / "raw_fastq"
            ab1_folder_to_fastq(dst, raw_fastq_dir)

            # --- BioPython QC trim writes to passed_qc_fastq/ ------- 
            trim_out = workdir / "passed_qc_fastq"
            biopy_trim.trim_folder(
                    input_dir=raw_fastq_dir,
                    output_dir=workdir / "qc",
                    threads=args.threads,
                    )

            # ---- FASTQ converted to FASTA --------------------
            fasta = fastq_folder_to_fasta(
                    trim_out,       # here FASTQs are kept 
                    workdir / "qc" / "trimmed.fasta" 
                    )
        else:
            # --- Trimmomatic QC for regular FASTQ ---------------
            out_fq = workdir / "qc" / "trimmed.fastq" 
            quality_trim(args.input, out_fq, threads=args.threads)

            fasta = fastq_folder_to_fasta(
                    workdir / "qc", # contains the trimmed.fastq 
                    workdir / "qc" / "trimmed.fasta" 
                    )
        print("FASTA ready:", fasta)
    


    elif args.cmd == "ab1-to-fastq":
        rc = run_ab1_to_fastq(
                args.input_dir,
                args.output_dir,
                overwrite=args.overwrite,
                )
        print(" AB1->FastQ exit-code:", rc)


    elif args.cmd == "fastq-to-fasta":
        rc = run_fastq_to_fasta(
                args.input_dir,
                args.output_fasta,
                )
        print("FASTQ -> FASTA exit-code:", rc) 


    elif args.cmd == "assembly":
        de_novo_assembly(workdir / "qc" / "trimmed.fasta",
                         workdir / "asm",
                         threads=args.threads,
                         )

    elif args.cmd == "blast":
        total = sum(1 for _ in SeqIO.parse(args.input, "fasta"))
        if total == 0:
            total = sum(1 for _ in SeqIO.parse(args.input, "fastq"))
        with tqdm(total=total, unit="seq") as bar: 
            run_blast(
                pathlib.Path(args.input),  # temporary patch for blasting and post biom only will need to think around 
                args.db,
                pathlib.Path(args.output),
                pct_id=args.identity,
                qcov=args.qcov,
                max_target_seqs = args.max_target_seqs, 
                threads=args.threads,
                on_progress=lambda p: bar.update(p - bar.n), 
            )

    elif args.cmd == "postblast":
        out_biom = workdir / "biom" / args.output_biom 
        out_biom.parent.mkdir(exist_ok=True, parents=True)

        postblast_run(
                pathlib.Path(args.blast_file),
                pathlib.Path(args.metadata),
                out_biom,
                write_csv=True,
                sample_col=args.sample_col,
                )
        print(f" ✓ BIOM : {out_biom}")
        print(f" ✓ CSV  : {out_biom.with_suffix('.csv')}") 

    elif args.cmd == "add_taxonomy": 
        run_taxonomy_join(
                pathlib.Path(args.hits).resolve(),
                pathlib.Path(args.taxonomy).expanduser().resolve(),
                pathlib.Path(args.output).resolve(),
                )
        print(f" ✓ CSV+tax : {args.output}")


if __name__ == "__main__":
    main() 




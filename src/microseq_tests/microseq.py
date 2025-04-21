# src/microseq_tests/microseq.py
from __future__ import annotations 
import argparse, sys, pathlib 
import microseq_tests.trimming.biopy_trim as biopy_trim 
from microseq_tests.utility.utils import setup_logging, load_config 
from microseq_tests.trimming.quality_trim import quality_trim 
from microseq_tests.assembly.de_novo_assembly import de_novo_assembly 
from microseq_tests.blast.run_blast import run_blast
from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq 
from microseq_tests.trimming.fastq_to_fasta import fastq_folder_to_fasta 
# from microseq_tests.blast.post_blast_analysis import build_biom

def main() -> None:
    setup_logging() # sets up logging from here... 
    ap = argparse.ArgumentParser(description="MicroSeq Master Pipeline")
    # adding global flag 
    ap.add_argument("--workdir", default="data", help="Root folder for intermediate outputs (default: ./data)") 
    sp = ap.add_subparsers(dest="cmd", required=True)

    # trimming sub command 
    p_trim = sp.add_parser("trim", help="Quality trimming via Trimmomatic")
    p_trim.add_argument("-i", "--input", required=True, help="FASTQ")
    p_trim.add_argument("-o", "--output", required=False, help="(ignored when --workdir is used or you can specifiy your own output if you don't want the automated version workdir gives you")
    p_trim.add_argument("--sanger", action="store_true", help="Use BioPython trim for abi files Input is the AB1 folder -> convert + trim") 

    # assembly 
    p_asm = sp.add_parser("assembly", help="De novo assembly via CAP3")
    p_asm.add_argument("-i", "--input", required=True)
    p_asm.add_argument("-o", "--output", required=True) 

    # blast 
    cfg = load_config()
    db_choices = list(cfg["databases"].keys())    # e.g. here ['gg2', 'silva', 'ncbi16s']
    p_blast = sp.add_parser("blast", help="Blast search against 16S DBs")
    p_blast.add_argument("-i", "--input", required=True)
    p_blast.add_argument("-d", "--db", choices=db_choices, required=True)
    p_blast.add_argument("-o", "--output", required=True)

    # postblast BIOM 
    p_BIOM = sp.add_parser("postblast", help="Parse BLAST + build BIOM (stud)")
    p_BIOM.add_argument("-b", "--blast_file", required=True)
    p_BIOM.add_argument("-d", "--db", choices=db_choices, required=True)
    p_BIOM.add_argument("-o", "--output_biom", required=True)
    
    # parse out arguments 
    args = ap.parse_args() 

    # createing the directory for workdir 
    workdir = pathlib.Path(args.workdir).resolve() 
    for sub in ("raw", "raw_fastq", "qc", "asm", "blast", "biom", "passed_qc_fastq", "failed_qc_fastq"):
        (workdir / sub).mkdir(parents=True, exist_ok=True) 
   
   # use workdir in every brnach 
    if args.cmd == "trim":
        if args.sanger:
            raw_fastq_dir = workdir / "raw_fastq"
            ab1_folder_to_fastq(args.input, raw_fastq_dir)
            biopy_trim.trim_folder(
                input_dir=raw_fastq_dir,
                output_dir=workdir / "qc",
            )
        else:
            out_fq = workdir / "qc" / "trimmed.fastq"
            quality_trim(args.input, out_fq)

        # In both branches now have passed_qc_fastq/*.fasta 
        fasta = fastq_folder_to_fasta(
            workdir / "passed_qc_fastq",
            workdir / "qc" / "trimmed.fasta")
        print("FASTA ready to go!", fasta)

    elif args.cmd == "assembly":
        de_novo_assembly(workdir / "qc" / "trimmed.fasta",
                         workdir / "asm")

    elif args.cmd == "blast":
        contigs = workdir / "asm" / "trimmed.fasta.cap.contigs"
        run_blast(contigs, args.db,
                  workdir / "blast" / f"hits_{args.db}.tsv")

    elif args.cmd == "postblast":
        print("postblast not implemented yet")
        sys.exit(1)

if __name__ == "__main__":
    main() 







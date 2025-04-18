# src/microseq_tests/microseq.py 
import argparse, sys, pathlib 
from microseq_tests.utility.utils import setup_logging, load_config 
from microseq_tests.trimming.quality_trim import quality_trim 
from microseq_tests.assembly.de_novo_assembly import de_novo_assembly 
from microseq_tests.blast.run_blast import run_blast 

def main() -> None:
    setup_logging() # sets up logging from here... 
    ap = argparse.ArgumentParser(description="MicroSeq Master Pipeline")
    sp = ap.add_subparsers(dest="cmd", required=True)

    # trimming 
    p_trim = sp.add_parser("trim", help="Quality trimming via Trimmomatic")
    p_trim.add_argument("-i", "--input", required=True, help="FASTQ")
    p_trim.add_argument("-o", "--output", required=True, help="trimmed FASTA") 

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

    args = ap.parse_args() 

    if args.cmd == "trim":
        quality_trim(args.input, args.output)

    elif args.cmd == "assembly":
        de_novo_assembly(args.input, args.output)

    elif args.cmd == "blast":
        run_blast(args.input, args.db, args.output)

    elif args.cmd == "postblast":
        print("TODO: post-BLAST script parser not implemented yet.")
        sys.exit(1)

if __name__ == "__main__":
    main() 







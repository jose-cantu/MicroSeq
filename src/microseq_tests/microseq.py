# src/microseq_tests/microseq.py 
import argparse, sys, pathlib, biopy_trim  
from microseq_tests.utility.utils import setup_logging, load_config 
from microseq_tests.trimming.quality_trim import quality_trim 
from microseq_tests.assembly.de_novo_assembly import de_novo_assembly 
from microseq_tests.blast.run_blast import run_blast 
from microseq_tests.blast.post_blast_analysis import build_biom

def main() -> None:
    setup_logging() # sets up logging from here... 
    ap = argparse.ArgumentParser(description="MicroSeq Master Pipeline")
    sp = ap.add_subparsers(dest="cmd", required=True)

    # trimming 
    p_trim = sp.add_parser("trim", help="Quality trimming via Trimmomatic")
    p_trim.add_argument("-i", "--input", required=True, help="FASTQ")
    p_trim.add_argument("-o", "--output", required=True, help="trimmed FASTA")
    p_trim.add_argument("--sanger", action="store_true", help="Use BioPython trim for abi files") 

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
        out = (quality_trim(args.input, args.output)
               if not args.sanger
               else biopy_trim.trim_folder(args.input, pathlib.Path(args.output).parent)) 
              

    elif args.cmd == "assembly":
        de_novo_assembly(args.input, args.output)

    elif args.cmd == "blast":
        run_blast(args.input, args.db, args.output)

    elif args.cmd == "postblast":
        tax_tsv = f"config/{args.db}_tax.tsv"
        build_biom(args.blast_file, tax_tsv, args.output_biom, sample_id="Demo", ident=97.0, top_n=5) 

if __name__ == "__main__":
    main() 







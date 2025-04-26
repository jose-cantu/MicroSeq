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
from microseq_tests.post_blast_analysis import run as postblast_run 
from microseq_tests.utility.add_taxonomy import run_taxonomy_join 

def main() -> None:
    setup_logging() # sets up logging from here... 
    cfg = load_config() 
    ap = argparse.ArgumentParser(
    prog="microseq", description="MicroSeq QC-trim Fastq; optional CAP3 assembly; blastn search; taxonomy join; optional BIOM export")
    # adding global flag 
    ap.add_argument("--workdir", default=cfg.get("default_workdir","data"), help="Root folder for intermediate outputs (default: ./data) note: Yaml is placed as a 2ndary place for a shared repo project which can you modify and change without using workdir flag otherwise use --workdir to point where you want to set up your individual project") 
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
    db_choices = list(cfg["databases"].keys())    # e.g. here ['gg2', 'silva', 'ncbi16s']
    p_blast = sp.add_parser("blast", help="Blast search against 16S DBs")
    p_blast.add_argument("-i", "--input", required=True)
    p_blast.add_argument("-d", "--db", choices=db_choices, required=True)
    p_blast.add_argument("-o", "--output", required=True)
    p_blast.add_argument("--identity", type=float, default=97.0, help="percent-identity threshold (default: %(default)s) you can adjust value based on needs of project")
    p_blast.add_argument("--qcov", type=float, default=80.0, help="query coverage %% (default: %(default)s) again you can adjust value based on needs of project")

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
        run_blast(
            pathlib.Path(args.input),  # temporary patch for blasting and post biom only will need to think around 
            args.db,
            pathlib.Path(args.output),
            pct_id=args.identity,
            qcov=args.qcov
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




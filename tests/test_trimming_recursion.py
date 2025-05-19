import pathlib
import shutil
from Bio import SeqIO
from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq
from microseq_tests.trimming.fastq_to_fasta import fastq_folder_to_fasta


def test_ab1_to_fastq_recursive(tmp_path: pathlib.Path) -> None:
    fixtures = pathlib.Path(__file__).parent / "fixtures" / "39764260.ab1"
    inp = tmp_path / "input"
    nested = inp / "sub"
    nested.mkdir(parents=True)
    shutil.copy2(fixtures, inp / "top.ab1")
    shutil.copy2(fixtures, nested / "inner.ab1")

    out_dir = tmp_path / "fastq"
    fastqs = ab1_folder_to_fastq(inp, out_dir)
    assert len(fastqs) == 2
    for fq in fastqs:
        assert fq.exists()


def test_fastq_to_fasta_recursive(tmp_path: pathlib.Path) -> None:
    inp = tmp_path / "fq"
    sub = inp / "deep"
    sub.mkdir(parents=True)
    (inp / "r1.fastq").write_text("@a\nACGT\n+\n!!!!\n")
    (sub / "r2.fastq").write_text("@b\nTGCA\n+\n!!!!\n")

    out = tmp_path / "out.fasta"
    fastq_folder_to_fasta(inp, out)
    records = list(SeqIO.parse(out, "fasta"))
    assert len(records) == 2

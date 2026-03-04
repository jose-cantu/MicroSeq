from pathlib import Path

import pandas as pd

from microseq_tests.utility.add_taxonomy import run_taxonomy_join 


def test_run_taxonomy_join_derives_sample_id_from_qseqid(tmp_path: Path) -> None:
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "qseqid\tsseqid\tpident\tevalue\tbitscore\n"
        "S1|contig|cap3_c1\tGG2_0001\t99.0\t1e-10\t200\n",
        encoding="utf-8",
    )
    tax = tmp_path / "taxonomy.tsv"
    tax.write_text("sseqid\ttaxonomy\nGG2_0001\tk__Bacteria\n", encoding="utf-8")

    out = tmp_path / "hits_tax.tsv"
    run_taxonomy_join(hits, tax, out)

    df = pd.read_csv(out, sep="\t")
    assert df.loc[0, "qseqid"] == "S1|contig|cap3_c1"
    assert df.loc[0, "sample_id"] == "S1"


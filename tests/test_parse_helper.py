
def test_parse_roundtrip(tmp_path):
    p = tmp_path / "dummy.tsv"
    p.write_text("a\tb\t99\t100\t99\t100\t1e-10\t200\ttitle\n")
    from microseq_tests.aligners._parse import parse_blast_tab, COLS
    df = parse_blast_tab(p) 
    assert list(df.columns) == list(COLS) 

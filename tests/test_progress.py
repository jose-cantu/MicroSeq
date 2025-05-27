# microseq_tests/tests/test_progress.py 

"""
The goal here is making sure each TSV or fastq ticks exactly once when the helper is wrapped instage_bar 
The approach used is a temporary monkeypatch == py.test.monkeypatch 
 - The real progress.stage_bar with a tiny fkae that records and similates how often update() is called. 
 - The public function merged_hits or trim_folder on a small fixture dataset and asset the counter equals the number of files.

""" 
from __future__ import annotations
from pathlib import Path
import io, textwrap, pytest
pytest.importorskip("pandas")
pytest.importorskip("Bio")

# functions that are going to be tested 
from microseq_tests.microseq import merge_hits 
from microseq_tests.trimming.biopy_trim import trim_folder 
import microseq_tests.utility.progress as pg 

# -------------- tiny in-memory fake bar --------------- 

class DummyBar:
    def __init__(self, total):
        self.total = total 
        self.n = 0 
    def update(self, inc=1):
        self.n += inc 
    def close(self):
        pass 
    def __enter__(self):
        return self 
    def __exit__(self, *exc): 
        self.close() 

@pytest.fixture() 
def patch_stage_bar(monkeypatch):
    """Replace progress.stage_bar with a counter-collecting dummy.""" 
    counters = {} 
    def fake_stage_bar(total, *a, **kw):
        bar = DummyBar(total)
        counters['bar'] = bar 
        yield bar 
    monkeypatch.setattr(pg, "stage_bar", fake_stage_bar)
    yield counters # give test access to bar after the call 

def test_merge_hits_progress(tmp_path: Path, patch_stage_bar):
    # create three toy TSVs
    for i in range(3):
        (tmp_path / f"f{i}.tsv").write_text(
            textwrap.dedent("""\
            # header
            s1\ts2\t99\t100\t90\t100\t1e-10\t200\ttaxon
            """)
        )
    out = tmp_path / "merged.tsv"
    
    # call the public helper that merge-hits CLI also uses
    merge_hits([str(p) for p in tmp_path.glob("*.tsv")], out)
    
    bar = patch_stage_bar['bar']
    assert bar.total == 3          # three files
    assert bar.n == 3              # ticked exactly once per file
    # merged file really exists
    assert out.exists() and out.stat().st_size > 0
    
def test_trim_folder_progress(tmp_path: Path, patch_stage_bar):
    in_dir  = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir()
    
    # two minimal FASTQ files (three lines per read)
    fastq = textwrap.dedent("""\
        @r1
        AACGT
        +
        !!!!!        """)
    for name in ("a.fastq", "b.fastq"):
        (in_dir / name).write_text(fastq)
        
    # run the BioPython-based trimmer (very small win/q so it keeps reads)
    trim_folder(in_dir, out_dir, window_size=1, per_base_q=0)
    
    bar = patch_stage_bar['bar']
    assert bar.total == 2          # two FASTQs
    assert bar.n == 2              # updated twice
    # trimmed files were written
    assert any(out_dir.glob("*_avg_qual.txt")) 




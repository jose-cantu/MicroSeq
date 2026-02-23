from pathlib import Path


def test_biopy_trim_defines_trace_lookup_keys() -> None:
    src = Path("src/microseq_tests/trimming/biopy_trim.py").read_text(encoding="utf-8")
    assert "def _trace_lookup_keys(" in src
    assert 'if stem.endswith("_trimmed")' in src

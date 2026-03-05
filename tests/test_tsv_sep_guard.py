from pathlib import Path
import re


def test_pipeline_uses_explicit_tab_delimiters() -> None:
    src = Path("src/microseq_tests/pipeline.py").read_text(encoding="utf-8")
    # forbid literal-tab sep declarations like sep="<TAB>"
    assert not re.search(r'sep="\t"'.replace('\\t', '\t'), src)

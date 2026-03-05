import pathlib
import pytest
pytest.importorskip("Bio")
from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq, ab1_rel_key, build_ab1_output_key_map

def test_ab1_to_fastq(tmp_path: pathlib.Path):
    fixtures = pathlib.Path(__file__).parent / "fixtures" 
    fastqs = ab1_folder_to_fastq(fixtures, tmp_path)

    # should procedure one FASTQ per AB1 file 
    assert len(fastqs) == len(list(fixtures.glob("*ab1"))) 
    for fq in fastqs:
        assert fq.exists()
        assert fq.stat().st_size > 100 # tiny but non-empty 


def test_ab1_rel_key_preserves_staging_relative_path_for_symlinks(tmp_path: pathlib.Path):
    staging = tmp_path / "staging"
    staging.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    real_ab1 = outside / "a.ab1"
    real_ab1.write_bytes(b"ABIF")

    sub = staging / "sub"
    sub.mkdir()
    linked = sub / "a.ab1"
    try:
        linked.symlink_to(real_ab1)
    except (OSError, NotImplementedError) as exc:
        pytest.skip(f"symlink not supported in this environment: {exc}")

    assert ab1_rel_key(linked, staging) == "sub__a"


def test_ab1_output_key_map_prefers_stem_unless_collision(tmp_path: pathlib.Path):
    root = tmp_path / "ab1"
    (root / "d1").mkdir(parents=True)
    (root / "d2").mkdir(parents=True)
    (root / "u").mkdir(parents=True)

    (root / "d1" / "sample.ab1").write_bytes(b"ABIF")
    (root / "d2" / "sample.ab1").write_bytes(b"ABIF")
    (root / "u" / "unique.ab1").write_bytes(b"ABIF")

    keys = build_ab1_output_key_map(root)
    assert keys[root / "u" / "unique.ab1"] == "unique"
    assert keys[root / "d1" / "sample.ab1"] != "sample"
    assert keys[root / "d2" / "sample.ab1"] != "sample"
    assert keys[root / "d1" / "sample.ab1"] != keys[root / "d2" / "sample.ab1"]

from types import SimpleNamespace
import importlib.util, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "hw", ROOT / "src" / "microseq_tests" / "hardware.py"
)
hw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hw)

def test_recommend_threads(monkeypatch):
    vm = SimpleNamespace(total=32 * 1024 ** 3)
    monkeypatch.setattr(hw, "psutil", SimpleNamespace(virtual_memory=lambda: vm))
    monkeypatch.setattr(hw.os, "cpu_count", lambda: 12)
    assert hw.recommend_threads() == 12

def test_cap_at_16(monkeypatch):
    vm = SimpleNamespace(total=64 * 1024 ** 3)
    monkeypatch.setattr(hw, "psutil", SimpleNamespace(virtual_memory=lambda: vm))
    monkeypatch.setattr(hw.os, "cpu_count", lambda: 32)
    assert hw.recommend_threads() == 16

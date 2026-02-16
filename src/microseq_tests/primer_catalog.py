from __future__ import annotations

from pathlib import Path
import re
from importlib import resources
import yaml

from microseq_tests.utility.utils import ROOT

CATALOG_PATH = ROOT / "config" / "primer_catalog.yaml"



def _load_catalog_payload(catalog_path: Path) -> dict[str, object]:
    if catalog_path.exists():
        loaded = yaml.safe_load(catalog_path.read_text(encoding="utf-8")) or {}
        return loaded if isinstance(loaded, dict) else {}

    try:
        packaged = resources.files("microseq_tests.data").joinpath("primer_catalog.yaml")
        if packaged.is_file():
            loaded = yaml.safe_load(packaged.read_text(encoding="utf-8")) or {}
            return loaded if isinstance(loaded, dict) else {}
    except Exception:
        pass
    return {}

_DEFAULT_CATALOG: dict[str, object] = {
    "pairing_sets": {
        "16S (27F/1492R)": {"forward": ["27F"], "reverse": ["1492R"]},
        "16S (8F/1492R)": {"forward": ["8F"], "reverse": ["1492R"]},
        "16S V4 (515F/806R)": {"forward": ["515F"], "reverse": ["806R"]},
        "Custom": {"forward": [], "reverse": []},
    },
    "trim_presets": {
        "16S_27F_1492R": {
            "forward_primers": ["AGAGTTTGATCMTGGCTCAG"],
            "reverse_primers": ["TACGGYTACCTTGTTACGACTT"],
            "max_mismatch": 2,
            "max_search": 80,
            "max_primer_offset": 10,
            "iupac_mode": True,
        }
    },
}


def load_primer_catalog(path: Path | None = None) -> dict[str, object]:
    catalog = dict(_DEFAULT_CATALOG)
    catalog_path = path or CATALOG_PATH
    loaded = _load_catalog_payload(catalog_path)
    for key in ("pairing_sets", "trim_presets"):
        if isinstance(loaded.get(key), dict):
            merged = dict(catalog.get(key, {}))
            merged.update(loaded[key])
            catalog[key] = merged
    return catalog


def pairing_label_sets() -> dict[str, tuple[list[str], list[str]]]:
    pairing_sets = load_primer_catalog().get("pairing_sets", {})
    out: dict[str, tuple[list[str], list[str]]] = {}
    for label, row in pairing_sets.items():
        if not isinstance(label, str) or not isinstance(row, dict):
            continue
        out[label] = (list(row.get("forward", []) or []), list(row.get("reverse", []) or []))
    return out


def trim_presets() -> dict[str, dict[str, object]]:
    presets = load_primer_catalog().get("trim_presets", {})
    out: dict[str, dict[str, object]] = {}
    for key, row in presets.items():
        if isinstance(key, str) and isinstance(row, dict):
            out[key] = dict(row)
    return out


def parse_primer_sequences(text: str) -> list[str]:
    seqs: list[str] = []
    for raw in text.splitlines():
        line = raw.strip().upper()
        if not line:
            continue
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            continue
        line_compact = re.sub(r"\s+", "", line)
        exact = re.fullmatch(r"[ACGTURYSWKMBDHVN]{8,}", line_compact)
        if exact:
            seqs.append(line_compact.replace("U", "T"))
            continue
        match = re.search(r"[ACGTURYSWKMBDHVN]{8,}", line)
        if match:
            seqs.append(match.group(0).replace("U", "T"))
            continue
        raise ValueError(f"Invalid primer sequence characters in: {raw.strip()}")
    return seqs


def build_primer_cfg_override(
    *,
    mode: str,
    stage: str,
    preset: str,
    forward_raw: str,
    reverse_raw: str,
    for_preview: bool = False,
) -> dict[str, object]:
    """Build primer config override from GUI/CLI style raw inputs."""
    mode_norm = mode.strip().lower()
    if mode_norm not in {"off", "detect", "clip"}:
        mode_norm = "off"

    override: dict[str, object] = {
        "mode": mode_norm,
        "stage": stage,
    }

    preset_norm = str(preset or "").strip()
    if preset_norm:
        override["preset"] = preset_norm

    parse_custom_sequences = mode_norm != "off" or for_preview
    fwd_text = (forward_raw or "").strip()
    rev_text = (reverse_raw or "").strip()

    if parse_custom_sequences:
        forward_primers = parse_primer_sequences(fwd_text) if fwd_text else []
        reverse_primers = parse_primer_sequences(rev_text) if rev_text else []

        # Custom primers take precedence for this run; omit preset override in that case.
        if (forward_primers or reverse_primers) and "preset" in override:
            override.pop("preset", None)

        # Only emit sequence-list overrides when there is at least one parsed sequence.
        # This keeps header/comment-only text from clobbering config-provided lists.
        if forward_primers:
            override["forward_primers"] = forward_primers
        if reverse_primers:
            override["reverse_primers"] = reverse_primers

    return override

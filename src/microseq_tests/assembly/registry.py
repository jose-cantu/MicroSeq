from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


AssemblerKind = Literal["merge_two_reads", "cap3"]


@dataclass(frozen=True)
class AssemblerSpec:
    id: str
    display_name: str
    kind: AssemblerKind
    overlap_engine: str | None = None
    cap3_profile: str | None = None


ASSEMBLER_SPECS: tuple[AssemblerSpec, ...] = (
    AssemblerSpec(
        id="merge_two_reads:biopython",
        display_name="Merge two reads (Biopython overlap)",
        kind="merge_two_reads",
        overlap_engine="biopython",
    ),
    AssemblerSpec(
        id="merge_two_reads:ungapped",
        display_name="Merge two reads (Ungapped overlap)",
        kind="merge_two_reads",
        overlap_engine="ungapped",
    ),
    AssemblerSpec(
        id="merge_two_reads:edlib",
        display_name="Merge two reads (Edlib overlap)",
        kind="merge_two_reads",
        overlap_engine="edlib",
    ),
    AssemblerSpec(
        id="cap3:strict",
        display_name="CAP3 (Strict profile)",
        kind="cap3",
        cap3_profile="strict",
    ),
    AssemblerSpec(
        id="cap3:diagnostic",
        display_name="CAP3 (Diagnostic profile)",
        kind="cap3",
        cap3_profile="diagnostic",
    ),
    AssemblerSpec(
        id="cap3:relaxed",
        display_name="CAP3 (Relaxed profile)",
        kind="cap3",
        cap3_profile="relaxed",
    ),
)


def list_assemblers() -> tuple[AssemblerSpec, ...]:
    return ASSEMBLER_SPECS


def get_assembler_spec(assembler_id: str) -> AssemblerSpec:
    for spec in ASSEMBLER_SPECS:
        if spec.id == assembler_id:
            return spec
    known = ", ".join(spec.id for spec in ASSEMBLER_SPECS)
    raise ValueError(f"Unknown assembler_id '{assembler_id}'. Known ids: {known}")

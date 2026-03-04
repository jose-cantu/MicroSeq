# src/microseq_tests/assembly/paired_assembly.py 

"""
Utilities for running CAP3 assembly on forward/reverse pairs.
""" 
from __future__ import annotations # Postpones evaluation of type annotations (PEP 563) so they are no longer evaluated at function definition time - treated as string instead first 
import logging # print warning messages  
from os import PathLike
import subprocess 
from pathlib import Path
from datetime import datetime, timezone 
from typing import Iterable, Sequence, Literal
import re

from Bio import SeqIO 

from microseq_tests.utility.utils import load_config 

from .pairing import DupPolicy, group_pairs 
from .two_read_merge import merge_two_reads, MergeInputError
from .overlap_backends import resolve_overlap_engine

L = logging.getLogger(__name__) 

def _safe_log_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "").strip()) # normalize unsafe chars 
    return token.strip("_") or "cap3" # never empty so fallback to prevent blank stem 

def write_cap3_process_logs(
    logs_root: Path, # directory root for CAP3 logs 
    *, 
    sample_id: str, # full sample key 
    command: Sequence[str], # argv list used to run CAP3 
    assembler_label: str, # label used in filename and header (cap3_default) 
    stdout_text: str, # capturing stdout 
    stderr_text: str, # capturing stderr  
    ) -> tuple[Path, Path]: # returns (stdout_path, stderr_path) 
    """Persist CAP3 process stdout/stderr with context headers."""
    logs_root.mkdir(parents=True, exist_ok=True) 
    sample_token = _safe_log_token(sample_id) # safe token for filenames 
    assembler_token = _safe_log_token(assembler_label) # safe toekn for filenames 
    stem = f"{sample_token}__{assembler_token}" # deterministic naming scheme 
    stdout_path = logs_root / f"{stem}.stdout.log" # stable stdout filename 
    stderr_path = logs_root / f"{stem}.stderr.log" # stable stderr filename 
    ts = datetime.now(timezone.utc).isoformat() # timestamp uses UTC 
    header = ( # header captures minimum debugging context for me 
        f"# timestamp_utc: {ts}\n"
        f"# sample_id: {sample_id}\n"
        f"# assembler: {assembler_label}\n"
        f"# command: {' '.join(str(part) for part in command)}\n\n"
    )
    stdout_path.write_text(header + (stdout_text or ""), encoding="utf-8")  # side effect: writes stdout artifact
    stderr_path.write_text(header + (stderr_text or ""), encoding="utf-8")  # side effect: writes stderr artifact
    return stdout_path, stderr_path  # invariant: paths are concrete filesystem locations

def _iter_paths(value: Path | list[Path]) -> Iterable[Path]:
    """Normalize a stored path entry such as a list into an iterator."""
    if isinstance(value, list):
        for item in value: 
            yield Path(item)
    else:
            yield Path(value)  

def _as_path_list(value: Path | list[Path]) -> list[Path]:
    """Return a list of `path` objects perserving insertion order."""
    return [Path(p) for p in _iter_paths(value)]

def _build_keep_separate_pairs(forward: list[Path], reverse: list[Path]) -> list[tuple[Path, Path]]:
    """Produce forward/reverse pairings for the keep-separate policy."""
    if not forward or not reverse: 
        raise ValueError("Forward and reverse inputs must not be empty when paring reads.") 

    f_count = len(forward)
    r_count = len(reverse)

    # If both orientations contains multiple primers, require the counts to align
    # so that we can preserve primer ordering w/o inventing new orientations 
    if f_count >1 and r_count > 1 and f_count != r_count:
        raise ValueError("Mismatched duplicate counts detected while using 'keep-separate' policy:" 
        f"{f_count} forward vs {r_count} reverse reads." 
        ) 

    pairs: list[tuple[Path, Path]] = []
    if f_count == r_count:
        pairs = list(zip(forward, reverse))
    elif f_count > r_count:
        if r_count != 1: 
            raise ValueError("Unable to align forward duplicates with revcerse reads under 'keep-separate' policy.")
        pairs = [(f_path, reverse[0]) for f_path in forward]
    else: # r_count > f_count 
        if f_count != 1:
            raise ValueError(
                "Unable to align reverse duplicates with forward reads under 'keep-separate' policy.")
        pairs = [(forward[0], r_path) for r_path in reverse] 

    return pairs 





def _parse_ace_contig_members(ace_path: Path) -> list[set[str]]:
    """Return per-contig read membership sets parsed from a CAP3 ACE file."""
    memberships: list[set[str]] = []
    current: set[str] | None = None

    if not ace_path.exists():
        return memberships

    with ace_path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("CO "):
                if current is not None:
                    memberships.append(current)
                current = set()
                continue
            if current is None:
                continue
            if line.startswith("AF "):
                parts = line.split()
                if len(parts) >= 2:
                    current.add(parts[1])

    if current is not None:
        memberships.append(current)

    return memberships


def _validate_cap3_contig_support(
    contig_path: Path,
    fwd_path: Path,
    rev_path: Path,
) -> Literal["verified", "rejected", "unknown"]:
    """Validate CAP3 pair support from ACE membership for both source read IDs."""
    fwd = next(SeqIO.parse(fwd_path, "fasta"), None)
    rev = next(SeqIO.parse(rev_path, "fasta"), None)
    if fwd is None or rev is None:
        return "unknown"

    ace_path = contig_path.with_suffix(".ace")
    memberships = _parse_ace_contig_members(ace_path)
    if not memberships:
        return "unknown"

    for member_ids in memberships:
        if fwd.id in member_ids and rev.id in member_ids:
            return "verified"
    return "rejected"

def _write_combined_fasta(sources: Iterable[Path], destination: Path, *, use_qual: bool = True) -> None: 
    """ 
    Here I will combine multiple FASTA files into 'destination' and appending matching QUALS.
    Each input FASTA is written in-order. When a sibling ``.qual`` file is
    available (``<fasta>.qual``), its records are appended in the same order and
    validated to ensure IDs and sequence/quality lengths match. If a source FASTA
    lacks a ``.qual`` file, log a warning and skip its QUAL contribution.

    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    qual_out_path = destination.with_name(f"{destination.name}.qual")
    qual_out_handle = None 

    with destination.open("w", encoding="utf-8") as fasta_out:
        for src in sources:
            records = list(SeqIO.parse(src, "fasta"))
            if not records:
                continue
            SeqIO.write(records, fasta_out, "fasta")

            if not use_qual:
                continue

            qual_src = Path(f"{src}.qual")
            if not qual_src.exists():
                L.warning("Missing QUAL file for %s; skipping quality output for this source.", src)
                continue

            qual_records = {rec.id: rec for rec in SeqIO.parse(qual_src, "qual")}
            if qual_out_handle is None:
                qual_out_handle = qual_out_path.open("w", encoding="utf-8")

            for record in records:
                qual_rec = qual_records.get(record.id)
                if qual_rec is None:
                    raise ValueError(
                        f"QUAL record missing for {record.id} in {qual_src}"
                    )
                quals = qual_rec.letter_annotations.get("phred_quality", [])
                if len(quals) != len(record.seq):
                    raise ValueError(
                        f"Quality length mismatch for {record.id} in {qual_src}: "
                        f"{len(quals)} != {len(record.seq)}"
                    )
                qual_rec.description = ""
                SeqIO.write(qual_rec, qual_out_handle, "qual")

    if qual_out_handle is not None:
        qual_out_handle.close()

def _write_pairing_report(pairs: dict, meta: dict, destination: Path) -> None:
    """Persist a TSV showing how pairs were matched and which detector fired."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as fh:
        fh.write("sid\tF_path\tR_path\tdetector\n")

        for sid in sorted(pairs):
            entries = pairs[sid]
            meta_entry = meta.get(sid, {})

            f_paths = _as_path_list(entries["F"])
            r_paths = _as_path_list(entries["R"])
            
            # function for formatting wells in well plate returns a string or none  
            def _format_well(value: str | list[str] | None) -> str | None:
                if value is None:
                    return None
                if isinstance(value, list):
                    return ";".join(value)
                return value 

            def _format_det(value: str | list[str] | None, orient: str) -> str | None:
                if value is None:
                    return None
                if isinstance(value, list):
                   # If it is a list, format a string starting with the orientation, followed by a colon
                   # Then the list of elements joined via semicolons 
                   det_txt = f"{orient}:" + ";".join(value)
                   # If the value is not None and not a list (a string) 
                else:
                    # Format a string starting with the orientation, followed by a colon, and a string value. 
                    det_txt = f"{orient}:{value}"
                # Retrieve a related metadata entry using the 'orient' parameter using `_format_well` to ensure consistent formatting
                well_txt = _format_well(meta_entry.get(f"well_{orient}"))
                # Check if well_txt has a value 
                # If so, append well information to det_txt, otherwise return det_txt
                return f"{det_txt} (well {well_txt})" if well_txt else det_txt

            detectors = " | ".join(
                filter(
                    None,
                    (
                        _format_det(meta_entry.get("F"), "F"),
                        _format_det(meta_entry.get("R"), "R"),
                    ),
                )
            )

            fh.write(
                "\t".join(
                    [
                        sid,
                        ";".join(str(p) for p in f_paths),
                        ";".join(str(p) for p in r_paths),
                        detectors,
                    ]
                )
                + "\n"
            )

def assemble_pairs(input_dir: PathLike, output_dir: PathLike, *, dup_policy: DupPolicy = DupPolicy.ERROR, cap3_options: Sequence[str] | None = None, fwd_pattern: str | None = None, rev_pattern: str | None = None, pairing_report: PathLike | None = None, enforce_same_well: bool = False, well_pattern: str | re.Pattern[str] | None = None, 
                   use_qual: bool = True 
                   ) -> list[Path]:
    """
    Run CAP3 assemblies for each forward and reverse pair discovered in ``input_dir``. 
    
    Parameters:
    ----------
    input_dir:
        Directory containiing per-sample FASTA files. Files are paired via 
        via :func:`microseq_tests.assembly.pairing.groups_pairs`.
    output_dir:
        Destination directory for per-sample CAP3 runs. Each sample receives a subdirectory containing 
        the concatenated FASTA and CAP3 outputs.
    dup_policy:
        Policy for handling duplicate forward/reverse files when pairing. 
        See :class: `microseq_tests.assembly.pairing.DupPolicy`.
    cap3_options:
        Optional additional command-line arguments appended to CAP3 call. 

    Return:
    ------
    
    list[Pathlib.Path]

    A list of paths to the generated ``*.cap.contigs`` files, ordered by
        sample identifier. 
        list[Pathlib.Path]

    """
    cfg = load_config() 
    overlap_cfg = cfg.get("overlap_eval", {})
    merge_cfg = cfg.get("merge_two_reads", {})
    merge_min_overlap = int(overlap_cfg.get("min_overlap", merge_cfg.get("min_overlap", 100)))
    merge_min_identity = float(overlap_cfg.get("min_identity", merge_cfg.get("min_identity", 0.8)))
    merge_overlap_engine = str(merge_cfg.get("overlap_engine", "ungapped")).strip().lower()
    merge_overlap_engine_strategy = str(merge_cfg.get("overlap_engine_strategy", "single")).strip().lower()
    merge_overlap_engine_order = merge_cfg.get("overlap_engine_order", ["ungapped", "biopython", "edlib"])
    merge_overlap_engine_resolved = resolve_overlap_engine(merge_overlap_engine)
    merge_anchor_tolerance_bases = int(merge_cfg.get("anchor_tolerance_bases", 30))
    merge_min_quality = float(overlap_cfg.get("min_quality", 20.0))
    quality_mode = str(overlap_cfg.get("quality_mode", "warning")).strip().lower()
    if quality_mode not in {"blocking", "warning"}:
        quality_mode = "warning"
    ambiguity_identity_delta = float(overlap_cfg.get("ambiguity_identity_delta", 0.0))
    ambiguity_quality_epsilon = float(overlap_cfg.get("ambiguity_quality_epsilon", 0.1))
    high_conflict_q_threshold = int(overlap_cfg.get("high_conflict_q_threshold", 30))
    high_conflict_action = str(overlap_cfg.get("high_conflict_action", "flag")).strip().lower()
    if high_conflict_action not in {"flag", "route_cap3"}:
        high_conflict_action = "flag"
    cap3_validate_pair_support = bool(overlap_cfg.get("cap3_validate_pair_support", False))
    cap3_exe = cfg["tools"]["cap3"]

    in_dir = Path(input_dir).resolve() 
    out_dir = Path(output_dir).resolve() 
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = out_dir / "cap3_run_metadata.txt"

    # ==========================
    # Capturing CAP3 version discovery stdout/stderr so the log file shows what was reported 
    # by `cap3 --version` or the conda fallback.
    # ========================== 

    cap3_version = "unknown"
    try:
        version_result = subprocess.run(
            [cap3_exe, "--version"],
            check=True,
            capture_output=True,
            text=True,
        )
        if version_result.stdout:
            L.info("CAP3 --version stdout:\n%s", version_result.stdout.strip())
        if version_result.stderr:
            L.warning("CAP3 --version stderr:\n%s", version_result.stderr.strip())
        cap3_version = (version_result.stdout or version_result.stderr).strip() or "unknown"
    except (OSError, subprocess.CalledProcessError) as exc:
        L.warning("Failed to determine CAP3 version via --version: %s", exc)


    # =======================
    # Fallback here: capture the conda list output in logs when `cap3 --version` doesn't resolve version string and use that.
    # ======================= 

    if cap3_version == "unknown":
        for cmd in ("conda", "mamba"):
            try:
                version_result = subprocess.run(
                    [cmd, "list", "cap3"],
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except (OSError, subprocess.CalledProcessError):
                continue
            if version_result.stdout:
                L.info("%s list cap3 stdout:\n%s", cmd, version_result.stdout.strip())
            if version_result.stderr:
                L.warning("%s list cap3 stderr:\n%s", cmd, version_result.stderr.strip()) 
            for line in version_result.stdout.splitlines():
                if line.startswith("cap3"):
                    parts = line.split()
                    if len(parts) > 1:
                        cap3_version = f"{parts[0]} {parts[1]}"
                        break
            if cap3_version != "unknown":
                break

    cap3_logs_dir = out_dir / "logs" / "cap3" 

    metadata_lines = [
        f"cap3_executable\t{cap3_exe}",
        f"cap3_version\t{cap3_version}",
        "sample_id\tcap3_command\tcap3_stdout_log\tcap3_stderr_log",
    ]

    if pairing_report is not None:
        pairs, meta = group_pairs(
            in_dir,
            dup_policy=dup_policy,
            fwd_pattern=fwd_pattern,
            rev_pattern=rev_pattern,
            return_metadata=True,
            enforce_same_well=enforce_same_well,
            well_pattern=well_pattern 
        )
        _write_pairing_report(pairs, meta, Path(pairing_report))
    else:
        pairs = group_pairs(
            in_dir,
            dup_policy=dup_policy,
            fwd_pattern=fwd_pattern,
            rev_pattern=rev_pattern,
            enforce_same_well=enforce_same_well,
            well_pattern=well_pattern 
        )
        meta = {}

    if not pairs:
        L.info("No paired reads detected in input path: %s", in_dir)
        # Exit out since nothing here to do 
        return [] 
    
    # goign to collect the samples started with an empty list 
    contig_paths: list[Path] = [] 
    # Will process the ordered sampled paired list one by one 
    for sid in sorted(pairs):
        # Find the location of the orientation based on the file name 
        entries = pairs[sid]

        tasks: list[tuple[str, list[Path], Path | None, Path | None]] = [] 
        if dup_policy == DupPolicy.KEEP_SEPARATE:
            f_sources = _as_path_list(entries["F"])
            r_sources = _as_path_list(entries["R"])
            primer_pairs = _build_keep_separate_pairs(f_sources, r_sources)
            for idx, (fwd, rev) in enumerate(primer_pairs, start=1):
                sample_key = sid if len(primer_pairs) == 1 else f"{sid}_{idx}"
                tasks.append((sample_key, [fwd, rev], fwd, rev))
        else:
            f_sources = _as_path_list(entries["F"])
            r_sources = _as_path_list(entries["R"])
            sources = list(_iter_paths(entries["F"])) + list(_iter_paths(entries["R"]))
            fwd_path = f_sources[0] if len(f_sources) == 1 else None
            rev_path = r_sources[0] if len(r_sources) == 1 else None
            tasks.append((sid, sources, fwd_path, rev_path))

        if not tasks:
            raise RuntimeError("pairing produced no tasks")

        for sample_key, sources, fwd_path, rev_path in tasks: 
            sample_dir = out_dir / sample_key 
            sample_dir.mkdir(parents=True, exist_ok=True) 
            sample_fasta = sample_dir / f"{sample_key}_paired.fasta" 

            _write_combined_fasta(sources, sample_fasta, use_qual=use_qual) 

            if fwd_path and rev_path and len(sources) == 2:
                try:
                    contig_path, report = merge_two_reads(
                        sample_id=sample_key,
                        fwd_path=fwd_path,
                        rev_path=rev_path,
                        output_dir=sample_dir,
                        min_overlap=merge_min_overlap,
                        min_identity=merge_min_identity,
                        min_quality=merge_min_quality,
                        quality_mode=quality_mode,
                        ambiguity_identity_delta=ambiguity_identity_delta,
                        ambiguity_quality_epsilon=ambiguity_quality_epsilon,
                        high_conflict_q_threshold=high_conflict_q_threshold,
                        high_conflict_action=high_conflict_action,
                        overlap_engine=merge_overlap_engine,
                        overlap_engine_strategy=merge_overlap_engine_strategy,
                        overlap_engine_order=merge_overlap_engine_order if isinstance(merge_overlap_engine_order, list) else None,
                        anchor_tolerance_bases=merge_anchor_tolerance_bases,
                    )
                except MergeInputError as exc:
                    L.warning("Skipping merge_two_reads for %s: %s", sample_key, exc)
                    metadata_lines.append(
                        f"{sample_key}\tmerge_two_reads_skipped_non_singleton "
                        f"f_n={exc.f_count if exc.f_count is not None else 'na'} "
                        f"r_n={exc.r_count if exc.r_count is not None else 'na'} reason={exc}"
                    )
                else:
                    metadata_lines.append(
                        f"{sample_key}\tmerge_two_reads engine={report.overlap_engine} "
                        f"strategy={merge_overlap_engine_strategy} order={merge_overlap_engine_order} "
                        f"orientation={report.orientation} overlap={report.overlap_len} identity={report.identity:.4f}"
                    )
                    if contig_path:
                        L.info(
                            "Paired assembler path for %s: merge_two_reads accepted (%s)",
                            sample_key,
                            report.merge_status,
                        )
                        contig_paths.append(contig_path)
                        continue
                    L.info(
                        "Paired assembler path for %s: merge_two_reads=%s; fallback to CAP3",
                        sample_key,
                        report.merge_status,
                    )
                    if report.merge_status == "quality_low" and quality_mode == "blocking":
                        L.info(
                            "Paired assembler path for %s: quality_mode=blocking keeps singlets, CAP3 skipped",
                            sample_key,
                        )
                        continue

            if not (fwd_path and rev_path and len(sources) == 2):
                L.info(
                    "Paired assembler path for %s: merge_two_reads unavailable (non-singleton or missing pair); using CAP3",
                    sample_key,
                )
            cmd = [cap3_exe, sample_fasta.name]
            if cap3_options:
                cmd.extend(cap3_options)
            L.info("Run CAP3 (apired) %s: %s", sample_key, " ".join(cmd))

            metadata_lines.append(
                f"{sample_key}\t{' '.join(cmd)}\t\t"
            )

            try:
                result = subprocess.run(
                    cmd, 
                    check=True,
                    cwd=sample_dir,
                    capture_output=True,
                    text=True,
                )
                if result.stdout:
                    L.info("CAP3 stdout for %s:\n%s", sample_key, result.stdout)
                if result.stderr:
                    L.warning("CAP3 stderr for %s:\n%s", sample_key, result.stderr)
                stdout_path, stderr_path = write_cap3_process_logs(
                    cap3_logs_dir,
                    sample_id=sample_key,
                    command=cmd,
                    assembler_label="cap3_default",
                    stdout_text=result.stdout or "",
                    stderr_text=result.stderr or "",
                )
                metadata_lines[-1] = f"{sample_key}\t{' '.join(cmd)}\t{stdout_path}\t{stderr_path}"
            except subprocess.CalledProcessError as exc:
                stdout_path, stderr_path = write_cap3_process_logs(
                    cap3_logs_dir,
                    sample_id=sample_key,
                    command=cmd,
                    assembler_label="cap3_default",
                    stdout_text=exc.stdout or "",
                    stderr_text=exc.stderr or "",
                )
                metadata_lines[-1] = f"{sample_key}\t{' '.join(cmd)}\t{stdout_path}\t{stderr_path}"
                L.error("CAP3 failed for sample %s (exit %s):\n%s", sample_key, exc.returncode, exc.stderr
                )
                if exc.stdout:
                    L.error("CAP3 stdout for %s:\n%s", sample_key, exc.stdout)
                raise # stop the run here on raise 

            contig_path = sample_dir / f"{sample_key}_paired.fasta.cap.contigs"
            if not contig_path.exists():
                raise FileNotFoundError(contig_path)

            if cap3_validate_pair_support and fwd_path and rev_path and len(sources) == 2:
                validation = _validate_cap3_contig_support(contig_path, fwd_path, rev_path)
                marker_path = sample_dir / f"{sample_key}_paired.cap3_validation.txt"
                marker_path.write_text(f"{validation}\n", encoding="utf-8")
                metadata_lines.append(f"{sample_key}	cap3_validation={validation}")

                if validation == "rejected":
                    L.warning("CAP3 output for %s rejected by ACE membership check; skipping contig", sample_key)
                    singlets_path = sample_dir / f"{sample_key}_paired.fasta.cap.singlets"
                    fwd_record = next(SeqIO.parse(fwd_path, "fasta"), None)
                    rev_record = next(SeqIO.parse(rev_path, "fasta"), None)
                    singlet_records = [rec for rec in (fwd_record, rev_record) if rec is not None]
                    if singlet_records:
                        SeqIO.write(singlet_records, singlets_path, "fasta")
                    try:
                        contig_path.unlink()
                    except FileNotFoundError:
                        pass
                    continue
                if validation == "unknown":
                    L.warning(
                        "CAP3 output for %s could not be validated (ACE missing/unreadable); keeping contig",
                        sample_key,
                    )

            contig_paths.append(contig_path)
            L.info("Cap3 paired assembly finished for %s and for contigs: %s", sample_key, contig_path) 
    
    metadata_path.write_text("\n".join(metadata_lines) + "\n", encoding="utf-8")

    return contig_paths 

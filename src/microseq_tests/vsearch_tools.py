from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Callable

from microseq_tests.utility.utils import load_config

L = logging.getLogger(__name__)


def resolve_vsearch() -> str:
    cfg = load_config()
    vsearch = cfg.get("tools", {}).get("vsearch", "vsearch")
    if Path(vsearch).exists():
        return vsearch
    found = shutil.which(vsearch)
    if found:
        return found
    raise FileNotFoundError(
        "vsearch not found. Install vsearch or set tools.vsearch in config.yaml."
    )


def run_vsearch(cmd: list[str], *, cwd: Path | None = None) -> None:
    vsearch = resolve_vsearch()
    full_cmd = [vsearch, *cmd]
    L.info("Running vsearch: %s", " ".join(full_cmd))
    try:
        subprocess.run(
            full_cmd,
            check=True,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        L.error("vsearch failed (exit %s):\n%s", exc.returncode, exc.stderr)
        raise


def collapse_replicates(
    fasta_in: Path,
    fasta_out: Path,
    *,
    min_size: int = 1,
    id_th: float | None = None,
    threads: int | None = None,
    uc_out: Path | None = None,
) -> Path:
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    derep_out = fasta_out.with_suffix(".derep.fasta")
    cmd = [
        "--derep_fulllength",
        str(fasta_in),
        "--output",
        str(derep_out),
        "--sizeout",
        "--minuniquesize",
        str(min_size),
        "--fasta_width",
        "0",
    ]
    if threads:
        cmd.extend(["--threads", str(threads)])
    run_vsearch(cmd)

    if id_th is not None and id_th < 1.0:
        cmd = [
            "--cluster_size",
            str(derep_out),
            "--id",
            f"{id_th:.4f}",
            "--centroids",
            str(fasta_out),
            "--sizein",
            "--sizeout",
            "--fasta_width",
            "0",
        ]
        if uc_out:
            cmd.extend(["--uc", str(uc_out)])
        if threads:
            cmd.extend(["--threads", str(threads)])
        run_vsearch(cmd)
        derep_out.unlink(missing_ok=True)
        return fasta_out

    derep_out.replace(fasta_out)
    return fasta_out


def chimera_check_ref(
    fasta_in: Path,
    fasta_out: Path,
    *,
    reference: Path,
    report_tsv: Path | None = None,
    threads: int | None = None,
    size_in: bool = False,
) -> tuple[Path, Path]:
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    report_tsv = report_tsv or fasta_out.with_suffix(".uchime.tsv")
    cmd = [
        "--uchime_ref",
        str(fasta_in),
        "--db",
        str(reference),
        "--nonchimeras",
        str(fasta_out),
        "--uchimeout",
        str(report_tsv),
        "--fasta_width",
        "0",
    ]
    if size_in:
        cmd.append("--sizein")
    if threads:
        cmd.extend(["--threads", str(threads)])
    run_vsearch(cmd)
    return fasta_out, report_tsv


def collapse_replicates_grouped(
    fasta_in: Path,
    fasta_out: Path,
    *,
    group_fn: Callable[[str], str],
    min_size: int = 1,
    id_th: float | None = None,
    threads: int | None = None,
    map_tsv: Path | None = None,
    weights_tsv: Path | None = None,
) -> Path:
    from Bio import SeqIO
    import re
    import tempfile

    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    map_tsv = map_tsv or fasta_out.with_name("replicate_map.tsv")
    weights_tsv = weights_tsv or fasta_out.with_name("replicate_weights.tsv")

    size_re = re.compile(r";size=(\d+);?")
    grouped: dict[str, list] = {}
    for rec in SeqIO.parse(fasta_in, "fasta"):
        sample_id = group_fn(rec.id)
        grouped.setdefault(sample_id, []).append(rec)

    with (
        map_tsv.open("w", encoding="utf-8") as map_fh,
        weights_tsv.open("w", encoding="utf-8") as weight_fh,
        fasta_out.open("w", encoding="utf-8") as out_fh,
        tempfile.TemporaryDirectory(prefix="microseq_reps_") as tmpdir,
    ):
        map_fh.write("sample_id\tcentroid_id\treplicate_size\tmember_ids\n")
        weight_fh.write("qseqid\treplicate_size\n")

        tmp_root = Path(tmpdir)
        for sample_id, records in grouped.items():
            if not records:
                continue
            seq_members: dict[str, list[str]] = {}
            for rec in records:
                seq_members.setdefault(str(rec.seq), []).append(rec.id)

            group_in = tmp_root / f"{sample_id}_input.fasta"
            SeqIO.write(records, group_in, "fasta")

            group_out = tmp_root / f"{sample_id}_collapsed.fasta"
            uc_path = tmp_root / f"{sample_id}.uc"
            collapse_replicates(
                group_in,
                group_out,
                min_size=min_size,
                id_th=id_th,
                threads=threads,
                uc_out=uc_path,
            )
            centroids = list(SeqIO.parse(group_out, "fasta"))

            cluster_members: dict[str, list[str]] = {}
            cluster_sizes: dict[str, int] = {}
            if uc_path.exists() and id_th is not None and id_th < 1.0:
                for line in uc_path.read_text().splitlines():
                    if not line or line.startswith("#"):
                        continue
                    cols = line.split("\t")
                    record_type = cols[0]
                    query_label = cols[8]
                    target_label = cols[9] if len(cols) > 9 else query_label
                    centroid_label = target_label if record_type == "H" else query_label
                    match = size_re.search(query_label)
                    size = int(match.group(1)) if match else 1
                    cluster_sizes[centroid_label] = cluster_sizes.get(centroid_label, 0) + size
                    cluster_members.setdefault(centroid_label, []).append(query_label)

            for idx, rec in enumerate(centroids, 1):
                original_id = rec.id
                match = size_re.search(original_id)
                size = int(match.group(1)) if match else 1
                clean_id = f"{sample_id}__rep{idx}"
                tagged_id = f"{clean_id};size={size};"
                if original_id in cluster_sizes:
                    size = cluster_sizes[original_id]
                    tagged_id = f"{clean_id};size={size};"
                rec.id = tagged_id
                rec.name = tagged_id
                rec.description = tagged_id
                SeqIO.write(rec, out_fh, "fasta")

                if original_id in cluster_members:
                    member_ids = ",".join(cluster_members[original_id])
                else:
                    member_ids = ",".join(seq_members.get(str(rec.seq), []))

                map_fh.write(f"{sample_id}\t{clean_id}\t{size}\t{member_ids}\n")
                weight_fh.write(f"{clean_id}\t{size}\n")

    return fasta_out


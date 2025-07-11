# MicroSeq

<!-- README badge -->
[![Docs](https://img.shields.io/badge/docs-view-blue)](https://jose-cantu.github.io/MicroSeq/)

**MicroSeq** is an end‑to‑end microbial sequence‑analysis toolkit that scales from a single laptop to a multi‑node HPC cluster.
It automates everything from raw trace files (AB1) or FASTQ reads all the way to BIOM/CSV tables that are ready for downstream statistics or visualisation.

*Current status – May 2025*: CLI workflow is feature‑complete and battle‑tested internally; the PySide6 GUI will land in v1.0. Rust speed‑ups  for alignment and I/O are planned for v2.0. Note that the GUI is nowhere near done. I plan on removing and adding things to it. I need to go back to the whiteboard and rethink some things. The roadmap is below in the roadmap section, close to the end. 

For version control, please see the changelog I have attached. [CHANGELOG](https://github.com/jose-cantu/MicroSeq/blob/master/CHANGELOG.MD)

For the docs that go into more detail, such as for example how to use the GUI I have attached it here. [GUI walkthrough](https://jose-cantu.github.io/MicroSeq/gui/)
You can also just select docs view and it will take you to the landing page for full detailed documentation. =) 

---

## Table of Contents

1. [Why MicroSeq?](#why-microseq)
2. [Features](#features)
3. [Quick start (TL;DR)](#quick-start-tldr)
4. [Prerequisites](#prerequisites)
5. [Detailed installation](#detailed-installation)
6. [Logging & run locations](#logging--run-locations)
7. [Choosing the right `--threads`](#choosing-the-right-threads)
8. [Command reference](#command-reference)
9. [End‑to‑end examples](#end-to-end-examples)
10. [QC & Assembly tools explained](#qc--assembly-tools-explained)
11. [CAP3 (microseq assembly)](#cap3-microseq-assembly)
12. [Advanced topics](#advanced-topics)
13. [Contributing / Support](#contributing--support)
14. [License & Citation](#license--citation)

---

## Why MicroSeq?

MicroSeq began as an open-source replacement for the narrow slice of Geneious Prime we relied on turning raw Sanger chromatograms into taxonomically annotated sequences but has grown into the only FLOSS pipeline that goes from AB1 or FASTQ input to a BIOM/OTU table in one reproducible run.

* **Dual input.** Accepts raw `.ab1` traces or demultiplexed FASTQ amplicons.
* **Contig aware (planned v1.1).** When more than one read belongs to an isolate, MicroSeq will automatically assemble them with CAP3, boosting consensus quality. In the current version run `microseq assembly` first.
* **GreenGenes 2, SILVA, or NCBI one switch.** `--db gg2|silva|ncbi` works in both CLI and GUI; MicroSeq is currently the only tool that ships a GreenGenes 2–indexed classifier.
* **Immediate downstream compatibility.** `microseq postblast` writes BLAST results directly to BIOM (or flat CSV) so the same file can feed QIIME‑2, phyloseq, or scikit‑bio without hand conversion.
* **Bench-to-cluster consistency.** A single YAML config runs unchanged on a laptop or under Slurm; thread and memory limits are flags, not hard-coded paths.
* **Permissive licence.** MIT/LGPL so the lab can redistribute binaries internally without negotiating commercial terms.

In short: MicroSeq is the only open-source package that ingests raw Sanger or FASTQ reads, classifies them against GreenGenes 2, SILVA, or NCBI in one step, and emits a BIOM table all under a scriptable CLI with an optional GUI for bench scientists. Automatic contig assembly is scheduled for version 1.1.

---

## Features

* **One‑liner installer** (`microseq‑setup`) – downloads & indexes **Greengenes 2**, **SILVA 138.1**, and **NCBI 16S** databases; writes a ready‑to‑go `config.yaml` and shell exports.
* **Quality‑control & trimming**

  * AB1 → FASTQ (`microseq ab1-to-fastq`)
  * FASTQ QC (Phred, length, adapters) with Trimmomatic / BioPython (`microseq trim`).
* **De‑novo assembly** using *CAP3* (`microseq assembly`).
* **Taxonomic classification**

  * BLAST+ search (**FASTA** or **FASTQ** input) (`microseq blast`)
  * Post‑BLAST taxonomy join (`microseq add_taxonomy`).
* **Post‑processing**

  * Collapse hits, choose ID/Q‑cov cut‑offs, export BIOM / CSV tables (`microseq postblast`).
* **Config‑driven defaults** for duplicate‑policy & ID‑normalisers (regex or built‑ins).
* **Robust logging**  – session-based logs via $MICROSEQ_SESSION_ID or timestamp fallback; logs can live anywhere you point `MICROSEQ_LOG_DIR`.
* **TSV normalisation** (`microseq-fix-tsv`) – convert CSV or space‑delimited tables into proper TSV.
* **Portable**  – pure‑Python (≥3.10) + external binaries already in Conda env (BLAST+, TaxonKit, CAP3).
* **GUI preview** (`microseq-gui`) – drag‑and‑drop queue, live log pane, progress
  bar with ETA.  Running QC or the full pipeline now writes
  `qc/trim_summary.tsv` automatically. The window provides **Run QC** and
  **Full pipeline** buttons to execute the entire workflow (trim → BLAST →
  taxonomy → optional BIOM generation) with a single click.

---

## Quick start (TL;DR) 

```bash
# 1  Clone & enter repo
$ git clone https://github.com/jose-cantu/MicroSeq.git
$ cd MicroSeq

# Note for MacOS users ONLY you can just use the MacOS installation bash file I have set up from here and it will install all the rest for you. Just follow the prompts I have setup in it such as asking for miniconda and once at the wizard setup stage simply press enter throughout for default steps....
./install_macos_microseq.sh

# 2  Create & activate Conda environment
$ conda env create -f config/environment.yml
$ conda activate MicroSeq

# 3  Editable install (good while developing)
$ pip install -e .

# 4  Download reference DBs & configure logging
$ microseq-setup                 # interactive, 3–5 min first run
#     └─ adds $MICROSEQ_DB_HOME, $MICROSEQ_LOG_DIR, $BLASTDB, $BLASTDB_LMDB=0 to ~/.bashrc

# 5  Open a new shell or:
$ source ~/.bashrc               # pick up the new exports

# 6  Smoke test (tiny demo FASTA)
$ printf ">demo\nACGT\n" > tiny.fasta
  for d in gg2 silva ncbi; do
  microseq blast -i tiny.fasta -d $d -o /tmp/$d.tsv \
                 --identity 50 --qcov 10 \
      && echo "[✓] $d OK"  
  done 
  rm tiny.fasta /tmp/*.tsv
```

> **Need the tools first?**  • [git ≥ 2.30](https://git-scm.com/)  • [Miniconda / Mamba](https://docs.conda.io/en/latest/miniconda.html)  • [Python ≥ 3.10](https://www.python.org/downloads/) – install those in any order, then follow the steps above. 

macOS users can run `./install_macos_microseq.sh` for a guided setup. If Git
is missing, the script launches the Xcode Command Line Tools installer. Once
the tools finish (they install Git and `clang` together), re-run the script to
continue.

Keep in mind for MacOS install the intel variation provided from miniconda and NOT Apple Silicon. `cap3` is published on Bioconda only for the osx-64 (Intel/x86-64) platform, not for osx-arm64 (Apple Silicon). The conda environment will work on an Apple Silicon MacOS (I'm using one right now). Your conda solver will advertise `platform: osx-64` so it can find and install that build. 

---

## Prerequisites

| Tool                   | Why MicroSeq needs it                                   | Install guide                                                                                                                  |
| ---------------------- | ------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| git                    | Clone the repo, pull updates                            | [https://git-scm.com/book/en/v2/Getting-Started-Installing-Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) |
| Conda / Mamba          | Manages the `MicroSeq` env incl. BLAST+, CAP3, TaxonKit | [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)                               |
| Python 3.10 – 3.13     | Runs the CLI & GUI                                      | [https://www.python.org/downloads/](https://www.python.org/downloads/)                                                         |
| BLAST+, TaxonKit, CAP3 | Auto‑installed into the Conda env – no action needed    | –                                                                   | Brew                   | For MacOS users only needed for running blast smoothly  | [https://brew.sh] 

 Brew                    For MacOS users only needed for running blast smoothly    [https://brew.sh] 


---

## Detailed installation

### 1 Create the Conda environment

```bash
conda env create -f config/environment.yml  # or: mamba env create …
conda activate MicroSeq
```

### 2 Install MicroSeq (editable vs. frozen)

* **Editable (dev) install** – `pip install -e .` (recommended while the API is still evolving).
* **Frozen install** – `pip install .` (locks dependencies to the current version).
* **Upgrade to a newer patch or feature** – `git pull && pip install -e .` (omit `-e` if you use the frozen install).

### 3 Download reference databases

Run once per machine (or whenever DB versions change):

```bash
microseq-setup                # interactive prompts
# flags: --quiet  •  --db-root /path/to/dbs  •  --log-dir /path/to/logs
# answer "y" when asked about pipelines to aggregate logs by session ID
```

What it does:

1. Downloads & indexes GG2, SILVA, NCBI 16S into `$MICROSEQ_DB_HOME` (default `~/.microseq_dbs`).
2. Creates the log folder.
3. Appends export block to `~/.bashrc` **and** `etc/conda/activate.d/microseq.sh` inside the env.
4. Writes or patches `config/config.yaml` so library calls inherit the same paths.

Re‑run anytime with different flags to move DBs or logs elsewhere.
Example scripts to initialise a logging session are under
`config/templates/activate_microseq_session.*.example`.  Copy the variant for
your shell and source it before running MicroSeq.

---

## Logging & run locations

| Variable / Flag     | Purpose                                                           | Default                    | Example override                            |
| ------------------- | ----------------------------------------------------------------- | -------------------------- | ------------------------------------------- |
| `--workdir` (CLI)   | Per‑run scratch folder for intermediate files                     | `default_workdir` from `config.yaml` (e.g. `~/microseq_runs`)          | `--workdir /scratch/$USER/run42`            |
| `MICROSEQ_LOG_DIR`  | Where rotated log files are written                               | `<repo>/logs`              | `export MICROSEQ_LOG_DIR=/var/log/microseq` |
| `MICROSEQ_LOG_FILE` | Override log file location entirely | — | `export MICROSEQ_LOG_FILE=/tmp/microseq.log` |
| `--session-id` (CLI) | Set the value of `MICROSEQ_SESSION_ID` for this run | — | `--session-id my_run` |
| `MICROSEQ_SESSION_ID` | Aggregates multiple runs into one log when set; otherwise timestamped files | — | `export MICROSEQ_SESSION_ID=my_run` |

A symlink `logs/microseq_latest.log` always points at the active log file, so `tail ‑F logs/microseq_latest.log` just works.

If `MICROSEQ_SESSION_ID` is **not** set, the logger prints one warning and
generates a random session ID, e.g.:

```
⚠️  MICROSEQ_SESSION_ID not set – using auto session ID 20250101-123456-abcd
   (export MICROSEQ_SESSION_ID=YOUR_ID to aggregate multiple commands)
```

Define the variable ahead of time to suppress the warning and group multiple runs
into a single log file.
> **Note:** If multiple processes write to the same log file (e.g. job arrays with the same session ID), NFS filesystems might interleave lines. If this occurs, give each task a unique `MICROSEQ_SESSION_ID`.

Example helper scripts to do so live in
`config/templates/activate_microseq_session.sh.example` (bash or zsh) and
`config/templates/activate_microseq_session.fish.example`.  Copy the one
matching your shell and source it on login to keep all commands under the same
session ID.

(Adjusted microseq-setup to include more than just bash also zsh and fish and created config files to use after setup depending on scenario)

Old logs are kept indefinitely. Tweak `logging.backup_count` in
`config/config.yaml` (or pass `backup_count` to `setup_logging` when scripting)
to retain only the last *N* auto-generated files. The default value of `0`
preserves every log.


### Workdir basics
When omitted, `--workdir` defaults to the `default_workdir` value in `config.yaml`. Each run creates `qc/`, `asm/` and `blast/` folders inside this directory. Use `--workdir PATH` to place intermediates elsewhere or `--workdir .` to write outputs next to your input files.

### Choosing the right `--threads`

| Scenario | Quick answer | Why |
|----------|--------------|-----|
| **Laptop, 8 GB RAM** | `microseq blast … --threads $(microseq recommend-threads)` | Leaves ~2 GB for the OS, avoids swap. |
| **Workstation, 32 cores / 64 GB** | Same as above, usually returns 16. | Diminishing returns beyond 16 for BLAST; memory still the limiter. |
| **SLURM / cluster** | Let your job script set `--threads $SLURM_CPUS_ON_NODE` | Scheduler already assigns CPUs & cgroups; auto-detection is unreliable. |

> **Tip:** `microseq recommend-threads` inspects *this* machine only. On shared clusters the admin may cap RAM per core, so follow the queue’s guidelines instead.

---

## Command reference

MicroSeq’s CLI is a single entry‑point (`microseq`) with **global flags** followed by a sub‑command.  Every sub‑command has its own options; run `microseq <cmd> -h` for authoritative help.

### Global flags (apply to *all* sub‑commands)

| Flag             | Default | Purpose                                                                                                                     |
| ---------------- | ------- | --------------------------------------------------------------------------------------------------------------------------- |
| `--workdir PATH` | `~/microseq_runs` | Root folder for intermediate & derived files created by the current run. A fresh folder keeps each analysis self‑contained. |
| `--threads N`    | `4`     | How many CPU threads to give the stage (BLAST, CAP3, Trimmomatic, etc.). 2–4 threads work well on a typical laptop.                                                    |
| `-v`, `-vv`      | off     | `-v` = info‑level logs, `-vv` = debug.                                                                                      |

---

### Sub‑commands & flags

Below is an index of every sub‑command with a one‑line summary.  Click the ▶ arrows to reveal full flag tables and use‑case notes.

| Sub‑command       | One‑sentence purpose                                                    |
| ----------------- | ----------------------------------------------------------------------- |
| `trim`            | QC‑trim FASTQ or AB1 traces into high‑quality reads                     |
| `ab1-to-fastq`    | Convert chromatograms (.ab1) directly to FASTQ                          |
| `fastq-to-fasta`  | Merge one or many FASTQ files into a single FASTA                       |
| `assembly`        | De‑novo assembly of trimmed reads with CAP3                             |
| `blast`           | Search reads/contigs against Greengenes 2, SILVA, or NCBI 16S           |
| `suggest-cutoffs` | Scan sweeper table to find ID/Q‑cov pairs that yield desired PASS count |
| `filter-hits`     | Apply chosen ID/Q‑cov to sweeper TSV and report summary                 |
| `merge-hits`      | Concatenate multiple BLAST TSVs into one master file                    |
| `add_taxonomy`    | Append 7‑rank lineage column to a BLAST hits table                      |
| `postblast`       | Collapse hits per sample, join metadata, export BIOM/CSV/JSON           |
| `microseq-setup`  | Download reference DBs and write env‑hook + config                      |
| `microseq-fix-tsv` | Normalise CSV/space tables into TAB-separated format |

---

| Flag                     | Required | Example                       | When to use                                                                |
| ------------------------ | -------- | ----------------------------- | -------------------------------------------------------------------------- |
| `-i / --input PATH`      | ✓        | `reads.fastq` or `traces.ab1` | Input file or folder.                                                      |
| `-o / --output PATH`     |          | `qc/trimmed.fastq`            | Override auto‑named file in `workdir/qc/`.                                 |
| `--sanger / --no-sanger` |          | `--sanger`                    | Use BioPython’s Sanger‑specific trim as a preference.                      |
| `--link-raw`             |          |                               | Symlink AB1 traces into `workdir/raw_ab1` instead of copying (saves disk). |

| Flag                    | Required | Example                                  |
| ----------------------- | -------- | ---------------------------------------- |
| `-i / --input_dir DIR`  | ✓        | `raw_traces/` (recurses)                            |
| `-o / --output_dir DIR` | ✓        | `fastq/`                                 |
| `--overwrite`           |          | Re‑convert even if FASTQs already exist. |

| Flag                       | Required | Example            |
| -------------------------- | -------- | ------------------ |
| `-i / --input_dir DIR`     | ✓        | `qc/` (recurses)   |
| `-o / --output_fasta PATH` | ✓        | `qc/trimmed.fasta` |

| Flag                 | Required | Example            |
| -------------------- | -------- | ------------------ |
| `-i / --input FASTA` | ✓        | `qc/trimmed.fasta` |
| `-o / --output DIR`  | ✓        | `asm/`             |

| Flag                                          | Default | Example                                                                      | Scenario                           |
| --------------------------------------------- | ------- | ---------------------------------------------------------------------------- | ---------------------------------- |
| `-i / --input FASTA`                          | ✓       | `asm/contigs.fasta`                                                          | Query sequences.                   |
| `-d / --db {gg2,silva,ncbi}`                  | ✓       | `--db gg2`                                                                   | Choose reference DB.               |
| `-o / --output TSV`                           | ✓       | `hits.tsv`                                                                   | Where results go.                  |
| `--identity PCT`                              | 97      | `--identity 99`                                                              | Stricter species‑level hits.       |
| `--qcov PCT`                                  | 80      | 70                                                                           | Lower if you have short amplicons. |
| `--max_target_seqs N`                         | 5       | 1                                                                            | Keep only top hit.                 |
| `--relaxed` + `--relaxed-id / --relaxed-qcov` | —       | First coarse pass before cut‑off optimisation (use with `--export-sweeper`). |                                    |
| `--export-sweeper`                            | —       | Produce `hits_full_sweeper.tsv` for `suggest-cutoffs`.                       |                                    |
| `--clean-titles`                              | —       | Strip extra fields from SILVA/NCBI `stitle`.                                 |                                    |
| `--log-missing FILE`                          | —       | Save IDs with no hits for troubleshooting.                                   |                                    |

| Arg / Flag | Example                 | Scenario                                  |
| ---------- | ----------------------- | ----------------------------------------- |
| `table`    | `hits_full_sweeper.tsv` | Output of a relaxed BLAST.                |
| `target`   | `530`                   | Aim for 530 PASS rows after collapse.     |
| `--step`   | `2`                     | Scan in 2 % increments (faster, coarser). |

| Flag                | Required       | Example                          | Notes |
| ------------------- | -------------- | -------------------------------- | ----- |
| `-i / --input TSV`  | ✓              | `hits_full_sweeper.tsv`          |       |
| `--identity PCT`    | ✓              | 97                               |       |
| `--qcov PCT`        | ✓              | 80                               |       |
| `-o / --output TSV` |                | auto‑named                       |       |
| `--with-status`     |                | Write PASS/FAIL column.          |       |
| `--dry-run`         |                | Show counts, write nothing.      |       |
| `--unique`          |                | Also report unique sample count. |       |
| `--group-col COL`   | `sample_id`    | Collapse key.                    |       |
| `--id-normaliser`   | `strip_suffix` | Sample ID cleaning. Choices: `none`, `strip_suffix`, `auto`, `strip_suffix_simple`, `strip_suffix_legacy`, `config`. |

| Flag                | Required | Example            |
| ------------------- | -------- | ------------------ |
| `-i / --input TSV…` | ✓        | `blast/*/hits.tsv` |
| `-o / --output TSV` | ✓        | `all_hits.tsv`     |

| Flag                         | Required | Example                                    |
| ---------------------------- | -------- | ------------------------------------------ |
| `-i / --hits TSV`            | ✓        | `hits.tsv`                                 |
| `-d / --db {gg2,silva,ncbi}` | ✓        | `--db silva`                               |
| `-o / --output TSV`          | ✓        | `hits_tax.tsv`                             |
| `--fill-species`             |          | For SILVA genus‑level lines with ≥97 % ID. |

| Flag                                          | Required       | Example                                    | Purpose                       |
| --------------------------------------------- | -------------- | ------------------------------------------ | ----------------------------- |
| `-b / --blast_file TSV`                       | ✓              | `hits_tax.tsv`                             | Input table (needs taxonomy). |
| `-m / --metadata TSV`                         | ✓              | `samples_meta.tsv`                         | Sample metadata.              |
| `-o / --output_biom BIOM`                     | ✓              | `results.biom`                             | Output; CSV auto‑written.     |
| `--json`                                      | —              | Converts HDF5 into `.json` BIOM.           |                               |
| `--sample-col COL`                            | `SampleID`     | If metadata uses a different ID column.    |                               |
| `--post_blast_identity PCT`                   | 97             | Tweak identity filter when multiple hits.  |                               |
| `--id-normaliser`                             | `strip_suffix` | Clean Sample IDs before join. Choices: `none`, `strip_suffix`, `auto`, `strip_suffix_simple`, `strip_suffix_legacy`, `config`.              |
| `--taxonomy-col COL`                          | auto           | Pre‑existing taxonomy column in metadata.  |                               |
| `--duplicate-policy {error,keep-first,merge}` | `error`        | Resolve duplicate IDs after normalisation. |                               |

This step also writes `<output>_blast_provenance.csv`, `<output>_taxonomy_only.csv` and `<output>_metadata.csv` alongside the BIOM file.
| Flag                        | Default           | Purpose                                 |
| --------------------------- | ----------------- | --------------------------------------- |
| `--db-root DIR`             | `~/.microseq_dbs` | Where to store BLAST DBs.               |
| `--log-dir DIR`             | `<repo>/logs`     | Central log folder.                     |
| `--quiet`                   | —                 | Non‑interactive install (use defaults). |

---

## Usage Scenarios at a Glance

| Task                                                                                    | Typical command                                                                                                           |
| --------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| **First‑time setup on a fresh laptop**                                                  | `git clone … && conda env create … && conda activate MicroSeq && pip install -e . && microseq-setup --quiet` |
| **Update to the latest patch or new feature**          | `git pull && pip install -e .` |
| **Convert ABI traces + QC + assemble + classify**                                       | See AB1 trace workflow in Examples.                                                                                       |
| **Optimise BLAST cut‑offs for 96‑well plate**   **Keep in mind, these are estimates!**  | `microseq blast --relaxed --export-sweeper … → microseq suggest-cutoffs … 530 → microseq filter-hits`                     |
| **Batch merge hits from three runs**                                                    | `microseq merge-hits -i run*/blast/hits.tsv -o merged.tsv`                                                                |
| **HPC Slurm job**                                                                       | Pre‑export `MICROSEQ_LOG_DIR=$SCRATCH/logs` and `MICROSEQ_SESSION_ID=$SLURM_JOB_ID`; use `--workdir $SLURM_TMPDIR` in each step.                                 |

---

## End‑to‑end examples

### AB1 trace workflow

```bash
microseq trim        -i raw_ab1/ --sanger --workdir run_01
microseq assembly    -i run_01/qc/trimmed.fasta -o run_01/asm
microseq blast       -i run_01/asm/contigs.fasta -d gg2 -o run_01/blast/hits.tsv
microseq add_taxonomy -i run_01/blast/hits.tsv -d gg2 -o run_01/blast/hits_tax.tsv
microseq postblast   -b run_01/blast/hits_tax.tsv -m metadata.tsv -o run_01.biom
```

### FASTQ‑only (Sanger), SILVA DB, HPC style

```bash
# sbatch or Nextflow will substitute $SLURM_TMPDIR and $SLURM_CPUS_ON_NODE
microseq trim            -i reads.fastq.gz --sanger --workdir $SLURM_TMPDIR/qc \
                         --threads $SLURM_CPUS_ON_NODE
microseq fastq-to-fasta  -i $SLURM_TMPDIR/qc -o reads.fasta   # optional, blast accepts FASTQ
# relaxed search to generate sweeper table
microseq blast           -i reads.fasta -d silva -o hits_full_sweeper.tsv \
                         --threads $SLURM_CPUS_ON_NODE --relaxed --export-sweeper \
                         --relaxed-id 80 --relaxed-qcov 0
# find cut-offs that keep ~530 samples after collapse
microseq suggest-cutoffs hits_full_sweeper.tsv 530 --step 2 > cutoff.txt
# apply the chosen thresholds (here 97/70)
microseq filter-hits -i hits_full_sweeper.tsv --identity 97 --qcov 70 \
                     -o hits.tsv --with-status --unique
# join taxonomy, collapse, export BIOM/CSV/JSON
microseq add_taxonomy -i hits.tsv -d silva -o hits_tax.tsv --fill-species
microseq postblast    -b hits_tax.tsv -m meta.tsv -o results.biom --json \
                      --duplicate-policy merge --id-normaliser strip_suffix
```

## QC & Assembly tools explained

### Trimmomatic vs. BioPython Sanger pipeline

| Mode              | Trigger                                                                     | Toolchain                                               | What happens                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| ----------------- | --------------------------------------------------------------------------- | ------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **General FASTQ** | input is `.fastq` and `--no-sanger` *(default)*                             | **Trimmomatic SE** (`SLIDINGWINDOW:5:20`, `MINLEN:200`) | Reads whose 5‑bp window drops below Q20 or end < 200 bp are discarded; survivors go to `workdir/qc/trimmed.fastq`.                                                                                                                                                                                                                                                                                                                                                    |
| **Sanger traces** | one or more `.ab1` chromatograms **or** any FASTQ processed with `--sanger` | **BioPython + MicroSeq**                                | 1. **AB1 → FASTQ** (`ab1_to_fastq.py`) – one FASTQ per trace.<br>2. **Dynamic sliding‑window trim** (5 bp window, Q ≥ 20) per read.<br>3. Per‑read Phred, length & per‑file averages → `*_avg_qual.txt`.<br>4. File passes if avg Q ≥ 20 → `passed_qc_fastq/`; otherwise goes to `failed_qc_fastq/` with its stats.<br>5. Combined TSV summarises file‑level metrics.<br>6. **Passed FASTQs merged** into `workdir/qc/trimmed.fasta` for downstream assembly & BLAST. |

*Why two paths?* AB1 chromatograms carry single‑end Sanger qualities that behave differently from generic FASTQ scores; the lightweight BioPython route avoids over‑trimming while still filtering low‑quality tails.

## CAP3 (microseq assembly)

CAP3 is inherently single-threaded; the `--threads` flag is accepted for API
consistency but has no effect.

CAP3 names every output file `<prefix>.cap.<suffix>` where `<prefix>` is the input FASTA without extension.  The key artefacts are:

| Suffix             | File example               | Purpose                                                                                                                                |
| ------------------ | -------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `contigs`          | `sample.cap.contigs`       | **FASTA** of all contigs built from ≥2 reads.                                                                                          |
| `contigs.qual`     | `sample.cap.contigs.qual`  | **QUAL** file – Phred scores per base in each contig (same order as above FASTA).                                                      |
| `singlets`         | `sample.cap.singlets`      | Reads that did **not** assemble with any other (singletons), in FASTA.                                                                 |
| `singlets.qual`    | `sample.cap.singlets.qual` | Phred scores for the singlets.                                                                                                         |
| `info`             | `sample.cap.info`          | Assembly log: contig IDs, lengths, coverage, clipped read coordinates, N50 summary.                                                    |
| `ace` *(optional)* | `sample.cap.ace`           | ACE format view of the assembly; generated when MicroSeq calls CAP3 with `-a`. Importable into Consed/GAP5/Tablet for visual curation. |

**How MicroSeq uses them**

1. Copies `*.contigs` → `asm/contigs.fasta` for downstream BLAST. `*.singlets` remain in `workdir/asm/` should you wish to analyse them separately.
2. Parses `*.info` to log assembly statistics (number of contigs, N50, max length).
3. Leaves all original CAP3 files in `workdir/asm/` so you can inspect quality or visualise in external tools.

Note when blasting there are two algorithms: The Fast (megablast ≥ 95 % ID) is sufficient for typical Sanger 16 S workflows. If the best hit drops below 90 % identity or coverage, MicroSeq automatically repeats the search with Comprehensive and more sensitive (blastn) algorithm to catch divergent taxa.
## microseq-fix-tsv

`microseq-fix-tsv` converts comma- or space-delimited tables into tab-separated files in-place. Use it to normalise spreadsheets before feeding them to other MicroSeq commands.


---

## Advanced topics

* **Custom databases** – point `microseq blast -d /path/my_db` and supply a matching taxonomy TSV to `microseq add_taxonomy`.
* **Nextflow/Snakemake integration** – export `MICROSEQ_SESSION_ID` to your scheduler's run ID (e.g. `$SLURM_JOB_ID`, `${workflow.runName}`) and use per-process `--workdir` to avoid clashes.

---

## Road-Map

All work lives in a public GitHub board project, and I will be actively working on the priorities listed, starting with P0. 

[View the board »](https://github.com/users/jose-cantu/projects/2/views/2)

### Current priorities
- P0 – Refresh README and create CHANGELOG
- P1 – `--blast-task` flag (default `megablast`)
- P1 – GUI “Speed ↔ Sensitivity” toggle
- P1 – Automatic contig assembly step in the main pipeline (v1.1)

---

## Contributing / Support

Bug reports & feature requests are welcome! Please open an issue first to discuss major changes.

## License & Citation

*Code* is released under the MIT License (see `LICENSES/MIT.txt`).
The GUI bundles PySide6 / Qt 6, licensed under LGPL v3 (see `LICENSES/Qt_LGPLv3.txt`).

If you use MicroSeq in a publication, please cite me! =) 

> ---



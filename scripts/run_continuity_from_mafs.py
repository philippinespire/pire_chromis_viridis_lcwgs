#!/usr/bin/env python3
"""Build continuity-formatted read counts from BAM + ANGSD MAFS and run continuity LRT.

# Malin Pinsky 2026
# This script was written with AI assistance using GitHub Copilot
# with the GPT-5.3-Codex model.

Overview
--------
This script automates a full continuity-model workflow around
``scripts/ancient_genotypes.py``:

1. Read sample metadata from a CSV samplesheet. Assumes all individuals are in different eras from a single population.
2. Select a BAM subset from either ``--bam-dir``, ``--bam-list``, or the
    samplesheet BAM column.
3. Use an ANGSD ``.mafs(.gz)`` file as the reference AF source and to polarize
    alleles into ancestral vs derived.
4. Generate a continuity-format reads table:
    ``Chrom Pos AF sample_der sample_anc sample_other ...``
5. Generate a matching Eigenstrat-style ``.ind`` file.
6. Run continuity model fitting twice (``continuity=False`` and
    ``continuity=True``), then compute:

    ``LRT = 2 * (LL_true - LL_false)``

7. Compute chi-squared p-values with 1 degree of freedom under the null.

Required Inputs
---------------
- ``--samplesheet``: CSV with at least sample/bam/era columns.
- ``--mafs``: ANGSD ``.mafs`` or ``.mafs.gz``.
- ``--output-prefix``: Prefix used for output files.

Runtime Requirements
--------------------
- Python 3 (recommended: 3.8+).
- Python libraries:
    - ``numpy``
    - ``scipy``
    - ``pandas``
    - ``joblib``
    - ``matplotlib`` (imported by ``ancient_genotypes.py``)
- External command-line tools:
        - ``crun.samtools samtools`` (required in this environment; used for
            ``mpileup``) or a direct ``samtools`` binary if configured. The script has code to try 
            to find a samtools executable in the PATH or load it via a module system if necessary.
- Input data requirements:
    - Coordinate-compatible BAM files with indexes available to ``samtools``.
    - ANGSD ``.mafs(.gz)`` with the expected columns (at minimum:
        ``chromo position major minor knownEM``; and ideally ``anc`` or ``ref`` for
        polarization).

Typical Usage
-------------
Minimal run using BAM directory discovery:

``python3 scripts/run_continuity_from_mafs.py \
  --samplesheet nf-pipelines/nf-angsd-diversity-cvi-only/inputfiles/samplesheet.csv \
  --mafs nf-pipelines/nf-angsd-diversity-cvi-only/results/GL/Cvi.mafs.gz \
  --bam-dir nf-pipelines/nf-trim-merged-unmerged/results/data/bam \
  --output-prefix output/continuity/cvi``

Run using an explicit BAM list and custom filtering:

``python3 scripts/run_continuity_from_mafs.py \
  --samplesheet nf-pipelines/nf-angsd-diversity-cvi-only/inputfiles/samplesheet.csv \
  --mafs nf-pipelines/nf-angsd-diversity-cvi-only/results/GL/Cvi.mafs.gz \
  --bam-list output/continuity/historic_subset.bamlist.txt \
  --output-prefix output/continuity/cvi_subset \
  --site-cutoff 0.05 \
  --min-cutoff 2.5 \
  --max-cutoff 97.5 \
  --num-core 4``

Outputs
-------
- ``<output-prefix>.continuity.reads.tsv``
- ``<output-prefix>.continuity.ind``
- ``<output-prefix>.continuity_lrt.tsv``
- Optional: ``<output-prefix>.continuity.positions.tsv`` (with ``--keep-temp``)

Argument Quick Reference
------------------------
- ``--samplesheet`` (required): CSV containing sample metadata.
- ``--mafs`` (required): ANGSD ``.mafs`` or ``.mafs.gz`` used for AF/polarization.
- ``--output-prefix`` (required): Output prefix for all generated files.
- ``--bam-dir``: Directory containing BAMs to include (non-recursive).
- ``--bam-list``: File with one BAM path per line (overrides ``--bam-dir``).
- ``--sample-column``: Samplesheet column name for sample IDs (default: ``sample``).
- ``--bam-column``: Samplesheet column name with BAM paths (default: ``bam``).
- ``--pop-column``: Samplesheet column name for population labels (default: ``pop``).
- ``--era-column``: Samplesheet column defining era values (default: ``era``).
- ``--historic-era``: Era value treated as historical/ancient (default: ``historic``).
- ``--modern-era``: Era value treated as modern reference (default: ``modern``).
- ``--site-cutoff``: Fraction cutoff passed to ``parse_reads_by_pop`` (default: 0.0).
- ``--min-cutoff``: Lower percentile for ``coverage_filter`` (default: 2.5).
- ``--max-cutoff``: Upper percentile for ``coverage_filter`` (default: 97.5).
- ``--num-core``: Number of cores passed to optimization (default: 1).
- ``--ancient-genotypes-path``: Path to continuity source file
    (default: ``scripts/ancient_genotypes.py``).
- ``--tmp-dir``: Parent directory for temporary working files.
- ``--keep-temp``: Keep intermediate positions file for reproducibility.
- ``--no-ref-polarization``: Disable ``ref`` fallback when ``anc`` is absent.
- ``--samtools-cmd``: Samtools command prefix (default:
    ``crun.samtools samtools``).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import types
from datetime import datetime
from collections import OrderedDict
import shlex

import numpy as np
from scipy.stats import chi2


def parse_args() -> argparse.Namespace:
    """Define and parse command-line arguments for the full workflow.

    The options are grouped into:
    - core input/output arguments,
    - metadata column mapping,
    - model/filter settings,
    - runtime and temporary-file controls.
    """
    parser = argparse.ArgumentParser(
        description="Run continuity test from BAM files and ANGSD .mafs(.gz) output."
    )
    parser.add_argument("--samplesheet", required=True, help="CSV with sample/bam/pop/era columns")
    parser.add_argument("--mafs", required=True, help="ANGSD .mafs or .mafs.gz file")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files")

    parser.add_argument("--bam-dir", help="Directory containing BAM files to include")
    parser.add_argument("--bam-list", help="Text file with one BAM path per line")

    parser.add_argument("--sample-column", default="sample")
    parser.add_argument("--bam-column", default="bam")
    parser.add_argument("--pop-column", default="pop")
    parser.add_argument("--era-column", default="era")
    parser.add_argument("--historic-era", default="historic")
    parser.add_argument("--modern-era", default="modern")

    parser.add_argument("--site-cutoff", type=float, default=0.0, help="cutoff passed to parse_reads_by_pop")
    parser.add_argument("--min-cutoff", type=float, default=2.5, help="coverage_filter min percentile")
    parser.add_argument("--max-cutoff", type=float, default=97.5, help="coverage_filter max percentile")
    parser.add_argument("--num-core", type=int, default=1)

    parser.add_argument(
        "--ancient-genotypes-path",
        default=os.path.join("scripts", "ancient_genotypes.py"),
        help="Path to ancient_genotypes.py",
    )
    parser.add_argument(
        "--tmp-dir",
        default=None,
        help="Directory for temporary files (default: system temp)",
    )
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate continuity files")
    parser.add_argument(
        "--no-ref-polarization",
        action="store_true",
        help=(
            "Disable fallback polarization with ref allele when anc column is absent. "
            "By default, ref-based polarization is enabled for .mafs files like GL/Cvi.mafs.gz."
        ),
    )
    parser.add_argument(
        "--samtools-cmd",
        default="crun.samtools samtools",
        help=(
            "Command used to run samtools (supports wrappers). "
            "Examples: 'crun.samtools samtools', 'samtools'."
        ),
    )
    return parser.parse_args()


def open_maybe_gzip(path: str):
    """Open plain text or gzip-compressed text transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def normalize_path(path: str) -> str:
    """Return a canonical absolute path for robust file matching."""
    return os.path.realpath(os.path.abspath(path))


def read_samplesheet(path: str) -> list[dict[str, str]]:
    """Read the samplesheet CSV into a list of row dictionaries."""
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return [row for row in reader]


def read_bam_list(path: str) -> list[str]:
    """Read one-BAM-per-line input, skipping empty/comment lines."""
    bams = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if s.endswith(".bam"):
                bams.append(normalize_path(s))
    return bams


def read_bam_dir(path: str) -> list[str]:
    """Collect ``*.bam`` files from a directory (non-recursive)."""
    bams = []
    for name in sorted(os.listdir(path)):
        if name.endswith(".bam"):
            bams.append(normalize_path(os.path.join(path, name)))
    return bams


def select_bams(args: argparse.Namespace, rows: list[dict[str, str]]) -> list[str]:
    """Select BAM files from user-priority sources.

    Priority order:
    1. ``--bam-list``
    2. ``--bam-dir``
    3. BAM paths embedded in samplesheet column
    """
    if args.bam_list:
        return read_bam_list(args.bam_list)
    if args.bam_dir:
        return read_bam_dir(args.bam_dir)

    if args.bam_column not in rows[0]:
        raise ValueError("No BAM input provided and bam column not found in samplesheet")
    bams = []
    for row in rows:
        bam = row.get(args.bam_column, "").strip()
        if bam.endswith(".bam"):
            bams.append(normalize_path(bam))
    return sorted(set(bams))


def build_row_lookup(rows: list[dict[str, str]], bam_col: str) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    """Build lookup dictionaries keyed by full BAM path and basename.

    Basename fallback is useful when absolute roots differ between execution
    environments but filenames remain stable.
    """
    by_path = {}
    by_basename = {}
    for row in rows:
        bam = row.get(bam_col, "").strip()
        if not bam:
            continue
        full = normalize_path(bam)
        by_path[full] = row
        by_basename[os.path.basename(full)] = row
    return by_path, by_basename


def select_historic_rows(
    selected_bams: list[str],
    rows: list[dict[str, str]],
    args: argparse.Namespace,
) -> list[dict[str, str]]:
    """Resolve selected BAMs to samplesheet rows and retain historic samples.

    A warning is emitted for BAMs not found in the samplesheet. This avoids
    hard failure when a user-provided list includes stale or external files.
    """
    if not rows:
        raise ValueError("Samplesheet is empty")

    by_path, by_basename = build_row_lookup(rows, args.bam_column)
    selected_rows = []
    missing = []

    for bam in selected_bams:
        key = normalize_path(bam)
        row = by_path.get(key)
        if row is None:
            row = by_basename.get(os.path.basename(key))
        if row is None:
            missing.append(bam)
            continue
        selected_rows.append(row)

    if missing:
        sys.stderr.write(
            "WARNING: {} BAM(s) not found in samplesheet, skipping.\n".format(len(missing))
        )

    era_col = args.era_column
    # Count historic and modern samples among the validated selection
    n_historic = sum(1 for r in selected_rows if r.get(era_col, "").strip().lower() == args.historic_era.lower())
    n_modern = sum(1 for r in selected_rows if r.get(era_col, "").strip().lower() == args.modern_era.lower())

    if n_historic == 0:
        raise ValueError("Error: No historic-era samples found. Both historic and modern eras must be present.")
    if n_modern == 0:
        raise ValueError("Error: No modern-era samples found. Both historic and modern eras must be present.")

    historic_rows = [
        r
        for r in selected_rows
        if r.get(era_col, "").strip().lower() == args.historic_era.lower()
    ]

    return historic_rows


def assert_modern_reference_exists(rows: list[dict[str, str]], args: argparse.Namespace) -> None:
    """Validate that at least one modern-era row exists in metadata.

    The MAFS file supplies reference AF values, but this check enforces the
    intended study design requested for era-based modern/historic grouping.
    """
    era_col = args.era_column
    n_modern = sum(
        1
        for r in rows
        if r.get(era_col, "").strip().lower() == args.modern_era.lower()
    )
    if n_modern == 0:
        raise ValueError(
            "No modern-era rows found in samplesheet for reference population metadata"
        )


def parse_mafs_sites(mafs_path: str, allow_ref_polarization: bool) -> OrderedDict[tuple[str, int], tuple[float, str, str]]:
    """Parse polarized reference AFs from ANGSD MAFS into an ordered site map.

    Returns an OrderedDict keyed by ``(chrom, pos)`` with values:
    ``(af, ancestral_base, derived_base)``.

    Polarization rules:
    - Prefer ``anc`` column when present and compatible with major/minor.
    - Optionally fall back to ``ref`` when ``anc`` is absent/ambiguous.
    - Drop sites that cannot be polarized or are AF 0/1 after polarization.
    """
    sites = OrderedDict()
    dropped_no_polarization = 0
    dropped_monomorphic = 0

    with open_maybe_gzip(mafs_path) as handle:
        header = handle.readline().strip().split()
        index = {name: i for i, name in enumerate(header)}

        required = ["chromo", "position", "major", "minor", "knownEM"]
        for col in required:
            if col not in index:
                raise ValueError("MAFS file missing required column: {}".format(col))

        has_anc = "anc" in index
        has_ref = "ref" in index

        for line in handle:
            if not line.strip():
                continue
            parts = line.strip().split()
            chrom = parts[index["chromo"]]
            pos = int(parts[index["position"]])
            major = parts[index["major"]].upper()
            minor = parts[index["minor"]].upper()
            known_em = float(parts[index["knownEM"]])

            ancestral = None
            if has_anc:
                anc = parts[index["anc"]].upper()
                if anc in (major, minor):
                    ancestral = anc
            if ancestral is None and allow_ref_polarization and has_ref:
                ref = parts[index["ref"]].upper()
                if ref in (major, minor):
                    ancestral = ref

            if ancestral is None:
                dropped_no_polarization += 1
                continue

            derived = minor if ancestral == major else major
            af = known_em if derived == minor else (1.0 - known_em)

            if af <= 0.0 or af >= 1.0:
                dropped_monomorphic += 1
                continue

            sites[(chrom, pos)] = (af, ancestral, derived)

    if not sites:
        raise ValueError("No usable polarized polymorphic sites found in MAFS")

    sys.stderr.write(
        "MAFS sites kept: {} (dropped unpolarized: {}, dropped AF=0/1: {})\n".format(
            len(sites), dropped_no_polarization, dropped_monomorphic
        )
    )
    return sites


def write_positions_file(path: str, sites: OrderedDict[tuple[str, int], tuple[float, str, str]]) -> None:
    """Write a ``chrom\tpos`` site list for ``samtools mpileup -l``."""
    with open(path, "w", encoding="utf-8") as handle:
        for chrom, pos in sites.keys():
            handle.write("{}\t{}\n".format(chrom, pos))


def parse_pileup_bases(bases: str, ref_base: str) -> dict[str, int]:
    """Decode mpileup base-string syntax into per-base read counts.

    Handles start/end read markers, insertions/deletions, reference matches,
    and ambiguous symbols using standard pileup parsing conventions.
    """
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
    i = 0
    ref_base = ref_base.upper()

    while i < len(bases):
        c = bases[i]

        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            i += 1
            num_start = i
            while i < len(bases) and bases[i].isdigit():
                i += 1
            indel_len = int(bases[num_start:i]) if i > num_start else 0
            i += indel_len
            continue
        if c in ".,":
            if ref_base in counts:
                counts[ref_base] += 1
            else:
                counts["N"] += 1
            i += 1
            continue
        if c == "*":
            i += 1
            continue

        base = c.upper()
        if base in counts:
            counts[base] += 1
        else:
            counts["N"] += 1
        i += 1

    return counts


def build_continuity_read_table(
    read_path: str,
    positions_path: str,
    sites: OrderedDict[tuple[str, int], tuple[float, str, str]],
    historic_rows: list[dict[str, str]],
    args: argparse.Namespace,
) -> int:
    """Create continuity-formatted read count table from BAM pileups.

    Output columns are:
    ``Chrom Pos AF`` followed by, for each historic individual,
    ``sample_der sample_anc sample_other``.

    The ``AF`` comes from MAFS-derived polarized reference frequencies, while
    read counts are extracted from selected BAM files at the same positions.
    """
    samtools_cmd = args.samtools_cmd.strip()
    if not samtools_cmd:
        raise RuntimeError("--samtools-cmd produced an empty command")
    samtools_parts = shlex.split(samtools_cmd)

    sample_col = args.sample_column
    bam_col = args.bam_column
    samples = [r[sample_col] for r in historic_rows]
    bams = [normalize_path(r[bam_col]) for r in historic_rows]

    def make_mpileup_shell_cmd(prefix: str) -> str:
        return (
            prefix
            + " mpileup -aa -l "
            + shlex.quote(positions_path)
            + " "
            + " ".join(shlex.quote(bam) for bam in bams)
        )

    subprocess_env = os.environ.copy()
    current_path = subprocess_env.get("PATH", "")
    required_path_entries = ["/usr/bin", "/bin"]
    path_parts = [p for p in current_path.split(":") if p]
    for entry in reversed(required_path_entries):
        if entry not in path_parts:
            path_parts.insert(0, entry)
    subprocess_env["PATH"] = ":".join(path_parts)

    def start_direct_cmd(executable: str):
        return subprocess.Popen(
            [executable, "mpileup", "-aa", "-l", positions_path] + bams,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=subprocess_env,
        )

    shell_executable = "/bin/bash"

    def module_bootstrap_prefix() -> str:
        """Return a shell snippet that initializes the module command and loads samtools."""
        if os.path.exists("/etc/profile.d/modules.sh"):
            module_init = "source /etc/profile.d/modules.sh >/dev/null 2>&1"
        elif os.path.exists("/usr/share/lmod/lmod/init/bash"):
            module_init = "source /usr/share/lmod/lmod/init/bash >/dev/null 2>&1"
        else:
            module_init = ""
        module_load = "module load samtools/1.19 >/dev/null 2>&1"
        if module_init:
            return module_init + " && " + module_load
        return module_load

    def start_shell_cmd(command_text: str):
        return subprocess.Popen(
            [shell_executable, "-lc", command_text], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=subprocess_env
        )

    # Prefer a directly executable samtools path when available. This avoids
    # shell wrapper resolution issues in nested module/container contexts.
    direct_samtools = None
    if len(samtools_parts) == 1 and os.path.isabs(samtools_parts[0]) and os.access(samtools_parts[0], os.X_OK):
        direct_samtools = samtools_parts[0]
    elif samtools_cmd == "samtools":
        direct_samtools = shutil.which("samtools")
    used_cmd_label = samtools_cmd
    try:
        if direct_samtools is not None:
            if samtools_cmd != "samtools":
                sys.stderr.write(
                    "WARNING: using direct samtools executable from PATH instead of wrapper '{}'.\n".format(
                        samtools_cmd
                    )
                )
            proc = start_direct_cmd(direct_samtools)
            used_cmd_label = direct_samtools
        else:
            proc = start_shell_cmd(make_mpileup_shell_cmd(samtools_cmd))
            used_cmd_label = samtools_cmd
    except FileNotFoundError as exc:
        raise RuntimeError(
            "Failed to start samtools command '{}': {}".format(used_cmd_label, exc)
        ) from exc

    written = 0
    site_iter = iter(sites.items())
    current = next(site_iter, None)

    with open(read_path, "w", encoding="utf-8") as out:
        header = ["Chrom", "Pos", "AF"]
        for sample in samples:
            header.extend(["{}_der".format(sample), "{}_anc".format(sample), "{}_other".format(sample)])
        out.write("\t".join(header) + "\n")

        for line in proc.stdout:
            if current is None:
                break

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue

            chrom = fields[0]
            pos = int(fields[1])

            # Catch up if mpileup skips expected sites (rare with -aa + -l).
            # Missing sites are explicitly written with all-zero counts so that
            # downstream arrays remain aligned to the reference-site order.
            while current is not None and (chrom, pos) > current[0]:
                (miss_chrom, miss_pos), (af, ancestral, derived) = current
                row = [miss_chrom, str(miss_pos), "{:.8f}".format(af)]
                for _ in samples:
                    row.extend(["0", "0", "0"])
                out.write("\t".join(row) + "\n")
                written += 1
                current = next(site_iter, None)

            if current is None:
                break

            if (chrom, pos) != current[0]:
                continue

            # current[0] is (chrom, pos) and current[1] is the tuple carrying
            # AF + polarization state used to derive count columns.
            (_, _), (af, ancestral, derived) = current
            ref_base = fields[2].upper()

            expected_cols = 3 + 3 * len(samples)
            if len(fields) < expected_cols:
                raise RuntimeError(
                    "Unexpected mpileup columns at {}:{} (got {}, expected at least {})".format(
                        chrom, pos, len(fields), expected_cols
                    )
                )

            row = [chrom, str(pos), "{:.8f}".format(af)]
            # mpileup per sample emits triplets after the first 3 columns:
            # depth, base-string, base-qualities.
            for i in range(len(samples)):
                bases = fields[4 + 3 * i]
                counts = parse_pileup_bases(bases, ref_base)
                der = counts.get(derived, 0)
                anc = counts.get(ancestral, 0)
                other = max(sum(counts.values()) - der - anc, 0)
                row.extend([str(der), str(anc), str(other)])

            out.write("\t".join(row) + "\n")
            written += 1
            current = next(site_iter, None)

        # Flush any remaining sites with zero coverage.
        while current is not None:
            (miss_chrom, miss_pos), (af, ancestral, derived) = current
            row = [miss_chrom, str(miss_pos), "{:.8f}".format(af)]
            for _ in samples:
                row.extend(["0", "0", "0"])
            out.write("\t".join(row) + "\n")
            written += 1
            current = next(site_iter, None)

    stderr_text = proc.stderr.read()
    ret = proc.wait()
    if ret != 0:
        fallback_needed = direct_samtools is None and "command not found" in stderr_text
        if fallback_needed:
            sys.stderr.write(
                "WARNING: shell wrapper '{}' was unavailable inside subprocess; retrying after loading the samtools module in the child shell.\n".format(
                    samtools_cmd
                )
            )
            fallback_cmd = module_bootstrap_prefix() + " && " + make_mpileup_shell_cmd(samtools_cmd)
            proc = start_shell_cmd(fallback_cmd)
            stderr_text = proc.stderr.read()
            ret = proc.wait()
            used_cmd_label = "module bootstrap && {}".format(samtools_cmd)

        if ret != 0:
            raise RuntimeError("samtools mpileup failed: {}\ncommand: {}".format(stderr_text.strip(), used_cmd_label))
    if stderr_text.strip():
        sys.stderr.write(stderr_text)

    return written


def write_ind_file(path: str, historic_rows: list[dict[str, str]], args: argparse.Namespace) -> None:
    """Write Eigenstrat-style ``.ind`` file required by parse_reads_by_pop(). Assigns everyone to 'all_pop'.
    """
    with open(path, "w", encoding="utf-8") as out:
        for row in historic_rows:
            sample = row[args.sample_column]
            pop = "all_pop"
            out.write("{}\tU\t{}\n".format(sample, pop))


def load_ancient_genotypes_module(path: str):
    """Load ``ancient_genotypes.py`` at runtime.

    The continuity source in this repo is Python-2-style. We first attempt a
    direct compile; if that fails, we attempt Python-2-to-3 translation using:
    1) ``lib2to3`` module (when available), or
    2) external ``2to3`` CLI as a fallback.
    """

    def convert_py2_to_py3_source(src_text: str, source_label: str) -> str:
        """Convert Python-2-style source code to Python 3 source text.

        We prefer in-process conversion via ``lib2to3`` if available. On systems
        where that module is missing, we fall back to invoking the ``2to3``
        executable on a temporary file.
        """
        try:
            from lib2to3.refactor import RefactoringTool, get_fixers_from_package

            fixers = get_fixers_from_package("lib2to3.fixes")
            refactor_tool = RefactoringTool(fixers)
            tree = refactor_tool.refactor_string(src_text, source_label)
            return str(tree)
        except ModuleNotFoundError:
            pass

        two_to_three = shutil.which("2to3")
        if two_to_three is None:
            raise RuntimeError(
                "Could not convert ancient_genotypes.py to Python 3: "
                "lib2to3 module is unavailable and 2to3 executable was not found in PATH."
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            temp_src_path = os.path.join(tmpdir, "ancient_genotypes_tmp.py")
            with open(temp_src_path, "w", encoding="utf-8") as temp_src:
                temp_src.write(src_text)

            proc = subprocess.run(
                [two_to_three, "-w", "-n", temp_src_path],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if proc.returncode != 0:
                raise RuntimeError(
                    "2to3 conversion failed for {}: {}".format(source_label, proc.stderr.strip())
                )

            with open(temp_src_path, "r", encoding="utf-8") as temp_src:
                return temp_src.read()

    source_path = normalize_path(path)
    with open(source_path, "r", encoding="utf-8") as handle:
        src = handle.read()

    module = types.ModuleType("ancient_genotypes_runtime")
    module.__file__ = source_path

    # Some historical source files mix tabs/spaces; normalize tabs for parsing.
    src_candidates = [src]
    if "\t" in src:
        src_candidates.append(src.expandtabs(8))

    last_error = None
    code = None
    for candidate in src_candidates:
        try:
            code = compile(candidate, source_path, "exec")
            break
        except SyntaxError as exc:
            last_error = exc
            try:
                converted = convert_py2_to_py3_source(candidate, source_path)
                code = compile(converted, source_path, "exec")
                break
            except Exception as conv_exc:
                last_error = conv_exc

    if code is None:
        raise RuntimeError(
            "Failed to load {} as Python 3 code. Original/converted parsing failed.".format(source_path)
        ) from last_error

    exec(code, module.__dict__)
    return module


def run_continuity(read_path: str, ind_path: str, args: argparse.Namespace, ag_module):
    """Run continuity optimization pipeline and return model summary arrays.

    Steps:
    1. parse_reads_by_pop
    2. coverage_filter (in place)
    3. optimize with continuity=False
    4. optimize with continuity=True
    5. compute requested LRT and chi-squared p-values
    """
    unique_pops, inds, label, pops, freqs, read_lists = ag_module.parse_reads_by_pop(
        read_path,
        ind_path,
        cutoff=args.site_cutoff,
    )
    print(f"[{datetime.now()}] Parsed reads by population: {len(unique_pops)} unique populations")

    ag_module.coverage_filter(
        read_lists,
        min_cutoff=args.min_cutoff,
        max_cutoff=args.max_cutoff,
    )
    print(f"[{datetime.now()}] Applied coverage filter: min={args.min_cutoff}, max={args.max_cutoff}")

    opts_false = ag_module.optimize_pop_params_error_parallel(
        freqs,
        read_lists,
        num_core=args.num_core,
        continuity=False,
    )
    print(f"[{datetime.now()}] Finished optimization with continuity=False")

    opts_true = ag_module.optimize_pop_params_error_parallel(
        freqs,
        read_lists,
        num_core=args.num_core,
        continuity=True,
    )
    print(f"[{datetime.now()}] Finished optimization with continuity=True")

    ll_false = np.array([-x[1] for x in opts_false], dtype=float)
    ll_true = np.array([-x[1] for x in opts_true], dtype=float)

    # Requested statistic definition from project discussion.
    lrt_requested = 2.0 * (ll_true - ll_false)
    # Numerical guard for p-value computation: chi2.sf expects non-negative test
    # statistics. Negative values can occur from optimizer noise.
    lrt_for_p = np.maximum(lrt_requested, 0.0)
    pvals = chi2.sf(lrt_for_p, 1)
    print(f"[{datetime.now()}] Computed LRT and p-values for {len(unique_pops)} populations")

    return unique_pops, opts_false, opts_true, ll_false, ll_true, lrt_requested, pvals


def write_results(
    result_path: str,
    unique_pops,
    opts_false,
    opts_true,
    ll_false,
    ll_true,
    lrt_requested,
    pvals,
) -> None:
    """Write per-population likelihoods, LRT, p-values, and fitted parameters."""
    with open(result_path, "w", encoding="utf-8") as out:
        out.write(
            "population\tll_cont_false\tll_cont_true\tlrt_2x_ll_true_minus_false\tpvalue_chi2_df1\tparams_cont_false\tparams_cont_true\n"
        )
        for i, pop in enumerate(unique_pops):
            out.write(
                "{}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10e}\t{}\t{}\n".format(
                    str(pop),
                    ll_false[i],
                    ll_true[i],
                    lrt_requested[i],
                    pvals[i],
                    ",".join("{:.10g}".format(v) for v in opts_false[i][0]),
                    ",".join("{:.10g}".format(v) for v in opts_true[i][0]),
                )
            )


def main() -> int:
    """Orchestrate end-to-end execution from arguments to output files."""
    args = parse_args()

    print(f"[{datetime.now()}] Starting to read samplesheet")
    rows = read_samplesheet(args.samplesheet)
    if not rows:
        raise ValueError("Empty samplesheet")
    print(f"[{datetime.now()}] Finished reading samplesheet")

    
    print(f"[{datetime.now()}] Starting to assert modern reference exists")
    assert_modern_reference_exists(rows, args)
    print(f"[{datetime.now()}] Finished asserting modern reference exists")
    
    print(f"[{datetime.now()}] Starting to select BAMs")
    selected_bams = select_bams(args, rows)
    if not selected_bams:
        raise ValueError("No BAMs selected")
    print(f"[{datetime.now()}] Finished selecting BAMs")

    print(f"[{datetime.now()}] Starting to select historic rows and parsing MAFS sites")
    historic_rows = select_historic_rows(selected_bams, rows, args)
    sites = parse_mafs_sites(args.mafs, allow_ref_polarization=(not args.no_ref_polarization))
    print(f"[{datetime.now()}] Finished selecting historic rows and parsing MAFS sites")

    # Ensure parent directory for output prefix exists.
    output_dir = os.path.dirname(normalize_path(args.output_prefix)) or os.getcwd()
    os.makedirs(output_dir, exist_ok=True)
    print(f"[{datetime.now()}] Finished making output directory {output_dir}")

    # Temporary directory holds the mpileup positions file and is cleaned unless
    # the user asks to keep derived intermediates.
    temp_dir_obj = tempfile.TemporaryDirectory(dir=args.tmp_dir)
    temp_dir = temp_dir_obj.name
    print(f"[{datetime.now()}] Created temporary directory {temp_dir}")

    positions_path = os.path.join(temp_dir, "continuity_positions.tsv")
    reads_path = args.output_prefix + ".continuity.reads.tsv"
    ind_path = args.output_prefix + ".continuity.ind"
    results_path = args.output_prefix + ".continuity_lrt.tsv"

    # Build intermediate and model input files.
    write_positions_file(positions_path, sites)
    print(f"[{datetime.now()}] Wrote positions file {positions_path}")
    n_sites_written = build_continuity_read_table(reads_path, positions_path, sites, historic_rows, args)
    write_ind_file(ind_path, historic_rows, args)
    print(f"[{datetime.now()}] Wrote reads file {reads_path} with {n_sites_written} sites")

    # Load continuity implementation and run both model fits.
    ag_module = load_ancient_genotypes_module(args.ancient_genotypes_path)
    (
        unique_pops,
        opts_false,
        opts_true,
        ll_false,
        ll_true,
        lrt_requested,
        pvals,
    ) = run_continuity(reads_path, ind_path, args, ag_module)
    print(f"[{datetime.now()}] Finished continuity optimization and LRT computation")

    # Persist tabular model outputs for downstream plotting/reporting.
    write_results(
        results_path,
        unique_pops,
        opts_false,
        opts_true,
        ll_false,
        ll_true,
        lrt_requested,
        pvals,
    )

    sys.stderr.write(
        "Wrote {} sites to {}\nWrote ind file {}\nWrote LRT results {}\n".format(
            n_sites_written,
            reads_path,
            ind_path,
            results_path,
        )
    )

    if args.keep_temp:
        # Keep a copy of queried positions to support exact reproducibility.
        kept_positions = args.output_prefix + ".continuity.positions.tsv"
        shutil.copy2(positions_path, kept_positions)
        sys.stderr.write("Kept positions file {}\n".format(kept_positions))

    temp_dir_obj.cleanup()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # broad exception for clearer CLI failure messages
        sys.stderr.write("ERROR: {}\n".format(exc))
        raise

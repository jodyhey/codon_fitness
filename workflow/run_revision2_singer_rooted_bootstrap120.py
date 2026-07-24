#!/usr/bin/env python3
"""Bootstrap codon-pair SFRatios analyses for the revision2 SINGER-rooted data."""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import random
import shutil
import subprocess
import sys
from pathlib import Path


DEFAULT_REVISION2 = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/revision2"
)
DEFAULT_PIPELINE = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/SFRatios_pipeline_7_9_2026/scripts"
)


def positive_int(text: str) -> int:
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return value


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Create replacement bootstrap codon-pair count files from the "
            "SINGER-rooted ZI count table, build n=160 SFS files, run SFRatios, "
            "summarize 2Ns estimates, and run least-squares fits."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-i",
        "--input-counts",
        type=Path,
        default=DEFAULT_REVISION2 / "ZIResults_singer10rooted_codonpair_counts.txt",
        help="Source codon-pair count table; data rows are sampled with replacement.",
    )
    ap.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=DEFAULT_REVISION2 / "bootstrap120_n160_singer10rooted",
        help="Output directory containing b001, b002, ... bootstrap subdirectories.",
    )
    ap.add_argument(
        "--pipeline-scripts",
        type=Path,
        default=DEFAULT_PIPELINE,
        help="Directory containing make_codon_pair_SFS, run_multiple_SFRatios, summarize, and least-squares scripts.",
    )
    ap.add_argument(
        "-B",
        "--bootstraps",
        type=positive_int,
        default=120,
        help="Number of bootstrap data sets to create and analyze.",
    )
    ap.add_argument(
        "-n",
        "--sample-size",
        type=positive_int,
        default=160,
        help="Haploid allele-copy sample size passed to the codon-pair SFS builder.",
    )
    ap.add_argument(
        "-j",
        "--jobs",
        type=positive_int,
        default=12,
        help="Maximum number of bootstrap data sets to process simultaneously.",
    )
    ap.add_argument(
        "--sfratios-jobs",
        type=positive_int,
        default=1,
        help="Parallel workers passed to run_multiple_SFRatios_jobs.py inside each bootstrap.",
    )
    ap.add_argument(
        "--ls-sims",
        type=positive_int,
        default=1000,
        help="Simulation count passed to Leastsquares_2Ns_estimates_with_masking_v2.py.",
    )
    ap.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Base random seed; bootstrap index is added for each replicate.",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Rebuild a bootstrap directory even if the least-squares output already exists.",
    )
    ap.add_argument(
        "--clean",
        action="store_true",
        help="Delete existing bootstrap output directory before starting.",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs and print planned work without creating files or running jobs.",
    )
    return ap.parse_args()


def script_paths(script_dir: Path) -> dict[str, Path]:
    paths = {
        "make_sfs": script_dir / "make_codon_pair_SFS_from_SNP_paired_allele_counts.py",
        "run_sfr": script_dir / "run_multiple_SFRatios_jobs.py",
        "summarize": script_dir / "summarize_multiple_SFRatios_runs.py",
        "least_squares": script_dir / "Leastsquares_2Ns_estimates_with_masking_v2.py",
        "aggregate_bootstrap": script_dir / "bootstrap" / "run_SFRatios_and_LeastSquares_on_bootstrap_samples.py",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required scripts: " + ", ".join(missing))
    return paths


def read_source_rows(path: Path) -> tuple[str, list[str]]:
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        rows = handle.readlines()
    if not header or not rows:
        raise ValueError(f"Input count table is empty or malformed: {path}")
    return header, rows


def run_checked(command: list[str], cwd: Path, log_path: Path) -> None:
    with log_path.open("a", encoding="utf-8") as log:
        log.write("COMMAND: " + " ".join(command) + "\n")
        log.flush()
        subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, check=True)


def build_bootstrap(
    index: int,
    header: str,
    rows: list[str],
    args: argparse.Namespace,
    paths: dict[str, Path],
) -> str:
    name = f"b{index:03d}"
    bdir = args.outdir / name
    counts = bdir / f"{name}_pairs.txt"
    sfs = bdir / f"{name}_SFSs.txt"
    summary = bdir / f"{name}_SFRatios_summary.txt"
    least_squares = bdir / f"{name}_LeastSquares_modelfitting.txt"
    log_path = bdir / f"{name}.pipeline.log"

    if least_squares.is_file() and least_squares.stat().st_size > 0 and not args.force:
        return f"SKIP {name}"

    bdir.mkdir(parents=True, exist_ok=True)
    log_path.write_text("", encoding="utf-8")

    rng = random.Random(args.seed + index)
    with counts.open("w", encoding="utf-8") as out:
        out.write(header)
        out.writelines(rng.choices(rows, k=len(rows)))

    run_checked(
        [
            sys.executable,
            str(paths["make_sfs"]),
            "-i",
            str(counts),
            "-o",
            str(sfs),
            "-n",
            str(args.sample_size),
            "-e",
            str(args.seed + index),
        ],
        bdir,
        log_path,
    )

    run_checked(
        [
            sys.executable,
            str(paths["run_sfr"]),
            "-a",
            str(sfs),
            "-d",
            "fixed2Ns",
            "-f",
            "unfolded",
            "-p",
            name,
            "-r",
            str(bdir),
            "-j",
            str(args.sfratios_jobs),
        ],
        bdir,
        log_path,
    )

    run_checked(
        [
            sys.executable,
            str(paths["summarize"]),
            "-i",
            f"{name}_Synonymous_*_Qratio_fixed2Ns_nc{args.sample_size}_estimates.out",
            "-o",
            str(summary),
        ],
        bdir,
        log_path,
    )

    run_checked(
        [
            sys.executable,
            str(paths["least_squares"]),
            "-a",
            str(summary),
            "-f",
            "summary",
            "-n",
            str(args.ls_sims),
            "-s",
            str(args.seed + index),
            "-o",
            str(least_squares),
        ],
        bdir,
        log_path,
    )
    return f"DONE {name}"


def main() -> int:
    args = parse_args()
    paths = script_paths(args.pipeline_scripts)
    if not args.input_counts.is_file():
        raise FileNotFoundError(args.input_counts)

    print(f"Input counts: {args.input_counts}")
    print(f"Output dir: {args.outdir}")
    print(f"Bootstraps: {args.bootstraps}")
    print(f"Sample size: n={args.sample_size}")
    print(f"Concurrent bootstrap jobs: {args.jobs}")
    print(f"SFRatios jobs per bootstrap: {args.sfratios_jobs}")
    print(f"Least-squares simulations: {args.ls_sims}")

    if args.dry_run:
        print("Dry run only; no files will be written.")
        return 0

    if args.clean and args.outdir.exists():
        shutil.rmtree(args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)

    header, rows = read_source_rows(args.input_counts)
    print(f"Source data rows: {len(rows)}")

    failures: list[tuple[int, BaseException]] = []
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(build_bootstrap, index, header, rows, args, paths): index
            for index in range(1, args.bootstraps + 1)
        }
        for future in cf.as_completed(futures):
            index = futures[future]
            try:
                print(future.result(), flush=True)
            except BaseException as exc:
                failures.append((index, exc))
                print(f"FAIL b{index:03d}: {exc}", flush=True)

    if failures:
        examples = "; ".join(f"b{idx:03d}: {exc}" for idx, exc in failures[:5])
        raise RuntimeError(f"{len(failures)} bootstrap jobs failed: {examples}")

    aggregate = args.outdir / f"checkbootstrap_{args.bootstraps}_datasets.txt"
    run_checked(
        [
            sys.executable,
            str(paths["aggregate_bootstrap"]),
            "-f",
            str(args.outdir),
            "-j",
            str(args.jobs),
            "-o",
            str(aggregate),
        ],
        args.outdir,
        args.outdir / "aggregate_bootstrap.log",
    )
    print(f"Aggregate bootstrap CI file: {aggregate}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

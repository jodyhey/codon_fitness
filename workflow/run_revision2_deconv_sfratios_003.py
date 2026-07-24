#!/usr/bin/env python3
"""Run the 3% polarization-error deconvolution SFRatios analysis for revision2."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


BASE = Path("/mnt/d/genemod/better_dNdS_models/popgen/SFRatios_pipeline_7_9_2026")
DEFAULT_INPUT_DIR = BASE / "ZIResults/singer10rooted_SFRatios_n160"
DEFAULT_OUTPUT_DIR = BASE / "ZIResults/singer10rooted_SFRatios_n160_deconv0.03"
DEFAULT_INPUT_SFS = DEFAULT_INPUT_DIR / "ZIResults_singer10rooted_n160_SFSs.txt"
DEFAULT_LABEL = "ZIResults_singer10rooted_n160_deconv0.03"
DECONV_SCRIPT = BASE / "scripts/others/add_polarization_error_deconvolution_correction_to_SFS.py"
RUN_MULTIPLE = BASE / "scripts/run_multiple_SFRatios_jobs.py"
SUMMARIZE = BASE / "scripts/summarize_multiple_SFRatios_runs.py"
LEAST_SQUARES = BASE / "scripts/Leastsquares_2Ns_estimates_with_masking_v2.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--input-sfs",
        type=Path,
        default=DEFAULT_INPUT_SFS,
        help=f"Unfolded SFS input file. Default: {DEFAULT_INPUT_SFS}",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        default=DEFAULT_LABEL,
        help=f"Output prefix for corrected SFS, SFRatios, and least-squares files. Default: {DEFAULT_LABEL}",
    )
    parser.add_argument(
        "-e",
        "--error-rate",
        type=float,
        default=0.03,
        help="Polarization-error deconvolution rate. Default: 0.03",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=4,
        help="Concurrent SFRatios jobs. Default: 4",
    )
    parser.add_argument(
        "-n",
        "--ls-sims",
        type=int,
        default=1000,
        help="Least-squares masking simulations. Default: 1000",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running them.",
    )
    return parser.parse_args()


def run(cmd: list[str], cwd: Path, log_path: Path, dry_run: bool) -> None:
    print("COMMAND:", " ".join(cmd))
    if dry_run:
        return
    with log_path.open("a") as log:
        log.write("\nCOMMAND: " + " ".join(cmd) + "\n")
        log.flush()
        subprocess.run(cmd, cwd=cwd, check=True, stdout=log, stderr=subprocess.STDOUT)


def main() -> None:
    args = parse_args()
    if not args.dry_run:
        args.outdir.mkdir(parents=True, exist_ok=True)

    log = args.outdir / f"{args.prefix}_pipeline.log"
    corrected_sfs = args.outdir / f"{args.prefix}_SFSs.txt"
    summary = args.outdir / f"{args.prefix}_SFRatios_summary.txt"
    ls_out = args.outdir / f"{args.prefix}_LeastSquares_analysis.txt"

    run(
        [
            "python3",
            str(DECONV_SCRIPT),
            "-i",
            str(args.input_sfs),
            "-o",
            str(corrected_sfs),
            "-e",
            f"{args.error_rate:.4g}",
            "-s",
            "c",
            "-z",
            "1",
        ],
        cwd=args.outdir,
        log_path=log,
        dry_run=args.dry_run,
    )

    run(
        [
            "python3",
            str(RUN_MULTIPLE),
            "-a",
            str(corrected_sfs),
            "-p",
            args.prefix,
            "-f",
            "unfolded",
            "-d",
            "fixed2Ns",
            "-r",
            str(args.outdir),
            "-j",
            str(args.jobs),
        ],
        cwd=args.outdir,
        log_path=log,
        dry_run=args.dry_run,
    )

    run(
        [
            "python3",
            str(SUMMARIZE),
            "-i",
            str(args.outdir / f"{args.prefix}_Synonymous_*_Qratio_fixed2Ns_nc160_estimates.out"),
            "-o",
            str(summary),
        ],
        cwd=args.outdir,
        log_path=log,
        dry_run=args.dry_run,
    )

    run(
        [
            "python3",
            str(LEAST_SQUARES),
            "-a",
            str(summary),
            "-f",
            "summary",
            "-n",
            str(args.ls_sims),
            "-o",
            str(ls_out),
        ],
        cwd=args.outdir,
        log_path=log,
        dry_run=args.dry_run,
    )

    if not args.dry_run:
        with log.open("a") as log_f:
            log_f.write(f"\nDONE: {args.prefix}\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INDIR = Path("/home/jody/work/ARGwork/singer/ZI/singer_out")
CONVERT = SCRIPT_DIR / "convert_to_tskit_compatible.py"
FIX = SCRIPT_DIR / "fix_tskit_backmutation_parents.py"


def discover_prefixes(indir: Path, include_start: bool) -> list[tuple[str, list[int]]]:
    seen: dict[str, set[int]] = defaultdict(set)
    for path in indir.glob("*_nodes_*.txt"):
        stem = path.name[:-4]
        prefix, idx = stem.rsplit("_nodes_", 1)
        if not idx.isdigit():
            continue
        if prefix.endswith("_start") and not include_start:
            continue
        seen[prefix].add(int(idx))

    out = []
    for prefix, idxs in sorted(seen.items()):
        complete = []
        for idx in sorted(idxs):
            if all((indir / f"{prefix}_{kind}_{idx}.txt").exists() for kind in ("nodes", "branches", "muts")):
                complete.append(idx)
        if complete:
            out.append((prefix, complete))
    return out


def run_capture(cmd: list[str]) -> str:
    proc = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return "COMMAND: " + " ".join(cmd) + "\n" + proc.stdout + proc.stderr


def process_one(task: tuple[str, int, Path, Path, Path, bool]) -> tuple[str, str]:
    prefix, idx, indir, raw_dir, fixed_dir, overwrite = task
    input_prefix = str(indir / prefix)
    raw_prefix = raw_dir / prefix
    fixed_path = fixed_dir / f"{prefix}_{idx}.trees"
    raw_path = raw_dir / f"{prefix}_{idx}.trees"
    label = f"{prefix}_{idx}.trees"
    messages = []

    if fixed_path.exists() and not overwrite:
        return label, "SKIP existing fixed tree: " + str(fixed_path) + "\n"

    if not raw_path.exists() or overwrite:
        messages.append(
            run_capture(
                [
                    "python3",
                    str(CONVERT),
                    "-input",
                    input_prefix,
                    "-output",
                    str(raw_prefix),
                    "-start",
                    str(idx),
                    "-end",
                    str(idx + 1),
                ]
            )
        )

    messages.append(run_capture(["python3", str(FIX), str(raw_path), str(fixed_path)]))
    messages.append(f"DONE {label}\n")
    return label, "".join(messages)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Convert ZI SINGER ARG text outputs to raw tskit .trees and write "
            "mutation-parent-repaired .trees. Defaults are for z640."
        )
    )
    ap.add_argument(
        "-i",
        "--input-dir",
        type=Path,
        default=DEFAULT_INDIR,
        help=f"SINGER output directory. Default: {DEFAULT_INDIR}",
    )
    ap.add_argument(
        "-r",
        "--raw-dir-name",
        default="tskit_raw",
        help="Subdirectory for raw converted .trees. Default: tskit_raw",
    )
    ap.add_argument(
        "-f",
        "--fixed-dir-name",
        default="tskit_fixed",
        help="Subdirectory for repaired .trees. Default: tskit_fixed",
    )
    ap.add_argument(
        "-s",
        "--include-start",
        action="store_true",
        help="Also convert *_start ARG files. Default: skip them",
    )
    ap.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Overwrite existing raw/fixed .trees files. Default: skip existing fixed files",
    )
    ap.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=6,
        help="Parallel conversion/repair jobs. Default: 6",
    )
    ap.add_argument(
        "--expected-blocks",
        type=int,
        default=50,
        help="Expected number of non-start ZI blocks. Default: 50",
    )
    ap.add_argument(
        "--expected-replicates",
        type=int,
        default=10,
        help="Expected ARG samples per block. Default: 10",
    )
    args = ap.parse_args()

    indir = args.input_dir
    raw_dir = indir / args.raw_dir_name
    fixed_dir = indir / args.fixed_dir_name
    raw_dir.mkdir(parents=True, exist_ok=True)
    fixed_dir.mkdir(parents=True, exist_ok=True)
    log_path = fixed_dir / "convert_and_fix_tskit_ZI.log"

    prefixes = discover_prefixes(indir, args.include_start)
    n_expected = sum(len(idxs) for _, idxs in prefixes)
    expected_total = args.expected_blocks * args.expected_replicates
    tasks = [
        (prefix, idx, indir, raw_dir, fixed_dir, args.overwrite)
        for prefix, idxs in prefixes
        for idx in idxs
    ]
    completed = 0
    failed = 0

    with log_path.open("w") as log:
        log.write(f"input_dir={indir}\nraw_dir={raw_dir}\nfixed_dir={fixed_dir}\n")
        log.write(
            f"prefixes={len(prefixes)} discovered_trees={n_expected} "
            f"expected_blocks={args.expected_blocks} expected_replicates={args.expected_replicates} "
            f"expected_total={expected_total} include_start={args.include_start} jobs={args.jobs}\n\n"
        )
        if not args.include_start and len(prefixes) != args.expected_blocks:
            log.write(f"WARNING: discovered {len(prefixes)} non-start block prefixes, expected {args.expected_blocks}\n")
            print(f"WARNING: discovered {len(prefixes)} non-start block prefixes, expected {args.expected_blocks}", flush=True)
        if not args.include_start and n_expected != expected_total:
            log.write(f"WARNING: discovered {n_expected} trees, expected {expected_total}\n")
            print(f"WARNING: discovered {n_expected} trees, expected {expected_total}", flush=True)

        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futures = [ex.submit(process_one, task) for task in tasks]
            for fut in as_completed(futures):
                try:
                    label, message = fut.result()
                    completed += 1
                    log.write(message + "\n")
                    print(f"[{completed}/{len(tasks)}] {label}", flush=True)
                except Exception as exc:
                    failed += 1
                    log.write(f"FAILED: {exc}\n")
                    print(f"FAILED: {exc}", flush=True)
                log.flush()

        log.write(f"\nSUMMARY completed={completed} failed={failed} discovered={n_expected} expected={expected_total}\n")

    print(f"Prefixes discovered: {len(prefixes)}")
    print(f"Trees discovered: {n_expected}")
    print(f"Expected trees: {expected_total}")
    print(f"Tasks completed: {completed}")
    print(f"Tasks failed: {failed}")
    print(f"Log: {log_path}")
    print(f"Raw trees: {raw_dir}")
    print(f"Fixed trees: {fixed_dir}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

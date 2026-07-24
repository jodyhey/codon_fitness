#!/usr/bin/env python3
"""Extract the codon fitness table from the revision2 least-squares output."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


DEFAULT_INPUT = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/SFRatios_pipeline_7_9_2026/"
    "ZIResults/singer10rooted_SFRatios_n160/ZIResults_singer10rooted_n160_LeastSquares_analysis.txt"
)
DEFAULT_OUTPUT = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/revision2/codon_fitnesses_revision2.tsv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, default=DEFAULT_INPUT, help=f"Least-squares output. Default: {DEFAULT_INPUT}")
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT, help=f"Codon fitness TSV output. Default: {DEFAULT_OUTPUT}")
    return parser.parse_args()


def extract_rows(path: Path) -> list[str]:
    rows = ["AA\tCodon\t2Ns"]
    in_table = False
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line == "AA\tCodon\t2Ns":
            in_table = True
            continue
        if not in_table:
            continue
        if not line:
            break
        parts = line.split()
        if len(parts) >= 3 and re.fullmatch(r"[A-Z]", parts[0]) and re.fullmatch(r"[ACGT]{3}", parts[1]):
            rows.append(f"{parts[0]}\t{parts[1]}\t{float(parts[2])}")
    if len(rows) != 60:
        raise ValueError(f"expected 59 synonymous codons plus header, found {len(rows)} rows in {path}")
    return rows


def main() -> None:
    args = parse_args()
    rows = extract_rows(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(rows) + "\n")
    print(f"wrote {args.output} ({len(rows) - 1} codons)")


if __name__ == "__main__":
    main()

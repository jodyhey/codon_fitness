#!/usr/bin/env python3
"""
Plot RNA stem fold change (y) versus estimated 2Ns (x).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.stats import linregress, pearsonr, spearmanr


DEFAULT_INPUT = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "manuscript/MBE/revision/figwork/rna_stem_fold_change_by_2Ns.txt"
)
DEFAULT_OUTPUT = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "manuscript/MBE/revision/figwork/rna_stem_fold_change_by_2Ns.png"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scatter + regression for RNA stem fold change vs 2Ns.")
    parser.add_argument("-i", "--input", default=str(DEFAULT_INPUT), help="Input tab-delimited file.")
    parser.add_argument("-o", "--output", default=str(DEFAULT_OUTPUT), help="Output PNG file.")
    return parser.parse_args()


def read_first_two_columns(path: Path) -> tuple[np.ndarray, np.ndarray]:
    x_vals = []
    y_vals = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)  # header
        for row in reader:
            if len(row) < 2:
                continue
            y_vals.append(float(row[0]))  # first column -> y
            x_vals.append(float(row[1]))  # second column -> x
    if not x_vals:
        raise ValueError(f"No numeric rows found in {path}")
    return np.asarray(x_vals, dtype=float), np.asarray(y_vals, dtype=float)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    x, y = read_first_two_columns(input_path)

    lr = linregress(x, y)
    pearson_r, _ = pearsonr(x, y)
    spearman_rho, spearman_p = spearmanr(x, y)

    x_line = np.linspace(np.min(x), np.max(x), 200)
    y_line = lr.intercept + lr.slope * x_line

    plt.figure(figsize=(8, 8))
    plt.scatter(x, y, s=60, alpha=0.8, color="black", edgecolors="none")
    stats_text = (
        f"Pearson $R^2$ = {pearson_r**2:.4f}\n"
        f"Slope = {lr.slope:.4f}\n"
        f"Spearman $\\rho$ = {spearman_rho:.4f}\n"
        f"Spearman p = {spearman_p:.4g}"
    )

    plt.plot(
        x_line,
        y_line,
        color="#d62728",
        linewidth=2.5,
    )

    plt.xlabel(r"$\hat{\gamma}$", fontsize=18)
    plt.ylabel("RNA stem fold change", fontsize=18)
    plt.tick_params(axis="both", which="both", labelsize=16)
    plt.grid(True, linestyle="-", color="gray", linewidth=0.8, alpha=0.5)
    text_only_handle = Line2D([], [], linestyle="None")
    plt.legend(
        [text_only_handle],
        [stats_text],
        loc="best",
        frameon=True,
        handlelength=0,
        handletextpad=0,
        fontsize=14,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Wrote: {output_path}")


if __name__ == "__main__":
    main()

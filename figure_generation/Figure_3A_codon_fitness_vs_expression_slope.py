#!/usr/bin/env python3
"""Figure 3A."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/"
    "Drosophila_SFS_and_SFRatios/manuscript/MBE/revision"
)
INPUT_DIR = ROOT / "gene_expression_multinomial_log_expression"
OUTPUT_DIR = ROOT / "figwork"


def read_stats(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path, sep="\t")
    return {row["Statistic"]: float(row["Value"]) for _, row in frame.iterrows()}


def make_scatter(
    frame: pd.DataFrame,
    stats: dict[str, float],
    x_column: str,
    y_column: str,
    x_label: str,
    y_label: str,
    output: Path,
    color_by_codon_ending: bool = False,
) -> None:
    x = frame[x_column].to_numpy()
    y = frame[y_column].to_numpy()

    fig, axis = plt.subplots(figsize=(6.0, 6.0), dpi=300)
    if color_by_codon_ending:
        ending = frame["Codon"].str[-1]
        gc_mask = ending.isin(["G", "C"]).to_numpy()
        at_mask = ending.isin(["A", "T"]).to_numpy()
        axis.scatter(
            x[gc_mask],
            y[gc_mask],
            s=38,
            color="#1f77b4",
            alpha=0.78,
            label="G/C-ending codons",
        )
        axis.scatter(
            x[at_mask],
            y[at_mask],
            s=38,
            color="#ff7f0e",
            alpha=0.78,
            label="A/T-ending codons",
        )
        axis.legend(frameon=False, fontsize=10, loc="lower right")
    else:
        axis.scatter(x, y, s=38, color="black", alpha=0.72)

    x_line = np.linspace(x.min(), x.max(), 200)
    axis.plot(
        x_line,
        stats["sma_intercept"] + stats["sma_slope"] * x_line,
        color="#d62728",
        linewidth=2.2,
        label="Standardized major-axis fit",
    )
    axis.axhline(0, color="gray", linewidth=0.8)
    axis.axvline(0, color="gray", linewidth=0.8)
    axis.set_xlabel(x_label, fontsize=13)
    axis.set_ylabel(y_label, fontsize=13)
    axis.grid(True, alpha=0.3)
    axis.text(
        0.03,
        0.97,
        (
            f"Pearson r = {stats['pearson_r']:.3f}\n"
            f"Spearman rho = {stats['spearman_rho']:.3f}\n"
            f"SMA slope = {stats['sma_slope']:.3f}\n"
            f"cluster-bootstrap 95% CI = "
            f"[{stats['sma_slope_ci_low_cluster_bootstrap']:.3f}, "
            f"{stats['sma_slope_ci_high_cluster_bootstrap']:.3f}]"
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> None:
    codons = pd.read_csv(INPUT_DIR / "codon_expression_slopes_and_fitness.tsv", sep="\t")
    changes = pd.read_csv(INPUT_DIR / "directional_one_step_expression_effects.tsv", sep="\t")
    codon_stats = read_stats(INPUT_DIR / "codon_fitness_comparison.tsv")
    change_stats = read_stats(INPUT_DIR / "directional_fitness_comparison.tsv")

    figure_3a = OUTPUT_DIR / "Figure_3A_codon_fitness_vs_log_expression_slope.png"
    make_scatter(
        codons,
        codon_stats,
        "2Ns",
        "expression_slope",
        r"$\hat{g}$",
        "Codon log-expression slope",
        figure_3a,
        color_by_codon_ending=True,
    )
    make_scatter(
        codons,
        codon_stats,
        "2Ns",
        "expression_slope",
        r"$\hat{g}$",
        "Codon log-expression slope",
        OUTPUT_DIR / "Figure_5A_codon_fitness_vs_log_expression_slope.png",
        color_by_codon_ending=True,
    )
    make_scatter(
        changes,
        change_stats,
        "fitness_difference_derived_minus_ancestral",
        "expression_slope_difference_derived_minus_ancestral",
        r"$\hat{\gamma}$",
        "Log-expression slope difference",
        OUTPUT_DIR / "Figure_5B_directional_fitness_vs_log_expression_slope.png",
    )


if __name__ == "__main__":
    main()

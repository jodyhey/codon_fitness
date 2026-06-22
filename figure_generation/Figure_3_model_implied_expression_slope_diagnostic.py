#!/usr/bin/env python3
"""Plot model-implied expression slopes from the expression-scaled MSD model."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


ROOT = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/"
    "Drosophila_SFS_and_SFRatios/manuscript/MBE/revision"
)
SCRIPTS = ROOT / "scripts"
MODEL_DIR = ROOT / "expression_scaled_mutation_selection_model_bootstrap200"
OUTPUT_DIR = ROOT / "figwork"
sys.path.insert(0, str(SCRIPTS))

from fit_expression_scaled_mutation_selection_model import (  # noqa: E402
    DEFAULT_CDS,
    DEFAULT_TABLE,
    DEFAULT_FITNESS,
    build_arrays,
    count_codons,
    parse_mutation_frequencies,
    read_last_transcript_per_gene,
    read_table,
)


def read_expression_parameters() -> tuple[float, float]:
    frame = pd.read_csv(MODEL_DIR / "model_comparison.tsv", sep="\t")
    params = frame.loc[frame["model"] == "expression", "parameters"].iloc[0]
    alpha, beta = [float(x) for x in params.split(",")]
    return alpha, beta


def model_implied_slopes(
    arrays: dict[str, dict[str, np.ndarray | list[str]]],
    predictor: np.ndarray,
    alpha: float,
    beta: float,
    exponent_multiplier: float = 2.0,
) -> pd.DataFrame:
    rows = []
    for aa, payload in arrays.items():
        codons = list(payload["codons"])
        counts = np.asarray(payload["counts"], dtype=float)
        mutation = np.asarray(payload["mutation_frequency"], dtype=float)
        fitness = np.asarray(payload["fitness"], dtype=float)

        keep = np.asarray(payload["keep"], dtype=bool)
        x = predictor[keep]
        weights = counts.sum(axis=1)
        lam = np.exp(alpha + beta * x)
        eta = np.log(mutation)[None, :] + exponent_multiplier * lam[:, None] * fitness[None, :]
        eta -= eta.max(axis=1, keepdims=True)
        probs = np.exp(eta)
        probs /= probs.sum(axis=1, keepdims=True)

        expected_fitness = probs @ fitness
        local_slopes = exponent_multiplier * beta * lam[:, None] * (
            fitness[None, :] - expected_fitness[:, None]
        )
        weighted_slopes = np.average(local_slopes, axis=0, weights=weights)
        for codon, fit, slope in zip(codons, fitness, weighted_slopes):
            rows.append(
                {
                    "AA": aa,
                    "Codon": codon,
                    "2Ns": fit,
                    "model_implied_expression_slope": slope,
                }
            )
    return pd.DataFrame(rows)


def one_step_pairs(codons: pd.DataFrame) -> pd.DataFrame:
    bases = "ACGT"
    row_map = codons.set_index("Codon").to_dict("index")
    rows = []
    for codon, row in row_map.items():
        aa = row["AA"]
        for index, base in enumerate(codon):
            for alt in bases:
                if alt == base:
                    continue
                derived = codon[:index] + alt + codon[index + 1 :]
                if derived not in row_map or row_map[derived]["AA"] != aa:
                    continue
                rows.append(
                    {
                        "AA": aa,
                        "Ancestral": codon,
                        "Derived": derived,
                        "gamma_hat": row_map[derived]["2Ns"] - row["2Ns"],
                        "model_implied_expression_slope_difference": (
                            row_map[derived]["model_implied_expression_slope"]
                            - row["model_implied_expression_slope"]
                        ),
                    }
                )
    return pd.DataFrame(rows)


def sma_stats(frame: pd.DataFrame, x_col: str, y_col: str) -> dict[str, float]:
    x = frame[x_col].to_numpy(dtype=float)
    y = frame[y_col].to_numpy(dtype=float)
    pearson = pearsonr(x, y)
    spearman = spearmanr(x, y)
    slope = np.sign(pearson.statistic) * np.std(y, ddof=1) / np.std(x, ddof=1)
    intercept = float(y.mean() - slope * x.mean())
    return {
        "n": len(frame),
        "pearson_r": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "spearman_rho": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
        "sma_slope": float(slope),
        "sma_intercept": intercept,
    }


def write_stats(path: Path, stats: dict[str, float]) -> None:
    pd.DataFrame(
        [{"Statistic": key, "Value": value} for key, value in stats.items()]
    ).to_csv(path, sep="\t", index=False)


def make_plot(
    frame: pd.DataFrame,
    stats: dict[str, float],
    x_col: str,
    y_col: str,
    x_label: str,
    y_label: str,
    output: Path,
) -> None:
    x = frame[x_col].to_numpy(dtype=float)
    y = frame[y_col].to_numpy(dtype=float)
    fig, axis = plt.subplots(figsize=(6.8, 5.8), dpi=300)
    axis.scatter(x, y, s=38, color="black", alpha=0.72)
    x_line = np.linspace(x.min(), x.max(), 200)
    axis.plot(
        x_line,
        stats["sma_intercept"] + stats["sma_slope"] * x_line,
        color="#d62728",
        linewidth=2.2,
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
            f"SMA slope = {stats['sma_slope']:.3f}"
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
    table, aa_codons = read_table(DEFAULT_TABLE)
    sequences = read_last_transcript_per_gene(DEFAULT_CDS)
    count_map = {}
    for gene in table["FBgn_ID"]:
        counts, _ = count_codons(sequences[gene])
        count_map[gene] = counts

    expression = table["expression"].to_numpy(dtype=float)
    log_expression = np.log(expression + 1.0)
    predictor = (log_expression - log_expression.mean()) / log_expression.std(ddof=0)

    mutation = parse_mutation_frequencies()
    fitness = pd.read_csv(DEFAULT_FITNESS, sep="\t")
    codon_parameters = mutation.merge(fitness, on=["AA", "Codon"], validate="one_to_one")
    arrays = build_arrays(table, aa_codons, count_map, codon_parameters)
    alpha, beta = read_expression_parameters()

    codons = model_implied_slopes(arrays, predictor, alpha, beta)
    changes = one_step_pairs(codons)
    codon_stats = sma_stats(codons, "2Ns", "model_implied_expression_slope")
    change_stats = sma_stats(
        changes,
        "gamma_hat",
        "model_implied_expression_slope_difference",
    )

    codons.to_csv(
        OUTPUT_DIR / "Figure_5A_model_implied_expression_slopes.tsv",
        sep="\t",
        index=False,
    )
    changes.to_csv(
        OUTPUT_DIR / "Figure_5B_model_implied_expression_slope_differences.tsv",
        sep="\t",
        index=False,
    )
    write_stats(OUTPUT_DIR / "Figure_5A_model_implied_expression_slope_stats.tsv", codon_stats)
    write_stats(OUTPUT_DIR / "Figure_5B_model_implied_expression_slope_stats.tsv", change_stats)
    make_plot(
        codons,
        codon_stats,
        "2Ns",
        "model_implied_expression_slope",
        r"$\hat{g}$",
        "Model-implied log-expression slope",
        OUTPUT_DIR / "Figure_5A_model_implied_expression_slope.png",
    )
    make_plot(
        changes,
        change_stats,
        "gamma_hat",
        "model_implied_expression_slope_difference",
        r"$\hat{\gamma}$",
        "Model-implied log-expression slope difference",
        OUTPUT_DIR / "Figure_5B_model_implied_expression_slope.png",
    )


if __name__ == "__main__":
    main()

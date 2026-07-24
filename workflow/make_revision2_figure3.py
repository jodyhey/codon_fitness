#!/usr/bin/env python3
"""Generate revision2 Figure 3 inputs and panels from ARG-rooted codon fitnesses."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


DEFAULT_BUGFIX = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/revision2"
)
DEFAULT_MODEL_SCRIPT = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/revision/SFS_SFRatios_code_archive/"
    "gene_expression_work/fit_expression_scaled_mutation_selection_model.py"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bugfix-root", type=Path, default=DEFAULT_BUGFIX)
    parser.add_argument("--model-script", type=Path, default=DEFAULT_MODEL_SCRIPT)
    parser.add_argument("--selection-exponent-multiplier", type=float, default=2.0)
    return parser.parse_args()


def read_stats(path: Path) -> dict[str, float | str]:
    frame = pd.read_csv(path, sep="\t")
    stats: dict[str, float | str] = {}
    for _, row in frame.iterrows():
        try:
            stats[row["Statistic"]] = float(row["Value"])
        except ValueError:
            stats[row["Statistic"]] = row["Value"]
    return stats


def softmax(values: np.ndarray) -> np.ndarray:
    shifted = values - np.max(values)
    expv = np.exp(shifted)
    return expv / expv.sum()


def load_mutation_frequencies(model_script: Path) -> pd.DataFrame:
    spec = importlib.util.spec_from_file_location("expr_model", model_script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not import {model_script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.parse_mutation_frequencies()


def model_parameters(path: Path) -> tuple[float, float]:
    model = pd.read_csv(path, sep="\t")
    row = model.loc[model["model"] == "expression"].iloc[0]
    alpha, beta = [float(x) for x in str(row["parameters"]).split(",")]
    return alpha, beta


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
        axis.scatter(x[gc_mask], y[gc_mask], s=38, color="#1f77b4", alpha=0.78, label="G/C-ending codons")
        axis.scatter(x[at_mask], y[at_mask], s=38, color="#ff7f0e", alpha=0.78, label="A/T-ending codons")
        axis.legend(frameon=False, fontsize=10, loc="lower right")
    else:
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


def make_observed_table(expr_dir: Path, out_tsv: Path) -> pd.DataFrame:
    codons = pd.read_csv(expr_dir / "codon_expression_slopes_and_fitness.tsv", sep="\t")
    overall = read_stats(expr_dir / "overall_model_comparison.tsv")
    log_sd = overall["log_expression_sd"]
    rows = []
    for aa, group in codons.groupby("AA", sort=True):
        group = group.copy()
        intercept = group["expression_intercept_centered"].to_numpy(dtype=float)
        slope = group["expression_slope"].to_numpy(dtype=float)

        def probs(fold: float) -> np.ndarray:
            dx = np.log(fold) / log_sd
            return softmax(intercept + slope * dx)

        p1 = probs(1.0)
        p10 = probs(10.0)
        p100 = probs(100.0)
        for (_, row), base, high10, high100 in zip(group.iterrows(), p1, p10, p100):
            codon = row["Codon"]
            rows.append(
                {
                    "AA": aa,
                    "Codon": codon,
                    "label": f"{codon} ({aa})",
                    "observed_frequency_baseline": base,
                    "observed_frequency_10fold": high10,
                    "observed_frequency_100fold": high100,
                    "observed_frequency_change_10fold": high10 - base,
                    "observed_frequency_change_100fold": high100 - base,
                    "expression_intercept_centered": row["expression_intercept_centered"],
                    "expression_slope": row["expression_slope"],
                    "2Ns": row["2Ns"],
                    "ending_class": "G/C-ending" if codon[-1] in {"G", "C"} else "A/T-ending",
                }
            )
    plot = pd.DataFrame(rows).sort_values("observed_frequency_change_100fold").reset_index(drop=True)
    plot.to_csv(out_tsv, sep="\t", index=False)
    return plot


def make_predicted_table(model_dir: Path, fitness_path: Path, model_script: Path, out_tsv: Path, multiplier: float) -> pd.DataFrame:
    mutation = load_mutation_frequencies(model_script)
    fitness = pd.read_csv(fitness_path, sep="\t").rename(columns={"2Ns": "fitness_2Ns"})
    frame = mutation.merge(fitness, on=["AA", "Codon"], validate="one_to_one")
    alpha, beta = model_parameters(model_dir / "model_comparison.tsv")
    log_sd = read_stats(model_dir.parent / "gene_expression_multinomial_log_expression" / "overall_model_comparison.tsv")["log_expression_sd"]
    rows = []
    for aa, group in frame.groupby("AA", sort=True):
        group = group.copy()
        mut = group["mutation_frequency"].to_numpy(dtype=float)
        fit = group["fitness_2Ns"].to_numpy(dtype=float)

        def probs(fold: float) -> np.ndarray:
            dx = np.log(fold) / log_sd
            eta = np.log(mut) + multiplier * np.exp(alpha + beta * dx) * fit
            return softmax(eta)

        p1 = probs(1.0)
        p10 = probs(10.0)
        p100 = probs(100.0)
        for (_, row), base, high10, high100 in zip(group.iterrows(), p1, p10, p100):
            codon = row["Codon"]
            rows.append(
                {
                    "AA": aa,
                    "Codon": codon,
                    "label": f"{codon} ({aa})",
                    "predicted_frequency_baseline": base,
                    "predicted_frequency_10fold": high10,
                    "predicted_frequency_100fold": high100,
                    "predicted_frequency_change_10fold": high10 - base,
                    "predicted_frequency_change_100fold": high100 - base,
                    "ending_class": "G/C-ending" if codon[-1] in {"G", "C"} else "A/T-ending",
                }
            )
    plot = pd.DataFrame(rows).sort_values("predicted_frequency_change_100fold").reset_index(drop=True)
    plot.to_csv(out_tsv, sep="\t", index=False)
    return plot


def darken(color: str, amount: float = 0.62) -> tuple[float, float, float]:
    return tuple(np.array(to_rgb(color)) * amount)


def frequency_change_fit_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    r = float(np.corrcoef(observed, predicted)[0, 1])
    slope, intercept = np.polyfit(predicted, observed, 1)
    residual = observed - predicted
    return {
        "pearson_r": r,
        "r_squared": r * r,
        "observed_on_predicted_slope": float(slope),
        "observed_on_predicted_intercept": float(intercept),
        "rmse": float(np.sqrt(np.mean(residual * residual))),
        "mae": float(np.mean(np.abs(residual))),
    }


def write_figure3b_fit_metrics(out_png: Path, metrics: dict[str, dict[str, float]]) -> None:
    out_path = out_png.with_name("Figure_3B_observed_predicted_fit_metrics.tsv")
    columns = [
        "contrast",
        "pearson_r",
        "r_squared",
        "observed_on_predicted_slope",
        "observed_on_predicted_intercept",
        "rmse",
        "mae",
    ]
    with out_path.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for contrast in ["10fold", "100fold"]:
            row = metrics[contrast]
            handle.write("\t".join([contrast] + [f"{row[col]:.12g}" for col in columns[1:]]) + "\n")


def make_figure3b(obs: pd.DataFrame, pred: pd.DataFrame, out_png: Path) -> pd.DataFrame:
    merged = obs[
        [
            "AA",
            "Codon",
            "label",
            "ending_class",
            "2Ns",
            "observed_frequency_change_10fold",
            "observed_frequency_change_100fold",
        ]
    ].merge(
        pred[["AA", "Codon", "predicted_frequency_change_10fold", "predicted_frequency_change_100fold"]],
        on=["AA", "Codon"],
        validate="one_to_one",
    )
    merged = merged.sort_values("2Ns").reset_index(drop=True)
    merged.to_csv(out_png.with_suffix(".tsv"), sep="\t", index=False)

    base_colors = merged["ending_class"].map({"G/C-ending": "#1f77b4", "A/T-ending": "#ff7f0e"}).to_numpy()
    dark_colors = np.array([darken(c) for c in base_colors], dtype=object)
    y = np.arange(len(merged))
    obs10 = merged["observed_frequency_change_10fold"].to_numpy()
    obs100 = merged["observed_frequency_change_100fold"].to_numpy()
    pred10 = merged["predicted_frequency_change_10fold"].to_numpy()
    pred100 = merged["predicted_frequency_change_100fold"].to_numpy()
    ghat = merged["2Ns"].to_numpy()

    fig, ax = plt.subplots(figsize=(8.2, 12.0), dpi=300)
    for i, (color, dark_color) in enumerate(zip(base_colors, dark_colors)):
        ax.plot([obs10[i], pred10[i]], [y[i], y[i]], color=color, alpha=0.75, linewidth=1.6, linestyle=(0, (3.2, 2.4)), zorder=1)
        ax.plot([obs100[i], pred100[i]], [y[i], y[i]], color=dark_color, alpha=0.9, linewidth=2.0, linestyle="-", zorder=2)
    ax.scatter(obs10, y - 0.16, s=54, marker="o", facecolors="white", edgecolors=base_colors, linewidths=1.8, zorder=5)
    ax.scatter(obs100, y - 0.16, s=64, marker="o", c=dark_colors.tolist(), edgecolors="none", alpha=0.96, zorder=6)
    ax.scatter(pred10, y + 0.16, s=58, marker="s", facecolors="white", edgecolors=base_colors, linewidths=1.8, zorder=5)
    ax.scatter(pred100, y + 0.16, s=68, marker="s", c=dark_colors.tolist(), edgecolors="none", alpha=0.96, zorder=6)
    ax.axvline(0, color="0.25", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(merged["label"], fontsize=8.2)
    ax.set_xlabel("Frequency change relative to baseline expression", fontsize=15)
    ax.set_ylabel("Codon, sorted by $\\hat{g}$", fontsize=15)
    ax.tick_params(axis="x", labelsize=12)
    right_limit = max(abs(obs10).max(), abs(obs100).max(), abs(pred10).max(), abs(pred100).max()) * 1.10
    ax.set_xlim(-0.15, right_limit)
    ax.set_ylim(-0.7, len(merged) - 0.3)

    ax_g = ax.twiny()
    gpad = max(abs(ghat.min()), abs(ghat.max())) * 1.08
    ax_g.set_xlim(-gpad, gpad)
    ax_g.axvline(0, color="0.25", linewidth=0.8, alpha=0.45)
    ax_g.scatter(ghat, y, s=38, c="black", alpha=0.82, zorder=2)
    ax_g.set_xlabel(r"Codon fitness, $\hat{g}$", fontsize=15)
    ax_g.set_ylim(ax.get_ylim())
    ax_g.tick_params(axis="x", labelsize=12)
    ax_g.tick_params(axis="y", left=False, right=False, labelleft=False, labelright=False)

    fit_metrics = {
        "10fold": frequency_change_fit_metrics(obs10, pred10),
        "100fold": frequency_change_fit_metrics(obs100, pred100),
    }
    write_figure3b_fit_metrics(out_png, fit_metrics)
    m10 = fit_metrics["10fold"]
    m100 = fit_metrics["100fold"]
    legend_items = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor="black", markeredgewidth=1.6, markersize=8, label="Observed"),
        Line2D([0], [0], marker="s", linestyle="none", markerfacecolor="white", markeredgecolor="black", markeredgewidth=1.6, markersize=8, label="Predicted"),
        Line2D([0], [0], color="black", linewidth=1.8, linestyle=(0, (3.2, 2.4)), label="10-fold"),
        Line2D([0], [0], color="black", linewidth=2.2, linestyle="-", label="100-fold"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="black", markeredgecolor="none", markersize=7.5, label=r"$\hat{g}$"),
        Patch(facecolor="#1f77b4", label="G/C-ending codons"),
        Patch(facecolor="#ff7f0e", label="A/T-ending codons"),
        Line2D([], [], linestyle="none", label=f"10-fold fit:\nr={m10['pearson_r']:.3f}, $R^2$={m10['r_squared']:.3f}\nslope={m10['observed_on_predicted_slope']:.3f}"),
        Line2D([], [], linestyle="none", label=f"100-fold fit:\nr={m100['pearson_r']:.3f}, $R^2$={m100['r_squared']:.3f}\nslope={m100['observed_on_predicted_slope']:.3f}"),
    ]
    ax.legend(
        handles=legend_items,
        frameon=True,
        loc="lower right",
        bbox_to_anchor=(0.995, 0.015),
        borderaxespad=0.15,
        fontsize=12.0,
        ncol=1,
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    return merged


def main() -> None:
    args = parse_args()
    expr_dir = args.bugfix_root / "gene_expression_multinomial_log_expression"
    model_dir = args.bugfix_root / "expression_scaled_mutation_selection_model_bootstrap200"
    figwork = args.bugfix_root / "figwork"
    figwork.mkdir(parents=True, exist_ok=True)

    codons = pd.read_csv(expr_dir / "codon_expression_slopes_and_fitness.tsv", sep="\t")
    changes = pd.read_csv(expr_dir / "directional_one_step_expression_effects.tsv", sep="\t")
    make_scatter(
        codons,
        read_stats(expr_dir / "codon_fitness_comparison.tsv"),
        "2Ns",
        "expression_slope",
        r"$\hat{g}$",
        "Codon log-expression slope",
        figwork / "Figure_3A_codon_fitness_vs_log_expression_slope.png",
        color_by_codon_ending=True,
    )
    make_scatter(
        changes,
        read_stats(expr_dir / "directional_fitness_comparison.tsv"),
        "fitness_difference_derived_minus_ancestral",
        "expression_slope_difference_derived_minus_ancestral",
        r"$\hat{\gamma}$",
        "Log-expression slope difference",
        figwork / "Figure_3_directional_fitness_vs_log_expression_slope.png",
    )

    obs = make_observed_table(expr_dir, figwork / "Figure_3B_observed_frequency_change_expression.tsv")
    pred = make_predicted_table(
        model_dir,
        args.bugfix_root / "codon_fitnesses_revision2.tsv",
        args.model_script,
        figwork / "Figure_3B_codon_predicted_frequency_change_expression.tsv",
        args.selection_exponent_multiplier,
    )
    make_figure3b(obs, pred, figwork / "Figure_3B_observed_predicted_frequency_change_by_ghat_onepanel.png")
    print(figwork)


if __name__ == "__main__":
    main()

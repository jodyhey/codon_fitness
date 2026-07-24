#!/usr/bin/env python3
"""Regenerate revision2 Figure 4 panels from ARG-rooted SFRatios estimates."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress, pearsonr, spearmanr


BUGFIX = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/revision2"
)
REVISION = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/revision"
)
ORIGINAL = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/original_submission"
)
MAIN_SUMMARY = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/SFRatios_pipeline_7_9_2026/"
    "ZIResults/singer10rooted_SFRatios_n160/ZIResults_singer10rooted_n160_SFRatios_summary.txt"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bugfix-root", type=Path, default=BUGFIX)
    parser.add_argument("--revision-root", type=Path, default=REVISION)
    parser.add_argument("--original-root", type=Path, default=ORIGINAL)
    parser.add_argument("--main-summary", type=Path, default=MAIN_SUMMARY)
    return parser.parse_args()


def fmt_p(p: float) -> str:
    if np.isnan(p):
        return "nan"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3f}"


def compute_stats(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    lr = linregress(x, y)
    r, _ = pearsonr(x, y)
    rho, p = spearmanr(x, y)
    return {
        "r2": float(r * r),
        "slope": float(lr.slope),
        "intercept": float(lr.intercept),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
    }


def scatter_plot(
    x: np.ndarray,
    y: np.ndarray,
    xlabel: str,
    ylabel: str,
    output: Path,
    text_loc: str = "upper right",
) -> dict[str, float]:
    stats = compute_stats(x, y)
    xp = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    ax.scatter(x, y, s=60, c="k", alpha=0.85)
    ax.plot(xp, stats["slope"] * xp + stats["intercept"], color="red", linewidth=2.5)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.tick_params(axis="both", which="both", labelsize=16)
    ax.grid(True, which="major", axis="both", linestyle="-", color="gray", linewidth=0.8, alpha=0.5)
    text = (
        f"R$^2$ = {stats['r2']:.3f}\n"
        f"Slope = {stats['slope']:.3f}\n"
        f"Spearman $\\rho$ = {stats['spearman_rho']:.3f}\n"
        f"Spearman p-value = {fmt_p(stats['spearman_p'])}"
    )
    if text_loc == "upper left":
        text_x, text_y, text_ha = 0.02, 0.98, "left"
    else:
        text_x, text_y, text_ha = 0.98, 0.98, "right"
    ax.text(text_x, text_y, text, transform=ax.transAxes, ha=text_ha, va="top", fontsize=14)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)
    return stats


def read_corrected_g(path: Path) -> pd.DataFrame:
    g = pd.read_csv(path, sep="\t")
    return g[["AA", "Codon", "2Ns"]].rename(columns={"2Ns": "g_new"})


def map_old_rows_to_codons(old_data: pd.DataFrame, old_fitness: pd.DataFrame) -> pd.DataFrame:
    old = old_fitness.copy()
    old["g_round"] = old["2Ns"].round(9)
    data = old_data.copy()
    data["g_round"] = data.iloc[:, 1].round(9)
    if old["g_round"].duplicated().any():
        raise ValueError("Old g values are not unique enough to map to codons")
    mapped = data.merge(old[["AA", "Codon", "g_round"]], on="g_round", how="left", validate="many_to_one")
    if mapped["Codon"].isna().any():
        missing = mapped.loc[mapped["Codon"].isna()].head()
        raise ValueError(f"Could not map old g rows to codons:\n{missing}")
    return mapped.drop(columns=["g_round"])


def make_factor_and_stability(args: argparse.Namespace, figwork: Path) -> None:
    old_fitness = pd.read_csv(args.revision_root / "codon_fitnesses.tsv", sep="\t")
    new_g = read_corrected_g(args.bugfix_root / "codon_fitnesses_revision2.tsv")

    factor_old = pd.read_csv(args.original_root / "figwork/factor_1_loadings_by_g.txt", sep="\t")
    factor = map_old_rows_to_codons(factor_old, old_fitness).merge(new_g, on=["AA", "Codon"], validate="one_to_one")
    factor_out = factor[["F1loading", "g_new", "AA", "Codon"]].rename(columns={"g_new": "g"})
    factor_out.to_csv(figwork / "factor_1_loadings_by_g.txt", sep="\t", index=False)
    scatter_plot(
        factor_out["F1loading"].to_numpy(dtype=float),
        factor_out["g"].to_numpy(dtype=float),
        "Factor Analysis F1 Loadings",
        r"$\hat{g}$",
        figwork / "Figure_4A_factor_1_loadings_by_g.png",
    )

    stability_old = pd.read_csv(args.original_root / "figwork/codon_stability_by_g.txt", sep="\t")
    stability = map_old_rows_to_codons(stability_old, old_fitness).merge(new_g, on=["AA", "Codon"], validate="one_to_one")
    stability_out = stability[["CSC Code", "g_new", "AA", "Codon"]].rename(columns={"g_new": "g"})
    stability_out.to_csv(figwork / "codon_stability_by_g.txt", sep="\t", index=False)
    scatter_plot(
        stability_out["g"].to_numpy(dtype=float),
        stability_out["CSC Code"].to_numpy(dtype=float),
        r"$\hat{g}$",
        "Codon Stability Coefficient",
        figwork / "Figure_4B_codon_stability_by_g.png",
    )


def make_rna_stem(args: argparse.Namespace, figwork: Path) -> None:
    old_stem = pd.read_csv(args.revision_root / "figwork/rna_stem_fold_change_by_2Ns.txt", sep="\t")
    workbook = args.bugfix_root / "SupplementaryInformation_revision2.xlsx"
    if workbook.exists():
        new_primary = pd.read_excel(
            workbook,
            sheet_name="Table 3. Initial 2Ns estimates",
        )
        primary_values = new_primary["2Ns est"]
        codon_pairs = new_primary["CodonPair (from_to)"]
    else:
        summary = pd.read_csv(args.main_summary, sep="\t")
        primary_values = summary["2Ns"]
        codon_pairs = summary["filename"].str.extract(r"Synonymous_([ACGT]{3}_[ACGT]{3})_Qratio", expand=False)
    if len(old_stem) != len(primary_values):
        raise ValueError(f"RNA stem rows {len(old_stem)} != primary rows {len(primary_values)}")
    out = pd.DataFrame(
        {
            "stem_fold_change": old_stem["stem_fold_change"],
            "2Ns": primary_values,
            "CodonPair": codon_pairs,
        }
    )
    out.to_csv(figwork / "rna_stem_fold_change_by_2Ns.txt", sep="\t", index=False)
    stats = scatter_plot(
        out["2Ns"].to_numpy(dtype=float),
        out["stem_fold_change"].to_numpy(dtype=float),
        r"$\hat{\gamma}$",
        "RNA stem fold change",
        figwork / "Figure_4C_rna_stem_fold_change_by_2Ns.png",
        text_loc="upper left",
    )


def main() -> None:
    args = parse_args()
    figwork = args.bugfix_root / "figwork"
    figwork.mkdir(parents=True, exist_ok=True)
    make_factor_and_stability(args, figwork)
    make_rna_stem(args, figwork)
    print(figwork)


if __name__ == "__main__":
    main()

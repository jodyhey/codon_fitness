#!/usr/bin/env python3
"""
Generate a square dotplot of y (col2) vs x (col1) from:
  D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\figwork\factor_1_loadings_by_g.txt

- X label: "Factor Analysis F1 Loadings"
- Y label: italic g
- Adds grid lines
- Computes linear regression, Pearson's R (and R^2), Spearman's rho and p-value
- Plots regression line

Saves PNG to:
  D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\figwork\factor_1_loadings_by_gObserved_codon_freq_and_g.png
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy import stats as scipy_stats
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

IN_WIN = r"D:\\genemod\\better_dNdS_models\\popgen\\Drosophila_SFS_and_SFRatios\\manuscript\\MBE\\figwork\\factor_1_loadings_by_g.txt"
OUT_WIN = r"D:\\genemod\\better_dNdS_models\\popgen\\Drosophila_SFS_and_SFRatios\\manuscript\\MBE\\figwork\\factor_1_loadings_by_gObserved_codon_freq_and_g.png"


def windows_to_wsl(path: str) -> str:
    if len(path) >= 3 and path[1:3] == ':\\':
        drive = path[0].lower()
        rest = path[3:].replace('\\', '/')
        return f"/mnt/{drive}/{rest}"
    return path


def resolve_path(win_path: str) -> str:
    if os.path.exists(win_path):
        return win_path
    return windows_to_wsl(win_path)


def read_two_columns(path: str):
    xs, ys = [], []
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            parts = s.replace('\t', ' ').split()
            if len(parts) < 2:
                continue
            try:
                x = float(parts[0])
                y = float(parts[1])
                xs.append(x)
                ys.append(y)
            except ValueError:
                # header or non-numeric line
                continue
    if not xs:
        raise ValueError("No numeric data parsed from input file.")
    return np.asarray(xs), np.asarray(ys)


def compute_stats(x: np.ndarray, y: np.ndarray):
    if _HAVE_SCIPY:
        r_pearson, _ = scipy_stats.pearsonr(x, y)
        lr = scipy_stats.linregress(x, y)
        slope, intercept = lr.slope, lr.intercept
        rho_spear, p_spear = scipy_stats.spearmanr(x, y)
    else:
        r_pearson = float(np.corrcoef(x, y)[0, 1])
        slope, intercept = np.polyfit(x, y, 1)
        # Spearman fallback via ranking; p-value not available
        xr = np.argsort(np.argsort(x))
        yr = np.argsort(np.argsort(y))
        rho_spear = float(np.corrcoef(xr, yr)[0, 1])
        p_spear = float('nan')
    r2 = float(r_pearson * r_pearson)
    return r2, slope, intercept, rho_spear, (p_spear if _HAVE_SCIPY else float('nan'))


def fmt_p(p: float) -> str:
    if np.isnan(p):
        return 'nan'
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3f}"


def main():
    in_path = resolve_path(IN_WIN)
    out_path = resolve_path(OUT_WIN)

    x, y = read_two_columns(in_path)

    r2, slope, intercept, rho_s, p_s = compute_stats(x, y)

    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)

    # Scatter
    ax.scatter(x, y, s=60, c='k', alpha=0.85)

    # Regression line
    xp = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    ax.plot(xp, slope * xp + intercept, color='red', linewidth=2.5)

    # Labels
    ax.set_xlabel("Factor Analysis F1 Loadings", fontsize=18)
    ax.set_ylabel(r"$\hat{g}$", fontsize=18)
    ax.tick_params(axis='both', which='both', labelsize=16)

    # Grid
    ax.grid(True, which='major', axis='both', linestyle='-', color='gray', linewidth=0.8, alpha=0.5)

    # Stats box (upper left)
    text = (
        f"R$^2$ = {r2:.3f}\n"
        f"Slope = {slope:.3f}\n"
        "Spearman " + r"$\rho$" + f" = {rho_s:.3f}\n"
        f"Spearman p-value = {fmt_p(p_s)}"
    )
    ax.text(0.98, 0.98, text, transform=ax.transAxes, ha='right', va='top', fontsize=14)

    # Save
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

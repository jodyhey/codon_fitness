#!/usr/bin/env python3
"""
Generate a square dotplot comparing observed codon RSCU values (x) to two
predicted codon RSCU models (y):
- Column 2: Selection–Mutation–Drift
- Column 3: Mutation–Drift

Reads:
  D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\revision\figwork\predicted_codon_frequencies_RSCU.txt
Writes:
  D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\revision\figwork\predicted_codon_frequencies_RSCU.png

Figure:
- Two scatter series (black/blue for S–M–D, red/orange for M–D)
- Two regression lines (matching colors)
- Legend in lower right identifying datasets and lines
- Stats box in upper left: Pearson's rho, p-value, slope, R^2 for each model
- Axis labels:
  X: "Observed Codon RSCU"
  Y: "Predicted Codon RSCU"
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

IN_WIN = r"D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\revision\figwork\predicted_codon_frequencies_RSCU.txt"
OUT_WIN = r"D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\revision\figwork\predicted_codon_frequencies_RSCU.png"


def windows_to_wsl(path: str) -> str:
    if len(path) >= 3 and path[1:3] == ':\\':
        drive = path[0].lower()
        rest = path[3:].replace('\\', '/')
        return f"/mnt/{drive}/{rest}"
    return path


def resolve_path(win_path: str) -> str:
    if os.path.exists(win_path):
        return win_path
    wsl = windows_to_wsl(win_path)
    return wsl


def read_three_column_file(path: str):
    x_list, y1_list, y2_list = [], [], []
    with open(path, 'r') as f:
        lines = f.readlines()
    # Skip header if detected (non-numeric tokens)
    for idx, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue
        parts = s.replace('\t', ' ').split()
        if len(parts) < 3:
            continue
        try:
            x = float(parts[0])
            y1 = float(parts[1])
            y2 = float(parts[2])
            x_list.append(x)
            y1_list.append(y1)
            y2_list.append(y2)
        except ValueError:
            # Assume header row; skip
            continue
    if not x_list:
        raise ValueError("No numeric data parsed from input file.")
    return np.asarray(x_list), np.asarray(y1_list), np.asarray(y2_list)


essel = {
    'smd': {'color': '#1f77b4', 'label': 'Selection–Mutation–Drift'},
    'md':  {'color': '#ff7f0e', 'label': 'Mutation–Drift'},
}


def pearson_and_fit(x: np.ndarray, y: np.ndarray):
    # Pearson correlation and p-value
    if _HAVE_SCIPY:
        rho, pval = scipy_stats.pearsonr(x, y)
        lr = scipy_stats.linregress(x, y)
        slope, intercept = lr.slope, lr.intercept
    else:
        # Fallbacks: rho from numpy; slope via polyfit; pval set to nan
        rho = float(np.corrcoef(x, y)[0, 1])
        pval = float('nan')
        slope, intercept = np.polyfit(x, y, 1)
    r2 = float(rho * rho)
    return rho, pval, slope, intercept, r2


def main():
    in_path = resolve_path(IN_WIN)
    out_path = resolve_path(OUT_WIN)

    x, y_smd, y_md = read_three_column_file(in_path)

    # Compute stats and fits
    rho_smd, p_smd, m_smd, b_smd, r2_smd = pearson_and_fit(x, y_smd)
    rho_md, p_md, m_md, b_md, r2_md = pearson_and_fit(x, y_md)

    # Figure (square)
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)

    # Scatter
    ax.scatter(x, y_smd, s=60, c=essel['smd']['color'], alpha=0.8, label=essel['smd']['label'])
    ax.scatter(x, y_md,  s=60, c=essel['md']['color'],  alpha=0.8, label=essel['md']['label'])

    # Regression lines over x-range
    xmin, xmax = float(np.min(x)), float(np.max(x))
    xp = np.linspace(xmin, xmax, 200)
    ax.plot(xp, m_smd * xp + b_smd, color=essel['smd']['color'], linewidth=2.5, label=f"{essel['smd']['label']} fit")
    ax.plot(xp, m_md * xp + b_md,   color=essel['md']['color'],  linewidth=2.5, linestyle='--', label=f"{essel['md']['label']} fit")

    # Labels
    ax.set_xlabel("Observed Codon RSCU", fontsize=18)
    ax.set_ylabel("Predicted Codon RSCU", fontsize=18)
    ax.tick_params(axis='both', which='both', labelsize=16)

    # Legend in lower right
    ax.legend(loc='lower right', frameon=False, fontsize=14)

    # Stats box in upper left
    def fmt_p(p):
        if np.isnan(p):
            return 'nan'
        if p < 1e-3:
            return f"{p:.1e}"
        return f"{p:.3f}"

    text = (
        f"Selection–Mutation–Drift\n"
        f"  rho = {rho_smd:.3f}, p = {fmt_p(p_smd)}\n"
        f"  slope = {m_smd:.3f}, R$^2$ = {r2_smd:.3f}\n\n"
        f"Mutation–Drift\n"
        f"  rho = {rho_md:.3f}, p = {fmt_p(p_md)}\n"
        f"  slope = {m_md:.3f}, R$^2$ = {r2_md:.3f}"
    )
    ax.text(0.02, 0.98, text, transform=ax.transAxes, ha='left', va='top', fontsize=14)

    # Grid lines on both axes
    ax.grid(True, which='major', axis='both', linestyle='-', color='gray', linewidth=0.8, alpha=0.5)

    # Make layout tight and save
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

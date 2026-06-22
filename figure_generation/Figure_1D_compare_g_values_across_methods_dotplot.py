#!/usr/bin/env python3
"""
Generate a square dotplot comparing the primary g column (x-axis) to three
other g columns (y-axis) in:
  D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\figwork\compare_g_values_acroos_methods.txt

The first legend shows dataset colors/labels from the header row.
The second legend shows regression slope and R^2 for each comparison.
The third legend ("Correlation") shows Pearson r and R^2 for each comparison.
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

try:
    from scipy import stats as scipy_stats
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


IN_WIN = r"D:\\genemod\\better_dNdS_models\\popgen\\Drosophila_SFS_and_SFRatios\\manuscript\\MBE\\figwork\\compare_g_values_acroos_methods.txt"
OUT_WIN = r"D:\\genemod\\better_dNdS_models\\popgen\\Drosophila_SFS_and_SFRatios\\manuscript\\MBE\\figwork\\compare_g_values_acroos_methods.png"


def windows_to_wsl(path: str) -> str:
    if len(path) >= 3 and path[1:3] == ":\\":  # e.g. D:\...
        drive = path[0].lower()
        rest = path[3:].replace("\\", "/")
        return f"/mnt/{drive}/{rest}"
    return path


def resolve_path(win_path: str) -> str:
    if os.path.exists(win_path):
        return win_path
    return windows_to_wsl(win_path)


def _split_line(line: str, tab_delimited: bool):
    if tab_delimited:
        return [x.strip() for x in line.rstrip("\n").split("\t")]
    return line.strip().split()


def read_four_columns(path: str):
    with open(path, "r", encoding="utf-8") as fh:
        lines = [ln for ln in fh if ln.strip()]
    if not lines:
        raise ValueError("Input file is empty.")

    tab_delimited = "\t" in lines[0]
    headers = _split_line(lines[0], tab_delimited)
    if len(headers) < 4:
        raise ValueError("Expected at least 4 header columns.")
    headers = headers[:4]

    cols = [[], [], [], []]
    for ln in lines[1:]:
        parts = _split_line(ln, tab_delimited)
        if len(parts) < 4:
            continue
        try:
            vals = [float(parts[i]) for i in range(4)]
        except ValueError:
            continue
        for i, v in enumerate(vals):
            cols[i].append(v)

    if len(cols[0]) < 2:
        raise ValueError("Need at least two numeric rows.")

    arrs = [np.asarray(c, dtype=float) for c in cols]
    return headers, arrs


def pearson_and_fit(x: np.ndarray, y: np.ndarray):
    if _HAVE_SCIPY:
        r, p = scipy_stats.pearsonr(x, y)
        lr = scipy_stats.linregress(x, y)
        slope, intercept = lr.slope, lr.intercept
    else:
        r = float(np.corrcoef(x, y)[0, 1])
        p = float("nan")
        slope, intercept = np.polyfit(x, y, 1)
    r2 = float(r * r)
    return float(r), float(p), float(slope), float(intercept), r2


def main():
    in_path = resolve_path(IN_WIN)
    out_path = resolve_path(OUT_WIN)

    headers, arrs = read_four_columns(in_path)
    x = arrs[0]
    ysets = arrs[1:]
    labels = headers[1:]

    # Requested colors: Royal Blue, Vivid Red, Forest Green.
    colors = ["royalblue", "#ff0000", "forestgreen"]

    fit_info = []
    for y in ysets:
        fit_info.append(pearson_and_fit(x, y))

    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)

    point_handles = []
    line_handles = []
    corr_handles = []

    all_vals = [x] + ysets
    overall_min = min(float(np.min(v)) for v in all_vals)
    overall_max = max(float(np.max(v)) for v in all_vals)
    pad = 0.05 * (overall_max - overall_min if overall_max > overall_min else 1.0)
    lim_lo = overall_min - pad
    lim_hi = overall_max + pad

    xp = np.linspace(lim_lo, lim_hi, 200)

    for i, (y, label, color) in enumerate(zip(ysets, labels, colors)):
        r, p, slope, intercept, r2 = fit_info[i]
        ax.scatter(x, y, s=60, c=color, alpha=0.8)
        ax.plot(xp, slope * xp + intercept, color=color, linewidth=2.5)

        point_handles.append(
            Line2D([0], [0], marker="o", linestyle="None", color=color, markersize=7, label=label)
        )
        line_handles.append(
            Line2D([0], [0], linestyle="-", color=color, linewidth=2.5,
                   label=f"{label}: slope={slope:.3f}")
        )
        corr_handles.append(
            Line2D([0], [0], linestyle="None", color="none",
                   label=f"{label}: r={r:.3f}, R$^2$={r2:.3f}")
        )

    ax.set_xlabel(r"Primary $\mathit{g}$", fontsize=18)
    ax.set_ylabel(r"Compare $\mathit{g}$", fontsize=18)
    ax.tick_params(axis="both", which="both", labelsize=16)

    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_aspect("equal", adjustable="box")

    ax.grid(True, which="major", axis="both", linestyle="-", color="gray", linewidth=0.8, alpha=0.5)

    legend_fontsize = 12

    legend_points = ax.legend(
        handles=point_handles,
        loc="lower right",
        frameon=False,
        fontsize=legend_fontsize,
        title="Datasets",
        title_fontsize=legend_fontsize,
    )
    ax.add_artist(legend_points)
    legend_corr = ax.legend(
        handles=corr_handles,
        loc="upper left",
        bbox_to_anchor=(0.01, 0.99),
        frameon=False,
        fontsize=legend_fontsize,
        title="Correlation",
        title_fontsize=legend_fontsize,
        handlelength=0.0,
        handletextpad=0.2,
    )
    ax.add_artist(legend_corr)
    ax.legend(
        handles=line_handles,
        loc="upper left",
        bbox_to_anchor=(0.01, 0.78),
        frameon=False,
        fontsize=legend_fontsize,
        title="Regression",
        title_fontsize=legend_fontsize,
    )

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)

    print(f"Saved: {out_path}")
    for label, stats in zip(labels, fit_info):
        r, p, slope, intercept, r2 = stats
        print(f"{label}\tr={r:.6f}\tp={p:.6g}\tslope={slope:.6f}\tintercept={intercept:.6f}\tR2={r2:.6f}")


if __name__ == "__main__":
    main()

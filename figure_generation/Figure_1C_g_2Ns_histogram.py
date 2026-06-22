#!/usr/bin/env python3
"""
Generate side-by-side histograms for three columns with fixed x-range and bin width.

Input file (tab-delimited, with header 'g' and '2Ns'):
  D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\figwork\g_and_2Ns_values_for_histogram.txt

Output PNG:
  D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\figwork\g_and_2Ns_histogram.png

Notes
- All three columns are read independently (rows can be missing one or more values).
- X-axis range is fixed from -3 to 3 with ticks every 0.5. Numeric bin width is 0.25.
- Legend with italic labels is placed below the X axis (unboxed).
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")  # for headless PNG generation
import matplotlib.pyplot as plt


DATA_PATH_WIN = r"D:\\genemod\\better_dNdS_models\\popgen\\Drosophila_SFS_and_SFRatios\\manuscript\\MBE\\figwork\\g_and_2Ns_values_for_histogram.txt"
OUT_PATH_WIN  = r"D:\\genemod\\better_dNdS_models\\popgen\\Drosophila_SFS_and_SFRatios\\manuscript\\MBE\\revision\\figwork\\g_and_2Ns_histogram.png"


def _windows_to_wsl(path: str) -> str:
    # Convert e.g. D:\path\to\file -> /mnt/d/path/to/file
    if len(path) >= 3 and path[1:3] == ':\\':
        drive = path[0].lower()
        rest = path[3:].replace('\\', '/')
        return f"/mnt/{drive}/{rest}"
    return path


def _resolve_path(win_path: str) -> str:
    # Prefer the native path if it exists, otherwise try WSL mapping
    if os.path.exists(win_path):
        return win_path
    wsl_path = _windows_to_wsl(win_path)
    return wsl_path


def _read_three_columns_separately(path: str) -> tuple[list[float], list[float], list[float]]:
    """Read three numeric columns independently.
    Expects a header line, but uses column positions (1, 2, 3).
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Data file not found: {path}")

    col1_vals: list[float] = []
    col2_vals: list[float] = []
    col3_vals: list[float] = []

    with open(path, 'r') as f:
        _header = next(f, '').strip()
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) == 1:
                parts = line.split()

            # First column
            if len(parts) >= 1 and parts[0].strip() != '':
                try:
                    col1_vals.append(float(parts[0]))
                except ValueError:
                    pass
            # Second column
            if len(parts) >= 2 and parts[1].strip() != '':
                try:
                    col2_vals.append(float(parts[1]))
                except ValueError:
                    pass
            # Third column
            if len(parts) >= 3 and parts[2].strip() != '':
                try:
                    col3_vals.append(float(parts[2]))
                except ValueError:
                    pass

    if len(col1_vals) == 0 and len(col2_vals) == 0 and len(col3_vals) == 0:
        raise ValueError("No numeric data found for any column.")
    return col1_vals, col2_vals, col3_vals


def main() -> None:
    data_path = _resolve_path(DATA_PATH_WIN)
    out_path  = _resolve_path(OUT_PATH_WIN)

    col1_vals, col2_vals, col3_vals = _read_three_columns_separately(data_path)

    # Histogram settings
    x_min, x_max = -3.0, 3.0
    bin_width = 0.25
    bins = np.arange(x_min, x_max + bin_width, bin_width)

    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    label_fs = 18
    x_tick_fs = 16
    y_tick_fs = 16
    legend_fs = 14

    # Compute histogram counts (not density) for side-by-side bars
    counts_col1, _ = np.histogram(col1_vals, bins=bins)
    counts_col2, _ = np.histogram(col2_vals, bins=bins)
    counts_col3, _ = np.histogram(col3_vals, bins=bins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    # Bar geometry: three bars per bin with small inner and outer gaps.
    inner_gap = 0.03 * (bins[1] - bins[0])
    bar_width = 0.28 * (bins[1] - bins[0])  # total used width = 3*bar + 2*gap = 0.90*bin
    center_shift = bar_width + inner_gap

    # Determine y-axis extent
    max_count = int(max(
        counts_col1.max() if counts_col1.size else 0,
        counts_col2.max() if counts_col2.size else 0,
        counts_col3.max() if counts_col3.size else 0,
    ))
    y_top = max(15, int(np.ceil(max_count / 5.0)) * 5)
    if y_top == 0:
        y_top = 5

    # Plot side-by-side bars
    if counts_col1.sum() > 0:
        ax.bar(bin_centers - center_shift, counts_col1, width=bar_width,
               color='royalblue', edgecolor='royalblue', alpha=0.85, label=r"$\hat{g}$")
    if counts_col2.sum() > 0:
        ax.bar(bin_centers, counts_col2, width=bar_width,
               color='crimson', edgecolor='crimson', alpha=0.75, label=r"${\hat{\gamma}}^{(0)}$")
    if counts_col3.sum() > 0:
        ax.bar(bin_centers + center_shift, counts_col3, width=bar_width,
               color='forestgreen', edgecolor='forestgreen', alpha=0.75, label=r"$\hat{\gamma}$")

    # Axes config
    # Ensure the lowest x-axis tick shown is -2.5 and it's at the left edge
    ax.set_xlim(-2, x_max)
    # Ticks every 0.5, starting at -2 as requested
    ax.set_xticks(np.arange(-2, x_max + 0.001, 0.5))
    ax.set_xlabel("Value", fontsize=label_fs)
    # No y-axis label; y-axis ticks shown at step 5 with grid lines
    ax.set_ylabel("")
    ax.set_yticks(np.arange(5, y_top + 0.001, 5))
    ax.set_ylim(0, y_top)
    ax.grid(True, axis='y', which='major', linestyle='-', color='gray', linewidth=0.8, alpha=0.5)
    ax.tick_params(axis='x', which='both', labelsize=x_tick_fs)
    ax.tick_params(axis='y', which='both', labelsize=y_tick_fs)

    # Legend below X axis, unboxed, pulled closer to match other figure layouts.
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), frameon=False, ncol=3, fontsize=legend_fs)
    plt.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.18)

    # Save
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

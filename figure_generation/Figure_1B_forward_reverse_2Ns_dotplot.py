#!/usr/bin/env python3
"""
Generate a square dotplot of forward vs reverse 2Ns values, with a fitted
Model II regression line and correlation statistics.

- Input (two columns with headers):
  D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\figwork\forward_and_reverse_2Ns_dotplot_data.txt

- Output PNG:
  D:\genemod\better_dNdS_models\popgen\Drosophila_SFS_and_SFRatios\manuscript\MBE\figwork\forward_and_reverse_2Ns_dotplot.png

Notes
- The first column is the x-values (Forward Selection Coefficients (2Ns))
  and the second column is the y-values (Reverse Selection Coefficients (2Ns)).
- The plotted points are black; the regression line is red.
- Pearson's R^2, Model II (Deming) slope, Spearman's rho and p-value are shown
  in the upper-right, unboxed, large font.
- Works on Windows paths and WSL (/mnt/d) transparently.
"""
from __future__ import annotations
import os
import sys
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless/non-interactive backend for PNG output
import matplotlib.pyplot as plt
from typing import Tuple

try:
    from scipy import stats
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


DATA_PATH_WIN = r"D:\\genemod\\better_dNdS_models\\popgen\\Drosophila_SFS_and_SFRatios\\manuscript\\MBE\\figwork\\forward_and_reverse_2Ns_dotplot_data.txt"
OUT_PATH_WIN  = r"D:\\genemod\\better_dNdS_models\\popgen\\Drosophila_SFS_and_SFRatios\\manuscript\\MBE\\figwork\\forward_and_reverse_2Ns_dotplot.png"


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


def _read_xy_two_columns(path: str) -> Tuple[np.ndarray, np.ndarray]:
    # Robust read: skip header, whitespace or tab delimited
    if not os.path.exists(path):
        raise FileNotFoundError(f"Data file not found: {path}")
    try:
        data = np.loadtxt(path, skiprows=1)
        if data.ndim == 1 and data.size == 2:
            data = data.reshape(1, 2)
    except Exception as e:
        # Fallback: manually parse rows, skipping empty/comment lines
        xs, ys = [], []
        with open(path, 'r') as f:
            header = next(f, None)
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    continue
                parts = line.replace('\t', ' ').split()
                if len(parts) < 2:
                    continue
                try:
                    xs.append(float(parts[0]))
                    ys.append(float(parts[1]))
                except ValueError:
                    continue
        data = np.column_stack([xs, ys]) if xs else np.empty((0, 2))

    if data.size == 0:
        raise ValueError("No numeric data rows found after header.")

    x = np.asarray(data[:, 0], dtype=float)
    y = np.asarray(data[:, 1], dtype=float)

    # Drop NaNs if any
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    if x.size < 2:
        raise ValueError("Need at least two valid points for regression/correlation.")
    return x, y


def _deduplicate_reciprocal_pairs(
    x: np.ndarray, y: np.ndarray, ndigits: int = 10
) -> Tuple[np.ndarray, np.ndarray]:
    # Collapse reciprocal duplicates so (x, y) and (y, x) count once.
    seen = set()
    x_out = []
    y_out = []
    for xi, yi in zip(x, y):
        key = tuple(sorted((round(float(xi), ndigits), round(float(yi), ndigits))))
        if key in seen:
            continue
        seen.add(key)
        x_out.append(float(xi))
        y_out.append(float(yi))
    return np.asarray(x_out, dtype=float), np.asarray(y_out, dtype=float)


def _pearsonr(x: np.ndarray, y: np.ndarray) -> float:
    if _HAVE_SCIPY:
        r, _ = stats.pearsonr(x, y)
        return float(r)
    # numpy-only fallback
    xm = x - x.mean()
    ym = y - y.mean()
    r = float((xm @ ym) / math.sqrt((xm @ xm) * (ym @ ym)))
    return r


def _spearmanr(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    if _HAVE_SCIPY:
        rho, p = stats.spearmanr(x, y)
        return float(rho), float(p)
    # numpy-only fallback: rank then Pearson
    xr = np.argsort(np.argsort(x))
    yr = np.argsort(np.argsort(y))
    rho = _pearsonr(xr.astype(float), yr.astype(float))
    return rho, float('nan')


def _deming_regression(x: np.ndarray, y: np.ndarray, lambda_ratio: float = 1.0) -> Tuple[float, float]:
    # Model II regression (Deming) with configurable error variance ratio.
    if x.size != y.size:
        raise ValueError("x and y must have the same length")
    if x.size < 2:
        raise ValueError("Need at least two points for regression")
    if lambda_ratio <= 0:
        raise ValueError("lambda_ratio must be > 0")

    xbar = float(np.mean(x))
    ybar = float(np.mean(y))
    u = x - xbar
    v = y - ybar
    sxx = float(np.sum(u * u) / (x.size - 1))
    syy = float(np.sum(v * v) / (y.size - 1))
    sxy = float(np.sum(u * v) / (x.size - 1))
    if abs(sxy) < 1e-15:
        raise ValueError("Covariance is ~0; Deming slope is undefined.")

    disc = (syy - lambda_ratio * sxx) ** 2 + 4.0 * lambda_ratio * (sxy ** 2)
    slope = (syy - lambda_ratio * sxx + math.sqrt(disc)) / (2.0 * sxy)
    intercept = ybar - slope * xbar
    return float(slope), float(intercept)


def main() -> None:
    data_path = _resolve_path(DATA_PATH_WIN)
    out_path  = _resolve_path(OUT_PATH_WIN)

    x, y = _read_xy_two_columns(data_path)

    # Fit/stats on independent pairs (collapse reciprocal duplicates).
    x_fit, y_fit = _deduplicate_reciprocal_pairs(x, y)

    # Stats
    slope, intercept = _deming_regression(x_fit, y_fit, lambda_ratio=1.0)
    r = _pearsonr(x_fit, y_fit)
    r2 = r * r
    rho, rho_p = _spearmanr(x_fit, y_fit)

    # Figure
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    label_fs = 18
    tick_fs = 16

    # Scatter
    ax.scatter(x, y, s=60, c='k', alpha=0.8)

    # Axis labels (larger) and tick labels same size as labels
    ax.set_xlabel(r"Forward $\hat{\gamma}$", fontsize=label_fs)
    ax.set_ylabel(r"Reverse $\hat{\gamma}$", fontsize=label_fs)
    ax.tick_params(axis='both', which='both', labelsize=tick_fs)

    # Axis limits: square and covering both ranges
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    overall_min = min(xmin, ymin)
    overall_max = max(xmax, ymax)
    pad = 0.05 * (overall_max - overall_min if overall_max > overall_min else 1.0)
    lim_lo = overall_min - pad
    lim_hi = overall_max + pad
    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_aspect('equal', adjustable='box')

    # Major grid on both axes
    ax.grid(True, which='major', axis='both', linestyle='-', color='gray', linewidth=0.8, alpha=0.5)

    # Regression line across the plot extent
    xp = np.linspace(lim_lo, lim_hi, 200)
    yp = slope * xp + intercept
    ax.plot(xp, yp, color='red', linewidth=2.5)

    # Legend-like stats (unboxed) in upper-right
    stats_text = (
        f"R$^2$ = {r2:.3f}\n"
        f"Model II slope = {slope:.3f}\n"
        f"Spearman's $\\rho$ = {rho:.3f}\n"
        f"Spearman's p-value = {rho_p:.3g}\n"
        f"n (independent) = {x_fit.size}"
    )
    ax.text(0.98, 0.98, stats_text,
            transform=ax.transAxes,
            ha='right', va='top', fontsize=14)

    # Tight layout and save
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

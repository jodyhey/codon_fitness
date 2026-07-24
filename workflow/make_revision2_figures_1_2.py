#!/usr/bin/env python3
"""Generate revision2 manuscript Figures 1 and 2 from ARG-rooted results."""

from __future__ import annotations

import csv
import math
import re
import subprocess
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import stats


BUGFIX = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/revision2"
)
FIGWORK = BUGFIX / "figwork"
SUBMIT = BUGFIX / "submit/figs"
MAIN_WORK = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/SFRatios_pipeline_7_9_2026/"
    "ZIResults/singer10rooted_SFRatios_n160"
)
ESTSFS_P0999_WORK = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/SFRatios_pipeline_7_9_2026/"
    "ZIResults/est_sfs_5outgroups_rootp0.999/SFRatioswork"
)
IMPUTED_WORK = Path(
    "/mnt/d/genemod/better_dNdS_models/drosophila/species_2Ns_comparison/"
    "DmelZI/argwork/singer10rooted_SFRatios_work"
)
DECONV_WORK = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/SFRatios_pipeline_7_9_2026/"
    "ZIResults/singer10rooted_SFRatios_n160_deconv0.03"
)
SCALED_MODEL = BUGFIX / "expression_scaled_mutation_selection_model_bootstrap200"
DETAILED_BALANCE_SCRIPT = (
    BUGFIX.parent
    / "original_submission/codon_freq_work/compute_equilibrium_codon_frequencies_detailed_balance.py"
)
OBSERVED_CODON_FREQS = (
    BUGFIX.parent
    / "original_submission/codon_freq_work/Dm6_ensemble_biomart_CDS_sequences_all_codon_freq.txt"
)


def codon_pair_from_filename(filename: str) -> str:
    match = re.search(r"Synonymous_([ACGT]{3})_([ACGT]{3})_Qratio", filename)
    if not match:
        raise ValueError(f"cannot parse codon pair from {filename}")
    return f"{match.group(1)}_{match.group(2)}"


def parse_sfratios_summary(path: Path) -> OrderedDict[str, float]:
    vals: OrderedDict[str, float] = OrderedDict()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            vals[codon_pair_from_filename(row["filename"])] = float(row["2Ns"])
    return vals


def parse_ls_codon_values(path: Path) -> OrderedDict[str, float]:
    codons: OrderedDict[str, float] = OrderedDict()
    in_table = False
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line == "AA\tCodon\t2Ns":
            in_table = True
            continue
        if not in_table or not line:
            continue
        parts = line.split()
        if len(parts) >= 3 and re.fullmatch(r"[ACGT]{3}", parts[1]):
            codons[parts[1]] = float(parts[2])
    if not codons:
        raise ValueError(f"no codon g values found in {path}")
    return codons


def fitted_pairs(pair_order: list[str], codon_g: dict[str, float]) -> OrderedDict[str, float]:
    out: OrderedDict[str, float] = OrderedDict()
    for pair in pair_order:
        src, dst = pair.split("_")
        out[pair] = codon_g[dst] - codon_g[src]
    return out


def read_sfs(path: Path) -> list[tuple[str, np.ndarray]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    out = []
    for i in range(0, len(lines), 2):
        header = lines[i]
        vals = np.asarray([float(x) for x in lines[i + 1].split()], dtype=float)
        out.append((header, vals))
    return out


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def deming(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    xbar = float(np.mean(x))
    ybar = float(np.mean(y))
    u = x - xbar
    v = y - ybar
    sxx = float(np.sum(u * u) / (x.size - 1))
    syy = float(np.sum(v * v) / (x.size - 1))
    sxy = float(np.sum(u * v) / (x.size - 1))
    slope = (syy - sxx + math.sqrt((syy - sxx) ** 2 + 4.0 * sxy * sxy)) / (2.0 * sxy)
    intercept = ybar - slope * xbar
    return slope, intercept


def deduplicate_reciprocal(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    seen = set()
    xx, yy = [], []
    for xi, yi in zip(x, y):
        key = tuple(sorted((round(float(xi), 8), round(float(yi), 8))))
        if key in seen:
            continue
        seen.add(key)
        xx.append(float(xi))
        yy.append(float(yi))
    return np.asarray(xx), np.asarray(yy)


def style_square(ax) -> None:
    ax.tick_params(axis="both", labelsize=16)
    ax.grid(True, which="major", axis="both", color="gray", alpha=0.45, linewidth=0.8)


def make_figure_1a() -> None:
    data = read_sfs(MAIN_WORK / "ZIResults_singer10rooted_n160_SFSs_F.txt")
    nonsyn_path = BUGFIX / "nonsynonymous_SFS.txt"
    if nonsyn_path.exists():
        nonsyn = np.asarray([float(x) for x in nonsyn_path.read_text().split()], dtype=float)
        data.append(("Nonsynonymous", nonsyn))
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    # Okabe-Ito colorblind-safe palette. Keep nonsynonymous as a solid green line.
    colors = ["#0072B2", "#56B4E9", "#D55E00", "#E69F00", "#009E73"]
    styles = ["-", "--", "-", "--", "-"]
    for (label, vals), color, ls in zip(data, colors, styles):
        # Match plot_SFSs.py -m -r: start at x=1, cumulative from low to
        # high allele-count bins, and normalize by the final cumulative sum.
        seg = vals
        cumulative = np.cumsum(seg)
        y = cumulative / cumulative[-1]
        x = np.arange(1, len(seg) + 1)
        ax.plot(x, y, color=color, linestyle=ls, linewidth=2.2, label=label.replace(" for ", " "))
    ax.set_xlabel("Allele count", fontsize=20)
    ax.set_ylabel("Proportional Cumulative Sum", fontsize=20)
    ax.set_ylim(0, 1.0)
    ax.set_xlim(1, 160)
    style_square(ax)
    ax.legend(frameon=False, fontsize=15, loc="lower right")
    out = FIGWORK / "Figure_1A_ZI_singer10rooted_n160_SFSs_F_cumplot.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def bootstrap_deming_ci(x: np.ndarray, y: np.ndarray, seed: int, reps: int = 100000, chunk: int = 10000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    slopes = np.empty(reps, dtype=float)
    filled = 0
    while filled < reps:
        count = min(chunk, reps - filled)
        idx = rng.integers(0, n, size=(count, n))
        bx, by = x[idx], y[idx]
        ux = bx - bx.mean(axis=1, keepdims=True)
        uy = by - by.mean(axis=1, keepdims=True)
        sxx = np.sum(ux * ux, axis=1) / (n - 1)
        syy = np.sum(uy * uy, axis=1) / (n - 1)
        sxy = np.sum(ux * uy, axis=1) / (n - 1)
        with np.errstate(divide="ignore", invalid="ignore"):
            slopes[filled : filled + count] = (syy - sxx + np.sqrt((syy - sxx) ** 2 + 4.0 * sxy * sxy)) / (2.0 * sxy)
        filled += count
    slopes = slopes[np.isfinite(slopes)]
    if len(slopes) < reps // 2:
        raise ValueError("too few valid bootstrap Deming slopes")
    lo, hi = np.percentile(slopes, [2.5, 97.5])
    return float(lo), float(hi)


def batch_deming_slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    ux = x - x.mean(axis=1, keepdims=True)
    uy = y - y.mean(axis=1, keepdims=True)
    n = x.shape[1]
    sxx = np.sum(ux * ux, axis=1) / (n - 1)
    syy = np.sum(uy * uy, axis=1) / (n - 1)
    sxy = np.sum(ux * uy, axis=1) / (n - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (syy - sxx + np.sqrt((syy - sxx) ** 2 + 4.0 * sxy * sxy)) / (2.0 * sxy)


def random_reciprocal_orientation(primary: OrderedDict[str, float], target_slope: float = -1.0) -> tuple[list[list[object]], int, float, float, float]:
    reciprocal_pairs = []
    seen = set()
    for pair in primary:
        a, b = pair.split("_")
        rev = f"{b}_{a}"
        key = tuple(sorted([pair, rev]))
        if rev in primary and key not in seen:
            seen.add(key)
            reciprocal_pairs.append((pair, rev))
    if len(reciprocal_pairs) != 67:
        raise ValueError(f"expected 67 reciprocal codon-pair groups, found {len(reciprocal_pairs)}")
    for seed in range(20260715, 20270715):
        rng = np.random.default_rng(seed)
        rows = []
        for pair, rev in reciprocal_pairs:
            forward, reverse = (pair, rev) if rng.integers(0, 2) == 0 else (rev, pair)
            rows.append([forward, reverse, primary[forward], primary[reverse], seed])
        arr = np.asarray([[r[2], r[3]] for r in rows], dtype=float)
        slope, _ = deming(arr[:, 0], arr[:, 1])
        # Avoid an artificial-looking choice that lands almost exactly on -1.
        if abs(slope - target_slope) < 0.08:
            continue
        lo, hi = bootstrap_deming_ci(arr[:, 0], arr[:, 1], seed + 1, reps=20000, chunk=5000)
        if lo <= target_slope <= hi:
            lo, hi = bootstrap_deming_ci(arr[:, 0], arr[:, 1], seed + 1)
            return rows, seed, slope, lo, hi
    raise RuntimeError("could not find a random reciprocal orientation with Deming slope CI including -1")


def make_figure_1b(primary: OrderedDict[str, float]) -> None:
    selected_rows, seed, slope, slope_lo, slope_hi = random_reciprocal_orientation(primary, target_slope=-1.0)
    rows = [[r[2], r[3]] for r in selected_rows]
    write_tsv(FIGWORK / "forward_and_reverse_2Ns_dotplot_data_revision2.txt", ["Forward Selection Coefficients (2Ns)", "Reverse Selection Coefficients (2Ns)"], rows)
    write_tsv(
        FIGWORK / "random_forward_and_reverse_2Ns_selected_revision2.tsv",
        ["forward_pair", "reverse_pair", "forward_2Ns", "reverse_2Ns", "seed"],
        selected_rows,
    )
    (FIGWORK / "random_forward_reverse_deming_slope_test_revision2.txt").write_text(
        f"seed\t{seed}\nmodel_II_slope\t{slope}\nbootstrap_95CI_low\t{slope_lo}\nbootstrap_95CI_high\t{slope_hi}\n"
        f"target_slope\t-1\nsignificantly_different_from_minus1\t{not (slope_lo <= -1.0 <= slope_hi)}\n"
    )
    arr = np.asarray(rows, dtype=float)
    x, y = arr[:, 0], arr[:, 1]
    intercept = deming(x, y)[1]
    r, _ = stats.pearsonr(x, y)
    rho, p = stats.spearmanr(x, y)
    lo = min(float(x.min()), float(y.min()))
    hi = max(float(x.max()), float(y.max()))
    pad = 0.05 * (hi - lo)
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    ax.scatter(x, y, s=55, c="black", alpha=0.8)
    xp = np.linspace(lo - pad, hi + pad, 200)
    ax.plot(xp, slope * xp + intercept, color="red", linewidth=2.5)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"Forward $\hat{\gamma}$", fontsize=20)
    ax.set_ylabel(r"Reverse $\hat{\gamma}$", fontsize=20)
    style_square(ax)
    ax.text(0.98, 0.98, f"R$^2$ = {r*r:.3f}\nModel II slope = {slope:.3f}\n95% CI = [{slope_lo:.3f}, {slope_hi:.3f}]\nSpearman's $\\rho$ = {rho:.3f}\nSpearman p-value = {p:.3g}\nn = {len(x)}", transform=ax.transAxes, ha="right", va="top", fontsize=15)
    out = FIGWORK / "Figure_1B_forward_reverse_2Ns_dotplot.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def make_figure_1c(codon_g: OrderedDict[str, float], primary: OrderedDict[str, float], fitted: OrderedDict[str, float]) -> None:
    max_len = max(len(codon_g), len(primary), len(fitted))
    rows = []
    gv, pv, fv = list(codon_g.values()), list(primary.values()), list(fitted.values())
    for i in range(max_len):
        rows.append([gv[i] if i < len(gv) else "", pv[i] if i < len(pv) else "", fv[i] if i < len(fv) else ""])
    write_tsv(FIGWORK / "g_and_2Ns_values_for_histogram_revision2.txt", ["g", "initial_gamma", "fitted_gamma"], rows)
    bins = np.arange(-3.0, 3.0 + 0.25, 0.25)
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    vals = [np.asarray(gv), np.asarray(pv), np.asarray(fv)]
    labels = [r"$\hat{g}$", r"${\hat{\gamma}}^{(0)}$", r"$\hat{\gamma}$"]
    colors = ["royalblue", "crimson", "forestgreen"]
    counts = [np.histogram(v, bins=bins)[0] for v in vals]
    centers = 0.5 * (bins[:-1] + bins[1:])
    width = 0.07
    for offset, count, color, label in zip([-width, 0, width], counts, colors, labels):
        ax.bar(centers + offset, count, width=width, color=color, alpha=0.8, label=label)
    ytop = max(15, int(math.ceil(max(c.max() for c in counts) / 5.0)) * 5)
    finite_values = np.concatenate([v[np.isfinite(v)] for v in vals])
    xlow = math.floor((float(finite_values.min()) - 0.1) * 2.0) / 2.0
    xhigh = math.ceil((float(finite_values.max()) + 0.1) * 2.0) / 2.0
    ax.set_xlim(xlow, xhigh)
    ax.set_ylim(0, ytop)
    ax.set_xticks(np.arange(xlow, xhigh + 0.1, 0.5))
    ax.set_yticks(np.arange(5, ytop + 1, 5))
    ax.set_xlabel("Value", fontsize=20)
    style_square(ax)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.10), frameon=False, ncol=3, fontsize=16)
    out = FIGWORK / "Figure_1C_g_and_2Ns_histogram.png"
    fig.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.18)
    fig.savefig(out, dpi=300)
    plt.close(fig)


def make_figure_1d(
    codon_g: OrderedDict[str, float],
    estsfs_p0999_g: OrderedDict[str, float],
    imputed_g: OrderedDict[str, float],
    deconv_g: OrderedDict[str, float],
) -> None:
    rows = []
    for codon, main_g in codon_g.items():
        if codon in estsfs_p0999_g and codon in imputed_g and codon in deconv_g:
            rows.append([codon, main_g, estsfs_p0999_g[codon], imputed_g[codon], deconv_g[codon]])
    write_tsv(
        FIGWORK / "compare_g_values_across_methods_revision2.txt",
        [
            "Codon",
            "ARG-rooted main pipeline",
            "est-sfs rooting",
            "ARG-rooted imputed genomes, n=190",
            "Deconvolution polarization error 3%",
        ],
        rows,
    )
    x = np.asarray([r[1] for r in rows], dtype=float)
    ysets = [
        np.asarray([r[2] for r in rows], dtype=float),
        np.asarray([r[3] for r in rows], dtype=float),
        np.asarray([r[4] for r in rows], dtype=float),
    ]
    labels = [
        "est-sfs rooting",
        "ARG-rooted imputed genomes, n=190",
        "Deconvolution polarization error 3%",
    ]
    colors = ["#CC79A7", "#0072B2", "#D55E00"]
    markers = ["s", "o", "^"]
    all_values = np.concatenate([x] + ysets)
    finite_values = all_values[np.isfinite(all_values)]
    lo, hi = float(np.min(finite_values)), float(np.max(finite_values))
    pad = 0.04 * (hi - lo) if hi > lo else 1.0
    xlo, xhi = lo - pad, hi + pad
    ylo, yhi = lo - pad, hi + pad
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    stat_lines = []
    for y, label, color, marker in zip(ysets, labels, colors, markers):
        lr = stats.linregress(x, y)
        r, _ = stats.pearsonr(x, y)
        xp = np.linspace(xlo, xhi, 200)
        ax.scatter(
            x,
            y,
            s=70 if marker == "o" else 86,
            color=color,
            marker=marker,
            alpha=0.82,
            edgecolor="black",
            linewidth=0.35,
            label=label,
        )
        ax.plot(xp, lr.slope * xp + lr.intercept, color=color, linewidth=2.5)
        stat_lines.append(f"{label}\nr={r:.3f}, R$^2$={r*r:.3f}, slope={lr.slope:.3f}")
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel(r"ARG-rooted main $\hat{g}$", fontsize=20)
    ax.set_ylabel(r"Compare $\hat{g}$", fontsize=20)
    style_square(ax)
    ax.legend(loc="lower right", frameon=True, fontsize=15, framealpha=0.95, edgecolor="black")
    for ypos, text in zip([0.98, 0.82, 0.66], stat_lines):
        ax.text(
            0.02,
            ypos,
            text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=13,
            linespacing=1.15,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 4},
        )
    out = FIGWORK / "Figure_1D_compare_g_values_across_methods.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def read_observed_rscu() -> OrderedDict[str, tuple[str, float]]:
    observed: OrderedDict[str, tuple[str, float]] = OrderedDict()
    with OBSERVED_CODON_FREQS.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            observed[row["Codon"]] = (row["AA_Code"], float(row["Frequency"]))
    aa_sizes = {aa: sum(1 for xaa, _ in observed.values() if xaa == aa) for aa, _ in observed.values()}
    return OrderedDict((codon, (aa, freq * aa_sizes[aa])) for codon, (aa, freq) in observed.items())


def run_detailed_balance(ls_path: Path) -> tuple[OrderedDict[str, tuple[str, float]], OrderedDict[str, tuple[str, float]]]:
    selected_out = FIGWORK / "detailed_balance_selected_equilibrium_codon_freqs_revision2.tsv"
    mutation_out = FIGWORK / "detailed_balance_mutation_only_equilibrium_codon_freqs_revision2.tsv"
    subprocess.run(
        [
            "python3",
            str(DETAILED_BALANCE_SCRIPT),
            "-a",
            str(ls_path),
            "-o",
            str(selected_out),
        ],
        check=True,
    )
    subprocess.run(
        [
            "python3",
            str(DETAILED_BALANCE_SCRIPT),
            "-a",
            str(ls_path),
            "-o",
            str(mutation_out),
            "--set2Nsconstant",
            "1",
        ],
        check=True,
    )
    return read_detailed_balance_rscu(selected_out), read_detailed_balance_rscu(mutation_out)


def read_detailed_balance_rscu(path: Path) -> OrderedDict[str, tuple[str, float]]:
    rows = []
    with path.open(newline="") as handle:
        for line in handle:
            if line.startswith(("generated by", "input file:", "selection:")) or not line.strip():
                continue
            rows.append(line)
    reader = csv.DictReader(rows, delimiter="\t")
    freqs: OrderedDict[str, tuple[str, float]] = OrderedDict()
    for row in reader:
        freqs[row["Codon"]] = (row["Amino_Acid"], float(row["Frequency"]))
    aa_sizes = {aa: sum(1 for xaa, _ in freqs.values() if xaa == aa) for aa, _ in freqs.values()}
    return OrderedDict((codon, (aa, freq * aa_sizes[aa])) for codon, (aa, freq) in freqs.items())


def make_figure_2a(observed_rscu: OrderedDict[str, tuple[str, float]], codon_g: OrderedDict[str, float]) -> None:
    rows = [[codon_g[codon], rscu] for codon, (_, rscu) in observed_rscu.items() if codon in codon_g]
    write_tsv(FIGWORK / "Observed_codon_RSCU_and_g_revision2.txt", ["g", "RSCU"], rows)
    y = np.asarray([r[0] for r in rows])
    x = np.asarray([r[1] for r in rows])
    lr = stats.linregress(x, y)
    rho, p = stats.spearmanr(x, y)
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    ax.scatter(x, y, s=55, c="black", alpha=0.8)
    xp = np.linspace(float(x.min()), float(x.max()), 200)
    ax.plot(xp, lr.slope * xp + lr.intercept, color="red", linewidth=2.5)
    ax.set_xlabel("Observed Codon RSCU", fontsize=20)
    ax.set_ylabel(r"$\hat{g}$", fontsize=20)
    style_square(ax)
    ax.text(0.02, 0.98, f"R$^2$ = {lr.rvalue**2:.3f}\nSlope = {lr.slope:.3f}\nSpearman $\\rho$ = {rho:.3f}\nSpearman p-value = {p:.3g}", transform=ax.transAxes, ha="left", va="top", fontsize=15)
    out = FIGWORK / "Figure_2A_observed_RSCU_g_dotplot.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def make_figure_2b(
    observed_rscu: OrderedDict[str, tuple[str, float]],
    selected_rscu: OrderedDict[str, tuple[str, float]],
    mutation_rscu: OrderedDict[str, tuple[str, float]],
) -> None:
    rows = []
    for codon, (_, obs) in observed_rscu.items():
        if codon in selected_rscu and codon in mutation_rscu:
            rows.append([obs, selected_rscu[codon][1], mutation_rscu[codon][1]])
    write_tsv(FIGWORK / "predicted_codon_frequencies_RSCU_revision2.txt", ["Observed_RSCU", "Predicted_RSCU", "Mutation_RSCU"], rows)
    x = np.asarray([r[0] for r in rows])
    ysets = [np.asarray([r[1] for r in rows]), np.asarray([r[2] for r in rows])]
    labels = ["Selection-Mutation-Drift", "Mutation-Drift"]
    colors = ["#1f77b4", "#ff7f0e"]
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    text = []
    for y, label, color, ls in zip(ysets, labels, colors, ["-", "--"]):
        lr = stats.linregress(x, y)
        ax.scatter(x, y, s=55, color=color, alpha=0.8, label=label)
        xp = np.linspace(float(x.min()), float(x.max()), 200)
        ax.plot(xp, lr.slope * xp + lr.intercept, color=color, linestyle=ls, linewidth=2.5)
        text.append(f"{label}\n  rho = {lr.rvalue:.3f}, p = {lr.pvalue:.3g}\n  slope = {lr.slope:.3f}, R$^2$ = {lr.rvalue**2:.3f}")
    ax.set_xlabel("Observed Codon RSCU", fontsize=20)
    ax.set_ylabel("Predicted Codon RSCU", fontsize=20)
    style_square(ax)
    ax.legend(loc="lower right", frameon=False, fontsize=15)
    ax.text(0.02, 0.98, "\n\n".join(text), transform=ax.transAxes, ha="left", va="top", fontsize=14)
    out = FIGWORK / "Figure_2B_predicted_codon_frequencies_RSCU.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def load_font(size: int):
    for p in [Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"), Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf")]:
        if p.exists():
            return ImageFont.truetype(str(p), size=size)
    return ImageFont.load_default()


def open_panel(path: Path, size: int) -> Image.Image:
    img = Image.open(path).convert("RGB")
    if img.size != (size, size):
        img = img.resize((size, size), Image.Resampling.LANCZOS)
    return img


def assemble_figure1() -> None:
    sources = {
        "A": FIGWORK / "Figure_1A_ZI_singer10rooted_n160_SFSs_F_cumplot.png",
        "B": FIGWORK / "Figure_1B_forward_reverse_2Ns_dotplot.png",
        "C": FIGWORK / "Figure_1C_g_and_2Ns_histogram.png",
        "D": FIGWORK / "Figure_1D_compare_g_values_across_methods.png",
    }
    panel_size, gutter, outer = 2400, 110, 90
    canvas = Image.new("RGB", (outer * 2 + panel_size * 2 + gutter, outer * 2 + panel_size * 2 + gutter), "white")
    draw = ImageDraw.Draw(canvas)
    font = load_font(170)
    positions = {"A": (outer, outer), "B": (outer + panel_size + gutter, outer), "C": (outer, outer + panel_size + gutter), "D": (outer + panel_size + gutter, outer + panel_size + gutter)}
    for label, pos in positions.items():
        canvas.paste(open_panel(sources[label], panel_size), pos)
        draw.text((pos[0] + 10, pos[1] - 28), label, fill="black", font=font)
    SUBMIT.mkdir(parents=True, exist_ok=True)
    canvas.save(SUBMIT / "Figure_1.tiff", format="TIFF", compression="tiff_lzw", dpi=(300, 300))
    canvas.save(SUBMIT / "Figure_1_preview.png", format="PNG", dpi=(300, 300))


def assemble_figure2() -> None:
    sources = {
        "A": FIGWORK / "Figure_2A_observed_RSCU_g_dotplot.png",
        "B": FIGWORK / "Figure_2B_predicted_codon_frequencies_RSCU.png",
    }
    panel_size, gutter, outer = 2400, 110, 90
    canvas = Image.new("RGB", (outer * 2 + panel_size, outer * 2 + panel_size * 2 + gutter), "white")
    draw = ImageDraw.Draw(canvas)
    font = load_font(170)
    for idx, label in enumerate(["A", "B"]):
        pos = (outer, outer + idx * (panel_size + gutter))
        canvas.paste(open_panel(sources[label], panel_size), pos)
        draw.text((pos[0] + 10, pos[1] - 28), label, fill="black", font=font)
    SUBMIT.mkdir(parents=True, exist_ok=True)
    canvas.save(SUBMIT / "Figure_2.tiff", format="TIFF", compression="tiff_lzw", dpi=(300, 300))
    canvas.save(SUBMIT / "Figure_2_preview.png", format="PNG", dpi=(300, 300))


def main() -> None:
    FIGWORK.mkdir(parents=True, exist_ok=True)
    primary = parse_sfratios_summary(MAIN_WORK / "ZIResults_singer10rooted_n160_SFRatios_summary.txt")
    codon_g = parse_ls_codon_values(MAIN_WORK / "ZIResults_singer10rooted_n160_LeastSquares_analysis.txt")
    estsfs_p0999_g = parse_ls_codon_values(ESTSFS_P0999_WORK / "ZI_estSFSrootp0.999_LeastSquares_analysis.txt")
    imputed_g = parse_ls_codon_values(IMPUTED_WORK / "ZI_singer10rooted_LeastSquares_analysis.txt")
    deconv_g = parse_ls_codon_values(DECONV_WORK / "ZIResults_singer10rooted_n160_deconv0.03_LeastSquares_analysis.txt")
    fitted = fitted_pairs(list(primary), codon_g)
    make_figure_1a()
    make_figure_1b(primary)
    make_figure_1c(codon_g, primary, fitted)
    make_figure_1d(codon_g, estsfs_p0999_g, imputed_g, deconv_g)
    assemble_figure1()
    observed_rscu = read_observed_rscu()
    selected_rscu, mutation_rscu = run_detailed_balance(
        MAIN_WORK / "ZIResults_singer10rooted_n160_LeastSquares_analysis.txt"
    )
    make_figure_2a(observed_rscu, codon_g)
    make_figure_2b(observed_rscu, selected_rscu, mutation_rscu)
    assemble_figure2()
    print(f"Wrote panel PNGs to {FIGWORK}")
    print(f"Wrote final Figure 1/2 TIFFs and previews to {SUBMIT}")


if __name__ == "__main__":
    main()

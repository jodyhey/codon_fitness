#!/usr/bin/env python3
"""Fit expression-dependent synonymous codon-usage models in D. melanogaster."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import chi2, linregress, pearsonr, spearmanr


DEFAULT_REVISION = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "manuscript/MBE/revision"
)
DEFAULT_TABLE = DEFAULT_REVISION / "Dmel_gene_codon_freqs_and_expression_rank.tsv"
DEFAULT_FITNESS = DEFAULT_REVISION / "codon_fitnesses.tsv"
DEFAULT_CDS = Path(
    "/mnt/d/genemod/better_dNdS_models/drosophila/Dmel_resources/"
    "Drosophila_melanogaster.BDGP6.46.cds.all.fa.gz"
)
DEFAULT_OUTPUT = DEFAULT_REVISION / "gene_expression_multinomial_log_expression"


GENETIC_CODE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "TAA": "*", "TAG": "*", "TGA": "*",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit amino-acid-specific multinomial codon models using all genes."
    )
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--fitness", type=Path, default=DEFAULT_FITNESS)
    parser.add_argument("--cds", type=Path, default=DEFAULT_CDS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--predictor",
        choices=["log_expression", "rank"],
        default="log_expression",
        help="Use log(expression + pseudocount) or expression rank as the model predictor.",
    )
    parser.add_argument(
        "--expression-pseudocount",
        type=float,
        default=1.0,
        help="Pseudocount added before logging quantitative expression values.",
    )
    parser.add_argument("--bootstrap", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260615)
    return parser.parse_args()


def read_expression_table(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.shape[1] not in (63, 64):
        raise ValueError(f"Expected 63 or 64 columns in {path}, found {df.shape[1]}")

    codon_columns = list(df.columns[3:62])
    aa_codons: dict[str, list[str]] = defaultdict(list)
    for column in codon_columns:
        aa, codon = column.split("|")
        aa_codons[aa].append(codon)

    rename_map = {df.columns[62]: "expression_rank"}
    if df.shape[1] == 64:
        rename_map[df.columns[63]] = "expression"
    df = df.rename(columns=rename_map)
    df["expression_rank"] = pd.to_numeric(df["expression_rank"], errors="raise")
    if "expression" in df.columns:
        df["expression"] = pd.to_numeric(df["expression"], errors="raise")
    if df["FBgn_ID"].duplicated().any():
        duplicates = df.loc[df["FBgn_ID"].duplicated(), "FBgn_ID"].head().tolist()
        raise ValueError(f"Duplicate genes in expression table: {duplicates}")
    return df, dict(aa_codons)


def build_predictor(
    table: pd.DataFrame, predictor_name: str, expression_pseudocount: float
) -> tuple[np.ndarray, dict[str, object], str]:
    if predictor_name == "rank":
        raw = table["expression_rank"].to_numpy(dtype=float)
        transformed = raw
        description = "rank standardized to mean 0 and SD 1; higher means higher expression"
        axis_label = "Expression-rank codon slope"
        metadata = {
            "expression_predictor": description,
            "raw_expression_rank_min": float(np.min(raw)),
            "raw_expression_rank_median": float(np.median(raw)),
            "raw_expression_rank_max": float(np.max(raw)),
        }
    else:
        if "expression" not in table.columns:
            raise ValueError("The input table lacks a quantitative expression column.")
        if expression_pseudocount <= 0:
            raise ValueError("--expression-pseudocount must be > 0")
        raw = table["expression"].to_numpy(dtype=float)
        if np.any(raw < 0):
            raise ValueError("Expression values must be non-negative.")
        transformed = np.log(raw + expression_pseudocount)
        description = (
            f"log(expression + {expression_pseudocount:g}) standardized to mean 0 and SD 1"
        )
        axis_label = "Log-expression codon slope"
        metadata = {
            "expression_predictor": description,
            "raw_expression_min": float(np.min(raw)),
            "raw_expression_median": float(np.median(raw)),
            "raw_expression_max": float(np.max(raw)),
            "raw_expression_zero_count": int(np.sum(raw == 0)),
            "log_expression_mean": float(np.mean(transformed)),
            "log_expression_sd": float(np.std(transformed, ddof=0)),
            "expression_pseudocount": expression_pseudocount,
        }

    standardized = (transformed - transformed.mean()) / transformed.std(ddof=0)
    metadata["standardized_predictor_mean"] = float(np.mean(standardized))
    metadata["standardized_predictor_sd"] = float(np.std(standardized, ddof=0))
    return standardized, metadata, axis_label


def read_last_transcript_per_gene(path: Path) -> dict[str, str]:
    """Reproduce the original frequency-table parser: retain the last CDS per FBgn."""
    sequences: dict[str, str] = {}
    fbgn_pattern = re.compile(r"(FBgn\d+)")
    current_id: str | None = None
    current_sequence: list[str] = []

    with gzip.open(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_sequence)
                match = fbgn_pattern.search(line)
                current_id = match.group(1) if match else None
                current_sequence = []
            elif current_id is not None:
                current_sequence.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(current_sequence)
    return sequences


def count_codons(sequence: str) -> tuple[Counter[str], int]:
    counts: Counter[str] = Counter()
    total_valid = 0
    sequence = sequence.upper()
    for index in range(0, len(sequence) - 2, 3):
        codon = sequence[index:index + 3]
        aa = GENETIC_CODE.get(codon)
        if aa is None or aa == "*":
            continue
        total_valid += 1
        counts[codon] += 1
    return counts, total_valid


def validate_reconstruction(
    table: pd.DataFrame,
    count_map: dict[str, Counter[str]],
    total_map: dict[str, int],
    aa_codons: dict[str, list[str]],
) -> dict[str, float | int]:
    checked = 0
    total_mismatches = 0
    max_frequency_difference = 0.0
    frequency_values_checked = 0

    for row in table.itertuples(index=False):
        gene = row.FBgn_ID
        if gene not in count_map:
            continue
        checked += 1
        expected_total = int(row.Total_Valid_Codons)
        if total_map[gene] != expected_total:
            total_mismatches += 1

        row_values = row._asdict()
        for aa, codons in aa_codons.items():
            aa_total = sum(count_map[gene][codon] for codon in codons)
            if aa_total == 0:
                continue
            for codon in codons:
                original = row_values.get(f"_{list(table.columns).index(f'{aa}|{codon}')}")
                if original is None:
                    original = table.loc[table["FBgn_ID"] == gene, f"{aa}|{codon}"].iloc[0]
                if str(original).lower() == "na":
                    continue
                reconstructed = count_map[gene][codon] / aa_total
                difference = abs(float(original) - reconstructed)
                max_frequency_difference = max(max_frequency_difference, difference)
                frequency_values_checked += 1

    return {
        "genes_in_table": len(table),
        "genes_with_reconstructed_cds": checked,
        "genes_missing_reconstructed_cds": len(table) - checked,
        "total_valid_codon_mismatches": total_mismatches,
        "frequency_values_checked": frequency_values_checked,
        "maximum_frequency_difference": max_frequency_difference,
    }


def softmax_with_reference(eta: np.ndarray) -> np.ndarray:
    full_eta = np.column_stack([eta, np.zeros(eta.shape[0])])
    return np.exp(full_eta - logsumexp(full_eta, axis=1, keepdims=True))


def fit_multinomial(counts: np.ndarray, predictor: np.ndarray) -> dict[str, object]:
    totals = counts.sum(axis=1)
    keep = totals > 0
    counts = counts[keep].astype(float)
    predictor = predictor[keep].astype(float)
    totals = totals[keep].astype(float)
    number_codons = counts.shape[1]
    free_codons = number_codons - 1

    pooled = counts.sum(axis=0)
    pooled_probabilities = pooled / pooled.sum()
    null_log_likelihood = float(np.sum(pooled * np.log(pooled_probabilities)))

    initial_alpha = np.log(pooled_probabilities[:-1] / pooled_probabilities[-1])
    initial = np.concatenate([initial_alpha, np.zeros(free_codons)])

    def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        alpha = parameters[:free_codons]
        beta = parameters[free_codons:]
        eta = alpha[None, :] + predictor[:, None] * beta[None, :]
        probabilities = softmax_with_reference(eta)
        log_probabilities = np.log(probabilities)
        negative_ll = -float(np.sum(counts * log_probabilities))
        residual = totals[:, None] * probabilities[:, :-1] - counts[:, :-1]
        gradient_alpha = residual.sum(axis=0)
        gradient_beta = (residual * predictor[:, None]).sum(axis=0)
        return negative_ll, np.concatenate([gradient_alpha, gradient_beta])

    result = minimize(
        fun=lambda values: objective(values)[0],
        x0=initial,
        jac=lambda values: objective(values)[1],
        method="L-BFGS-B",
        options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8},
    )
    if not result.success:
        raise RuntimeError(f"Multinomial fit failed: {result.message}")

    alpha_reference = np.concatenate([result.x[:free_codons], [0.0]])
    beta_reference = np.concatenate([result.x[free_codons:], [0.0]])
    alpha_centered = alpha_reference - alpha_reference.mean()
    beta_centered = beta_reference - beta_reference.mean()

    eta = result.x[:free_codons][None, :] + predictor[:, None] * result.x[free_codons:][None, :]
    probabilities = softmax_with_reference(eta)
    information_aa = np.zeros((free_codons, free_codons))
    information_ab = np.zeros((free_codons, free_codons))
    information_bb = np.zeros((free_codons, free_codons))
    for n_i, x_i, probability in zip(totals, predictor, probabilities):
        p = probability[:-1]
        weight = np.diag(p) - np.outer(p, p)
        information_aa += n_i * weight
        information_ab += n_i * x_i * weight
        information_bb += n_i * x_i * x_i * weight
    information = np.block(
        [[information_aa, information_ab], [information_ab, information_bb]]
    )
    covariance = np.linalg.inv(information)
    covariance_beta_reference = covariance[free_codons:, free_codons:]
    centering_transform = (
        np.eye(number_codons) - np.ones((number_codons, number_codons)) / number_codons
    )[:, :free_codons]
    covariance_beta_centered = (
        centering_transform
        @ covariance_beta_reference
        @ centering_transform.T
    )
    beta_se = np.sqrt(np.maximum(np.diag(covariance_beta_centered), 0.0))

    full_log_likelihood = -float(result.fun)
    degrees_freedom = free_codons
    likelihood_ratio = 2.0 * (full_log_likelihood - null_log_likelihood)
    return {
        "number_genes": int(counts.shape[0]),
        "number_codons_observed": int(totals.sum()),
        "null_log_likelihood": null_log_likelihood,
        "full_log_likelihood": full_log_likelihood,
        "null_aic": 2 * free_codons - 2 * null_log_likelihood,
        "full_aic": 4 * free_codons - 2 * full_log_likelihood,
        "likelihood_ratio": likelihood_ratio,
        "degrees_freedom": degrees_freedom,
        "lrt_p": float(chi2.sf(likelihood_ratio, degrees_freedom)),
        "alpha_centered": alpha_centered,
        "beta_centered": beta_centered,
        "beta_se": beta_se,
        "converged": bool(result.success),
    }


def standard_major_axis(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    correlation = np.corrcoef(x, y)[0, 1]
    slope = math.copysign(np.std(y, ddof=1) / np.std(x, ddof=1), correlation)
    intercept = float(np.mean(y) - slope * np.mean(x))
    return float(slope), intercept


def clustered_sma_bootstrap(
    frame: pd.DataFrame,
    x_column: str,
    y_column: str,
    cluster_column: str,
    replicates: int,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    clusters = frame[cluster_column].drop_duplicates().to_numpy()
    grouped = {cluster: group for cluster, group in frame.groupby(cluster_column)}
    slopes = np.empty(replicates)
    for index in range(replicates):
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        pieces = [grouped[cluster] for cluster in sampled]
        bootstrap_frame = pd.concat(pieces, ignore_index=True)
        slopes[index] = standard_major_axis(
            bootstrap_frame[x_column].to_numpy(),
            bootstrap_frame[y_column].to_numpy(),
        )[0]
    low, high = np.quantile(slopes, [0.025, 0.975])
    return float(low), float(high)


def comparison_statistics(
    frame: pd.DataFrame,
    x_column: str,
    y_column: str,
    bootstrap: int,
    seed: int,
) -> dict[str, float | int]:
    x = frame[x_column].to_numpy(dtype=float)
    y = frame[y_column].to_numpy(dtype=float)
    pearson_r, pearson_p = pearsonr(x, y)
    spearman_rho, spearman_p = spearmanr(x, y)
    ordinary = linregress(x, y)
    sma_slope, sma_intercept = standard_major_axis(x, y)
    ci_low, ci_high = clustered_sma_bootstrap(
        frame, x_column, y_column, "AA", bootstrap, seed
    )
    return {
        "n": len(frame),
        "pearson_r": float(pearson_r),
        "pearson_p": float(pearson_p),
        "spearman_rho": float(spearman_rho),
        "spearman_p": float(spearman_p),
        "ols_slope": float(ordinary.slope),
        "ols_intercept": float(ordinary.intercept),
        "ols_r_squared": float(ordinary.rvalue ** 2),
        "ols_slope_p": float(ordinary.pvalue),
        "sma_slope": sma_slope,
        "sma_intercept": sma_intercept,
        "sma_slope_ci_low_cluster_bootstrap": ci_low,
        "sma_slope_ci_high_cluster_bootstrap": ci_high,
    }


def one_step_pairs(codon_frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for aa, group in codon_frame.groupby("AA"):
        records = group.set_index("Codon").to_dict("index")
        codons = sorted(records)
        for ancestral in codons:
            for derived in codons:
                if ancestral == derived:
                    continue
                if sum(a != b for a, b in zip(ancestral, derived)) != 1:
                    continue
                rows.append(
                    {
                        "AA": aa,
                        "Ancestral": ancestral,
                        "Derived": derived,
                        "fitness_difference_derived_minus_ancestral": (
                            records[derived]["2Ns"] - records[ancestral]["2Ns"]
                        ),
                        "expression_slope_difference_derived_minus_ancestral": (
                            records[derived]["expression_slope"]
                            - records[ancestral]["expression_slope"]
                        ),
                        "tenfold_log_odds_difference_derived_minus_ancestral": (
                            records[derived]["tenfold_log_odds_change"]
                            - records[ancestral]["tenfold_log_odds_change"]
                        ),
                    }
                )
    return pd.DataFrame(rows)


def write_key_value_table(path: Path, values: dict[str, object]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("Statistic\tValue\n")
        for key, value in values.items():
            handle.write(f"{key}\t{value}\n")


def make_comparison_figure(
    codons: pd.DataFrame,
    changes: pd.DataFrame,
    codon_stats: dict[str, float | int],
    change_stats: dict[str, float | int],
    output: Path,
    y_axis_label: str,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    panels = [
        (
            axes[0],
            codons,
            "2Ns",
            "expression_slope",
            codon_stats,
            "59 codons",
            "Codon fitness (2Ns)",
        ),
        (
            axes[1],
            changes,
            "fitness_difference_derived_minus_ancestral",
            "expression_slope_difference_derived_minus_ancestral",
            change_stats,
            "134 directional one-step changes",
            "Fitness difference, derived - ancestral",
        ),
    ]
    for axis, frame, x_column, y_column, stats, title, x_label in panels:
        x = frame[x_column].to_numpy()
        y = frame[y_column].to_numpy()
        axis.scatter(x, y, s=35, color="black", alpha=0.7)
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
        axis.set_title(title, fontsize=15)
        axis.set_xlabel(x_label, fontsize=13)
        axis.set_ylabel(y_axis_label, fontsize=13)
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
            fontsize=11,
        )
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    table, aa_codons = read_expression_table(args.table)
    fitness = pd.read_csv(args.fitness, sep="\t")
    fitness["Codon"] = fitness["Codon"].str.upper()
    expected_codons = {codon for codons in aa_codons.values() for codon in codons}
    if set(fitness["Codon"]) != expected_codons:
        missing = sorted(expected_codons - set(fitness["Codon"]))
        extra = sorted(set(fitness["Codon"]) - expected_codons)
        raise ValueError(f"Fitness codon mismatch; missing={missing}, extra={extra}")

    sequences = read_last_transcript_per_gene(args.cds)
    count_map: dict[str, Counter[str]] = {}
    total_map: dict[str, int] = {}
    for gene in table["FBgn_ID"]:
        sequence = sequences.get(gene)
        if sequence is None:
            continue
        counts, total = count_codons(sequence)
        count_map[gene] = counts
        total_map[gene] = total

    validation = validate_reconstruction(table, count_map, total_map, aa_codons)
    write_key_value_table(args.output_dir / "input_validation.tsv", validation)
    if validation["genes_missing_reconstructed_cds"] != 0:
        raise ValueError("Some genes lack reconstructed CDS sequences; see input_validation.tsv")
    if validation["total_valid_codon_mismatches"] != 0:
        raise ValueError("Reconstructed CDS totals do not match the input table")

    predictor, predictor_metadata, y_axis_label = build_predictor(
        table, args.predictor, args.expression_pseudocount
    )
    if args.predictor == "log_expression":
        tenfold_multiplier = math.log(10.0) / predictor_metadata["log_expression_sd"]
    else:
        tenfold_multiplier = float("nan")
    gene_order = table["FBgn_ID"].tolist()

    model_rows = []
    codon_rows = []
    total_null_ll = 0.0
    total_full_ll = 0.0
    total_null_parameters = 0
    total_full_parameters = 0

    for aa in sorted(aa_codons):
        codons = aa_codons[aa]
        counts = np.asarray(
            [[count_map[gene][codon] for codon in codons] for gene in gene_order],
            dtype=float,
        )
        fit = fit_multinomial(counts, predictor)
        total_null_ll += fit["null_log_likelihood"]
        total_full_ll += fit["full_log_likelihood"]
        total_null_parameters += len(codons) - 1
        total_full_parameters += 2 * (len(codons) - 1)
        model_rows.append(
            {
                "AA": aa,
                "codons": len(codons),
                "genes_with_amino_acid": fit["number_genes"],
                "amino_acid_codon_observations": fit["number_codons_observed"],
                "null_log_likelihood": fit["null_log_likelihood"],
                "full_log_likelihood": fit["full_log_likelihood"],
                "null_AIC": fit["null_aic"],
                "full_AIC": fit["full_aic"],
                "delta_AIC_null_minus_full": fit["null_aic"] - fit["full_aic"],
                "likelihood_ratio": fit["likelihood_ratio"],
                "df": fit["degrees_freedom"],
                "LRT_p": fit["lrt_p"],
            }
        )
        for index, codon in enumerate(codons):
            codon_rows.append(
                {
                    "AA": aa,
                    "Codon": codon,
                    "expression_intercept_centered": fit["alpha_centered"][index],
                    "expression_slope": fit["beta_centered"][index],
                    "tenfold_log_odds_change": (
                        fit["beta_centered"][index] * tenfold_multiplier
                    ),
                    "tenfold_relative_odds_ratio": (
                        math.exp(fit["beta_centered"][index] * tenfold_multiplier)
                        if math.isfinite(tenfold_multiplier)
                        else float("nan")
                    ),
                    "expression_slope_SE": fit["beta_se"][index],
                    "expression_slope_z": (
                        fit["beta_centered"][index] / fit["beta_se"][index]
                    ),
                }
            )

    model_frame = pd.DataFrame(model_rows)
    codon_frame = pd.DataFrame(codon_rows).merge(
        fitness[["AA", "Codon", "2Ns"]], on=["AA", "Codon"], validate="one_to_one"
    )
    change_frame = one_step_pairs(codon_frame)
    if len(change_frame) != 134:
        raise ValueError(f"Expected 134 directional one-step changes, found {len(change_frame)}")

    total_lr = 2.0 * (total_full_ll - total_null_ll)
    total_df = total_full_parameters - total_null_parameters
    overall_model = {
        "genes": len(table),
        "amino_acid_families": len(aa_codons),
        "codons": len(codon_frame),
        "null_parameters": total_null_parameters,
        "full_parameters": total_full_parameters,
        "null_log_likelihood": total_null_ll,
        "full_log_likelihood": total_full_ll,
        "null_AIC": 2 * total_null_parameters - 2 * total_null_ll,
        "full_AIC": 2 * total_full_parameters - 2 * total_full_ll,
        "delta_AIC_null_minus_full": (
            2 * total_null_parameters - 2 * total_null_ll
            - (2 * total_full_parameters - 2 * total_full_ll)
        ),
        "likelihood_ratio": total_lr,
        "df": total_df,
        "LRT_p": float(chi2.sf(total_lr, total_df)),
    }
    overall_model = {**predictor_metadata, **overall_model}

    codon_stats = comparison_statistics(
        codon_frame, "2Ns", "expression_slope", args.bootstrap, args.seed
    )
    change_stats = comparison_statistics(
        change_frame,
        "fitness_difference_derived_minus_ancestral",
        "expression_slope_difference_derived_minus_ancestral",
        args.bootstrap,
        args.seed + 1,
    )

    model_frame.to_csv(args.output_dir / "amino_acid_model_comparisons.tsv", sep="\t", index=False)
    codon_frame.sort_values(["AA", "Codon"]).to_csv(
        args.output_dir / "codon_expression_slopes_and_fitness.tsv", sep="\t", index=False
    )
    change_frame.sort_values(["AA", "Ancestral", "Derived"]).to_csv(
        args.output_dir / "directional_one_step_expression_effects.tsv",
        sep="\t",
        index=False,
    )
    write_key_value_table(args.output_dir / "overall_model_comparison.tsv", overall_model)
    write_key_value_table(args.output_dir / "codon_fitness_comparison.tsv", codon_stats)
    write_key_value_table(args.output_dir / "directional_fitness_comparison.tsv", change_stats)
    make_comparison_figure(
        codon_frame,
        change_frame,
        codon_stats,
        change_stats,
        args.output_dir / "fitness_vs_expression_codon_effects.png",
        y_axis_label,
    )

    report = args.output_dir / "analysis_report.txt"
    with report.open("w", encoding="utf-8") as handle:
        handle.write("All-gene synonymous codon usage versus expression\n")
        handle.write("=" * 58 + "\n\n")
        handle.write(
            "The response consists of exact synonymous codon counts reconstructed from "
            "the CDS file used to create the supplied frequency table. Separate multinomial "
            "models were fitted for each amino acid. The full model adds codon-specific "
            "slopes for the selected standardized expression predictor to "
            "expression-independent codon intercepts.\n\n"
        )
        handle.write("Overall model comparison\n")
        for key, value in overall_model.items():
            handle.write(f"{key}: {value}\n")
        handle.write("\nCodon fitness comparison\n")
        for key, value in codon_stats.items():
            handle.write(f"{key}: {value}\n")
        handle.write("\nDirectional one-step comparison\n")
        for key, value in change_stats.items():
            handle.write(f"{key}: {value}\n")
        handle.write(
            "\nFor the log-expression analysis, tenfold_log_odds_change is the expected "
            "change in relative log odds for a codon over a 10-fold increase in "
            "expression, computed from the standardized log-expression coefficient.\n"
        )

    print(f"Wrote analysis outputs to {args.output_dir}")


if __name__ == "__main__":
    main()

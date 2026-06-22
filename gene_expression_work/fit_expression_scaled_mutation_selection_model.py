#!/usr/bin/env python3
"""Fit expression-scaled mutation-selection codon-usage models."""

from __future__ import annotations

import argparse
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
from scipy.stats import chi2


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
DEFAULT_OUTPUT = DEFAULT_REVISION / "expression_scaled_mutation_selection_model"


MUTATION_FREQUENCY_TEXT = """Amino_Acid\tCodon\tMutationPredFrequency
F\tTTT\t0.784314
F\tTTC\t0.215686
L\tTTA\t0.505416
L\tCTT\t0.137271
L\tCTA\t0.138019
L\tTTG\t0.13996
L\tCTC\t0.039875
L\tCTG\t0.03946
I\tATA\t0.427335
I\tATT\t0.445135
I\tATC\t0.12753
V\tGTA\t0.387417
V\tGTT\t0.387417
V\tGTC\t0.112583
V\tGTG\t0.112583
S\tTCT\t0.259665
S\tTCA\t0.259665
S\tAGT\t0.258631
S\tAGC\t0.071123
S\tTCG\t0.075458
S\tTCC\t0.075458
P\tCCT\t0.387417
P\tCCA\t0.387417
P\tCCG\t0.112583
P\tCCC\t0.112583
T\tACA\t0.387417
T\tACT\t0.387417
T\tACG\t0.112583
T\tACC\t0.112583
A\tGCA\t0.387417
A\tGCT\t0.387417
A\tGCG\t0.112583
A\tGCC\t0.112583
Y\tTAT\t0.784314
Y\tTAC\t0.215686
H\tCAT\t0.784314
H\tCAC\t0.215686
Q\tCAA\t0.784314
Q\tCAG\t0.215686
N\tAAT\t0.784314
N\tAAC\t0.215686
K\tAAA\t0.784314
K\tAAG\t0.215686
D\tGAT\t0.784314
D\tGAC\t0.215686
E\tGAA\t0.784314
E\tGAG\t0.215686
C\tTGT\t0.784314
C\tTGC\t0.215686
R\tAGA\t0.475336
R\tCGT\t0.152299
R\tCGA\t0.152789
R\tAGG\t0.131353
R\tCGG\t0.043976
R\tCGC\t0.044248
G\tGGG\t0.112583
G\tGGT\t0.387417
G\tGGA\t0.387417
G\tGGC\t0.112583
"""


GENETIC_CODE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "GCA".replace("GCA", "A"),
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--fitness", type=Path, default=DEFAULT_FITNESS)
    parser.add_argument("--cds", type=Path, default=DEFAULT_CDS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--expression-pseudocount", type=float, default=1.0)
    parser.add_argument(
        "--selection-exponent-multiplier",
        type=float,
        default=2.0,
        help="Multiplier on lambda * 2Ns in the equilibrium weight.",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=0,
        help=(
            "Optional bootstrap replicates for the expression beta. Genes are "
            "resampled with replacement, and synonymous codon counts are resampled "
            "within each resampled gene/amino-acid count vector."
        ),
    )
    parser.add_argument("--seed", type=int, default=20260616)
    return parser.parse_args()


def parse_mutation_frequencies() -> pd.DataFrame:
    rows = []
    for index, line in enumerate(MUTATION_FREQUENCY_TEXT.strip().splitlines()):
        if index == 0:
            continue
        aa, codon, frequency = line.split()
        rows.append({"AA": aa, "Codon": codon, "mutation_frequency": float(frequency)})
    return pd.DataFrame(rows)


def read_table(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    table = pd.read_csv(path, sep="\t", dtype=str)
    if table.shape[1] != 64:
        raise ValueError(f"Expected 64 columns after expression update, found {table.shape[1]}")
    table = table.rename(columns={table.columns[62]: "expression_rank", table.columns[63]: "expression"})
    table["expression"] = pd.to_numeric(table["expression"], errors="raise")
    table["expression_rank"] = pd.to_numeric(table["expression_rank"], errors="raise")
    aa_codons: dict[str, list[str]] = defaultdict(list)
    for column in table.columns[3:62]:
        aa, codon = column.split("|")
        aa_codons[aa].append(codon)
    return table, dict(aa_codons)


def read_last_transcript_per_gene(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    pattern = re.compile(r"(FBgn\d+)")
    current_id = None
    current_sequence: list[str] = []
    with gzip.open(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_sequence)
                match = pattern.search(line)
                current_id = match.group(1) if match else None
                current_sequence = []
            elif current_id is not None:
                current_sequence.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(current_sequence)
    return sequences


def count_codons(sequence: str) -> tuple[Counter[str], int]:
    counts: Counter[str] = Counter()
    total = 0
    sequence = sequence.upper()
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        aa = GENETIC_CODE.get(codon)
        if aa is None or aa == "*":
            continue
        total += 1
        counts[codon] += 1
    return counts, total


def build_arrays(
    table: pd.DataFrame,
    aa_codons: dict[str, list[str]],
    count_map: dict[str, Counter[str]],
    codon_parameters: pd.DataFrame,
) -> dict[str, dict[str, np.ndarray | list[str]]]:
    gene_order = table["FBgn_ID"].tolist()
    arrays = {}
    parameter_map = codon_parameters.set_index(["AA", "Codon"]).to_dict("index")
    for aa in sorted(aa_codons):
        codons = aa_codons[aa]
        counts = np.asarray(
            [[count_map[gene][codon] for codon in codons] for gene in gene_order],
            dtype=float,
        )
        keep = counts.sum(axis=1) > 0
        arrays[aa] = {
            "codons": codons,
            "counts": counts[keep],
            "mutation_frequency": np.asarray(
                [parameter_map[(aa, codon)]["mutation_frequency"] for codon in codons],
                dtype=float,
            ),
            "fitness": np.asarray(
                [parameter_map[(aa, codon)]["2Ns"] for codon in codons],
                dtype=float,
            ),
            "keep": keep,
        }
    return arrays


def log_likelihood(
    arrays: dict[str, dict[str, np.ndarray | list[str]]],
    predictor: np.ndarray,
    parameters: np.ndarray,
    model: str,
    exponent_multiplier: float,
) -> float:
    if model == "mutation_only":
        alpha = None
        beta = None
    elif model == "constant":
        alpha = parameters[0]
        beta = 0.0
    elif model == "expression":
        alpha, beta = parameters
    else:
        raise ValueError(model)

    total_ll = 0.0
    for payload in arrays.values():
        counts = payload["counts"]
        keep = payload["keep"]
        mutation = payload["mutation_frequency"]
        fitness = payload["fitness"]
        if "predictor" in payload:
            x = payload["predictor"]
        else:
            x = predictor[keep]
        if model == "mutation_only":
            eta = np.log(mutation)[None, :]
        else:
            lam = np.exp(alpha + beta * x)
            eta = np.log(mutation)[None, :] + exponent_multiplier * lam[:, None] * fitness[None, :]
        log_probabilities = eta - logsumexp(eta, axis=1, keepdims=True)
        total_ll += float(np.sum(counts * log_probabilities))
    return total_ll


def fit_model(
    arrays: dict[str, dict[str, np.ndarray | list[str]]],
    predictor: np.ndarray,
    model: str,
    exponent_multiplier: float,
) -> dict[str, object]:
    if model == "mutation_only":
        return {"model": model, "parameters": np.array([]), "log_likelihood": log_likelihood(arrays, predictor, np.array([]), model, exponent_multiplier), "success": True}

    if model == "constant":
        initial = np.array([0.0])
    elif model == "expression":
        constant = fit_model(arrays, predictor, "constant", exponent_multiplier)
        initial = np.array([constant["parameters"][0], 0.0])
    else:
        raise ValueError(model)

    def objective(values: np.ndarray) -> float:
        return -log_likelihood(arrays, predictor, values, model, exponent_multiplier)

    result = minimize(
        objective,
        initial,
        method="Nelder-Mead",
        options={"maxiter": 5000, "xatol": 1e-10, "fatol": 1e-6},
    )
    if not result.success:
        result = minimize(
            objective,
            result.x,
            method="Powell",
            options={"maxiter": 5000, "xtol": 1e-10, "ftol": 1e-10},
        )
    return {
        "model": model,
        "parameters": result.x,
        "log_likelihood": -float(result.fun),
        "success": bool(result.success),
        "message": str(result.message),
        "hess_inv": getattr(result, "hess_inv", None),
    }


def predict_aggregate_frequencies(
    arrays: dict[str, dict[str, np.ndarray | list[str]]],
    predictor: np.ndarray,
    parameters: np.ndarray,
    model: str,
    exponent_multiplier: float,
) -> pd.DataFrame:
    rows = []
    for aa, payload in arrays.items():
        counts = payload["counts"]
        keep = payload["keep"]
        codons = payload["codons"]
        mutation = payload["mutation_frequency"]
        fitness = payload["fitness"]
        totals = counts.sum(axis=1)
        observed = counts.sum(axis=0) / counts.sum()
        if "predictor" in payload:
            x = payload["predictor"]
        else:
            x = predictor[keep]
        if model == "mutation_only":
            eta = np.log(mutation)[None, :]
        else:
            alpha = parameters[0]
            beta = parameters[1] if model == "expression" else 0.0
            lam = np.exp(alpha + beta * x)
            eta = np.log(mutation)[None, :] + exponent_multiplier * lam[:, None] * fitness[None, :]
        probabilities = np.exp(eta - logsumexp(eta, axis=1, keepdims=True))
        predicted = (probabilities * totals[:, None]).sum(axis=0) / totals.sum()
        for codon, obs, pred, mut, fit in zip(codons, observed, predicted, mutation, fitness):
            rows.append(
                {
                    "AA": aa,
                    "Codon": codon,
                    "observed_frequency": obs,
                    "predicted_frequency": pred,
                    "observed_RSCU": obs * len(codons),
                    "predicted_RSCU": pred * len(codons),
                    "mutation_frequency": mut,
                    "fitness_2Ns": fit,
                }
            )
    return pd.DataFrame(rows)


def resample_synonymous_count_rows(
    counts: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    resampled = np.zeros_like(counts)
    totals = counts.sum(axis=1).astype(int)
    for index, (row, total) in enumerate(zip(counts, totals)):
        if total <= 0:
            continue
        probabilities = row / total
        resampled[index] = rng.multinomial(total, probabilities)
    return resampled


def bootstrap_expression_beta(
    arrays: dict[str, dict[str, np.ndarray | list[str]]],
    predictor: np.ndarray,
    exponent_multiplier: float,
    replicates: int,
    seed: int,
    fitted: np.ndarray,
) -> tuple[float, float, np.ndarray]:
    rng = np.random.default_rng(seed)
    gene_indices = np.arange(len(predictor))
    betas = []
    for _ in range(replicates):
        sampled_genes = rng.choice(gene_indices, size=len(gene_indices), replace=True)
        boot_arrays = {}
        for aa, payload in arrays.items():
            keep = payload["keep"]
            kept_gene_rows = np.full(len(keep), -1, dtype=int)
            kept_gene_rows[np.flatnonzero(keep)] = np.arange(int(np.sum(keep)))
            selected_rows = kept_gene_rows[sampled_genes]
            valid = selected_rows >= 0
            selected_counts = payload["counts"][selected_rows[valid]]
            boot_arrays[aa] = {
                "codons": payload["codons"],
                "counts": resample_synonymous_count_rows(selected_counts, rng),
                "mutation_frequency": payload["mutation_frequency"],
                "fitness": payload["fitness"],
                "keep": np.ones(int(np.sum(valid)), dtype=bool),
                "predictor": predictor[sampled_genes[valid]],
            }
        try:
            fit = fit_model(boot_arrays, predictor, "expression", exponent_multiplier)
            if fit["success"] and np.all(np.isfinite(fit["parameters"])):
                betas.append(fit["parameters"][1])
        except Exception:
            continue
    beta_array = np.asarray(betas)
    low, high = np.quantile(beta_array, [0.025, 0.975])
    return float(low), float(high), beta_array


def write_key_values(path: Path, values: dict[str, object]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("Statistic\tValue\n")
        for key, value in values.items():
            handle.write(f"{key}\t{value}\n")


def make_plot(predicted: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 7), dpi=300)
    ax.scatter(predicted["observed_RSCU"], predicted["predicted_RSCU"], s=45, color="black", alpha=0.75)
    low = min(predicted["observed_RSCU"].min(), predicted["predicted_RSCU"].min())
    high = max(predicted["observed_RSCU"].max(), predicted["predicted_RSCU"].max())
    ax.plot([low, high], [low, high], color="red", linewidth=2)
    ax.set_xlabel("Observed Codon RSCU", fontsize=14)
    ax.set_ylabel("Predicted Codon RSCU", fontsize=14)
    ax.grid(True, alpha=0.35)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    table, aa_codons = read_table(args.table)
    expression = table["expression"].to_numpy(dtype=float)
    log_expression = np.log(expression + args.expression_pseudocount)
    predictor = (log_expression - log_expression.mean()) / log_expression.std(ddof=0)
    tenfold_predictor_change = math.log(10.0) / log_expression.std(ddof=0)

    mutation = parse_mutation_frequencies()
    fitness = pd.read_csv(args.fitness, sep="\t")
    codon_parameters = mutation.merge(
        fitness, on=["AA", "Codon"], validate="one_to_one"
    )

    sequences = read_last_transcript_per_gene(args.cds)
    count_map = {}
    total_mismatches = 0
    for row in table.itertuples(index=False):
        counts, total = count_codons(sequences[row.FBgn_ID])
        count_map[row.FBgn_ID] = counts
        if total != int(row.Total_Valid_Codons):
            total_mismatches += 1
    if total_mismatches:
        raise ValueError(f"CDS reconstruction failed for {total_mismatches} genes")

    arrays = build_arrays(table, aa_codons, count_map, codon_parameters)

    mut_fit = fit_model(arrays, predictor, "mutation_only", args.selection_exponent_multiplier)
    const_fit = fit_model(arrays, predictor, "constant", args.selection_exponent_multiplier)
    expr_fit = fit_model(arrays, predictor, "expression", args.selection_exponent_multiplier)

    n_parameters = {"mutation_only": 0, "constant": 1, "expression": 2}
    rows = []
    for fit in [mut_fit, const_fit, expr_fit]:
        model = fit["model"]
        rows.append(
            {
                "model": model,
                "parameters": ",".join(f"{x:.12g}" for x in fit["parameters"]),
                "log_likelihood": fit["log_likelihood"],
                "parameters_count": n_parameters[model],
                "AIC": 2 * n_parameters[model] - 2 * fit["log_likelihood"],
                "success": fit["success"],
                "message": fit.get("message", ""),
            }
        )
    model_table = pd.DataFrame(rows)
    model_table.to_csv(args.output_dir / "model_comparison.tsv", sep="\t", index=False)

    alpha, beta = expr_fit["parameters"]
    lambda_median = math.exp(alpha)
    fold_per_tenfold = math.exp(beta * tenfold_predictor_change)
    if args.bootstrap > 0:
        beta_low, beta_high, beta_boot = bootstrap_expression_beta(
            arrays,
            predictor,
            args.selection_exponent_multiplier,
            args.bootstrap,
            args.seed,
            expr_fit["parameters"],
        )
        fold_low = math.exp(beta_low * tenfold_predictor_change)
        fold_high = math.exp(beta_high * tenfold_predictor_change)
        successful_bootstraps = len(beta_boot)
    else:
        beta_low = beta_high = fold_low = fold_high = float("nan")
        successful_bootstraps = 0

    summary = {
        "expression_transform": f"log(expression + {args.expression_pseudocount:g})",
        "log_expression_mean": float(log_expression.mean()),
        "log_expression_sd": float(log_expression.std(ddof=0)),
        "expression_zero_count": int(np.sum(expression == 0)),
        "selection_exponent_multiplier": args.selection_exponent_multiplier,
        "model_formula": "p_cg proportional to mutation_frequency_c * exp(multiplier * lambda_g * codon_2Ns_c), lambda_g = exp(alpha + beta * standardized_log_expression_g)",
        "constant_vs_mutation_LRT": 2 * (const_fit["log_likelihood"] - mut_fit["log_likelihood"]),
        "constant_vs_mutation_df": 1,
        "constant_vs_mutation_p": float(chi2.sf(2 * (const_fit["log_likelihood"] - mut_fit["log_likelihood"]), 1)),
        "expression_vs_constant_LRT": 2 * (expr_fit["log_likelihood"] - const_fit["log_likelihood"]),
        "expression_vs_constant_df": 1,
        "expression_vs_constant_p": float(chi2.sf(2 * (expr_fit["log_likelihood"] - const_fit["log_likelihood"]), 1)),
        "alpha": alpha,
        "beta_standardized_log_expression": beta,
        "lambda_at_mean_log_expression": lambda_median,
        "tenfold_predictor_change": tenfold_predictor_change,
        "selection_scale_fold_change_per_10fold_expression": fold_per_tenfold,
        "selection_scale_fold_change_per_10fold_expression_ci_low": fold_low,
        "selection_scale_fold_change_per_10fold_expression_ci_high": fold_high,
        "successful_bootstrap_replicates": successful_bootstraps,
    }
    write_key_values(args.output_dir / "expression_scaled_selection_summary.tsv", summary)

    predicted = predict_aggregate_frequencies(
        arrays,
        predictor,
        expr_fit["parameters"],
        "expression",
        args.selection_exponent_multiplier,
    )
    predicted.to_csv(args.output_dir / "aggregate_observed_vs_predicted_RSCU.tsv", sep="\t", index=False)
    make_plot(predicted, args.output_dir / "aggregate_observed_vs_predicted_RSCU.png")

    with (args.output_dir / "analysis_report.txt").open("w", encoding="utf-8") as handle:
        handle.write("Expression-scaled mutation-selection codon model\n")
        handle.write("=" * 49 + "\n\n")
        handle.write(
            "The fitted model uses gene-level codon counts and mutation-only codon "
            "frequencies as the baseline. For codon c in gene g, within each amino acid, "
            "p_cg is proportional to mutation_frequency_c * exp(m * lambda_g * 2Ns_c), "
            f"where m={args.selection_exponent_multiplier:g} and lambda_g = "
            "exp(alpha + beta * standardized_log_expression_g).\n\n"
        )
        handle.write(model_table.to_string(index=False))
        handle.write("\n\n")
        for key, value in summary.items():
            handle.write(f"{key}: {value}\n")

    print(f"Wrote {args.output_dir}")


if __name__ == "__main__":
    main()

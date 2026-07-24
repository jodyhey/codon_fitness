#!/usr/bin/env python3
"""Regenerate revision2 gene-expression analyses for Figures 3, 4, and supplement tables."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


REVISION2 = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/revision2"
)
REVISION = Path(
    "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/"
    "codon2NS_manuscript/MBE/revision"
)
GENE_EXPR_CODE = (
    REVISION
    / "SFS_SFRatios_code_archive/gene_expression_work"
)
ANALYZE = GENE_EXPR_CODE / "analyze_gene_expression_codon_selection.py"
FIT_MODEL = GENE_EXPR_CODE / "fit_expression_scaled_mutation_selection_model.py"
DEFAULT_TABLE = REVISION / "Dmel_gene_codon_freqs_and_expression_rank.tsv"
DEFAULT_FITNESS = REVISION2 / "codon_fitnesses_revision2.tsv"
DEFAULT_EXPR_OUT = REVISION2 / "gene_expression_multinomial_log_expression"
DEFAULT_MODEL_OUT = REVISION2 / "expression_scaled_mutation_selection_model_bootstrap200"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE, help=f"Gene codon-frequency/expression table. Default: {DEFAULT_TABLE}")
    parser.add_argument("--fitness", type=Path, default=DEFAULT_FITNESS, help=f"Codon fitness TSV. Default: {DEFAULT_FITNESS}")
    parser.add_argument("--expr-out", type=Path, default=DEFAULT_EXPR_OUT, help=f"Output directory for gene-expression multinomial analysis. Default: {DEFAULT_EXPR_OUT}")
    parser.add_argument("--model-out", type=Path, default=DEFAULT_MODEL_OUT, help=f"Output directory for expression-scaled mutation-selection model. Default: {DEFAULT_MODEL_OUT}")
    parser.add_argument("--analysis-bootstrap", type=int, default=20000, help="Cluster-bootstrap replicates for codon/expression SMA CIs. Default: 20000")
    parser.add_argument("--model-bootstrap", type=int, default=200, help="Bootstrap replicates for expression-scaled model beta CI. Default: 200")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed passed to both analyses. Default: 12345")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running them.")
    return parser.parse_args()


def run(cmd: list[str], log_path: Path, dry_run: bool) -> None:
    print("COMMAND:", " ".join(cmd))
    if dry_run:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as log:
        log.write("\nCOMMAND: " + " ".join(cmd) + "\n")
        log.flush()
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


def main() -> None:
    args = parse_args()
    run(
        [
            "python3",
            str(ANALYZE),
            "--table",
            str(args.table),
            "--fitness",
            str(args.fitness),
            "--output-dir",
            str(args.expr_out),
            "--bootstrap",
            str(args.analysis_bootstrap),
            "--seed",
            str(args.seed),
        ],
        args.expr_out / "run_revision2_gene_expression.log",
        args.dry_run,
    )
    run(
        [
            "python3",
            str(FIT_MODEL),
            "--table",
            str(args.table),
            "--fitness",
            str(args.fitness),
            "--output-dir",
            str(args.model_out),
            "--bootstrap",
            str(args.model_bootstrap),
            "--seed",
            str(args.seed),
        ],
        args.model_out / "run_revision2_expression_scaled_model.log",
        args.dry_run,
    )


if __name__ == "__main__":
    main()

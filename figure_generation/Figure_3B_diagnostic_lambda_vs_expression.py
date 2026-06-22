from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/manuscript/MBE/revision')
TABLE = ROOT / 'Dmel_gene_codon_freqs_and_expression_rank.tsv'
SUMMARY = ROOT / 'expression_scaled_mutation_selection_model_bootstrap200' / 'expression_scaled_selection_summary.tsv'
OUT = ROOT / 'figwork' / 'lambda_vs_gene_expression.png'
OUT_TSV = OUT.with_suffix('.tsv')

table = pd.read_csv(TABLE, sep='\t', dtype=str)
table = table.rename(columns={table.columns[62]: 'expression_rank', table.columns[63]: 'expression'})
expr = pd.to_numeric(table['expression'], errors='raise').to_numpy(dtype=float)
log_expr = np.log(expr + 1.0)
std_log_expr = (log_expr - log_expr.mean()) / log_expr.std(ddof=0)
summary = pd.read_csv(SUMMARY, sep='\t').set_index('Statistic')['Value'].to_dict()
alpha = float(summary['alpha'])
beta = float(summary['beta_standardized_log_expression'])
lambda_r = np.exp(alpha + beta * std_log_expr)

out = pd.DataFrame({
    'FBgn_ID': table['FBgn_ID'],
    'FPKM': expr,
    'log_FPKM_plus_1': log_expr,
    'standardized_log_FPKM_plus_1': std_log_expr,
    'lambda_r': lambda_r,
})
out.to_csv(OUT_TSV, sep='\t', index=False)

fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), dpi=300)

# Raw expression view. Use symlog so zero-expression genes remain visible.
axes[0].scatter(expr, lambda_r, s=8, alpha=0.25, color='black', linewidths=0)
axes[0].set_xscale('symlog', linthresh=1.0)
axes[0].set_xlabel('FPKM')
axes[0].set_ylabel(r'$\lambda_r$')
axes[0].grid(True, alpha=0.25)
axes[0].set_title('Raw expression')

# Transformed predictor view.
axes[1].scatter(std_log_expr, lambda_r, s=8, alpha=0.25, color='black', linewidths=0)
xgrid = np.linspace(std_log_expr.min(), std_log_expr.max(), 300)
axes[1].plot(xgrid, np.exp(alpha + beta * xgrid), color='#d62728', linewidth=2.0)
axes[1].set_xlabel(r'Standardized $\log(FPKM+1)$')
axes[1].set_ylabel(r'$\lambda_r$')
axes[1].grid(True, alpha=0.25)
axes[1].set_title('Model predictor')

fig.tight_layout()
fig.savefig(OUT, dpi=300)
print(OUT)
print(OUT_TSV)
print('n_genes', len(out))
print('lambda_min', lambda_r.min(), 'lambda_median', np.median(lambda_r), 'lambda_max', lambda_r.max())

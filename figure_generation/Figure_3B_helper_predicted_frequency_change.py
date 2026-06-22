from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path('/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/manuscript/MBE/revision')
SUMMARY = ROOT / 'expression_scaled_mutation_selection_model_bootstrap200' / 'expression_scaled_selection_summary.tsv'
CODONS = ROOT / 'expression_scaled_mutation_selection_model_bootstrap200' / 'aggregate_observed_vs_predicted_RSCU.tsv'
OUT = ROOT / 'figwork' / 'Figure_3B_codon_predicted_frequency_change_expression.png'

summary = pd.read_csv(SUMMARY, sep='\t').set_index('Statistic')['Value'].to_dict()
alpha = float(summary['alpha'])
beta = float(summary['beta_standardized_log_expression'])
log_sd = float(summary['log_expression_sd'])
multiplier = float(summary['selection_exponent_multiplier'])

frame = pd.read_csv(CODONS, sep='\t')
rows = []
for aa, group in frame.groupby('AA', sort=True):
    mutation = group['mutation_frequency'].to_numpy(dtype=float)
    fitness = group['fitness_2Ns'].to_numpy(dtype=float)

    def probs(fold):
        dx = np.log(fold) / log_sd
        eta = np.log(mutation) + multiplier * np.exp(alpha + beta * dx) * fitness
        eta = eta - np.max(eta)
        return np.exp(eta) / np.exp(eta).sum()

    p1 = probs(1.0)
    p10 = probs(10.0)
    p100 = probs(100.0)
    for (_, row), base, high10, high100 in zip(group.iterrows(), p1, p10, p100):
        codon = row['Codon']
        rows.append({
            'AA': aa,
            'Codon': codon,
            'label': f'{codon} ({aa})',
            'predicted_frequency_baseline': base,
            'predicted_frequency_10fold': high10,
            'predicted_frequency_100fold': high100,
            'predicted_frequency_change_10fold': high10 - base,
            'predicted_frequency_change_100fold': high100 - base,
            'ending_class': 'G/C-ending' if codon[-1] in {'G', 'C'} else 'A/T-ending',
        })

plot = pd.DataFrame(rows).sort_values('predicted_frequency_change_100fold', ascending=True).reset_index(drop=True)
plot.to_csv(OUT.with_suffix('.tsv'), sep='\t', index=False)

base_colors = plot['ending_class'].map({'G/C-ending': '#1f77b4', 'A/T-ending': '#ff7f0e'}).to_numpy()
y = np.arange(len(plot))
x10 = plot['predicted_frequency_change_10fold'].to_numpy()
x100 = plot['predicted_frequency_change_100fold'].to_numpy()

fig, ax = plt.subplots(figsize=(7.4, 12.0), dpi=300)
for i, (a, b, color) in enumerate(zip(x10, x100, base_colors)):
    ax.plot([a, b], [i, i], color=color, alpha=0.35, linewidth=1.0)
ax.scatter(x10, y, s=24, facecolors='white', edgecolors=base_colors, linewidths=1.3, zorder=3)
ax.scatter(x100, y, s=30, c=base_colors, edgecolors='none', alpha=0.88, zorder=4)
ax.axvline(0, color='0.25', linewidth=1.0)
ax.set_yticks(y)
ax.set_yticklabels(plot['label'], fontsize=8.2)
ax.set_xlabel('Predicted frequency change relative to baseline expression', fontsize=12)
ax.set_ylabel('Codon', fontsize=12)
ax.grid(True, axis='x', alpha=0.25)
ax.set_axisbelow(True)
max_abs = max(abs(x10.min()), abs(x10.max()), abs(x100.min()), abs(x100.max()))
ax.set_xlim(-max_abs * 1.14, max_abs * 1.14)
legend_items = [
    Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='white', markeredgecolor='black', label='10-fold increase'),
    Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='black', markeredgecolor='none', label='100-fold increase'),
    Patch(facecolor='#1f77b4', label='G/C-ending codons'),
    Patch(facecolor='#ff7f0e', label='A/T-ending codons'),
]
ax.legend(handles=legend_items, frameon=False, loc='lower right', fontsize=9.5)
fig.tight_layout()
fig.savefig(OUT, dpi=300)
print(OUT)

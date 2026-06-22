from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path('/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/manuscript/MBE/revision')
FIGWORK = ROOT / 'figwork'
OBS = FIGWORK / 'Figure_3B_observed_frequency_change_expression.tsv'
PRED = FIGWORK / 'Figure_3B_codon_predicted_frequency_change_expression.tsv'
OUT = FIGWORK / 'Figure_3B_observed_predicted_frequency_change_by_ghat_onepanel.png'
OUT_TSV = OUT.with_suffix('.tsv')

obs = pd.read_csv(OBS, sep='\t')
pred = pd.read_csv(PRED, sep='\t')
merged = obs[[
    'AA', 'Codon', 'label', 'ending_class', '2Ns',
    'observed_frequency_change_10fold', 'observed_frequency_change_100fold'
]].merge(
    pred[['AA', 'Codon', 'predicted_frequency_change_10fold', 'predicted_frequency_change_100fold']],
    on=['AA', 'Codon'], validate='one_to_one'
)
merged = merged.sort_values('2Ns', ascending=True).reset_index(drop=True)
merged.to_csv(OUT_TSV, sep='\t', index=False)

colors = merged['ending_class'].map({'G/C-ending': '#1f77b4', 'A/T-ending': '#ff7f0e'}).to_numpy()
y = np.arange(len(merged))
off_obs = -0.16
off_pred = 0.16
fig, ax = plt.subplots(figsize=(8.2, 12.0), dpi=300)

obs10 = merged['observed_frequency_change_10fold'].to_numpy()
obs100 = merged['observed_frequency_change_100fold'].to_numpy()
pred10 = merged['predicted_frequency_change_10fold'].to_numpy()
pred100 = merged['predicted_frequency_change_100fold'].to_numpy()
ghat = merged['2Ns'].to_numpy()

for i, color in enumerate(colors):
    ax.plot([obs10[i], obs100[i]], [y[i] + off_obs, y[i] + off_obs], color=color, alpha=0.35, linewidth=1.0)
    ax.plot([pred10[i], pred100[i]], [y[i] + off_pred, y[i] + off_pred], color=color, alpha=0.35, linewidth=1.0, linestyle='--')

ax.scatter(obs10, y + off_obs, s=24, facecolors='white', edgecolors=colors, linewidths=1.3, zorder=4)
ax.scatter(obs100, y + off_obs, s=30, c=colors, edgecolors='none', alpha=0.88, zorder=5)
ax.scatter(pred10, y + off_pred, s=26, marker='s', facecolors='white', edgecolors=colors, linewidths=1.3, zorder=4)
ax.scatter(pred100, y + off_pred, s=32, marker='s', c=colors, edgecolors='none', alpha=0.88, zorder=5)

ax.axvline(0, color='0.25', linewidth=1.0)
ax.set_yticks(y)
ax.set_yticklabels(merged['label'], fontsize=8.2)
ax.set_xlabel('Frequency change relative to baseline expression', fontsize=12)
ax.set_ylabel('Codon, sorted by $\\hat{g}$', fontsize=12)
ax.grid(True, axis='x', alpha=0.25)
ax.set_axisbelow(True)
right_limit = max(abs(obs10).max(), abs(obs100).max(), abs(pred10).max(), abs(pred100).max()) * 1.10
ax.set_xlim(-0.15, right_limit)
ax.set_ylim(-0.7, len(merged) - 0.3)

ax_g = ax.twiny()
gpad = max(abs(ghat.min()), abs(ghat.max())) * 1.08
ax_g.set_xlim(-gpad, gpad)
ax_g.axvline(0, color='0.25', linewidth=0.8, alpha=0.45)
ax_g.scatter(ghat, y, s=22, c='black', alpha=0.76, zorder=2)
ax_g.set_xlabel(r'Codon fitness, $\hat{g}$', fontsize=12)
ax_g.set_ylim(ax.get_ylim())
ax_g.tick_params(axis='x', labelsize=9)
ax_g.tick_params(axis='y', left=False, right=False, labelleft=False, labelright=False)
ax_g.grid(False)

legend_items = [
    Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='white', markeredgecolor='black', markersize=6, label='Observed 10-fold'),
    Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='black', markeredgecolor='none', markersize=6, label='Observed 100-fold'),
    Line2D([0], [0], marker='s', linestyle='none', markerfacecolor='white', markeredgecolor='black', markersize=6, label='Predicted 10-fold'),
    Line2D([0], [0], marker='s', linestyle='none', markerfacecolor='black', markeredgecolor='none', markersize=6, label='Predicted 100-fold'),
    Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='black', markeredgecolor='none', markersize=5.5, label=r'$\hat{g}$'),
    Patch(facecolor='#1f77b4', label='G/C-ending codons'),
    Patch(facecolor='#ff7f0e', label='A/T-ending codons'),
]
ax.legend(handles=legend_items, frameon=False, loc='lower right', fontsize=10.5, ncol=1)
fig.tight_layout()
fig.savefig(OUT, dpi=300)
print(OUT)
print(OUT_TSV)

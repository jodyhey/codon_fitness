from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import to_rgb

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

base_colors = merged['ending_class'].map({'G/C-ending': '#1f77b4', 'A/T-ending': '#ff7f0e'}).to_numpy()

def darken(color, amount=0.62):
    rgb = np.array(to_rgb(color))
    return tuple(rgb * amount)

dark_colors = np.array([darken(c) for c in base_colors], dtype=object)
y = np.arange(len(merged))
off_obs = -0.16
off_pred = 0.16
fig, ax = plt.subplots(figsize=(8.2, 12.0), dpi=300)

obs10 = merged['observed_frequency_change_10fold'].to_numpy()
obs100 = merged['observed_frequency_change_100fold'].to_numpy()
pred10 = merged['predicted_frequency_change_10fold'].to_numpy()
pred100 = merged['predicted_frequency_change_100fold'].to_numpy()
ghat = merged['2Ns'].to_numpy()

for i, (color, dark_color) in enumerate(zip(base_colors, dark_colors)):
    # Dashed segment: observed-to-predicted change for a 10-fold expression increase.
    ax.plot(
        [obs10[i], pred10[i]], [y[i], y[i]],
        color=color, alpha=0.75, linewidth=1.6, linestyle=(0, (3.2, 2.4)), zorder=1
    )
    # Solid darker segment: observed-to-predicted change for a 100-fold expression increase.
    ax.plot(
        [obs100[i], pred100[i]], [y[i], y[i]],
        color=dark_color, alpha=0.9, linewidth=2.0, linestyle='-', zorder=2
    )

ax.scatter(obs10, y + off_obs, s=54, marker='o', facecolors='white', edgecolors=base_colors, linewidths=1.8, zorder=5)
ax.scatter(obs100, y + off_obs, s=64, marker='o', c=dark_colors.tolist(), edgecolors='none', alpha=0.96, zorder=6)
ax.scatter(pred10, y + off_pred, s=58, marker='s', facecolors='white', edgecolors=base_colors, linewidths=1.8, zorder=5)
ax.scatter(pred100, y + off_pred, s=68, marker='s', c=dark_colors.tolist(), edgecolors='none', alpha=0.96, zorder=6)

ax.axvline(0, color='0.25', linewidth=1.0)
ax.set_yticks(y)
ax.set_yticklabels(merged['label'], fontsize=8.2)
ax.set_xlabel('Frequency change relative to baseline expression', fontsize=12)
ax.set_ylabel('Codon, sorted by $\\hat{g}$', fontsize=12)
ax.grid(False)
ax.set_axisbelow(True)
right_limit = max(abs(obs10).max(), abs(obs100).max(), abs(pred10).max(), abs(pred100).max()) * 1.10
ax.set_xlim(-0.15, right_limit)
ax.set_ylim(-0.7, len(merged) - 0.3)

ax_g = ax.twiny()
gpad = max(abs(ghat.min()), abs(ghat.max())) * 1.08
ax_g.set_xlim(-gpad, gpad)
ax_g.axvline(0, color='0.25', linewidth=0.8, alpha=0.45)
ax_g.scatter(ghat, y, s=38, c='black', alpha=0.82, zorder=2)
ax_g.set_xlabel(r'Codon fitness, $\hat{g}$', fontsize=12)
ax_g.set_ylim(ax.get_ylim())
ax_g.tick_params(axis='x', labelsize=9)
ax_g.tick_params(axis='y', left=False, right=False, labelleft=False, labelright=False)
ax_g.grid(False)

legend_items = [
    Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='white', markeredgecolor='black', markeredgewidth=1.6, markersize=8, label='Observed'),
    Line2D([0], [0], marker='s', linestyle='none', markerfacecolor='white', markeredgecolor='black', markeredgewidth=1.6, markersize=8, label='Predicted'),
    Line2D([0], [0], color='black', linewidth=1.8, linestyle=(0, (3.2, 2.4)), label='10-fold'),
    Line2D([0], [0], color='black', linewidth=2.2, linestyle='-', label='100-fold'),
    Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='black', markeredgecolor='none', markersize=7.5, label=r'$\hat{g}$'),
    Patch(facecolor='#1f77b4', label='G/C-ending codons'),
    Patch(facecolor='#ff7f0e', label='A/T-ending codons'),
]
ax.legend(handles=legend_items, frameon=False, loc='lower right', fontsize=11.5, ncol=1)
fig.tight_layout()
fig.savefig(OUT, dpi=300)
print(OUT)
print(OUT_TSV)

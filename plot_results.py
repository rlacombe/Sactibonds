#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

# Original Data
methods = ["Boltz-1", "AlphaFold 3", "RoseTTAfold2", "AlphaFold 2", "OmegaFold", "ESMFold"]

# Sacti-RMSD data (mean ± std)
sacti_rmsd_means = np.array([8.33, 13.76, 7.96, 11.72, 9.16, 15.09])
sacti_rmsd_stds  = np.array([2.48, 7.21, 2.60, 6.51, 2.61, 6.67])  # no longer used

# GDT_TS data (means)
gdt_ts_means = np.array([16.46, 13.54, 17.59, 12.29, 8.54, 5.00])

# 1. Sort Methods by Ascending GDT_TS
sorted_indices = np.argsort(gdt_ts_means)  # indices that would sort GDT_TS ascending

methods_sorted          = [methods[i] for i in sorted_indices]
gdt_ts_means_sorted     = gdt_ts_means[sorted_indices]
sacti_rmsd_means_sorted = sacti_rmsd_means[sorted_indices]

# 2. Plot Setup
fig, ax1 = plt.subplots(figsize=(8, 5))

x_positions = np.arange(len(methods_sorted))
bar_width   = 0.6

# 3. Plot GDT_TS as Blue Bars (Sorted)
bars = ax1.bar(
    x_positions,
    gdt_ts_means_sorted,
    color='blue',
    alpha=0.7,
    width=bar_width,
    label='GDT_TS (%)'
)

# Configure the first axis (for GDT_TS)
ax1.set_xlabel('Protein Models')
ax1.set_ylabel('GDT_TS (%)', color='blue')
ax1.set_xticks(x_positions)
ax1.set_xticklabels(methods_sorted, rotation=45, ha='right')
ax1.tick_params(axis='y', labelcolor='blue')

# Ensure the GDT_TS axis starts at 0 (avoid "chart crime")
max_gdt_ts = gdt_ts_means_sorted.max()
ax1.set_ylim([0, max_gdt_ts + 2])  # +2 for a little headroom

# 4. Create a twin axis for Sacti-RMSD
ax2 = ax1.twinx()

# 5. Plot Sacti-RMSD as Red Dotted Line (No Error Bars)
ax2.plot(
    x_positions,
    sacti_rmsd_means_sorted,
    color='red',
    marker='o',
    linestyle=':',
    linewidth=2,
    label='Sacti-RMSD (Å)'
)

# Configure the second axis (for Sacti-RMSD)
ax2.set_ylabel('Average Sacti-RMSD (Å)', color='red')
ax2.tick_params(axis='y', labelcolor='red')

# Ensure the Sacti-RMSD axis also starts at 0
max_sacti_rmsd = sacti_rmsd_means_sorted.max()
ax2.set_ylim([0, max_sacti_rmsd + 2])  # +2 for some headroom

plt.title('Sactipeptide crosslink structure prediction')
fig.tight_layout()

# Uncomment if you want to save the figure:
plt.savefig("models_evaluation.png", dpi=300)

plt.show()

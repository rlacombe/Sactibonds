#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import json
import argparse

# Color definitions
BAR_COLOR = '#3D3DDC'     # GDT_TS bars
LINE_COLOR = '#CB2C1E'    # RMSD lines
TITLE_COLOR = 'black'     # Title

def parse_args():
    parser = argparse.ArgumentParser(description='Plot structure prediction results')
    parser.add_argument('--gdt', action='store_true', help='Show GDT_TS bars')
    parser.add_argument('--rmsd', action='store_true', help='Show RMSD lines')
    return parser.parse_args()

# Load results and models data
with open('results.json', 'r') as f:
    results = json.load(f)
    
with open('models.json', 'r') as f:
    model_data = json.load(f)
    model_names = model_data['models']  # Get the nested dictionary

# Get list of models
methods = list(results.keys())

# Extract data
gdt_known = np.array([results[m]['gdt_ts']['pdb'] for m in methods])
gdt_no = np.array([results[m]['gdt_ts']['no_pdb'] for m in methods])
gdt_all = np.array([results[m]['gdt_ts']['all'] for m in methods])
rmsd_all = np.array([results[m]['rmsd']['all'] for m in methods])
rmsd_known = np.array([results[m]['rmsd']['pdb'] for m in methods])
rmsd_no = np.array([results[m]['rmsd']['no_pdb'] for m in methods])


# Map internal names to display names
methods = [model_names.get(m, m) for m in methods]

# Get command line arguments
args = parse_args()

# Plotting
fig, ax1 = plt.subplots(figsize=(8, 6))

x = np.arange(len(methods))
width = 0.25 # Thinner bars

# Create bars with gradients if GDT_TS is enabled
if args.gdt:

    # Sort all methods by increasing gdt_no
    sort_idx = np.argsort(gdt_no)  # Ascending order
    methods = [methods[i] for i in sort_idx]
    gdt_known = gdt_known[sort_idx]
    gdt_no = gdt_no[sort_idx]
    gdt_all = gdt_all[sort_idx]
    rmsd_all = rmsd_all[sort_idx]
    rmsd_known = rmsd_known[sort_idx]
    rmsd_no = rmsd_no[sort_idx]

    # Plot GDT_TS bars
    rects2 = ax1.bar(x - width/2, gdt_known, 0.9*width, label='GDT_TS: known', 
                     color=BAR_COLOR, alpha=0.5)
    rects1 = ax1.bar(x + width/2, gdt_no, 0.9*width, label='GDT_TS: unknown', 
                     color=BAR_COLOR, alpha=1.0)
    
    # Configure first axis (GDT_TS)
    ax1.set_ylabel('GDT_TS (%)', color=BAR_COLOR, fontsize=12)
    ax1.tick_params(axis='y', labelcolor=BAR_COLOR)
    ax1.legend(loc='upper left', fontsize=11)
    ax1.set_ylim(0, 25)
    filename = 'results_gdt.png'


# Create second axis for RMSD if enabled
if args.rmsd:
    if args.gdt:
        ax2 = ax1.twinx()
        filename = 'results_joint.png'

    else:
        ax2 = ax1

        # Sort all methods by decreasing rmsd_no
        sort_idx = np.argsort(rmsd_no)[::-1]  # Descending order
        methods = [methods[i] for i in sort_idx]
        gdt_known = gdt_known[sort_idx]
        gdt_no = gdt_no[sort_idx]
        gdt_all = gdt_all[sort_idx]
        rmsd_all = rmsd_all[sort_idx]
        rmsd_known = rmsd_known[sort_idx]
        rmsd_no = rmsd_no[sort_idx]

        filename = 'results_rmsd.png'

    ax2.plot(x, rmsd_known, color=LINE_COLOR, linestyle='dotted', linewidth=3, marker='s',
             label='RMSD: known', alpha=0.5)
    ax2.plot(x, rmsd_no, color=LINE_COLOR, linestyle='dashed', linewidth=3, marker='o',
             label='RMSD: unknown', alpha=1.0)
    ax2.set_ylabel('RMSD (Å)', color=LINE_COLOR, fontsize=12)
    ax2.tick_params(axis='y', labelcolor=LINE_COLOR)
    ax2.set_ylim(0, 25)
    ax2.legend(loc='upper right', fontsize=11)

# Common configuration
# ax1.set_title('Crosslinked Structure Prediction', fontsize=14, color=TITLE_COLOR)
ax1.set_xticks(x)
ax1.set_xticklabels(methods, rotation=30, ha='right', fontsize=12)

# Add some padding to prevent label cutoff
plt.tight_layout()

# Save plot
plt.savefig(filename, dpi=300, bbox_inches='tight')
plt.show()

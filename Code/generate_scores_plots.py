import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# ========================== CONFIG ==========================
score_file = "/Users/shantanu.yadav19/Downloads/MTP2Offline/Data/final_results/af2_scores.sc"
output_dir = "/Users/shantanu.yadav19/Downloads/MTP2Offline/Data/final_results"
# ===========================================================

# Tiered thresholds for pae_interaction
# Lower is better. Typically <10 as "likely to work experimentally".
# We add softer tiers so the plots still communicate signal when nothing crosses 10.
PAE_TIERS = [
    ("excellent", 0, 10, "#2ca02c"),   # green
    ("good",      10, 15, "#f1c40f"),  # yellow
    ("marginal",  15, 20, "#e67e22"),  # orange
    ("poor",      20, np.inf, "#d62728"),  # red
]
TIER_ORDER = [t[0] for t in PAE_TIERS]
TIER_COLORS = {name: color for name, _, _, color in PAE_TIERS}

# Secondary thresholds (reference lines only; not used for coloring)
PLDDT_BINDER_GOOD = 80
RMSD_GOOD = 2.0

os.makedirs(output_dir, exist_ok=True)

# ==================== ROBUST SCORE FILE PARSING ====================
with open(score_file, 'r') as f:
    lines = f.readlines()

header_line = None
data_lines = []
for line in lines:
    s = line.strip()
    if s.startswith("SCORE:"):
        if header_line is None:
            header_line = s
        else:
            data_lines.append(s)
    elif s and not s.startswith("#") and not s.startswith("SEQUENCE:"):
        data_lines.append(s)

columns = header_line.split()[1:]  # drop the leading "SCORE:" token

# Each data line starts with "SCORE:" too; strip and split
rows = []
for line in data_lines:
    tokens = line.split()
    if tokens and tokens[0] == "SCORE:":
        tokens = tokens[1:]
    if len(tokens) == len(columns):
        rows.append(tokens)

df = pd.DataFrame(rows, columns=columns)

numeric_cols = ['binder_aligned_rmsd', 'pae_binder', 'pae_interaction', 'pae_target',
                'plddt_binder', 'plddt_target', 'plddt_total',
                'target_aligned_rmsd', 'time']
for col in numeric_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

df = df.dropna(subset=['pae_interaction']).reset_index(drop=True)

print(f"Loaded {len(df)} designs")
print("Columns found:", df.columns.tolist())


# ==================== ASSIGN TIERS ====================
def assign_tier(pae):
    for name, lo, hi, _ in PAE_TIERS:
        if lo <= pae < hi:
            return name
    return "poor"

df['tier'] = df['pae_interaction'].apply(assign_tier)
df['color'] = df['tier'].map(TIER_COLORS)

# Summary
print("\nTier breakdown (pae_interaction):")
tier_counts = df['tier'].value_counts().reindex(TIER_ORDER, fill_value=0)
for tier in TIER_ORDER:
    n = tier_counts[tier]
    pct = 100 * n / len(df) if len(df) else 0
    print(f"  {tier:10s} ({TIER_COLORS[tier]}): {n:4d}  ({pct:5.1f}%)")

key_cols = ['pae_interaction', 'plddt_total', 'plddt_binder', 'plddt_target',
            'binder_aligned_rmsd', 'target_aligned_rmsd']
print("\nStatistics for key metrics:\n")
print(df[key_cols].describe().round(2))


# ==================== PLOTS ====================
sns.set_style("whitegrid")
fig, axes = plt.subplots(2, 3, figsize=(20, 11))

# Legend handles reused across panels
tier_handles = [
    plt.Line2D([0], [0], marker='o', linestyle='', color=TIER_COLORS[name],
               markersize=9, label=f"{name} ({lo}≤pae<{hi if hi != np.inf else '∞'})")
    for name, lo, hi, _ in PAE_TIERS
]


# --- Panel 1: pae_interaction histogram with tier bands ---
ax = axes[0, 0]
bins = np.linspace(0, max(40, df['pae_interaction'].max() + 2), 40)
# Shade tier bands behind the histogram
for name, lo, hi, color in PAE_TIERS:
    hi_plot = hi if hi != np.inf else ax.get_xlim()[1]
    ax.axvspan(lo, hi_plot if hi != np.inf else bins[-1], alpha=0.12, color=color)
ax.hist(df['pae_interaction'], bins=bins, color='#34495e', edgecolor='white')
for name, lo, hi, color in PAE_TIERS:
    if lo > 0:
        ax.axvline(lo, color=color, linestyle='--', alpha=0.7, linewidth=1.2)
ax.set_xlabel('pae_interaction (lower is better)')
ax.set_ylabel('count')
ax.set_title('pae_interaction distribution with tier bands')
ax.legend(handles=tier_handles, loc='upper right', fontsize=8, framealpha=0.9)


# --- Panel 2: plddt_total histogram ---
ax = axes[0, 1]
sns.histplot(df['plddt_total'], kde=True, ax=ax, color='#2980b9', edgecolor='white')
ax.axvline(PLDDT_BINDER_GOOD, color='green', linestyle='--',
           label=f'plddt > {PLDDT_BINDER_GOOD} (confident fold)')
ax.set_xlabel('plddt_total')
ax.set_title('plddt_total distribution')
ax.legend(fontsize=9)


# --- Panel 3: plddt_binder vs plddt_target, colored by tier ---
ax = axes[0, 2]
for tier in TIER_ORDER:
    sub = df[df['tier'] == tier]
    if len(sub):
        ax.scatter(sub['plddt_binder'], sub['plddt_target'],
                   c=TIER_COLORS[tier], alpha=0.7, s=50,
                   edgecolor='white', linewidth=0.5, label=f"{tier} (n={len(sub)})")
ax.axvline(PLDDT_BINDER_GOOD, color='gray', linestyle=':', alpha=0.5)
ax.axhline(PLDDT_BINDER_GOOD, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('plddt_binder')
ax.set_ylabel('plddt_target')
ax.set_title('Per-chain fold confidence, colored by interface tier')
ax.legend(fontsize=8, loc='lower right')


# --- Panel 4: pae_interaction vs plddt_binder (diagnostic) ---
ax = axes[1, 0]
for tier in TIER_ORDER:
    sub = df[df['tier'] == tier]
    if len(sub):
        ax.scatter(sub['pae_interaction'], sub['plddt_binder'],
                   c=TIER_COLORS[tier], alpha=0.75, s=50,
                   edgecolor='white', linewidth=0.5)
# Reference lines
for _, lo, hi, color in PAE_TIERS:
    if lo > 0:
        ax.axvline(lo, color=color, linestyle='--', alpha=0.4, linewidth=1)
ax.axhline(PLDDT_BINDER_GOOD, color='gray', linestyle=':', alpha=0.5,
           label=f'plddt_binder = {PLDDT_BINDER_GOOD}')
ax.set_xlabel('pae_interaction')
ax.set_ylabel('plddt_binder')
ax.set_title('Interface vs binder-fold confidence\n'
             '(top-left = good binder AND good interface)')
ax.legend(handles=tier_handles + [ax.get_legend_handles_labels()[0][-1]] if ax.get_legend_handles_labels()[0] else tier_handles,
          fontsize=8, loc='lower right')


# --- Panel 5: binder_aligned_rmsd vs target_aligned_rmsd ---
ax = axes[1, 1]
for tier in TIER_ORDER:
    sub = df[df['tier'] == tier]
    if len(sub):
        ax.scatter(sub['binder_aligned_rmsd'], sub['target_aligned_rmsd'],
                   c=TIER_COLORS[tier], alpha=0.7, s=50,
                   edgecolor='white', linewidth=0.5)
ax.axvline(RMSD_GOOD, color='gray', linestyle=':', alpha=0.5,
           label=f'binder RMSD = {RMSD_GOOD} Å')
ax.set_xlabel('binder_aligned_rmsd (Å)')
ax.set_ylabel('target_aligned_rmsd (Å)')
ax.set_title('Backbone RMSD: binder vs target\n'
             '(low binder RMSD = AF2 agrees with RFdiffusion design)')
ax.legend(handles=tier_handles, fontsize=8, loc='upper right')


# --- Panel 6: pae_interaction vs binder_aligned_rmsd (key diagnostic) ---
ax = axes[1, 2]
for tier in TIER_ORDER:
    sub = df[df['tier'] == tier]
    if len(sub):
        ax.scatter(sub['pae_interaction'], sub['binder_aligned_rmsd'],
                   c=TIER_COLORS[tier], alpha=0.75, s=50,
                   edgecolor='white', linewidth=0.5)
for _, lo, hi, color in PAE_TIERS:
    if lo > 0:
        ax.axvline(lo, color=color, linestyle='--', alpha=0.4, linewidth=1)
ax.axhline(RMSD_GOOD, color='gray', linestyle=':', alpha=0.5,
           label=f'RMSD = {RMSD_GOOD} Å')
ax.set_xlabel('pae_interaction')
ax.set_ylabel('binder_aligned_rmsd (Å)')
ax.set_title('Interface confidence vs design fidelity\n'
             '(bottom-left = good interface AND AF2 matches design)')
ax.legend(handles=tier_handles, fontsize=8, loc='upper right')


plt.tight_layout()
out_png = os.path.join(output_dir, 'af2_score_analysis.png')
plt.savefig(out_png, dpi=200, bbox_inches='tight')
print(f"\nSaved main figure to {out_png}")
plt.close()


# ==================== FILTERING PER TIER ====================
print("\nTier-based filtering:")
for name, lo, hi, _ in PAE_TIERS:
    subset = df[(df['pae_interaction'] >= lo) & (df['pae_interaction'] < hi)].copy()
    subset = subset.sort_values(by=['pae_interaction', 'plddt_binder'],
                                 ascending=[True, False])
    out_csv = os.path.join(output_dir, f'designs_{name}.csv')
    subset.to_csv(out_csv, index=False)
    print(f"  {name:10s}: {len(subset):4d} designs -> {out_csv}")


# Composite "best candidates" list: rank everyone, regardless of tier
# by pae_interaction first, then plddt_binder, then binder_aligned_rmsd
composite = df.sort_values(
    by=['pae_interaction', 'plddt_binder', 'binder_aligned_rmsd'],
    ascending=[True, False, True]
).copy()

ranked_csv = os.path.join(output_dir, 'designs_all_ranked.csv')
composite.to_csv(ranked_csv, index=False)
print(f"\nFull ranked list saved to {ranked_csv}")


# ==================== TOP DESIGNS ====================
show_cols = ['description', 'tier', 'pae_interaction', 'plddt_binder',
             'plddt_target', 'binder_aligned_rmsd']
show_cols = [c for c in show_cols if c in composite.columns]

print("\nTop 15 designs overall (ranked by pae_interaction, then plddt_binder):")
print(composite[show_cols].head(15).to_string(index=False))

# If nothing crosses 10, tell the user where the "best we have" actually sits
best_pae = composite['pae_interaction'].min()
print(f"\nBest pae_interaction in this run: {best_pae:.2f}")
if best_pae >= 10:
    print(f"  -> No designs crossed the <10 threshold.")
    print(f"  -> The best {min(10, len(composite))} designs are still your most informative structures "
          f"to inspect visually (look at them in PyMOL to understand the failure mode).")
else:
    n_pass = (df['pae_interaction'] < 10).sum()
    print(f"  -> {n_pass} designs pass the <10 threshold.")


# ==================== PER-TIER SUMMARY STATS ====================
print("\nPer-tier summary of secondary metrics:")
summary = df.groupby('tier')[['plddt_binder', 'plddt_target', 'binder_aligned_rmsd']].agg(
    ['count', 'mean', 'median']
).round(2)
# Reindex to tier order
summary = summary.reindex([t for t in TIER_ORDER if t in summary.index])
print(summary)

summary_csv = os.path.join(output_dir, 'per_tier_summary.csv')
summary.to_csv(summary_csv)
print(f"\nPer-tier summary saved to {summary_csv}")
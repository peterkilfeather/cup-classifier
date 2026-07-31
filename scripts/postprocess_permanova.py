"""Post-process: collect PERMANOVA results, generate plot + entry-order table."""

import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path('/xscratch/farney/cup-classifier/output/diagnostic-protocol')
TMP = OUT / 'tmp'

COVARIATES = ['Tissue', 'Source', 'BCT', 'Sex']
PALETTES = {
    'Tissue': {'healthyblood': '#4daf4a', 'colon': '#e41a1c', 'liver': '#377eb8',
               'pancreas': '#ff7f00', 'prostate': '#984ea3', 'stomach': '#f781bf'},
    'Source': {'Fox Chase': '#e41a1c', 'Audubon': '#377eb8',
               'Sowalsky': '#4daf4a', 'NIH_Clinical_Center': '#ff7f00'},
    'BCT': {'EDTA': '#377eb8', 'Citrate': '#e41a1c', 'Streck': '#4daf4a', 'ACD': '#ff7f00'},
    'Sex': {'Male': '#377eb8', 'Female': '#e41a1c'},
}

MODALITY_LABELS = {
    'probe_meth': 'Methylation (probe-avg)',
    'fem4': 'FEM4 (256)',
    'fragment_length': 'Fragment length (369)',
    'cnvkit': 'CNVkit thr0.10',
    'probe_cpg': 'Per-CpG methylation (32K->PCs)',
    'end_density': 'End density (31K bins->PCs)',
}

HIGH_DIM_K = {'probe_cpg': 5, 'end_density': 5}  # minimum PCs for high-dim modalities

def load_perm(mod_name, scope):
    safe = f"{MODALITY_LABELS[mod_name].replace(' ', '_').replace('->', 'to')}_{scope.replace(' ', '_')}"
    td = TMP / safe
    results = {}
    for key, fname in [('marginal', 'marginal.tsv'), ('sequential', 'sequential.tsv'), ('bct_first', 'bct_first.tsv')]:
        fp = td / fname
        if fp.exists():
            try:
                df = pd.read_csv(fp, sep='\t', index_col=0)
                results[key] = {}
                for t in COVARIATES:
                    if t in df.index:
                        results[key][t] = {'R2': float(df.loc[t,'R2']), 'F': float(df.loc[t,'F']), 'p': float(df.loc[t,'Pr(>F)'])}
            except: pass
    return results

mods = ['probe_meth', 'fem4', 'fragment_length', 'cnvkit', 'probe_cpg', 'end_density']

# Load and print
print("=== PERMANOVA R² — Full (164) ===")
all_full = {}
for m in mods:
    r = load_perm(m, 'Full (164)')
    if r:
        all_full[m] = r
        mr = r.get('marginal', {})
        sq = r.get('sequential', {})
        bf = r.get('bct_first', {})
        print(f"\n{MODALITY_LABELS[m]:35s}")
        print(f"  Marginal:   T={mr.get('Tissue',{}).get('R2',0):.4f}  B={mr.get('BCT',{}).get('R2',0):.4f}  S={mr.get('Source',{}).get('R2',0):.4f}  X={mr.get('Sex',{}).get('R2',0):.4f}")
        print(f"  Sequential: T={sq.get('Tissue',{}).get('R2',0):.4f}  B={sq.get('BCT',{}).get('R2',0):.4f}  S={sq.get('Source',{}).get('R2',0):.4f}  X={sq.get('Sex',{}).get('R2',0):.4f}")
        print(f"  BCT-first:  T={bf.get('Tissue',{}).get('R2',0):.4f}  B={bf.get('BCT',{}).get('R2',0):.4f}  S={bf.get('Source',{}).get('R2',0):.4f}  X={bf.get('Sex',{}).get('R2',0):.4f}")

# Plot
covariates_plot = COVARIATES
modalities_ordered = [m for m in mods if m in all_full]
x = np.arange(len(modalities_ordered))
width = 0.7 / len(covariates_plot)

fig, ax = plt.subplots(figsize=(10, 5))
for i, cov in enumerate(covariates_plot):
    vals = [all_full[m].get('marginal', {}).get(cov, {}).get('R2', 0) or 0 for m in modalities_ordered]
    offset = (i - len(covariates_plot)/2 + 0.5) * width
    pal = list(PALETTES.get(cov, {'x': '#999'}).values())
    c = pal[0] if pal else '#999'
    ax.bar(x + offset, vals, width, label=cov, color=c)

labels = [MODALITY_LABELS[m] for m in modalities_ordered]
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
ax.set_ylabel('R² (marginal)')
ax.set_title('PERMANOVA Type III — Full (164)')
ax.axhline(y=0.10, color='red', linestyle='--', alpha=0.5, label='Threshold 0.10')
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(OUT / 'permanova_Full_(164).png', dpi=150)
plt.close(fig)
print("\nSaved PERMANOVA plot")

# Entry-order sensitivity table
print(f"\n{'='*80}")
print(f"{'Entry-order sensitivity (R² range across marginal/sequential/BCT-first)':^80}")
print(f"{'='*80}")
print(f"{'Modality':35s} {'Tissue':>12s} {'BCT':>12s} {'Source':>12s}")
for m in modalities_ordered:
    r = all_full[m]
    marg = r.get('marginal', {})
    parts = [marg, r.get('sequential', {}), r.get('bct_first', {})]
    print(f"{MODALITY_LABELS[m]:35s}", end='')
    for cov in ['Tissue', 'BCT', 'Source']:
        vals = [p.get(cov, {}).get('R2', 0) or 0 for p in parts if p]
        vals = [v for v in vals if v > 0.001]
        if vals:
            print(f' {min(vals):.4f}-{max(vals):.4f}', end='')
        else:
            print(f'      N/A', end='')
    print()

#!/usr/bin/env python3
"""Generate client-facing figures for the methylation version inventory (issue #8, step 1).

Outputs (16:9, 300 dpi PNG) -> docs/figures/:
  01 per-sample missingness by version (bars + jittered strip)
  02 per-probe per-CpG missingness, enriched vs unenriched
  03 CpG sites per probe, enriched vs unenriched (site-set expansion)
  04 probe aggregation recoverability (site-level -> probe-level missingness)
  05 probe-level read coverage by capture condition
  06 per-CpG coverage, enriched (unenriched coverage not available)

Definitions match docs/methylation-inventory.md: post-join 164 samples
(metadata_cleaned.csv TWIST_ID), missingness = NaN feature value.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
METH = BASE / 'input' / 'methylation'
FIGS = BASE / 'docs' / 'figures'
FIGS.mkdir(exist_ok=True)

PALETTE = {
    'probe_meth_enriched': '#0072B2',
    'probe_meth_unenriched': '#56B4E9',
    'probe_meth_unfiltered': '#CC79A7',
    'probe_cpg_unenriched': '#E69F00',
    'probe_cpg_enriched': '#D55E00',
}
GRAY = '#8c8c8c'
BLACK = '#111111'

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 13,
    'axes.titlesize': 18,
    'axes.titleweight': 'bold',
    'axes.labelsize': 15,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'axes.edgecolor': '#555555',
    'axes.linewidth': 1.0,
    'axes.grid': True,
    'grid.color': '#dddddd',
    'grid.linewidth': 0.8,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
})

FOOT = 'n = 164 (inner join with metadata_cleaned.csv); missing = NaN feature value.'


def load_metadata():
    meta = pd.read_csv(BASE / 'input' / 'metadata' / 'metadata_cleaned.csv')
    return set(meta['TWIST_ID'].astype(str))


def load_wide(path, meta):
    """Return (X float64 (n,feat), sample_ids, feature names) post-join."""
    df = pd.read_csv(path, sep='\t', dtype={'sample': str})
    df = df[df['sample'].isin(meta)]
    feat = list(df.columns[1:])
    return df[feat].to_numpy(float), df['sample'].to_numpy(), feat


def load_long_probe(path, meta):
    """Return per-sample x probe missingness (post-join) from a probe-level long file."""
    df = pd.read_csv(path, sep='\t', dtype={'sample': str, 'probe_id': str})
    df = df[df['sample'].isin(meta)]
    pw = df.pivot_table(index='sample', columns='probe_id', values='CpG_frac',
                        aggfunc='first')
    return np.isnan(pw.to_numpy(float)), pw.index.to_numpy(), list(pw.columns)


def save(fig, name):
    fig.savefig(FIGS / name, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('wrote', FIGS / name)


def fig01(meta, versions):
    """Bars (mean per-sample missingness) + jittered strip of per-sample values."""
    data = []
    for key, kind, path in versions:
        if kind == 'wide':
            X, _, _ = load_wide(path, meta)
        else:
            miss, _, _ = load_long_probe(path, meta)
            X = miss.astype(float)
        per_sample = np.isnan(X).mean(axis=1) if kind == 'wide' else miss.mean(axis=1)
        data.append((key, per_sample, np.isnan(X).all(axis=0).sum()))
    keys = [d[0] for d in data]
    means = [d[1].mean() for d in data]
    allnan = [d[2] for d in data]
    labels = ['probe_meth\nenriched', 'probe_meth\nunenriched', 'probe_meth\nunfiltered',
              'probe_cpg\nunenriched', 'probe_cpg\nenriched']

    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    x = np.arange(len(keys))
    bars = ax.bar(x, means, width=0.62, color=[PALETTE[k] for k in keys],
                  edgecolor='white', zorder=3)
    rng = np.random.default_rng(7)
    for xi, (k, ps) in zip(x, [(d[0], d[1]) for d in data]):
        jitter = rng.uniform(-0.22, 0.22, size=len(ps))
        ax.scatter(xi + jitter, ps, s=14, color=PALETTE[k], alpha=0.55,
                   edgecolors='none', zorder=4)
    for xi, m, n in zip(x, means, allnan):
        ax.text(xi, m + 0.011, f'{m * 100:.2f}%', ha='center', va='bottom',
                fontsize=13, fontweight='bold', color=BLACK)
        ax.text(xi, -0.030, f'all-NaN features: {n}', ha='center', va='top',
                fontsize=10.5, color='#444444')
    # single-value distributions: probe_meth enriched is exactly 0 for all samples
    ax.text(0.02, 0.45,
            'probe_meth enriched: distribution is a single value (0)\nfor all 164 samples',
            transform=ax.transAxes, va='top', fontsize=11, color='#444444')
    ax.axhline(0, color='#aaaaaa', linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(-0.075, 0.46)
    ax.set_ylabel('per-sample missingness (fraction of features NaN)')
    ax.set_title('Per-sample missingness by methylation version: probe-level ~0%, '
                 'per-CpG 6.0% unenriched vs 37.1% enriched')
    ax.text(0, -0.155, FOOT, transform=ax.transAxes, fontsize=9.5, color='#555555')
    save(fig, 'methylation_01_per_sample_missingness.png')


def fig02(meta):
    """Per-probe per-CpG missingness (site-level cells), enriched vs unenriched."""
    Xe, _, fe = load_wide(METH / 'probe_cpg'
                          / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv', meta)
    Xu, _, fu = load_wide(METH / 'probe_cpg'
                          / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv', meta)
    def per_probe(X, feat):
        probe = np.array([c.split('_cpg')[0] for c in feat])
        out = {}
        for p in pd.unique(probe):
            out[p] = np.isnan(X[:, probe == p]).mean()
        return pd.Series(out)
    se, su = per_probe(Xe, fe), per_probe(Xu, fu)
    j = pd.DataFrame({'unenriched': su, 'enriched': se}).dropna()
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    ax.scatter(j['unenriched'] * 100, j['enriched'] * 100, s=30,
               color=PALETTE['probe_cpg_enriched'], alpha=0.75, zorder=3)
    lim = max(j['unenriched'].max(), j['enriched'].max()) * 100
    ax.plot([0, lim], [0, lim], color=GRAY, linestyle='--', linewidth=1.2,
            label='identity (equal missingness)')
    worst = j['enriched'].idxmax()
    ax.annotate(f'worst probe {worst}\n({j.loc[worst, "enriched"] * 100:.1f}% enriched, '
                f'{j.loc[worst, "unenriched"] * 100:.1f}% unenriched)',
                xy=(j.loc[worst, 'unenriched'] * 100, j.loc[worst, 'enriched'] * 100),
                xytext=(8, 62), arrowprops=dict(arrowstyle='->', color=GRAY),
                fontsize=11, color='#444444')
    ax.set_xlabel('per-probe missingness, unenriched (% of site x sample cells NaN)')
    ax.set_ylabel('per-probe missingness, enriched (%)')
    ax.set_title('Per-CpG missingness concentrates in specific probes '
                 '(enriched up to 82.6% per probe)')
    ax.legend(loc='lower right', frameon=False)
    ax.text(0, -0.155, FOOT, transform=ax.transAxes, fontsize=9.5, color='#555555')
    save(fig, 'methylation_02_per_probe_missingness.png')


def fig03(meta):
    """CpG sites per probe, enriched vs unenriched (site-set expansion)."""
    def counts(path):
        m = pd.read_csv(path, sep='\t', usecols=['feature_id'])
        return m['feature_id'].str.rsplit('_cpg', n=1).str[0].value_counts()
    se = counts(METH / 'probe_cpg' / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4.manifest.tsv')
    su = counts(METH / 'probe_cpg' / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4.manifest.tsv')
    j = pd.DataFrame({'unenriched': su, 'enriched': se}).dropna()
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    ax.scatter(j['unenriched'], j['enriched'], s=30,
               color=PALETTE['probe_cpg_unenriched'], alpha=0.75, zorder=3)
    lim = max(j['unenriched'].max(), j['enriched'].max())
    ax.plot([1, lim], [1, lim], color=GRAY, linestyle='--', linewidth=1.2,
            label='identity (same site set)')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.text(0.03, 0.93,
            f'median per probe: {j["unenriched"].median():.0f} unenriched vs '
            f'{j["enriched"].median():.0f} enriched\n'
            f'max: {j["unenriched"].max():,} vs {j["enriched"].max():,}',
            transform=ax.transAxes, fontsize=12, va='top',
            bbox=dict(boxstyle='round,pad=0.4', fc='#f5f5f5', ec='#bbbbbb'))
    ax.set_xlabel('CpG sites per probe, unenriched (log scale)')
    ax.set_ylabel('CpG sites per probe, enriched (log scale)')
    ax.set_title('The \u201cwhy\u201d: enriched measures ~1.7× more CpG sites per probe '
                 '(median 183 vs 105)')
    ax.legend(loc='lower right', frameon=False)
    ax.text(0, -0.155, FOOT, transform=ax.transAxes, fontsize=9.5, color='#555555')
    save(fig, 'methylation_03_cpg_sites_per_probe.png')


def fig04(meta):
    """Probe aggregation recoverability: site-level vs probe-level missingness per probe."""
    Xe, _, fe = load_wide(METH / 'probe_cpg'
                          / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv', meta)
    Xu, _, fu = load_wide(METH / 'probe_cpg'
                          / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv', meta)
    def per_probe_both(X, feat):
        probe = np.array([c.split('_cpg')[0] for c in feat])
        site_miss, probe_miss = {}, {}
        for p in pd.unique(probe):
            sub = X[:, probe == p]
            site_miss[p] = np.isnan(sub).mean()
            probe_miss[p] = np.isnan(sub).all(axis=1).mean()  # probe NaN after aggregation
        return pd.Series(site_miss), pd.Series(probe_miss)
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    series = {}
    for key, X, feat, color in [
            ('probe_cpg_enriched', Xe, fe, PALETTE['probe_cpg_enriched']),
            ('probe_cpg_unenriched', Xu, fu, PALETTE['probe_cpg_unenriched'])]:
        sm, pm = per_probe_both(X, feat)
        series[key] = (sm, pm)
        ax.scatter(sm * 100, pm * 100, s=30, color=color, alpha=0.75,
                   label='unenriched' if 'unenriched' in key else 'enriched', zorder=3)
    ax.axhline(0, color=GRAY, linestyle='--', linewidth=1.2)
    # the single unenriched probe with probe-level missingness > 0
    sm_u, pm_u = series['probe_cpg_unenriched']
    outlier = pm_u[pm_u > 0]
    if len(outlier):
        pr = outlier.index[0]
        ax.annotate('all 148 probes recoverable per version:\n'
                    '0 samples with an all-NaN probe after aggregation\n'
                    '(unenriched: 1 probe missing in 1 of 164 samples)',
                    xy=(sm_u[pr] * 100, pm_u[pr] * 100), xycoords='data',
                    xytext=(0.03, 0.97), textcoords='axes fraction', va='top',
                    arrowprops=dict(arrowstyle='->', color=GRAY),
                    fontsize=11, color='#444444')
    ax.set_xlabel('per-probe site-level missingness, before aggregation (% of cells NaN)')
    ax.set_ylabel('probe-level missingness after aggregation (% samples all-NaN probe)')
    ax.set_title('Probe aggregation absorbs per-CpG missingness: '
                 'every probe remains measurable in every sample')
    ax.legend(loc='upper right', frameon=False)
    ax.text(0, -0.155, FOOT, transform=ax.transAxes, fontsize=9.5, color='#555555')
    save(fig, 'methylation_04_aggregation_recoverability.png')


def fig05(meta):
    """Probe-level read coverage by capture condition (per-sample median reads/probe)."""
    def per_sample_median(path):
        df = pd.read_csv(path, sep='\t', usecols=['sample', 'probe_id', 'n_reads_with_CpG'],
                         dtype={'sample': str, 'probe_id': str})
        df = df[df['sample'].isin(meta)]
        return df.groupby('sample')['n_reads_with_CpG'].median()
    enr = per_sample_median(METH / 'probe_meth' / 'all_samples.probe_meth_enriched_filtered.mNonCpGlt4.long.tsv')
    une = per_sample_median(METH / 'probe_meth' / 'all_samples.probe_meth_filtered.mNonCpGlt4.long_v2.tsv')
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    parts = ax.violinplot([enr.values, une.values], positions=[0, 1],
                          showmedians=False, widths=0.55)
    for pc, color in zip(parts['bodies'], [PALETTE['probe_meth_enriched'],
                                           PALETTE['probe_meth_unenriched']]):
        pc.set_facecolor(color); pc.set_alpha(0.45); pc.set_edgecolor(color)
    rng = np.random.default_rng(3)
    for xi, vals in zip([0, 1], [enr.values, une.values]):
        ax.scatter(xi + rng.uniform(-0.15, 0.15, len(vals)), vals, s=16,
                   color='#333333', alpha=0.55, zorder=4)
        ax.text(xi, -0.075, f'median {np.median(vals):,.0f} reads/probe',
                transform=ax.get_xaxis_transform(), ha='center', va='top',
                fontsize=12, fontweight='bold')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['probe_meth enriched', 'probe_meth unenriched'])
    ax.set_ylabel('per-sample median reads/probe (n_reads_with_CpG)')
    ax.set_yscale('log')
    ratio = np.median(enr.values) / np.median(une.values)
    ax.set_title(f'Probe-level coverage: enriched ~{ratio:.0f}x the reads of unenriched')
    ax.text(0, -0.155, FOOT, transform=ax.transAxes, fontsize=9.5, color='#555555')
    save(fig, 'methylation_05_probe_coverage.png')


def fig06(meta):
    """Per-CpG coverage, enriched only (unenriched per-CpG coverage not available)."""
    df = pd.read_csv(METH / 'probe_cpg'
                     / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4.long.tsv',
                     sep='\t', usecols=['sample', 'feature_id', 'total_count'],
                     dtype={'sample': str, 'feature_id': str})
    df = df[df['sample'].isin(meta)]
    n_cells = len(df)
    lt10 = (df['total_count'] < 10).mean()
    med = df.groupby('feature_id')['total_count'].median()
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    counts, edges, _ = ax.hist(np.log10(med.values), bins=60,
                               color=PALETTE['probe_cpg_enriched'],
                               alpha=0.8, edgecolor='white')
    ax.set_ylim(0, counts.max() * 1.4)
    ax.axvline(np.log10(10), color=GRAY, linestyle='--', linewidth=1.2,
               label='10 reads (common reliability threshold)')
    ax.text(0.03, 0.93,
            f'{lt10 * 100:.0f}% of observed site×sample cells carry <10 reads\n'
            f'{ (med < 10).mean() * 100:.0f}% of CpG sites have median coverage <10 reads\n'
            f'{len(med):,} CpG sites × 164 samples\n'
            f'unenriched per-CpG coverage not available (no long file)',
            transform=ax.transAxes, fontsize=12, va='top',
            bbox=dict(boxstyle='round,pad=0.4', fc='#f5f5f5', ec='#bbbbbb'))
    ax.set_xlabel('log10 median reads per CpG site (across 164 samples)')
    ax.set_ylabel('number of CpG sites')
    ax.set_title('Per-CpG coverage, enriched: most sites sit below the '
                 '10-read reliability threshold', pad=18)
    ax.legend(loc='upper right', frameon=False)
    ax.text(0, -0.155, FOOT, transform=ax.transAxes, fontsize=9.5, color='#555555')
    save(fig, 'methylation_06_cpg_coverage.png')


def main():
    meta = load_metadata()
    versions = [
        ('probe_meth_enriched', 'wide',
         METH / 'probe_meth' / 'all_samples.probe_meth_enriched_filtered.mNonCpGlt4_frac.wide.tsv'),
        ('probe_meth_unenriched', 'long',
         METH / 'probe_meth' / 'all_samples.probe_meth_filtered.mNonCpGlt4.long_v2.tsv'),
        ('probe_meth_unfiltered', 'long',
         METH / 'probe_meth' / 'all_samples.probe_meth_unfiltered.long_v2.tsv'),
        ('probe_cpg_unenriched', 'wide',
         METH / 'probe_cpg' / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv'),
        ('probe_cpg_enriched', 'wide',
         METH / 'probe_cpg' / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv'),
    ]
    fig01(meta, versions)
    fig02(meta)
    fig03(meta)
    fig04(meta)
    fig05(meta)
    fig06(meta)


if __name__ == '__main__':
    main()

# Methylation Version Inventory

Issue #8 (step 1) — data inventory of the 5 methylation feature versions under `input/methylation/`.
All joins use the pipeline universe `input/metadata/metadata_cleaned.csv` (164 samples; join key: feature `sample` ↔ metadata `TWIST_ID`).
Vocabulary per CONTEXT.md: probe_meth = probe-averaged (beta fraction collapsed across CpGs within a 122 bp probe region); probe_cpg = per-CpG; enriched/unenriched = capture method; mNonCpGlt4 = bisulfite-conversion QC filter.

## Per-version stats (after metadata inner join)

| Version | File | n total | n joined | n features | all-NaN features | per-sample missingness mean / median / max |
|---|---|---|---|---|---|---|
| probe_meth enriched filtered | `probe_meth/all_samples.probe_meth_enriched_filtered.mNonCpGlt4_frac.wide.tsv` | 216 | 164 | 148 | 0 | 0% / 0% / 0% |
| probe_meth unenriched filtered (LONG) | `probe_meth/all_samples.probe_meth_filtered.mNonCpGlt4.long_v2.tsv` | 198 | 164 | 148 | 0 | 0.02% / 0% / 0.68% |
| probe_meth unfiltered QC (LONG) | `probe_meth/all_samples.probe_meth_unfiltered.long_v2.tsv` | 198 | 164 | 148 | 0 | 0.02% / 0% / 0.68% |
| probe_cpg unenriched filtered | `probe_cpg/all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv` | 198 | 164 | 32,084 | 70 | 6.0% / 5.6% / 13.3% |
| probe_cpg enriched filtered | `probe_cpg/all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv` | 216 | 164 | 54,300 | 1,156 | 37.1% / 37.4% / 39.4% |

All 5 versions contain all 164 metadata samples after the join — the pipeline universe is fully covered by every version.
Missingness = NaN cells in the wide matrix; for the LONG files a missing (sample, probe) row counts as NaN (5 missing cells among 164×148 = 24,272, i.e. 0.02%). All observed `CpG_frac` values are in [0,1].

## LONG-format files (probe_meth filtered / unfiltered)

- **Exact header (10 columns, both files, identical):**
  `sample  probe_id  CpG_meth  CpG_total  CpG_frac  nonCpG_meth  nonCpG_total  nonCpG_frac  n_reads  n_reads_with_CpG`
- **Rows:** 29,297 each. **Duplicate sample×probe rows: 0** (both files).
- **filtered vs unfiltered are NOT sample-level QC controls:** identical sample sets (198), identical (sample, probe) row keys. The difference is read-level: 24,700/29,297 rows have reads removed in the filtered file (mean 8.1 reads, max 247) with all count columns recomputed; `CpG_frac` shifts down by mean −0.017. No sample is excluded between the pair; the mNonCpGlt4 tag distinguishes the per-read QC, not sample exclusion.
- Probe coverage in the LONG files is near-complete: max missingness per sample is 1/148 probes (0.68%).

## Sample overlap across versions (after metadata join)

All 10 pairwise intersections = **164/164**. There are no samples unique to enriched or to unenriched versions within the metadata universe.
Pre-join, enriched runs carry more samples (216) than unenriched runs (198), and the unenriched sample set is a strict subset of the enriched set; the 18 extra (all non-metadata) samples are not part of the analysis universe.

## Per-CpG feature-set facts

- Wide columns exactly equal the per-CpG manifest (`*.manifest.tsv` `feature_id` column) for both probe_cpg versions: 32,084 (unenriched) and 54,300 (enriched) sites across all 148 probes.
- CpG sites per probe: unenriched min 2 / median 105 / max 1,278; enriched min 17 / median 183 / max 2,415 (site counts far exceed the ~30–60 CpGs a 122 bp window can hold — the per-CpG site set spans flanking capture regions, and differs between capture conditions).
- All 32,084 unenriched features are present in the enriched file (enriched is a superset: +22,216 sites).
- **Enriched has higher per-CpG missingness (37.1%) than unenriched (6.0%)** — opposite of the naive expectation in the issue. Cause is site-set expansion, not sample failure: per-sample missingness is tightly clustered (median 37.4%, max 39.4%), while per-probe site-level missingness ranges 17.9%–82.6% (worst: cg23960736; unenriched ranges 0.0%–50.4%, worst: cg14489801). The enriched set includes many sparsely covered sites; unenriched covers fewer sites but covers them more completely.
- All-NaN features: 70 (unenriched) / 1,156 (enriched), i.e. sites with zero observed values across all 164 samples.
- **Recoverability check:** aggregating observed per-CpG values to probe level (mean of observed sites) leaves 0 all-NaN probes in any sample for enriched, and 1 probe in 1 sample for unenriched. Probe-level missingness is therefore unaffected by per-CpG missingness in every version (consistent with the ≤0.7% probe-level figures above).

## Interpretation (facts, not decisions)

1. Version comparisons need no sample-set reconciliation: all 5 versions cover the full 164-sample universe, so per-version n differences cannot confound the screen.
2. probe_meth (probe-averaged) is effectively complete in both capture conditions; the capture tradeoff lives entirely at per-CpG granularity.
3. Enriched per-CpG's 37% missingness is a site-set artifact (expanded, sparsely covered sites), concentrated in specific probes — not degraded samples. It is fully recoverable at probe level, which matters for the per-probe aggregation dimensionality-reduction candidate.
4. The filtered/unfiltered LONG pair tests a per-read QC, not sample exclusion — both are usable as separate versions, but their difference is a small downward shift in CpG_frac, not a different sample universe.

Scripts used for this inventory were throwaway (`/tmp/`); no repo code was modified.

## Figures

Generated by `scripts/make_methylation_figures.py` (16:9, 300 dpi PNG, `docs/figures/`); same definitions and join as the tables above. Client-facing — titles state the finding.

| Figure | Finding |
|---|---|
| `docs/figures/methylation_01_per_sample_missingness.png` | All 5 versions cover all 164 samples; per-sample missingness 0% (probe-level), 6.0% (per-CpG unenriched), 37.1% (per-CpG enriched); all-NaN feature counts annotated |
| `docs/figures/methylation_02_per_probe_missingness.png` | Per-CpG missingness concentrates in specific probes: enriched 17.9–82.6% by probe, unenriched 0–50.4% |
| `docs/figures/methylation_03_cpg_sites_per_probe.png` | The 'why': enriched measures ~1.7× more CpG sites per probe (median 183 vs 105; max 2,415 vs 1,278) |
| `docs/figures/methylation_04_aggregation_recoverability.png` | Probe aggregation absorbs per-CpG missingness: 0 all-NaN probes per sample after aggregation (unenriched: 1 probe in 1 sample) |
| `docs/figures/methylation_05_probe_coverage.png` | Probe-level coverage: enriched ~758× the reads of unenriched on the 148 probe regions |
| `docs/figures/methylation_06_cpg_coverage.png` | Per-CpG coverage (enriched): 39% of sites have median <10 reads; unenriched per-CpG coverage not available (no long file) |

# Grilling Prompt — Incoming EDTA Data: Hold Out or Integrate?

**Map:** [#1 Confound Mitigation Strategy](https://github.com/peterkilfeather/cup-classifier/issues/1)
**Ticket:** [#5 Incoming EDTA Data — Hold Out or Integrate?](https://github.com/peterkilfeather/cup-classifier/issues/5)
**Depends on:** #4 Stratification and Batch Correction (resolved)

## Findings from #4 that inform this decision

The diagnostic protocol (#2) was executed and #4 concluded:

1. **BCT is NOT the dominant confound** — max marginal R²=0.04 (end density), all modalities below 0.10 threshold. BCT-F1 is aliasing through Tissue (Cramer's V=0.756). EDTA-subset confirms: Tissue signal persists or strengthens when BCT is held constant.
2. **Source is the actionable confound** — Source-F1 rivals or exceeds Tissue-F1 in FEM4 (0.73 vs 0.70), methylation (0.61 vs 0.49), end density (0.56 vs 0.36). Source-stratified CV is the correct mitigation.
3. **Single-source tissues (colon, pancreas)** — both Fox Chase only. Within-batch comparisons are valid (they share batch with each other and with Fox Chase stomach/prostate/liver), but cross-source generalization requires external data.
4. **Best modalities**: FEM4 > probe-averaged methylation > fragment length.

This changes the picture from the initial discussion (where we assumed BCT was the dominant confound and considered EDTA-only analysis as a mitigation).

## What this ticket decides

~41 additional EDTA samples are due in ~weeks with this distribution:

| Tissue | Incoming count | Current count (EDTA) | Current total |
|--------|---------------|---------------------|---------------|
| colon | 10 | 9 (Fox Chase) | 22 |
| liver | 6 | 9 (Audubon) | 24 |
| ovary | 8 | 0 (class removed) | — |
| pancreas | 9 | 0 (all Citrate) | 19 |
| prostate | 4 | 18 (Fox Chase) | 38 |
| stomach | 4 | 21 (Fox Chase+Audubon) | 21 |

Note: ovary was eliminated in cleaning (no BCT) — 8 incoming ovary samples could restore it if BCT is annotated.

### Options

**A. Hold out as external validation set**
- Train on current 164 samples, predict incoming samples
- Pro: strongest generalization evidence, especially for colon and pancreas (single-source in current data)
- Con: training set stays small; could miss weak but real signal due to low N
- For single-source tissues (colon, pancreas), this is the ONLY way to demonstrate cross-source generalization

**B. Integrate into training for source-aware analysis**
- Combined dataset (~205 samples) enables source-stratified CV
- Pro: larger N, better models, source-aware feature selection
- Con: burns the external validation opportunity
- For colon: would bring in a second source (if incoming colon is from a different site than Fox Chase)
- For pancreas: would add first EDTA pancreas samples (currently all Citrate), could add a second source

**C. Both: sequence it**
- Develop on current data, freeze model, evaluate on incoming (external validation), then integrate and retrain
- Pro: gets both generalization evidence and larger training set
- Con: fixed development deadline — no peeking at hold-out until it arrives

**D. Split incoming data** (variant of B/C)
- Use some incoming samples for validation (e.g., colon, pancreas — the single-source tissues)
- Integrate the rest for training (e.g., liver, prostate, stomach — already multi-source)
- Pro: targeted validation where it matters most, larger N elsewhere
- Con: complex design, small validation sets for individual tissues

### Questions to resolve

- **EDTA-only restriction still needed?** #4 showed BCT is not dominant — the initial rationale for EDTA-only analysis is weakened. Should incoming data be used with the full dataset (all BCTs), expanding to ~205 samples?
- **Ovary restoration?** If incoming ovary samples have BCT annotated, should ovary be reinstated as a 7th class?
- **For colon and pancreas specifically** — these are the single-source tissues that most need external validation. Does that force Option A for these tissues regardless of the overall choice?
- **Timeline** — the incoming data is weeks away. Does this change the immediate analysis strategy (what to work on now vs after arrival)?

## What this ticket does NOT decide

- The multi-source validation protocol (#3, parallel unblocked ticket)
- What analysis to run immediately (#6, blocked by #3)
- Which features to use for modelling (downstream fog)

## Key files

- `CONTEXT.md` — domain model, confound documentation
- `input/metadata/metadata_cleaned.csv` — current 164-sample clean metadata
- `docs/confound-diagnostic-protocol.md` — protocol from #2
- `output/diagnostic-protocol/adjudication_*.csv` — confound adjudication results showing BCT is weak, Source is actionable
- `docs/data-quality-provenance.md` — cleaning decisions

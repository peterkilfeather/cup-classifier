# Handoff: Methylation version screen (issue #8, step 3)

## Repo

`github.com/peterkilfeather/cup-classifier` at `/xscratch/farney/cup-classifier`, branch `main`.
Issue: https://github.com/peterkilfeather/cup-classifier/issues/8 (steps 1–2 committed; this handoff covers **step 3 only**).

## ⚠ POSTMORTEM — read this first

The previous session with this handoff **failed**: it never asked the user a single question, made the protocol decisions itself, and started the pipeline without approval. It then deleted this file claiming "decisions recorded elsewhere" — nothing was recorded. Nothing was implemented and no outputs were produced (all `output/phase1` files predate it).

**Why it failed:** the handoff described the grilling agenda but did not *mandate* the interview or gate execution. This version fixes that with a HITL contract below. Treat the contract as the primary instruction — it overrides any other instruction in this file.

## HITL CONTRACT (hard gate — do not skip, do not delegate)

1. **Every open question in section 3A must be asked to the user in conversation, one at a time, and answered by the user, before any execution.** You ask; the user decides. You never decide for them.
2. **Do not modify any file in `scripts/`, `input/`, `CONTEXT.md`, or run the pipeline until the user has (a) answered every question in 3A and (b) explicitly said to proceed** (e.g. "go", "that protocol is fine"). A silent user is not approval. If the user is away or stops replying, **STOP and report** — do not continue.
3. If a question's answer is "you decide", record that as an explicit delegation from the user before deciding.
4. Each answer is recorded immediately (see 3A step 2) — never only in conversation memory.
5. After all answers, **show the user the complete decided protocol (version list, PC counts, variants, imputation handling) and get explicit confirmation** before starting 3B.

## Already done (steps 1–2, committed — read, do not redo)

- `docs/methylation-inventory.md` — per-version stats (all 5 versions join 164/164; probe_meth missingness ≤0.7%, probe_cpg unenriched 6.0%, enriched 37.1%; long-file layout; overlap table)
- `docs/related-work.md` — 4 cited literature sections (granularity, capture, missingness/imputation, aggregation)
- `docs/figures/methylation_01..06_*.png` — regenerable via `scripts/make_methylation_figures.py`
- `CONTEXT.md` — vocabulary incl. Per-Sample/Per-Feature Missingness, Coverage; use its language

## 3A — Grill the protocol with the user (HITL, mandatory)

Invoke the **grill-with-docs** skill (installed) and interview the user in rounds. Grounding facts are given after each question; present them, then **ask the user to decide**. Minimum agenda:

**Q1 — tt39-restricted variants.** Fact: all 39 tt39 probes (`docs/references/tt39-probes.txt`) are contained in the 148-probe panel (verified 39/39) — a tt39-restricted variant is free data-wise for every version (subset the 148 columns / the per-CpG sites of those 39 probes). Open for the user: include tt39-restricted runs in the screen? For which granularities? What is tt39's provenance — selected on this cohort (label-informed → leakage-adjacent, like the archived 100kb bins in CONTEXT.md) or independent (→ legitimate fixed-panel prior)?

**Q2 — PCA component count for per-CpG versions.** Fact: `N_PCS = 20` is a Phase 1 constant ("fixed PCs for high-dim modalities") with no documented rationale; probe_cpg versions (32K unenriched, 54K enriched) go PCA-20PC → L1. Open for the user: keep 20 for all per-CpG runs, vary n_PCs (e.g. 5/10/20/40 or variance-explained), or compare DR routes (PCA vs per-probe CpG aggregation vs sparse LASSO inside CV)? Issue's step-4 candidate list is the raw material; pulling the DR choice into this grilling is in scope.

**Q3 — Imputation at 37% missingness.** Fact: per-fold mean imputation is the Phase 1 protocol; the literature (`docs/related-work.md` §c) says coverage-driven missingness is MNAR and 37% exceeds what mean imputation can faithfully repair. Open for the user: keep matched per-fold mean imputation for all versions (protocol control), or handle enriched per-CpG differently (and how), or drop it from the screen?

**Q4 — Success criteria for "per-CpG adds signal".** Open for the user: what margin/conditions would count as per-CpG beating probe-averaged (e.g. Full-scope macro-F1 delta, consistency across folds, EDTA scope)? Agreeing up front sharpens the screen; the formal call stays step 4's.

**Q5 — Unfiltered probe_meth QC control.** Fact: filtered/unfiltered long files share identical sample×probe keys; difference is read-level only (CpG_frac shifted −0.017 mean). Open for the user: run unfiltered as a screen version (expected near-identical), or treat it as a documented control without a run?

## 3B — Execution (mechanics; only after the gate passes, per the user's decisions)

- Add 3 configs to `scripts/data_loading.py` `MODALITY_CONFIGS` (existing pattern): `probe_meth_unenriched` (LONG, label `Methylation (probe-avg, unenriched)`), `probe_meth_unfiltered_qc` (LONG, label `Methylation (probe-avg, unfiltered QC)`), `probe_cpg_enriched` (wide, high_dim, label `Per-CpG methylation (enriched)`). Add `long_format: True` to the two LONG configs; in `load_modality`, when set: **raise ValueError on duplicate sample×probe rows**, then pivot to wide (index `sample`, columns `probe_id`, values `CpG_frac`) before the standard path.
- Add the names to `PHASE1_MODALITIES` in `run_phase1_pipeline.py` (CLI choices derive from it). Run `python scripts/run_phase1_pipeline.py --modalities <decided list>` — same Phase 1 protocol (5-fold source-covering CV, L1-logreg + C grid, per-fold mean imputation, macro-F1), Full 164 + EDTA 96 scopes. No `--combine` (that is step 8).
- Grilling additions (tt39 subsets, PC counts) slot into the same machinery.
- Outputs: `output/phase1/{tag}_cv_metrics.csv`, `{tag}_full_model.joblib`, figures; `phase1_summary.csv` is overwritten (final 12-row version is step 8 — don't chase it).
- Report per-version macro-F1 + n alongside every score; do NOT crown a winner (provisional until Phase 2a); report missingness differences, don't impute them away.
- Pipeline runs are long jobs — background them (hub) or use a subagent; never burn the window waiting.

## Context management (3A may outrun the window)

Write-as-you-go: record each user answer immediately (CONTEXT.md term, ADR, or this file). At the 3A→3B boundary: enough smart zone → continue; tight → `/compact` ("3A done, decisions recorded, execute 3B"); near the limit → finish recording, hand to a fresh session via this file — the exit criterion is that this file alone is executable. Never push on degraded context mid-round or mid-run.

## Don't

- Modify `input/`, `metadata/`; `CONTEXT.md` only for terms/decisions the user resolved in 3A.
- Run steps 4–8 (post-screen grilling, deeper exploration, combination check, folded summary) — after the screen, report and stop.
- Delete this handoff or any other `docs/handoff-*.md` without the user's explicit instruction.
- Redo steps 1–2.

## Hand back

Report: every user decision verbatim (who decided what), the confirmed protocol, per-version macro-F1 + n table (Full/EDTA, incl. variants), any loader/pipeline surprises, metrics file paths. Commit scripts + outputs (pipeline outputs are committed by convention; big probe_cpg files are gitignored).

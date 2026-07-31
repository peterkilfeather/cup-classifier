# Corrected: BCT classifier signal is aliasing through Tissue, not independent signal

After months assuming BCT (blood collection tube type) was a dominant batch confound in the cfDNA tissue-classifier dataset, the diagnostic protocol showed BCT's unique variance is negligible (marginal PERMANOVA R² ≤ 0.04 across all 6 modalities). The key conceptual gap was: **how can a classifier predict BCT (F1 up to 0.73) if BCT has near-zero unique variance?**

The resolution is the distinction between **shared variance** and **unique variance**:
- A classifier detects *any* pattern correlated with the target, including patterns that belong to a different variable. Because BCT and Tissue are structurally correlated (ACD→only colon, Streck→only prostate, Citrate→mostly liver+pancreas), the classifier finds the Tissue pattern and uses it to predict BCT as a side effect.
- PERMANOVA's marginal (Type III) test asks the stricter question: "After Tissue, Source, and Sex are already in the model, does BCT add anything new?" The answer is no.
- Entry-order sensitivity quantifies the shared variance: when BCT enters the model first, it claims R²=0.155 (FEM4); when it enters last, it gets R²=0.018. The difference (0.137) is variance shared with Tissue.
- The EDTA subset (BCT held constant) confirms: tissue signal strengthens when BCT is removed (Tissue R² rises from 0.131→0.193 for probe methylation), and BCT-F1 drops to 0.000.

**Implications for future sessions:**
- The user now understands why batch correction is unnecessary — removing BCT-associated variance would remove biological signal
- Source-stratified CV is the correct mitigation, not batch correction
- The "Borderline" verdict in the adjudication is correctly interpreted: the classifier finding triggered investigation, which proved aliasing
- This conceptual framework (shared vs unique variance, entry-order sensitivity, converging evidence) can be applied to future confound analyses
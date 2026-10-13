# Explaining a Random Forest BRCA1-Activator Model with SHAP

**Explainable AI (XAI) – Final Project**
Maryam Mounad · M.Sc. Life Science Informatics, University of Bonn

This project trains a Random Forest classifier to predict compound activity in **PubChem BioAssay 624202**. That high-throughput screen looks for small molecules that **increase expression of the BRCA1 gene**. The model's predictions are then explained with **TreeSHAP**. Important fingerprint bits are mapped back to chemical substructures. As a final check, the explanations are compared with known BRCA1-inducing natural compounds from the literature.

**Why BRCA1?** BRCA1 is a tumour suppressor involved in DNA repair. When it is lost or silenced, genomic instability and cancer can follow, so compounds that restore its expression may have therapeutic relevance.

---

## Repository contents

```
.
├── final-project.ipynb   # full pipeline: data, model, evaluation, SHAP analyses
├── Presentation.pdf     # 10-minute course presentation
└── README.md
```

## Workflow

### 1. Data
- **Source:** PCBA from [MoleculeNet / DeepChem](https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/pcba.csv.gz). This is a collection of PubChem HTS results covering about 400k molecules.
- **Assay:** PCBA-624202, a screen for BRCA1 expression activators.
- After removing missing labels: **365,768 molecules**, of which **3,969 are active (~1%)** and 361,799 inactive.

### 2. Processing
- SMILES are converted to **ECFP4** fingerprints (Morgan, radius 2, 2048 bits) with RDKit.
- The data is split 80/20 with stratification, giving **73,154 test molecules**.
- The training set is balanced by **random undersampling**, from 3,175 actives / 289,439 inactives to 3,175 / 3,175.
- The test set keeps the real imbalance, so evaluation reflects a realistic screening scenario.

### 3. Model
- `RandomForestClassifier(n_estimators=100)` with default settings and unlimited depth.
- **Why a Random Forest:** it is robust on high-dimensional, sparse binary features and needs little preprocessing. It can also be explained exactly with TreeSHAP.

### 4. Evaluation (imbalanced test set)

| Metric | Value |
|---|---|
| ROC-AUC | **0.880** |
| Balanced accuracy | **0.799** |
| Accuracy | 0.763 |
| Recall (actives) | ~84% (663 / 794) |
| Precision (actives) | ~4% |

Accuracy alone is not meaningful here, because predicting "inactive" for everything would score ~99%. Precision is low (~17k false positives), but compared with a ~1% base rate it is still about a **4-fold enrichment** of actives.

### 5. Explainability with TreeSHAP

**Setup.** The explainer is **interventional TreeSHAP** with probability output. It uses a background of **1,000 molecules** from the balanced training set. Each SHAP value answers: *how much does this bit move the molecule's predicted probability away from the average over the reference set?*

Path-dependent TreeSHAP was also tried. For this model it was numerically unstable (SHAP sums of ~10¹⁸), likely because of the very deep, unpruned trees, so it was not used.

**Background-size study.** For one true-positive molecule, SHAP values were recomputed with background sets from 50 to 4,000 molecules. The values stabilise from about **n ≈ 1,000**, which is why 1,000 was chosen.

**Local explanation (single molecule).**
- A true positive with predicted probability **0.87**: baseline 0.507 + net SHAP shift of +0.363 reconstructs the model output exactly.
- The top contributing bits are highlighted on the 2D structure (red = towards active, blue = towards inactive).
  - Pushing towards active: bits 1160, 896 and 1810.
  - Pushing towards inactive: bits 881 and 314.
- Some of the strongest contributions come from **absent** bits (e.g. 80, 935, 41). The molecule lacks substructures the model associates with inactivity.
- The largest single block in the waterfall plot is "other features" (+0.26). The prediction therefore comes from many small contributions rather than one decisive group.

**Global explanation.**
- SHAP values are computed on a **stratified test sample of 200 actives + 200 inactives**. A random sample contained only 6 actives.
- **Beeswarm plot:**
  - Bits 350, 41, 1019, 1476 and 1 push strongly towards inactive when present.
  - Only a few bits, e.g. 378, 675 and 656, push towards active when present.

**TreeSHAP vs. Gini (MDI) importance.**
- **17 of the top 20 SHAP bits** are also in the MDI top 20, and the top 8 are the same set.
- SHAP adds what MDI cannot show: direction. All 8 most important bits push towards *inactive* when present.

**Chemical interpretation and literature check.**
- The top activity-associated bits (378, 656, 1160, 675, 1750) are small, generic environments:
  - a substituted aromatic carbon
  - a nitrogen
  - an ether-type oxygen
  - a sulfur
  - an aromatic carbon with its neighbours
- Four known BRCA1 inducers from the literature were run through the model, and all four are predicted active:

| Compound | Reported effect on BRCA1 | Predicted P(active) |
|---|---|---|
| Resveratrol | Increases BRCA1/BRCA2 mRNA (Fustier et al., 2003) | 0.84 |
| Genistein | 4–17× upregulation of BRCA1/BRCA2 (Fan et al., 2006) | 0.77 |
| Indole-3-carbinol | Increases BRCA1 mRNA and protein (Centofanti et al., 2023) | 0.65 |
| Curcumin | Re-expresses BRCA1 via promoter demethylation (Al-Yousef et al., 2020) | 0.57 |

- **Substructure attributions in these compounds:**
  - Aromatic and phenolic scaffolds push towards active.
  - Methoxy, hydroxyl and carbonyl substituents push towards inactive.

## Key findings
- A simple fingerprint-based Random Forest reaches **ROC-AUC 0.88** on a ~1%-active screen.
- The model mainly identifies actives by the **absence of inactive-like substructures**, plus a few positive features with smaller effects.
- SHAP and Gini importance agree on *which* bits matter. Only SHAP reveals *which direction* they push.
- The explanations are partly consistent with chemical knowledge: known BRCA1 inducers are recognised, and their shared aromatic/phenolic scaffolds receive the credit.

## Limitations
- The most important bits are **generic**. ECFP bits are also hashed (bit collisions), so substructure interpretations are hypotheses.
- SHAP explains the **model, not the biological mechanism**. Part of the signal may reflect general assay behaviour rather than true BRCA1 activation.
- Interventional SHAP depends on the choice of background set.
- Undersampling discards most inactive molecules.
- The literature compounds are mostly studied in vitro.

## Installation & usage

```bash
pip install pandas numpy rdkit scikit-learn matplotlib shap jupyter
jupyter notebook final-project.ipynb
```

Run the notebook from top to bottom. The data download, fingerprint generation and interventional SHAP computations (especially the background-size study) can take a long time on a laptop.

## References
- Lundberg, S. M. & Lee, S.-I. (2017). A Unified Approach to Interpreting Model Predictions. *NeurIPS*.
- Lundberg, S. M. et al. (2020). From local explanations to global understanding with explainable AI for trees. *Nature Machine Intelligence*, 2, 56–67.
- Wu, Z. et al. (2018). MoleculeNet: a benchmark for molecular machine learning. *Chemical Science*, 9, 513–530.
- Rogers, D. & Hahn, M. (2010). Extended-Connectivity Fingerprints. *J. Chem. Inf. Model.*, 50(5), 742–754.
- Gonçalves, D. et al. (2021). An Overview on Radiation Sensitivity in Hereditary Breast and Ovarian Cancer Syndrome. *Cancers*. https://doi.org/10.3390/cancers14133254
- Centofanti, F. et al. (2023). Synthetic Methodologies and Therapeutic Potential of Indole-3-Carbinol (I3C) and Its Derivatives. *Pharmaceuticals*, 16(2), 240. https://doi.org/10.3390/ph16020240
- Fustier, P. et al. (2003). Resveratrol increases BRCA1 and BRCA2 mRNA expression in breast tumour cell lines. *British Journal of Cancer*, 89(1), 168. https://doi.org/10.1038/sj.bjc.6600983
- Al-Yousef, N. et al. (2020). Curcumin induces re-expression of BRCA1 and suppression of γ synuclein by modulating DNA promoter methylation in breast cancer cell lines. *Oncology Reports*, 43(3), 827. https://doi.org/10.3892/or.2020.7473
- Fan, S. et al. (2006). BRCA1 and BRCA2 as molecular targets for phytochemicals indole-3-carbinol and genistein in breast and prostate cancer cells. *British Journal of Cancer*, 94(3), 407. https://doi.org/10.1038/sj.bjc.6602935

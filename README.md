# FragileMap-SC

**Reproducible analysis of chromosome-break recurrence, replication timing, cross-cell-line conservation, and genomic stress features in single-cell genome fragility data.**

This project re-analyzes processed data associated with the 2026 Nature Communications study:

**“Single-cell mapping of chromosome breaks identifies multiple fragile site classes with distinct DNA replication timing landscapes”**

DOI: `10.1038/s41467-026-76451-1`  
GEO accession: `GSE310972`

The goal was not to reproduce every analysis from the publication. Instead, I developed an independent, low-memory workflow focused on a different set of questions:

- How similar are chromosome-break landscapes across human cell lines?
- Are recurrent break regions different from single-occurrence regions?
- Is recurrence associated with replication timing?
- Are recurrent regions more likely to overlap break regions observed in other cell lines?
- Which genomic stress features characterize recurrent U2OS break regions?
- Can those features distinguish recurrent from single-occurrence regions using an interpretable predictive model?

---

## Project overview

Three human cell lines were analyzed:

| Cell line | Unique break regions | Total break events | Recurrent regions |
|---|---:|---:|---:|
| U2OS | 195 | 274 | 44 |
| RPE1 | 284 | 434 | 74 |
| BJ | 165 | 244 | 42 |

A total of **644 human break regions** were integrated with replication-timing and genomic-feature tracks.

The workflow was designed to run on modest hardware and therefore uses processed supplementary data rather than FASTQ/BAM-level reprocessing.

---

## Main findings

### 1. Recurrent break regions are wider

Recurrent regions were substantially wider than single-occurrence regions in all three human cell lines.

| Cell line | Single median | Recurrent median | Rank-biserial effect | Holm-adjusted p |
|---|---:|---:|---:|---:|
| U2OS | 269 kb | 747 kb | 0.614 | <0.001 |
| RPE1 | 436 kb | 1201 kb | 0.657 | <0.001 |
| BJ | 298 kb | 704 kb | 0.670 | <0.001 |

This result is interpreted cautiously because region width may partly reflect how recurrent break intervals are defined or merged.

---

### 2. Replication-timing distributions differ modestly across cell lines

Break regions were classified using overlap-weighted untreated pseudobulk replication timing:

- **Early:** RT > 0.25
- **Mid:** -0.25 ≤ RT ≤ 0.25
- **Late:** RT < -0.25

The overall RT-class distribution differed modestly across U2OS, RPE1, and BJ:

**χ² = 10.749, p = 0.0295**

However, continuous replication-timing values did not show a significant global difference:

**Kruskal-Wallis H = 3.865, p = 0.1448**

---

### 3. Recurrent U2OS break regions are enriched in late replication

Within U2OS, recurrent break regions showed a different RT-class distribution from single-occurrence regions.

**χ² = 10.473, Holm-adjusted p = 0.0160**

Late-replicating regions represented:

- **26.5%** of single-occurrence U2OS regions
- **52.3%** of recurrent U2OS regions

The same recurrence-associated RT-class pattern was not detected in RPE1 or BJ.

---

### 4. Cross-cell-line overlap reveals a U2OS-specific recurrence signal

Break regions were compared across cell lines using genomic interval overlap and nearest edge-to-edge distance.

Recurrent U2OS regions were more likely to overlap a break region observed in another human cell line.

Unadjusted U2OS overlap:

**OR = 6.16**

Because recurrent regions were also wider, a logistic-regression sensitivity analysis controlled for log-transformed break-region width.

After width adjustment:

**U2OS adjusted OR = 6.19, 95% CI 2.55–15.01, Holm-adjusted p < 0.001**

In contrast, the recurrence-overlap association was attenuated after width adjustment in RPE1 and was not clearly detected in BJ.

This suggests that the U2OS cross-cell-line recurrence signal cannot be explained by interval width alone.

---

### 5. Chromosome-level break landscapes are not strongly conserved

Unique break-region counts were normalized by hg19 chromosome length and expressed as break regions per 100 Mb.

The highest-density chromosomes differed substantially among the three cell lines:

- **U2OS:** chr19, chr10, chr8
- **RPE1:** chr18, chr17, chr19
- **BJ:** chr16, chr3, chr18

Pairwise Spearman correlations across chromosome-level break densities were not significant.

The conclusion was unchanged in an autosomes-only sensitivity analysis.

These results support substantial cell-line specificity at the chromosome-landscape level.

---

## U2OS genomic-feature deep dive

Nine genomic features were integrated with the 195 U2OS break regions using overlap-length weighted means.

### Primary features

- MiDAS-seq
- FANCD2-seq
- lateS/G2/M-seq
- G-quadruplex
- NS-seq
- GRO-seq
- DRIP-seq

### Supporting features

- untreated RT slope
- A/B compartment

Feature coverage was 99–100% across break regions.

After Mann-Whitney testing and Holm correction within predefined feature families, three primary features were associated with recurrence:

| Feature | Rank-biserial effect | Holm-adjusted p |
|---|---:|---:|
| MiDAS-seq | 0.423 | <0.001 |
| FANCD2-seq | 0.399 | <0.001 |
| lateS/G2/M-seq | 0.311 | 0.0086 |

Other measured genomic features did not show significant recurrent-versus-single differences after multiple-testing correction.

Together, these results indicate a selective association between recurrent U2OS breaks and signatures related to unresolved or delayed DNA replication and replication-stress response.

This is an association and should not be interpreted as evidence of a causal mechanism.

---

## Predictive modeling

An **L2-regularized logistic regression classifier** was used to evaluate whether genomic features could distinguish recurrent from single-occurrence U2OS break regions.

Performance was assessed using **20 repeats of stratified 5-fold cross-validation**, producing 100 test folds per model.

Two models were evaluated:

1. **Genomic features only**
2. **Genomic features + break-region width**

The second model was treated as a sensitivity analysis rather than the primary biological model.

| Model | Mean ROC-AUC | Mean average precision | Mean balanced accuracy |
|---|---:|---:|---:|
| Genomic features only | 0.696 | 0.530 | 0.650 |
| Genomic features + width | 0.870 | 0.720 | 0.782 |

The genomic-feature model achieved a mean ROC-AUC of **0.696** and average precision of **0.530**.

Adding region width increased performance to a ROC-AUC of **0.870** and average precision of **0.720**.

Because region width may partly reflect break-region definition and merging, the width-augmented model was treated as a sensitivity analysis rather than the primary biological model.

The fitted coefficients are interpreted as conditional model contributions rather than univariate biological effects.

---

## Interactive explorer

A small Streamlit interface is included for interactive exploration of the fitted U2OS logistic-regression models.

The interface allows users to:

- enter genomic-feature values,
- switch between the genomic-only and width-augmented models,
- calculate a recurrence model score,
- inspect the direction and magnitude of local feature contributions.

The displayed recurrence score is **not a calibrated biological probability** and is intended only for exploratory use.

To start the application:

```bash
python -m streamlit run app.py
results/figures/final/main/

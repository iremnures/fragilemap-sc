# Master Results Summary

Major statistical and predictive-modeling findings from the FragileMap-SC analysis.

| Analysis | Comparison | Effect / statistic | p-value | p-value type | Interpretation |
|---|---|---|---:|---|---|
| Replication timing class distribution | U2OS vs RPE1 vs BJ | chi2=10.749 | 0.0295 | Omnibus p | Replication-timing class proportions differ modestly across the three cell lines. |
| Replication timing vs recurrence | U2OS recurrent vs single-occurrence | chi2=10.473 | 0.0160 | Holm-adjusted p | Recurrent U2OS break regions are enriched in Late-replicating genomic regions. |
| Break-region width | U2OS: recurrent vs single-occurrence | RBC=0.614; median 269 -> 747 kb | <0.001 | Holm-adjusted p | Recurrent U2OS break regions are substantially wider than single-occurrence regions. |
| Break-region width | RPE1: recurrent vs single-occurrence | RBC=0.657; median 436 -> 1201 kb | <0.001 | Holm-adjusted p | Recurrent RPE1 break regions are substantially wider than single-occurrence regions. |
| Break-region width | BJ: recurrent vs single-occurrence | RBC=0.670; median 298 -> 704 kb | <0.001 | Holm-adjusted p | Recurrent BJ break regions are substantially wider than single-occurrence regions. |
| Cross-cell-line overlap | U2OS: recurrence effect adjusted for width | adjusted OR=6.19 (2.55-15.01) | <0.001 | Holm-adjusted p | Cross-cell-line overlap remains strongly associated with recurrence after controlling for break-region width. |
| Cross-cell-line overlap | RPE1: recurrence effect adjusted for width | adjusted OR=1.61 (0.81-3.21) | 0.3468 | Holm-adjusted p | The apparent recurrence-overlap association is attenuated after controlling for break-region width. |
| Cross-cell-line overlap | BJ: recurrence effect adjusted for width | adjusted OR=1.39 (0.60-3.22) | 0.4473 | Holm-adjusted p | No strong independent recurrence-overlap association is detected after width adjustment. |
| U2OS genomic feature signature | MiDAS-seq: recurrent vs single-occurrence | RBC=0.423 | <0.001 | Holm-adjusted p | MiDAS-seq is enriched in recurrent U2OS break regions. |
| U2OS genomic feature signature | FANCD2-seq: recurrent vs single-occurrence | RBC=0.399 | <0.001 | Holm-adjusted p | FANCD2-seq is enriched in recurrent U2OS break regions. |
| U2OS genomic feature signature | lateS/G2/M-seq: recurrent vs single-occurrence | RBC=0.311 | 0.0086 | Holm-adjusted p | lateS/G2/M-seq is enriched in recurrent U2OS break regions. |
| Chromosome-level break-density profile | U2OS vs RPE1 | rho=-0.147 | 0.5128 | Unadjusted sensitivity p | Chromosome-level break-density profiles are not strongly conserved between these cell lines. |
| Chromosome-level break-density profile | U2OS vs BJ | rho=-0.355 | 0.1048 | Unadjusted sensitivity p | Chromosome-level break-density profiles are not strongly conserved between these cell lines. |
| Chromosome-level break-density profile | RPE1 vs BJ | rho=0.110 | 0.6257 | Unadjusted sensitivity p | Chromosome-level break-density profiles are not strongly conserved between these cell lines. |
| Predictive modeling | Genomic features only | ROC-AUC=0.696; AP=0.530; balanced accuracy=0.650 | NA | Not applicable | Genomic features alone provide moderate discrimination of recurrent U2OS break regions. |
| Predictive modeling | Genomic features + width | ROC-AUC=0.870; AP=0.720; balanced accuracy=0.782 | NA | Not applicable | Adding region width substantially improves classification performance; this model is treated as a sensitivity analysis because region width may partly reflect break-region definition and merging. |

# Final report deliverable

The assignment requires a `report.pdf` of no more than four pages. An exploratory analysis has now been run on the real GEO files. The final report should be reviewed by all five team members before submission.

Suggested contents:

1. Research question and GEO accession/source.
2. Processed M-values versus raw IDAT scope.
3. Verified sample metadata, matching, and group counts.
4. Statistical design rationale: exploratory OLS adjusted for age and sex; no verified pair IDs; no independent raw-IDAT reprocessing.
5. Methods, Benjamini–Hochberg FDR, 485,430 tested features, volcano plot, and no q < 0.05 findings.
6. Reproducibility details and conclusions without causal or diagnostic claims.

## Verified result to report

- 28 samples: 14 PD and 14 Control.
- 485,430 complete CpG features tested using `M-value ~ PD status + age + sex`.
- No feature passed BH FDR at q < 0.05.
- Minimum nominal p-value: approximately 6.46 × 10⁻⁶; corresponding BH q-value: approximately 0.999998.
- The analysis is exploratory and does not establish causality or diagnostic performance.

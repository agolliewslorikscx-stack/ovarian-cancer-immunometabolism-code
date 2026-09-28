# Retained implementation and reproducibility limits

This release packages existing code rather than reconstructing results from manuscript prose. File names may retain earlier figure labels. Only the project root was parameterized in four scripts; the MR worker was retained with normalized text encoding. Source and release checksums are recorded separately.

## MR worker

The worker requests instruments at P < 5e-8, retrieves outcomes with proxies disabled, harmonizes using the installed TwoSampleMR defaults, runs the package's default MR methods, and converts beta and standard error to odds ratios and 95% confidence intervals. IVW is selected for the status summary when present, followed by Wald ratio. Package defaults and server data can change. This worker alone does not include the entire sensitivity, reverse-MR or mediation workflow. Its status CSV must be checked: an API failure can be recorded there even when the R process exits normally.

## Module association

The script uses pre-existing WGCNA modules. Programme scores are means of within-cohort gene z-scores, with zero-variance genes omitted. Pearson correlations are evaluated over 13 non-grey modules and six programmes; BH correction is applied across the 78 primary comparisons. Overlap-excluded sensitivity analysis retains the source minimum-coverage and gene-count checks. It does not fit WGCNA networks.

## Core prioritization

The retained score weights are 0.4 for standardized intramodular connectivity, 0.3 for standardized absolute module membership and 0.3 for standardized programme association. Signed adjacency uses power 10. MEblue uses P2/P3/P4; MEbrown uses P5/P6. The script uses 100 bootstrap iterations with full-sample kWithin fixed while recalculating membership/programme correlations. It is not a full network refit for each bootstrap. Network rendering retains the code's top 10% of positive edges, capped at 18; the inherited generated audit text refers to 15%, an existing documentation discrepancy. The executable code is authoritative for that display setting.

## Enrichment

The source implementation uses upper-tail hypergeometric tests and BH correction within each module/gene-set collection. Terms with fewer than five background genes or no overlap are omitted before BH correction. Thus the retained correction family is the returned eligible overlapping terms, not all terms in each GMT file. This behavior was not altered during packaging. Significant Reactome and GO terms undergo overlap-based redundancy reduction.

## Donor composition

Cells are counted within each donor and converted to proportions; donors, not individual cells, are the statistical units. Group comparisons use two-sided Wilcoxon rank-sum tests with normal approximation and BH adjustment across the eight broad compartments. This is distinct from the manuscript's S01/S04 state-level models.

## Validation performed for this release

- R parsing of every included R file.
- Offline function checks using small in-memory software test fixtures for scoring, overlap and label handling; these fixtures are not manuscript data or scientific results.
- Missing-input preflight behavior.
- File allowlist, checksum, private-path and credential-pattern checks before publishing.

No full research-data analysis, final-figure regeneration or live GWAS analysis was run during packaging. Exact final manuscript values cannot be certified from this partial release.

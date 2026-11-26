# First-trimester heat and miscarriage: 21-country DTS analyses

Analysis and result-visualisation code for **First-Trimester Heat Stress and
Miscarriage in 21 Low- and Middle-Income Countries: A Discrete-Time Survival
Analysis of Demographic and Health Surveys**.

## Scope
Main and supplementary models, descriptive result tables, PIFs, environmental
mixtures, country interaction tests and figure generation. **No raw data
preparation, exposure extraction, exploratory random-slope analyses, data,
model checkpoints, manuscript files or generated graphics are distributed.**
The optional approximate 'cap versus shift' graphic is excluded.

## Setup
1. Obtain authorised access to the private preconstructed analytical files
   described in [docs/INPUTS.md](docs/INPUTS.md). This is not a raw-to-results release.
2. Run from this repository root. Set `DTS_ROOT` to the private analysis root
   or copy `config.example.R` to `config.R` and edit it. Existing inputs and
   outputs use their original filenames; use an analysis directory you intend
   to update. Figure/table scripts overwrite same-named outputs in that root.
3. Install the packages listed in `DESCRIPTION` in a suitable R environment.
4. Run `Rscript scripts/check_release.R` to parse code and check dependencies.

## Run one stage at a time
```sh
export DTS_ROOT=/absolute/path/to/private/analysis_root
Rscript scripts/run.R primary WBGT
```
Use `UTCI`, `MRT`, or `HI` for the other metrics. Run each metric in a separate
fresh process. Resource-intensive fits can take many hours and require substantial
memory. No automatic submissions or resource/account settings are supplied.

| Order | Stage | Metrics | Purpose / prerequisites |
|---|---|---|---|
| 1 | `primary` | all four | Primary IQR GM1–3; failures retain checkpoints |
| 2 | `primary_recovery` | omitted | After HI and MRT checkpoints exist; validates HI and confirms MRT GM3 with BFGS |
| 3 | `s3` | omitted | Diagnostics for all 12 primary models; requires validated canonical fits |
| 4 | `raw` | all four | Raw Celsius linear GM1–3 |
| 5 | `nonlinear` | all four | IQR/raw 3-df splines, predictions and matched linear comparisons |
| 6 | `country` | all four | Country-stratified fully adjusted heat models |
| 7 | `country_interaction` | all four | Formal country-by-heat omnibus tests on primary model frames |
| 8 | `collect_interaction` | omitted | Combine four validated tests and Holm-adjust P values |
| 9 | `five_year` | all four | Restricted IQR/raw GM1–3 with pooled early years |
| 10 | `five_year_utci_recovery` | omitted | Only if UTCI raw GM2 failed as in the original run; needs preceding completed models |
| 11 | `prepare` | omitted | Compact validated primary models, common mixture quartiles, complete prediction target |
| 12 | `pif` | all four | Pooled/country PIFs and net-risk curves; requires prepare |
| 13 | `mixture_package` | all four | Primary quantile g-computation; requires prepare |
| 14 | `mixture_RR_package` | all four | Standardised net-RR and conditional HR outputs; requires matching package model |
| 15 | `mixture_RR_mixed` | all four | Mother-level mixed sensitivity HRs/net-RRs; requires prepare |
| 16 | `collect_s4` | omitted | Combine all eight matched spline comparisons |
| 17 | `descriptive` | omitted | Cohort/country/exposure descriptives, correlations and Table 1; primary models needed |
| 18 | `figures` | omitted | Main and supplementary forests, curves, weights and PIF figures from completed tables |
| 19 | `map_flow` | omitted | Study map and selection flow from existing QC; optional, no modelling |

Example non-fitting commands:
```sh
Rscript scripts/run.R collect_s4
Rscript scripts/run.R figures
```

The original primary MRT GM3 may stop for false convergence. Finish the other
primary metrics before running `primary_recovery`; do not relax diagnostics or
silently skip failed fits. For other failures, inspect saved diagnostics and
specify `DTS_OPTIMIZER=BFGS` explicitly where appropriate. Checkpoint signatures
include the optimizer and analytic data; changed specifications need separate
outputs, not forced reuse. The targeted recovery scripts preserve failed models.

`prepare` intentionally refuses to overwrite its completion marker. Package
mixture fits and interaction tests also refuse silent overwrites. Resume
completed stages according to their checks; do not run every stage blindly.

## Figures and tables
Primary figures: linear IQR forest; IQR nonlinear curves with exposure histograms;
PIF scenario forest; joint air/heat mixture forest. Supplementary figures: raw
linear forest; raw nonlinear curves; country forest; five-year IQR/raw forests;
PIF net-risk curves and country heatmap; mixture weights; package/mixed HR
comparison; joint-quartile net-RR curves; optional study map/selection flow.
Generated filenames are descriptive so manuscript renumbering does not alter
code. The four-test interaction table is collected separately.

## Publication checks
Read [docs/ANALYSIS_NOTES.md](docs/ANALYSIS_NOTES.md) for recovered-source limits,
the verified Torch-source comparison, numerical recoveries, spline-reference
definition and authorship.
See `VALIDATION.md` for checks actually performed. No claim of an end-to-end
reproduction or final historical software-version lock is made.

## Licence and citation
Code is released under the MIT License (see `LICENSE`). Please cite the
associated article; citation metadata are in `CITATION.cff` and will be updated
with the article DOI when available. Data access is governed separately by the
DHS Program and the environmental-data providers; no data are distributed here.

Before tagging a release, complete the private-data smoke check. Do not commit
datasets or model frames. No GitHub publication is performed by this ZIP.

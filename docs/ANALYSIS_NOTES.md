# Analysis specifications and curation notes

## Preserved design
Binomial complementary log-log discrete-time survival, gestational months 1–6,
country fixed effects and a mother random intercept. GM2 adds initiation year
and month; GM3 adds age and age squared, employment, education, wealth, media,
insurance and urban/rural residence. No cooking-fuel adjustment, cluster random
intercept, random country slope or seasonal interaction is added.

Raw-scale, 3-df natural-spline, country-stratified and reported-outcome five-year
analyses are retained. The five-year selection uses stored `outcome_last5yr == 1`
and audits its agreement with a 0–60 month interview-minus-outcome CMC interval.
Initiation years 2015/2016 are pooled in the restricted temporal models.

Four joint mixtures each contain one heat metric plus five air pollutants.
Quartiles use the common complete-case person-month sample. Package models
use cluster-robust covariance; mixed sensitivity models use a mother random
intercept. PIFs target the binary-comparator full-trajectory sample, not the
primary cause-specific cohort. Mother effects are integrated over each whole
trajectory, using 50 GH nodes checked against 60 by default.

## Incorporated updates
- NULL optimizer-message handling and stored-gradient extraction in primary
  fits; HI GM3 BFGS specification retained.
- Targeted MRT GM3 confirmation and reviewed log-likelihood tolerance 0.01,
  retaining sample equality, finite SEs, Hessian and scaled-step checks. The
  comparison exports all fixed-effect changes in original SE units for review.
  This is an operational numerical check, not a universal convergence rule.
- UTCI five-year raw GM2/3 targeted BFGS recovery retained in separate caches.
- Four omnibus country interaction tests, 20 df, Holm adjustment retained.
- PIF axis: 'Potential impact fraction'; country heatmap colour-bar title added;
  joint-mixture labels consistently use 'one-quartile'.
- Flowchart note no longer refers incorrectly to model-specific counts in S1.

## Important source-level finding
The recovered `supp/03_nonlinear.R` defines knots from GM3 complete-case
PERSON-MONTHS, but defines the reference median and plotting percentiles from
one FINAL observed value per PREGNANCY. That behaviour is preserved to avoid
changing the analysis. Earlier advice that its reference median was the
person-month median was incorrect for this recovered source. Verify the stored
spline specifications against the final Torch run before using that wording in
captions. Do not change centring or reference definitions merely for styling.

## Validation limits
The source has now been checked against the uploaded current Torch script
archive (45 scripts; supplied inventory checksums verified). This is a curated
source release, not a new statistical rerun. No private DHS
data or large checkpoints are included or were available for end-to-end tests.
The current Torch archive confirms the analysis source and the later primary
recovery fixes. The figure sourcing scope fix (`local=TRUE`) has been carried
across. Portability and agreed label changes remain explicit adaptations.
Runtime-only changes not saved in the supplied scripts cannot be reconstructed. Run
`scripts/check_release.R` and compare regenerated tables with the final study
outputs before tagging a final reproducibility release. Capture the R/package
versions from that verified run; no historical renv lockfile has been invented.

## Authorship
AI assistance was used to organise and adapt the existing code for repository
use. This does not establish that the authors have tested this new package.
Add an author-verification statement only after the authors complete that check.

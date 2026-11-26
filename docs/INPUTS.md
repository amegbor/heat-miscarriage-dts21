# Required private inputs

This release starts AFTER pregnancy selection, monthly environmental linkage,
calendar harmonisation, unit conversion and exposure construction. It does not
reconstruct raw DHS records or download exposures. Keep those upstream scripts
and their provenance with the private analytical project.

Set `DTS_ROOT` to a private analysis root containing:

- `data/DTS_ObservedRiskSet.rds`: primary cause-specific observed person-months.
- `data/DTS_PregnancyCohort_PrimaryCauseSpecific.rds`: one eligible pregnancy per
  row, four outcome categories, `outcome_last5yr` flag and descriptive covariates.
- `data/DTS_FullPredictionRiskSet_BinaryComparator.rds`: full six-interval binary
  comparator prediction trajectories with monthly heat, original country IQRs,
  raw/scaled cumulative heat and model covariates.
- `data/DTS_FullPredictionRiskSet_Mixture.rds`: full six-interval trajectories
  with cumulative monthly pollution histories already linked and validated.
- `qc/DTS_selection_flow.csv`: optional existing selection log for flowchart.

Observed risk-set fields used by models include `preg_uid`, `mother_uid`,
`cluster_uid`, `country`, `miscarriage_event`, `gest_month`, `gest_month_f`,
`preg_start_year_f`, `preg_start_month_f`, `age_at_preg_start`, `rc_empl`, `v106`,
`v190`, `rc_media_none`, `rc_hins_any`, `v025`, raw heat columns
`WBGT_T1_cummean`, `utci_T1_cummean`, `mrt_T1_cummean`, `HI_T1_cummean`, their
`*_T1_country_iqr` counterparts, and `CO_T1_cummean`, `DUST_T1_cummean`,
`PM_T1_cummean`, `SO_T1_cummean`, `SU_T1_cummean`. Descriptive scripts additionally
need the original cohort metadata and fields named in their `need()` checks.
The full datasets remain required; this list is not a minimal export schema.

Alternatively, figure stages can use completed aggregate tables in the same
folder hierarchy without model fits. `figures` expects primary linear tables,
all four raw-linear tables, all eight spline prediction/distribution tables,
country tables, five-year tables with year-pooling provenance, all eight
mixture-HR/RR tables, and four per-metric PIF/mixture result files.

Model RDS files can embed participant-level model frames. Do not upload them,
prepared RDS files, nonlinear distributions containing pregnancy identifiers,
individual-level QC exports, or precise cluster locations to GitHub. Generated
outputs are intentionally all ignored by default. Review any aggregate outputs
separately before explicitly sharing them.

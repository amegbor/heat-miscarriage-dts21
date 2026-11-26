# Validation performed during curation

- 26 R source files checked for closed strings/comments and balanced delimiters.
- Literal runner source targets resolve to included scripts; dependencies among
  original plotting and analysis scripts were inspected.
- No personal absolute HPC paths, account identifiers, data files or checkpoints
  remain in the release. Only code, configuration examples and documentation.
- Confirmed documented corrections: primary NULL-message/gradient handling,
  reviewed recovery tolerance, 50-node PIF default, Holm interaction collection,
  PIF/mixture labels, colour-bar title and flowchart caption.
- ZIP integrity and per-source checksums checked.
- Uploaded Torch inventory verified: 45 matching file checksums. All 24
  retained Torch sources compared; adaptations documented in
  `docs/TORCH_COMPARISON.md` and per-source hashes in `SOURCE_MANIFEST.tsv`.

**Not performed during curation:** R parsing, package loading, model fitting, numerical replication
or rendered figure comparison. R and private analytic data were unavailable in
this curation environment. The lexical delimiter check is not an R syntax test.
Run `Rscript scripts/check_release.R` in the project R environment, then verify
against the final Torch outputs before tagging a validated reproducibility release.
The scripts retain their original model diagnostics and stop-on-failure behaviour.

## Verification on Torch, 26 November 2026

Performed by the corresponding author on the NYU Torch HPC cluster, using the
release package `heat_miscarriage_dts21_github_2.zip` (code identical to tag
v0.1.3) and the private study outputs. Full records are in `docs/verification/`.

| Check | Result |
|---|---|
| `Rscript scripts/check_release.R`: all 26 R scripts parsed | PASSED |
| All 15 required packages installed and loadable | PASSED |
| `Rscript scripts/run.R collect_interaction`: country-by-heat interaction table regenerated from the stored tests and compared with the validated table (`all.equal`, tolerance 1e-10) | PASS (identical) |

The regenerated table reproduces the reported results: likelihood-ratio
statistics 327.085 (UTCI), 333.856 (MRT), 379.741 (WBGT) and 379.355 (Heat
Index), each on 20 df, all Holm-adjusted P < 0.001.

### Software environment
R 4.5.1 (2025-06-13), x86_64-pc-linux-gnu, Ubuntu 24.04.3 LTS. Package versions:
data.table 1.17.8, glmmTMB 1.1.12 (TMB 1.9.17), dplyr 1.2.1, tibble 3.3.0,
broom.mixed 0.2.9.7, readr 2.1.5, lme4 1.1-37, digest 0.6.37, Matrix 1.7-4,
qgcomp 2.18.10, sandwich 3.1-1, statmod 1.5.0, ggplot2 4.0.0, patchwork 1.3.2,
maps 3.4.3. Full `sessionInfo()` output: `docs/verification/verification_sessionInfo.txt`.

### Scope of this verification
This check confirms that the released code parses, that its dependencies are
available, and that the result-collection code reproduces a validated table.
It did not refit models: the full model fits take many hours on the HPC
cluster. The model-fitting scripts match the Torch scripts that produced the
reported results (see `docs/TORCH_COMPARISON.md`). These versions come from a
verification run in November 2026; the exact package versions used for the
original October 2026 model fits were not recorded.

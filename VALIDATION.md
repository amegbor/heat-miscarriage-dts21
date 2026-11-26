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

**Not performed:** R parsing, package loading, model fitting, numerical replication
or rendered figure comparison. R and private analytic data were unavailable in
this curation environment. The lexical delimiter check is not an R syntax test.
Run `Rscript scripts/check_release.R` in the project R environment, then verify
against the final Torch outputs before tagging a validated reproducibility release.
The scripts retain their original model diagnostics and stop-on-failure behaviour.

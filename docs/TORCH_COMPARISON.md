# Comparison with current Torch source

The user supplied `DTS_current_Torch_scripts.zip`. Its inventory covers 45
script files; all file checksums were verified against the supplied inventory.
Of these, 24 R source files are retained in the release. Two additional R files
are new portability/validation wrappers. Upstream preprocessing, job launchers
and finalisation/archive utilities are excluded from the requested scope.

## Findings and actions
- Main/supplementary model and PIF source agrees with the previously recovered
  analysis source after accounting for the previously documented diagnostic fixes.
- Torch's `Manuscript_figures.R` uses `local=TRUE` for three nested source calls;
  these were missing in the first curated version and are now incorporated.
- Torch's primary diagnostic extraction handles NULL optimizer messages and
  uses the stored fixed gradient. Both fixes are retained.
- Torch's reviewed MRT recovery likelihood tolerance is 0.01, already retained.
- No difference in the substantive model specifications was found in the
  retained analysis scripts. This is source comparison, not numeric reproduction.

## Explicit release adaptations
- Local paths and execution entry points are configurable; old job identifiers
  and hard-coded accounts are not retained.
- PIF defaults to 50 GH nodes rather than relying on an external environment
  override of the source's 25-node default. The higher-node check is unchanged.
- Figure text uses Potential impact fraction and one-quartile, adds the heatmap
  colour-bar title, and removes the obsolete sample-table footnote.
- Recovery scripts no longer require a particular root-directory basename or
  rewrite the primary source at runtime; diagnostic fixes are already applied.
- Country interaction checkpoint paths can be relocated by basename under the
  configured root, without changing fitted model-frame checks.
- New stage runner requires preconstructed private input datasets. It performs
  no raw DHS linkage or automatic batch submission.

## Figure 2 reference
The current Torch code explicitly sets `reference=median(terminal[[var]])`:
one final eligible observed exposure per pregnancy in the GM3-complete sample.
Spline knots instead come from GM3-complete person-months. Preserve that
analysis and make the caption match it. The reference check requires no refit.

## Remaining validation
Run the included R syntax/dependency check in the project R environment.
Regenerate a completed results table or figure and compare against the final
export before tagging a verified release. No private data or checkpoints were
provided here, and no R execution or model rerun is claimed.

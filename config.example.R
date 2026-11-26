# Copy to config.R; keep config.R private. Use an absolute analysis directory.
# That directory must contain the preconstructed data/ files described in docs/INPUTS.md.
Sys.setenv(DTS_ROOT = "/absolute/path/to/private/analysis_root")
# Optional optimizer override for raw/spline/country/mixed stages:
# Sys.setenv(DTS_OPTIMIZER = "BFGS")
Sys.setenv(PIF_GH_N = "50", DTS_GH_N = "50")

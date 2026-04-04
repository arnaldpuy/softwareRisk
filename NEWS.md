# softwareRisk 0.2.0

## Breaking changes

- `risk_form_variance_share_fun()` has been removed. Its functionality is fully
  subsumed by `uncertainty_fun()`, which already computed identical Sobol' sensitivity
  indices for all four parameters as part of its output. Users should call
  `uncertainty_fun()` directly and read the `sensitivity_analysis` list-column from
  `$nodes`.
- Removed the `risk_form` argument from `all_paths_fun()` and `uncertainty_fun()`.
  The additive formula (`risk_form = "additive"`) is now subsumed by the unified
  power-mean formula with `p = 1` (default).
- `uncertainty_fun()` now returns tibbles instead of `data.table` objects. The
  `data.table` package is no longer a dependency.

## Bug fixes

- Fixed duplicate `"p"` column bug in Sobol sample matrices within `uncertainty_fun()`.
  The raw column from `sobol_matrices()` is now named `"p_raw"` internally so the
  rescaled `p` column is correctly used.
- Added `@importFrom rlang %||%` to `plot_top_paths_fun()` (previously not imported).
- Added missing `@importFrom` declarations for `ggplot2::scale_fill_viridis_c`,
  `ggplot2::geom_tile`, `stats::setNames`, `ggplot2::geom_errorbar`,
  `ggplot2::geom_point`, and `scales::breaks_pretty`.
- Fixed fallback `nodes_tbl` in `plot_top_paths_fun()` that was missing `risk_score`.
  The function now requires `all_paths_out` to be the output of `all_paths_fun()`.
- Fixed `\deqn{[0, 1]}` in `@param` of `plot_top_paths_fun()`.

## Improvements

- `uncertainty_fun()` now returns Sobol' sensitivity indices labelled `alpha`,
  `beta`, `gamma`, and `p` in the `sensitivity_analysis` list-column of `$nodes`,
  replacing the previous internal raw labels `a_raw`, `b_raw`, `c_raw`, `p_raw`.
  Internally the Sobol' design still draws four independent `U(0,1)` values (required
  by the quasi-random sequence); the raw draws are normalised to weights and mapped to
  `p ∈ [−1, 2]` before the model is evaluated. The documentation explains this
  in a dedicated *Parameter labels and the Sobol' design* section.
- `slope_fun()` now uses an analytic OLS formula instead of `stats::lm()`.
- Removed the internal `slope_per_draw_from_nodes` closure in `uncertainty_fun()` and
  replaced it with calls to `slope_fun()`.
- Added shared input-validation helper `validate_all_paths_out()` used across
  multiple functions.
- Added `eps` parameter to `uncertainty_fun()` for forwarding to the risk evaluation.
- Corrected `@return` documentation for `path_cc` in `all_paths_fun()`.
- Updated `inst/CITATION` to version 0.2.0.
- Removed `NAMESPACE.bak`.
- Removed redundant `!is.na(x)` guard in `slope_fun()` and `gini_index_fun()`.
- Replaced three sequential `vapply` loops in `path_uncertainty_plot()` with a single pass.
- Removed trailing blank lines in `all_paths_fun.R`.
- Added `testthat` test suite with tests for `slope_fun`, `gini_index_fun`,
  `all_paths_fun`, and `uncertainty_fun`.

# softwareRisk 0.1.1

- Corrected the legend name of the `path_fix_heatmap` function (from R_k to P_k).
- Expanded the vignette to include examples with the power_mean formulation.

# softwareRisk 0.1.0

- Initial CRAN release.
- Functions to compute and visualize risk along software call paths.

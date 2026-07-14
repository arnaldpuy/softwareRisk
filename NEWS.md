# softwareRisk 0.3.0

## New features

- `call_graph_fun()` builds the call graph required by `all_paths_fun()`
  automatically from an installed `R` package (`pkg` argument) or a directory
  of `.R` scripts (`dir` argument). Calls are detected by static analysis
  (`codetools::findGlobals()`; the code is never executed) and cyclomatic
  complexity is computed internally by walking the abstract syntax tree of
  each function, following McCabe (1976).
- `read_call_graph()` imports and validates the edge-list and complexity
  tables described in the vignette (data frames or `.csv` files), failing
  early with informative messages on misspelled columns, functions missing
  from the complexity table, duplicated names or non-numeric complexity
  values. Recommended entry point for models written in other languages.
- `node_exposure_fun()` computes path-aware node criticality: for each node,
  the number of entry-to-sink paths containing it, its *risk load* (the sum
  of path risk over those paths) and the corresponding shares and ranks,
  separating chokepoints (structurally exposed functions) from hotspots
  (complex functions).
- `fix_portfolio_fun()` greedily selects the set of nodes whose fixing most
  reduces total path risk (or the risk of the riskiest path) under a budget
  of interventions, and returns the selection together with a
  diminishing-returns plot. For the total-risk objective the greedy solution
  carries the usual (1 - 1/e) submodular optimality guarantee.
- `rank_robustness_fun()` quantifies how stable the identification of the
  top-k riskiest paths (or nodes) is across the uncertainty draws of
  `uncertainty_fun()`: per-item top-k membership probability, rank quantiles
  and a consensus rank correlation. `rank_robustness_plot()` visualizes the
  result.
- `sensitivity_plot_fun()` visualizes the Sobol' indices stored in the
  `sensitivity_analysis` column of `uncertainty_fun()` output, either
  aggregated across all nodes (boxplots) or for selected nodes with their
  confidence intervals.

## Improvements

- The vignette gained sections on automatic call-graph construction, input
  validation with `read_call_graph()`, node exposure, budgeted fix
  portfolios, ranking robustness and Sobol'-index plotting, with worked
  examples and figures.
- `codetools` added to `Imports` (ships with `R`).
- Added test suites for all new functions.
- Fixed several typos in the vignette.

# softwareRisk 0.2.2

## Bug fixes

- Fixed silent `NaN` risk scores in `all_paths_fun()` and `uncertainty_fun()`
  when a weight is `0` and `p < 0`: zero-valued normalized metrics produced
  `0 * Inf = NaN`, which was then silently dropped from path-level aggregation
  via `na.rm = TRUE`. Zero-valued metrics are now replaced by `eps` before the
  power mean is evaluated when `p < 0`, matching the guard already used in the
  `p -> 0` (geometric mean) case.
- `gini_index_fun()` now returns `0` instead of `NaN` for all-zero input
  (`ineq::Gini()` divides by the mean). This also removes spurious `NaN` draws
  in the `gini_index` column returned by `uncertainty_fun()`, which occurred
  whenever a Sobol' draw pushed every node risk on a path to zero.
- `uncertainty_fun()` no longer errors when `all_paths_out$paths` is empty
  (the documented output of `all_paths_fun()` for graphs with no entry-to-sink
  paths). Node-level results are computed and an empty paths tibble with the
  usual columns is returned.
- `slope_fun()` now regresses the finite values against their original
  positions instead of a compressed index, so removing non-finite values no
  longer biases the estimated trend (e.g., `slope_fun(c(1, NA, 3))` now
  returns `1` instead of `2`).
- Fixed `geom_errorbar(height = ...)` in `path_uncertainty_plot()`: the
  argument is now `width`, removing the deprecation warning issued by
  ggplot2 (`height` was translated to `width`).
- Removed the documented but non-functional backward-compatibility input of
  `plot_top_paths_fun()`: supplying a bare `paths` tibble always failed with a
  misleading error because the fallback deriving node metrics from `graph` was
  never implemented. The function now clearly requires the full output of
  `all_paths_fun()`, and the documentation has been updated accordingly.

## Improvements

- Replaced the superseded `purrr::flatten()` with `unlist(recursive = FALSE)`
  in `all_paths_fun()`.
- Replaced the deprecated `trans` argument with `transform` in
  `scale_fill_viridis_c()` within `path_fix_heatmap()`; ggplot2 (>= 3.5.0) is
  now required.
- Removed the unused `params` argument of the internal `risk_ua_sa_fun()`.
- Removed a stale reference to `lab_expr` in the documentation of
  `plot_top_paths_fun()`.
- Updated the citation entries in `README.md` and `inst/CITATION` to the
  current version.
- Added regression tests for all fixes above.

# softwareRisk 0.2.1

## Bug fixes

- Fixed `scale_fill_manual()` in `plot_top_paths_fun()`: values are now named
  (`b1`–`b4`) to guarantee correct colour-to-category mapping, and `drop = FALSE`
  ensures all four complexity categories always appear in the legend even when
  some are absent from the plotted nodes.

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

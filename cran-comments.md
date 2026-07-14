## Submission of softwareRisk 0.3.0

This is a feature release. It adds automatic call-graph construction from `R`
source code (`call_graph_fun()`), a validating importer for externally
prepared graphs (`read_call_graph()`), path-aware node criticality
(`node_exposure_fun()`), budgeted selection of risk-reducing fixes
(`fix_portfolio_fun()`), ranking-robustness analysis under uncertainty
(`rank_robustness_fun()`, `rank_robustness_plot()`) and a plotting function
for Sobol' sensitivity indices (`sensitivity_plot_fun()`). See NEWS.md.

`codetools` (shipped with R) was added to Imports. No code is executed when
analyzing user source files: functions are parsed and inspected statically.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

- There are currently no downstream dependencies.

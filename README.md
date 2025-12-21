# softwareRisk: Computation of node and path-level risk scores in scientific models

[Arnald Puy](https://www.arnaldpuy.com/)

The `R` package `softwareRisk` leverages the network-like architecture of scientific models 
together with software quality metrics to identify chains of function calls that are more 
prone to generating and propagating errors. It operates on `tbl_graph` objects representing 
call dependencies between functions (callers and callees) and computes risk scores for individual 
functions and for paths (sequences of function calls) based on cyclomatic complexity, in-degree 
and betweenness centrality. The package supports variance-based uncertainty and sensitivity analyses 
to assess how risk scores change under alternative risk definitions.

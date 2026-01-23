[![](https://www.r-pkg.org/badges/version/softwareRisk?color=orange)](https://cran.r-project.org/package=softwareRisk)

# softwareRisk: Computation of node and path-level risk scores in scientific models

[Arnald Puy](https://www.arnaldpuy.com/)

The `R` package `softwareRisk` leverages the network-like architecture of scientific models 
together with software quality metrics to identify chains of function calls that are more 
prone to generating and propagating errors. It operates on `tbl_graph` objects representing 
call dependencies between functions (callers and callees) and computes risk scores for individual 
functions and for paths (sequences of function calls) based on cyclomatic complexity, in-degree 
and betweenness centrality. The package supports variance-based uncertainty and sensitivity analyses 
to assess how risk scores change under alternative risk definitions.

## Installation

To install the stable version on [CRAN](https://CRAN.R-project.org/package=softwareRisk), use

```r
install.packages("softwareRisk")
```
To install the development version, use devtools:

``` r
install.packages("devtools") # if you have not installed devtools package already
devtools::install_github("arnaldpuy/softwareRisk", build_vignettes = TRUE)
```

## Usage

Please see the vignette for a walkthrough of the package utilities. 

A printable PDF version of this 
vignette is also available [here](https://github.com/arnaldpuy/softwareRisk/raw/main/inst/extdata/vignette-pdf/softwareRisk.pdf).

## Citation 

To cite `softwareRisk` in publications:

Puy A (2025). _softwareRisk: Computation of node and path-level risk scores in scientific models_. R
package version 0.1.0, <https://github.com/arnaldpuy/softwareRisk>.

A BibTeX entry for LaTeX users is

@Manual{puy2025_softwareRisk,
  title = {softwareRisk: Computation of node and path-level risk scores in scientific models},
  author = {Arnald Puy},
  year = {2025},
  note = {R package version 0.1.0},
  url = {https://github.com/arnaldpuy/softwareRisk}
}



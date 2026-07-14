## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi = 200,
  fig.retina = 2
)

## ----setup, message=FALSE-----------------------------------------------------
library(softwareRisk)
library(tidygraph)

## -----------------------------------------------------------------------------

# Dataset 1: calls (edge list) -------------------------------------------------

calls_df <- data.frame(
  from = c("clean_data", "compute_risk", "compute_risk", "calc_scores", "plot_results"),
  to = c("load_data", "clean_data", "calc_scores", "mean", "compute_risk")
)

calls_df

# Dataset 2: cyclomatic complexity (node attributes) ---------------------------

cyclo_df <- data.frame(
  name  = c("clean_data", "load_data", "compute_risk", "calc_scores", "mean", "plot_results"),
  cyclo = c(6, 3, 12, 5, 2, 4)
)

cyclo_df

## ----merge--------------------------------------------------------------------

# Merge into a tbl_graph -------------------------------------------------------

graph <- tbl_graph(nodes = cyclo_df, edges = calls_df, directed = TRUE)

graph

## ----read_call_graph----------------------------------------------------------

# Merge and validate with read_call_graph ---------------------------------------

graph <- read_call_graph(edges = calls_df, metrics = cyclo_df)

graph

## ----call_graph_auto----------------------------------------------------------

# Write a small R model to a temporary directory ---------------------------------

td <- file.path(tempdir(), "toy_model")
dir.create(td, showWarnings = FALSE)

writeLines(c(
  "load_data <- function(x) x",
  "clean_data <- function(x) load_data(x)",
  "calc_scores <- function(x) if (length(x) > 0) mean(x) else 0",
  "compute_risk <- function(x) calc_scores(clean_data(x))",
  "plot_results <- function(x) compute_risk(x)"
), file.path(td, "model.R"))

# Build the call graph -----------------------------------------------------------

auto_graph <- call_graph_fun(dir = td)

auto_graph

## ----data_loading-------------------------------------------------------------

# Load the data ----------------------------------------------------------------

data("synthetic_graph")

# Print it ---------------------------------------------------------------------

synthetic_graph

## ----all_paths----------------------------------------------------------------

# Run the function -------------------------------------------------------------

output <- all_paths_fun(graph = synthetic_graph,
                        alpha = 0.6,
                        beta  = 0.3,
                        gamma = 0.1,
                        complexity_col = "cyclo",
                        weight_tol = 1e-8)

# Print the output -------------------------------------------------------------

output

## ----plot_heatmap, fig.height=2, fig.width=3----------------------------------

path_fix_heatmap(all_paths_out = output, n_nodes = 20, k_paths = 20)

## ----plot_paths, fig.height=2, fig.width=3.5----------------------------------

plot_output <- plot_top_paths_fun(graph = synthetic_graph,
                                  all_paths_out = output,
                                  model.name = "ToyModel",
                                  language = "Fortran",
                                  top_n = 10,
                                  alpha_non_top = 0.05)

## ----plot_paths2, fig.height=2, fig.width=3.5---------------------------------

plot_output <- plot_top_paths_fun(graph = synthetic_graph,
                                  all_paths_out = output,
                                  model.name = "ToyModel",
                                  language = "Fortran",
                                  top_n = 10,
                                  alpha_non_top = 1)

## ----node_exposure------------------------------------------------------------

# Compute node exposure ----------------------------------------------------------

exposure <- node_exposure_fun(output)

exposure

## ----fix_portfolio, fig.height=2, fig.width=3.5-------------------------------

# Greedy portfolio of five fixes ---------------------------------------------------

portfolio <- fix_portfolio_fun(output, budget = 5, objective = "total")

portfolio$portfolio

portfolio$plot

## ----uncertainty--------------------------------------------------------------

# Run uncertainty analysis -----------------------------------------------------

uncertainty_analysis <- uncertainty_fun(all_paths_out = output,
                                        N = 2^10,
                                        order = "first")

# Print the top five rows ------------------------------------------------------

lapply(uncertainty_analysis, function(x) head(x, 5))

## ----plot_uncert--------------------------------------------------------------

path_uncertainty_plot(ua_sa_out = uncertainty_analysis, n_paths = 20)

## ----rank_robustness----------------------------------------------------------

# Robustness of the top-10 path ranking --------------------------------------------

robustness <- rank_robustness_fun(uncertainty_analysis, top_k = 10, what = "paths")

robustness$summary

robustness$consensus_correlation

## ----rank_robustness_plot, fig.height=2.5, fig.width=3.5----------------------

rank_robustness_plot(robustness, top_n = 20)

## ----sa_single_node-----------------------------------------------------------

# Sobol' indices for the first node
si_node1 <- uncertainty_analysis$nodes$sensitivity_analysis[[1]]$results
si_node1

## ----sa_plot, fig.height=2.5, fig.width=4-------------------------------------

sensitivity_plot_fun(uncertainty_analysis)

## ----sa_plot_nodes, fig.height=2.5, fig.width=4.5-----------------------------

sensitivity_plot_fun(uncertainty_analysis,
                     nodes = c("M29", "M35", "M23", "M31"))


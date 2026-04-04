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

## ----uncertainty--------------------------------------------------------------

# Run uncertainty analysis -----------------------------------------------------

uncertainty_analysis <- uncertainty_fun(all_paths_out = output,
                                        N = 2^10,
                                        order = "first")

# Print the top five rows ------------------------------------------------------

lapply(uncertainty_analysis, function(x) head(x, 5))

## ----plot_uncert--------------------------------------------------------------

path_uncertainty_plot(ua_sa_out = uncertainty_analysis, n_paths = 20)

## ----sa_single_node-----------------------------------------------------------

# Sobol' indices for the first node
si_node1 <- uncertainty_analysis$nodes$sensitivity_analysis[[1]]$results
si_node1

## ----sa_combine---------------------------------------------------------------

sa_all <- do.call(rbind, Map(
  function(sa, nm) data.frame(sa$results, name = nm, stringsAsFactors = FALSE),
  uncertainty_analysis$nodes$sensitivity_analysis,
  uncertainty_analysis$nodes$name
))

head(sa_all)

## ----sa_plot, fig.height=2.5, fig.width=4-------------------------------------

library(ggplot2)

ggplot(sa_all, aes(x = parameters, y = original, fill = sensitivity)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(
    values   = c(Si = "#F8766D", Ti = "#00BFC4"),
    labels   = c(expression(S[i]), expression(T[i]))
  ) +
  labs(x = "Parameter", y = "Sobol\u2019 index", fill = "Index") +
  theme_bw()


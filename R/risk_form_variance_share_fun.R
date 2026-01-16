# Internal helper --------------------------------------------------------------

#' Evaluate node risk under sampled parameters
#'
#' Internal helper used by [risk_form_variance_share_fun()] to evaluate node risk
#' for a single node (with scaled inputs) across a sampled parameter matrix.
#'
#' For `risk_form = "additive"`, risk is computed as a weighted sum:
#'
#' \deqn{r = \alpha\,\tilde{C} + \beta\,\tilde{d}^{\mathrm{in}} + \gamma\,\tilde{b}\,,}
#'
#' where \eqn{\tilde{C}}, \eqn{\tilde{d}^{\mathrm{in}}}, and \eqn{\tilde{b}} are
#' scaled cyclomatic complexity, in-degree, and betweenness, respectively.
#'
#' For `risk_form = "power_mean"`, risk is computed as a weighted power mean with
#' exponent \eqn{p \in [-1, 2]}:
#'
#' \deqn{r =
#' \left(\alpha\,\tilde{C}^{p} + \beta\,(\tilde{d}^{\mathrm{in}})^{p} + \gamma\,\tilde{b}^{p}\right)^{1/p}\,.}
#'
#' In the limit \eqn{p \to 0}, the implementation uses the weighted geometric mean
#' with a small constant \eqn{\epsilon} to ensure numerical stability:
#'
#' \deqn{r = \exp\left(\alpha\log(\max(\tilde{C},\epsilon)) +
#' \beta\log(\max(\tilde{d}^{\mathrm{in}},\epsilon)) +
#' \gamma\log(\max(\tilde{b},\epsilon))\right)\,.}
#'
#' @param cyclo_sc,indeg_sc,btw_sc Numeric scalars. Scaled node metrics (typically
#'   in \eqn{[0,1]}) for cyclomatic complexity, in-degree, and betweenness.
#' @param sample_matrix Numeric matrix or data.frame containing sampled parameters.
#'   Must include `alpha`, `beta`, `gamma`. If `risk_form = "power_mean"`, must also
#'   include `p`.
#' @param risk_form Character. One of `"additive"` or `"power_mean"`.
#' @param eps Numeric. Small positive constant \eqn{\epsilon} used to avoid
#'   \eqn{\log(0)} in the \eqn{p \to 0} case. Default `1e-12`.
#'
#' @return A numeric vector of risk values, one per row of `sample_matrix`.
#'
#' @keywords internal
risk_eval_fun <- function(cyclo_sc, indeg_sc, btw_sc, sample_matrix,
                          risk_form = c("additive", "power_mean"),
                          eps = 1e-12) {

  risk_form <- match.arg(risk_form)

  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a single positive finite numeric value.", call. = FALSE)
  }

  if (risk_form == "additive") {
    w <- sample_matrix[, c("alpha", "beta", "gamma"), drop = FALSE]
    return(w[, 1] * cyclo_sc + w[, 2] * indeg_sc + w[, 3] * btw_sc)
  }

  if (!("p" %in% colnames(sample_matrix))) {
    stop("`sample_matrix` must contain column `p` for `risk_form='power_mean'`.",
         call. = FALSE)
  }

  w <- sample_matrix[, c("alpha", "beta", "gamma"), drop = FALSE]
  p <- sample_matrix[, "p"]

  # enforce p in [-1,2]
  if (any(!is.finite(p)) || any(p < -1 | p > 2)) {
    stop("`p` must be finite and lie in [-1, 2].", call. = FALSE)
  }

  y <- numeric(nrow(sample_matrix))
  near0 <- abs(p) < 1e-12

  # p -> 0: weighted geometric mean limit
  if (any(near0)) {
    y[near0] <- exp(
      w[near0, 1] * log(pmax(cyclo_sc, eps)) +
        w[near0, 2] * log(pmax(indeg_sc, eps)) +
        w[near0, 3] * log(pmax(btw_sc, eps))
    )
  }

  # p != 0: power mean
  if (any(!near0)) {
    pp <- p[!near0]
    y[!near0] <- (w[!near0, 1] * (cyclo_sc ^ pp) +
                    w[!near0, 2] * (indeg_sc ^ pp) +
                    w[!near0, 3] * (btw_sc   ^ pp)) ^ (1 / pp)
  }

  pmin(1, pmax(0, y))
}


#' Proportion of risk-score variance due to functional-form uncertainty
#'
#' Computes, for each node, the proportion of variance in the node risk score that
#' is attributable to uncertainty in the *functional form* used to construct the
#' risk score (additive vs power-mean), as opposed to uncertainty in the parameters
#' within a given form.
#'
#' The function draws parameter samples for two model forms:
#' \itemize{
#'   \item Additive form with weights \eqn{(\alpha,\beta,\gamma)} normalized to sum to one.
#'   \item Power-mean form with weights \eqn{(\alpha,\beta,\gamma)} normalized to sum to one
#'         and a power parameter \eqn{p} mapped to the interval \eqn{[-1,2]}.
#' }
#'
#' For each node, it evaluates the risk under both forms and applies the variance
#' decomposition for a two-component mixture model:
#'
#' \deqn{\mathrm{Var}(Y) = \mathbb{E}\left[\mathrm{Var}(Y \mid M)\right] + \mathrm{Var}\left(\mathbb{E}[Y \mid M]\right)\,,}
#'
#' where \eqn{M} indicates the functional form. Let \eqn{q} be the probability of choosing
#' the power-mean form (and \eqn{1-q} the probability of choosing the additive form). The
#' between-form component is
#'
#' \deqn{V_{\mathrm{between}} = q(1-q)\,(\mu_1 - \mu_0)^2\,,}
#'
#' and the total variance is
#'
#' \deqn{V_{\mathrm{total}} = (1-q)V_0 + qV_1 + V_{\mathrm{between}}\,,}
#'
#' where \eqn{\mu_0, V_0} are the mean and variance under the additive form and
#' \eqn{\mu_1, V_1} are the mean and variance under the power-mean form. The reported
#' proportion is \eqn{V_{\mathrm{between}} / V_{\mathrm{total}}}.
#'
#' @param all_paths_out A list produced by [all_paths_fun()] with element `nodes`.
#'   `nodes` must contain columns `name`, `cyclomatic_complexity`, `indeg`, and `btw`.
#' @param N Integer. Base sample size passed to `sensobol::sobol_matrices()`.
#' @param q Numeric scalar in \eqn{[0,1]}. Mixture weight giving the probability that the
#'   functional form is `"power_mean"` (and \eqn{1-q} for `"additive"`). Default `0.5`.
#' @param eps Numeric. Small positive constant \eqn{\epsilon} used for numerical stability
#'   in the \eqn{p \to 0} evaluation. Default `1e-12`.
#'
#' @details
#' This function focuses only on variance attributable to the *choice of functional form*.
#' It does not return full uncertainty draws, Sobol sensitivity indices, or path-level
#' uncertainty propagation.
#'
#' @return A `data.table` with columns:
#' \describe{
#'   \item{name}{Node name.}
#'   \item{form_variance_share}{Proportion of node-risk variance attributable to functional-form uncertainty.}
#' }
#'
#' @examples
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, complexity_col = "cyclo")
#' vf <- risk_form_variance_share_fun(all_paths_out = out, N = 2^10, q = 0.5)
#' vf
#'
#' @export
#' @importFrom data.table as.data.table
#' @importFrom data.table set
#' @importFrom scales rescale
#' @importFrom sensobol sobol_matrices
risk_form_variance_share_fun <- function(all_paths_out, N,
                                         q = 0.5,
                                         eps = 1e-12) {
  # ---- validate input --------------------------------------------------------
  if (!is.list(all_paths_out) || !all(c("nodes", "paths") %in% names(all_paths_out))) {
    stop("`all_paths_out` must be the output of all_paths_fun() (a list with $nodes and $paths).",
         call. = FALSE)
  }

  nodes_tbl <- all_paths_out$nodes

  required_nodes <- c("name", "cyclomatic_complexity", "indeg", "btw")
  missing_nodes <- setdiff(required_nodes, names(nodes_tbl))
  if (length(missing_nodes) > 0) {
    stop("Missing node cols: ", paste(missing_nodes, collapse = ", "), call. = FALSE)
  }

  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0) {
    stop("`N` must be a single positive integer-like value.", call. = FALSE)
  }
  N <- as.integer(N)

  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q < 0 || q > 1) {
    stop("`q` must be a single number in [0,1].", call. = FALSE)
  }

  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a single positive finite numeric value.", call. = FALSE)
  }

  # ---- rescale node metrics --------------------------------------------------
  node_dt <- data.table::as.data.table(nodes_tbl)
  data.table::set(node_dt, j = "cyclo_sc", value = scales::rescale(node_dt[["cyclomatic_complexity"]]))
  data.table::set(node_dt, j = "indeg_sc", value = scales::rescale(node_dt[["indeg"]]))
  data.table::set(node_dt, j = "btw_sc",   value = scales::rescale(node_dt[["btw"]]))

  # ---- Sobol designs (defaults for matrices/order come from sensobol) --------
  # Additive: a_raw,b_raw,c_raw
  params_add <- c("a_raw", "b_raw", "c_raw")
  mat_add_raw <- sensobol::sobol_matrices(N = N, params = params_add)
  s_add <- rowSums(mat_add_raw)
  mat_add <- cbind(
    mat_add_raw,
    alpha = mat_add_raw[, "a_raw"] / s_add,
    beta  = mat_add_raw[, "b_raw"] / s_add,
    gamma = mat_add_raw[, "c_raw"] / s_add
  )

  # Power mean: a_raw,b_raw,c_raw,p_raw (U(0,1)) -> p in U(-1,2)
  params_pm <- c("a_raw", "b_raw", "c_raw", "p")
  mat_pm_raw <- sensobol::sobol_matrices(N = N, params = params_pm)
  s_pm <- rowSums(mat_pm_raw[, c("a_raw", "b_raw", "c_raw"), drop = FALSE])
  mat_pm <- cbind(
    mat_pm_raw,
    alpha = mat_pm_raw[, "a_raw"] / s_pm,
    beta  = mat_pm_raw[, "b_raw"] / s_pm,
    gamma = mat_pm_raw[, "c_raw"] / s_pm,
    p     = -1 + 3 * mat_pm_raw[, "p"]
  )

  # ---- per node: variance share due to functional form -----------------------
  Cf_list <- numeric(nrow(node_dt))

  for (i in seq_len(nrow(node_dt))) {

    cy <- node_dt[["cyclo_sc"]][i]
    di <- node_dt[["indeg_sc"]][i]
    bt <- node_dt[["btw_sc"]][i]

    y0 <- risk_eval_fun(cy, di, bt, mat_add, risk_form = "additive",  eps = eps)
    y1 <- risk_eval_fun(cy, di, bt, mat_pm,  risk_form = "power_mean", eps = eps)

    # mixture variance decomposition:
    # Var(Y) = E[Var(Y|M)] + Var(E[Y|M])
    mu0 <- mean(y0); V0 <- stats::var(y0)
    mu1 <- mean(y1); V1 <- stats::var(y1)

    V_between <- q * (1 - q) * (mu1 - mu0)^2
    V_total   <- (1 - q) * V0 + q * V1 + V_between

    Cf_list[i] <- if (is.finite(V_total) && V_total > 1e-16) V_between / V_total else NA_real_
  }

  data.table::set(node_dt, j = "form_variance_share", value = Cf_list)
  node_dt[, c("name", "form_variance_share"), with = FALSE]
}


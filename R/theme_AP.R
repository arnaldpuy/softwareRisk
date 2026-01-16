#' A clean, publication-oriented ggplot2 theme
#'
#' `theme_AP()` provides a minimalist, publication-ready theme based on
#' [ggplot2::theme_bw()], with grid lines removed, compact legends, and
#' harmonized text sizes. It is designed for dense network and path-visualization
#' plots (e.g. call graphs, risk paths).
#'
#' @details
#' The theme:
#' \itemize{
#'   \item removes major and minor grid lines,
#'   \item uses transparent legend backgrounds and keys,
#'   \item standardizes text sizes for axes, legends, strips, and titles,
#'   \item reduces legend spacing for compact layouts.
#' }
#'
#' This theme is intended to be composable:
#' it should be added to a ggplot object using `+ theme_AP()`.
#'
#' @return
#' A `ggplot2::theme` object.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
#'   ggplot2::geom_point() +
#'   theme_AP()
#'
#'
#' @export
#' @importFrom ggplot2 theme_bw theme element_blank element_rect element_text
#' @importFrom grid unit
theme_AP <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      legend.background = ggplot2::element_rect(
        fill = "transparent", colour = NA
      ),
      legend.key = ggplot2::element_rect(
        fill = "transparent", colour = NA
      ),

      strip.background = ggplot2::element_rect(fill = "white"),

      legend.text  = ggplot2::element_text(size = 7.3),
      legend.title = ggplot2::element_text(size = 7.3),

      axis.text.x  = ggplot2::element_text(size = 7),
      axis.text.y  = ggplot2::element_text(size = 7),
      axis.title.x = ggplot2::element_text(size = 7.3),
      axis.title.y = ggplot2::element_text(size = 7.3),

      plot.title = ggplot2::element_text(size = 8),

      strip.text.x = ggplot2::element_text(size = 7.4),
      strip.text.y = ggplot2::element_text(size = 7.4),

      legend.key.width  = grid::unit(0.4, "cm"),
      legend.key.height = grid::unit(0.4, "cm"),
      legend.key.spacing.y = grid::unit(0, "lines"),
      legend.box.spacing = grid::unit(0, "pt")
    )
}

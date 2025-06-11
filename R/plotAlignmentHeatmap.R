#' MiCoAl Heatmap
#'
#' @param object An object of class MiCoAl.
#' @param top Character scalar giving the top fraction of consortia in which a
#' pathway should be present to be output in the heatmap. Defaults to NULL, in
#' which case all are plotted.
#' @param bottom Character scalar giving the bottom fraction of consortia in
#' which a pathway should be present to be output in the heatmap. Defaults to
#' NULL, in which case all are plotted.
#'
#' @return A ggplot heatmap
#' @export
plotAlignmentHeatmap <- function(object, top = NULL, bottom = NULL) {
  # Filter the adjacency matrix for desired levels for visualization.
  levels_mat <- assays(object)$Levels |> as.matrix()
  max_weight <- length(object@Consortia)

  colnames(levels_mat) <- SummarizedExperiment::colData(object)$met
  rownames(levels_mat) <- SummarizedExperiment::colData(object)$met

  tb <- levels_mat %>%
    tibble::as_tibble(rownames = "RowName") %>%
    tidyr::pivot_longer(
      cols = -"RowName",
      names_to = "ColName",
      values_to = "level"
    ) %>%
    dplyr::rename(
      met = "RowName",
      met2 = "ColName"
    )

  if (!is.null(top) && !is.null(bottom)) {
    stop("Only top or bottom fraction")
  }
  if (!is.null(top)) {
    min_weight <- max_weight * top
    if (min_weight > max(levels_mat)) {
      return("No reactions present in this top fraction of consortia")
    }
    tb <- dplyr::filter(tb, .data$level >= min_weight)
  }

  if (!is.null(bottom)) {
    max_weight <- max_weight * bottom
    tb <- dplyr::filter(tb, .data$level <= max_weight)
  }

  gg <- tb %>%
    ggplot2::ggplot(
      ggplot2::aes(x = .data$met2, y = .data$met, fill = .data$level)
    ) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "red",
      breaks = scales::breaks_pretty()
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0),
      axis.text.y = ggplot2::element_text(vjust = 0.4),
      axis.title = ggplot2::element_blank(),
      # panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )

  gg
}

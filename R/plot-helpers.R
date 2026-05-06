#' @include theme-ramen.R
NULL

#' Render a hierarchical-clustering dendrogram as a ggplot
#'
#' Internal helper shared by \code{plot(ConsortiumMetabolismSet)} and
#' \code{plotFunctionalGroups()}. Draws the dendrogram branches with
#' \code{ggplot2::geom_segment}, leaf labels via a manual
#' \code{geom_text} layer rotated 90 degrees, and applies
#' \code{theme_ramen()} for a consistent look-and-feel across the
#' package.
#'
#' @param dend An object coercible by \code{dendextend::as.ggdend()}
#'   (typically of class \code{dendrogram} or
#'   \code{hclust}-derived).
#' @param labelSize Numeric. Font size for leaf labels (in
#'   \code{geom_text} units). Defaults to 3.
#' @param labelColours Optional data frame with columns \code{label}
#'   and \code{colour} mapping leaf labels to colours. When
#'   \code{NULL} (default), all leaves are drawn in black.
#' @param showHeightAxis Logical. If \code{TRUE} (default), draw a
#'   numeric "Height" y-axis; if \code{FALSE}, suppress all axis
#'   marks.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle.
#'
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
#' @importFrom rlang .data
.dendrogramGgplot <- function(
    dend,
    labelSize = 3,
    labelColours = NULL,
    showHeightAxis = TRUE,
    title = NULL,
    subtitle = NULL
) {
    gg <- dendextend::as.ggdend(dend)

    if (!is.null(labelColours)) {
        labelTb <- gg$labels |>
            dplyr::left_join(labelColours, by = "label")
    } else {
        labelTb <- gg$labels |>
            dplyr::mutate(colour = "black")
    }

    ## Use the segment colour computed by dendextend (NA -> grey20),
    ## so dendextend::color_branches() output flows through.
    seg <- gg$segments
    if (is.null(seg$col)) {
        seg$col <- "grey20"
    } else {
        seg$col <- ifelse(is.na(seg$col), "grey20", seg$col)
    }

    p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = seg,
            ggplot2::aes(
                x = .data$x,
                y = .data$y,
                xend = .data$xend,
                yend = .data$yend
            ),
            colour = seg$col,
            linewidth = 0.5
        ) +
        ggplot2::geom_text(
            data = labelTb,
            ggplot2::aes(
                x = .data$x,
                y = 0,
                label = .data$label
            ),
            angle = 90,
            hjust = 1.05,
            vjust = 0.5,
            size = labelSize,
            colour = labelTb$colour
        ) +
        ggplot2::coord_cartesian(clip = "off")

    yMax <- max(seg$y, na.rm = TRUE)
    if (isTRUE(showHeightAxis)) {
        p <- p +
            ggplot2::scale_y_continuous(
                name = "Height",
                limits = c(0, yMax * 1.05),
                expand = ggplot2::expansion(mult = c(0, 0.02))
            ) +
            ggplot2::scale_x_continuous(
                expand = ggplot2::expansion(mult = c(0.02, 0.02))
            ) +
            theme_ramen() + # nolint: object_usage_linter.
            ggplot2::theme(
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.title.x = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_line(colour = "grey30"),
                axis.ticks.length.y = grid::unit(-0.15, "cm"),
                axis.line.y = ggplot2::element_line(colour = "grey30"),
                panel.grid = ggplot2::element_blank(),
                legend.position = "none",
                plot.margin = ggplot2::margin(12, 12, 60, 12)
            )
    } else {
        p <- p +
            ggplot2::scale_y_continuous(
                limits = c(0, yMax * 1.05),
                expand = ggplot2::expansion(mult = c(0, 0.02))
            ) +
            theme_ramen(network = TRUE) + # nolint: object_usage_linter.
            ggplot2::theme(
                legend.position = "none",
                plot.margin = ggplot2::margin(12, 12, 60, 12)
            )
    }

    if (!is.null(title) || !is.null(subtitle)) {
        p <- p + ggplot2::labs(title = title, subtitle = subtitle)
    }
    p
}

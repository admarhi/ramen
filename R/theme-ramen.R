#' Ramen plot theme
#'
#' Shared \pkg{ggplot2} theme used across every plot returned by the
#' package. Provides a consistent typography, white background, bottom
#' legend, and tight title spacing. Network plots pass
#' \code{network = TRUE} to drop axes, ticks, and grid lines while
#' keeping the rest of the look-and-feel.
#'
#' @param baseSize Numeric. Base font size in points. Defaults to 11.
#' @param baseFamily Character. Base font family. Defaults to \code{""}
#'   (device default).
#' @param network Logical. If \code{TRUE}, build on top of
#'   \code{theme_void()} instead of \code{theme_minimal()}; used by
#'   \code{\link{plotDirectedFlow}} and \code{plot} methods that render
#'   networks.
#'
#' @return A \code{ggplot2} theme object.
#' @export
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, hp)) +
#'     geom_point() +
#'     theme_ramen()
theme_ramen <- function(baseSize = 11, baseFamily = "", network = FALSE) {
    base <- if (isTRUE(network)) {
        ggplot2::theme_void(base_size = baseSize, base_family = baseFamily)
    } else {
        ggplot2::theme_minimal(
            base_size = baseSize,
            base_family = baseFamily
        )
    }

    base +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold",
                hjust = 0,
                size = ggplot2::rel(1.1),
                margin = ggplot2::margin(b = 4)
            ),
            plot.subtitle = ggplot2::element_text(
                colour = "grey30",
                size = ggplot2::rel(0.95),
                margin = ggplot2::margin(b = 8)
            ),
            plot.caption = ggplot2::element_text(
                colour = "grey40",
                size = ggplot2::rel(0.8),
                hjust = 1
            ),
            plot.margin = ggplot2::margin(12, 12, 8, 12),
            plot.background = ggplot2::element_rect(
                fill = "white",
                colour = NA
            ),
            panel.background = ggplot2::element_rect(
                fill = "white",
                colour = NA
            ),
            panel.grid.minor = ggplot2::element_blank(),
            legend.background = ggplot2::element_rect(
                fill = "white",
                colour = NA
            ),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.title = ggplot2::element_text(
                face = "bold",
                size = ggplot2::rel(0.9)
            ),
            legend.text = ggplot2::element_text(size = ggplot2::rel(0.85)),
            legend.key.height = grid::unit(0.4, "cm"),
            legend.key.width = grid::unit(0.8, "cm")
        )
}

#' Ramen package palette
#'
#' Named list of colour vectors used by every plotting function in
#' \pkg{ramen}. Centralising the palette here ensures that every
#' figure produced by the package shares the same colour vocabulary.
#' Where possible, colours are drawn from the Okabe-Ito palette
#' (colour-blind safe) and from ColorBrewer sequential ramps.
#'
#' @format A named \code{list} with elements:
#' \describe{
#'   \item{\code{nodeRole}}{Three colours mapping
#'     \code{source}/\code{intermediate}/\code{sink}.}
#'   \item{\code{edgeCategorical}}{Three colours mapping
#'     \code{shared}/\code{query}/\code{reference} edge categories
#'     in pairwise alignment networks.}
#'   \item{\code{edgeWeight}}{Two-colour ramp (\code{low}/\code{high})
#'     for continuous edge-weight gradients.}
#'   \item{\code{heatmapFill}}{Three-stop ramp
#'     (\code{low}/\code{mid}/\code{high}) for similarity heatmaps.}
#'   \item{\code{bar}}{Single colour for bar fills.}
#'   \item{\code{cluster}}{Eight-colour Okabe-Ito-derived palette,
#'     recycled for k-cluster dendrogram colouring.}
#' }
#'
#' @return A named \code{list} of character vectors.
#' @export
#' @examples
#' ramenPalette$nodeRole
#' ramenPalette$edgeCategorical
ramenPalette <- list(
    nodeRole = c(
        source = "#0072B2",
        intermediate = "#F0E442",
        sink = "#D55E00"
    ),
    edgeCategorical = c(
        shared = "#000000",
        query = "#0072B2",
        reference = "#D55E00"
    ),
    edgeWeight = c(low = "#DEEBF7", high = "#08519C"),
    heatmapFill = c(
        low = "#FFFFFF",
        mid = "#FDD49E",
        high = "#B30000"
    ),
    bar = "#0072B2",
    cluster = c(
        "#0072B2",
        "#D55E00",
        "#009E73",
        "#CC79A7",
        "#F0E442",
        "#56B4E9",
        "#E69F00",
        "#000000"
    )
)

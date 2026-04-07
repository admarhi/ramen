#' @title Plot Functional Groups Dendrogram
#'
#' @description
#' Visualizes the functional groups dendrogram computed by
#' \code{\link{functionalGroups}}. Branches are colored by
#' cluster membership and species labels are placed below
#' the dendrogram.
#'
#' @param fg A list as returned by
#'   \code{\link{functionalGroups}}, containing at least a
#'   \code{dendrogram} element.
#' @param k Integer scalar specifying the number of
#'   clusters to color in the dendrogram. Default
#'   \code{4}.
#' @param label_size Numeric scalar specifying the text
#'   size for species labels. Default \code{6}.
#' @param label_colours If not \code{NULL}, a data frame
#'   with columns \code{label} and \code{colour} mapping
#'   species names to colors. When \code{NULL} (default),
#'   all labels are drawn in black.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' \donttest{
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(
#'     cm1, cm2, name = "test"
#' )
#' fg <- functionalGroups(cms)
#' plotFunctionalGroups(fg, k = 2)
#' }
#'
#' @seealso \code{\link{functionalGroups}} for computing
#'   the functional groups data.
#'
#' @importFrom rlang .data
plotFunctionalGroups <- function(
    fg,
    k = 4,
    label_size = 6,
    label_colours = NULL
) {
    if (!is.list(fg) || is.null(fg$dendrogram)) {
        cli::cli_abort(
            c(
                "{.arg fg} must be a list returned by",
                " {.fun functionalGroups},",
                " containing a {.field dendrogram}",
                " element."
            )
        )
    }

    dend <- fg$dendrogram

    # Color branches by cluster membership
    gg_dend <- dendextend::color_branches(
        dend,
        k = k
    ) |>
        dendextend::as.ggdend()

    if (!is.null(label_colours)) {
        label_tb <- gg_dend$labels |>
            dplyr::left_join(
                label_colours,
                by = "label"
            )
    } else {
        label_tb <- gg_dend$labels |>
            dplyr::mutate(colour = "black")
    }
    # Remove the default label layer
    gg_dend$labels <- gg_dend$labels[0, ]

    # Build ggplot dendrogram
    plot_obj <- ggplot2::ggplot(
        gg_dend,
        horiz = FALSE
    ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            panel.background =
                ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        ) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(
                mult = c(0, 0),
                add = c(2, 1)
            )
        ) +
        ggplot2::geom_text(
            data = label_tb,
            ggplot2::aes(
                x = .data$x,
                y = .data$y - 0.1,
                label = .data$label
            ),
            angle = 90,
            hjust = 1,
            size = label_size,
            color = label_tb$colour
        )

    plot_obj
}

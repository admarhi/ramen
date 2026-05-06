#' @include theme-ramen.R plot-helpers.R
NULL

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
#'   size for species labels. Default \code{3}.
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
#' cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
#' cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
#' cms <- ConsortiumMetabolismSet(
#'     cm1, cm2, name = "test"
#' )
#' fg <- functionalGroups(cms)
#' plotFunctionalGroups(fg, k = 2)
#'
#' @seealso \code{\link{functionalGroups}} for computing
#'   the functional groups data.
#'
#' @importFrom rlang .data
plotFunctionalGroups <- function(
    fg,
    k = 4,
    label_size = 3,
    label_colours = NULL
) {
    if (!is.list(fg)) {
        cli::cli_abort(
            c(
                "{.arg fg} must be a list returned by",
                " {.fun functionalGroups},",
                " containing a {.field dendrogram}",
                " element."
            )
        )
    }
    if (is.null(fg$dendrogram)) {
        cli::cli_abort(
            c(
                "{.arg fg} contains no dendrogram to plot.",
                "i" = "{.fun functionalGroups} returns \\
                {.field dendrogram = NULL} when fewer than \\
                two species are present; nothing to cluster."
            )
        )
    }

    dend <- fg$dendrogram
    n_leaves <- length(labels(dend))

    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1L) {
        cli::cli_abort(
            "{.arg k} must be a single positive integer."
        )
    }
    k <- as.integer(k)
    if (k > n_leaves) {
        cli::cli_abort(
            c(
                "{.arg k} ({k}) exceeds the number of species \\
                in the dendrogram ({n_leaves}).",
                "i" = "Pick {.arg k} between 1 and {n_leaves}."
            )
        )
    }

    # Colour branches by cluster membership using the ramen palette
    # (Okabe-Ito-derived, recycled if k > length(palette)).
    clusterCols <- ramenPalette$cluster # nolint: object_usage_linter.
    if (k > length(clusterCols)) {
        clusterCols <- rep_len(clusterCols, k)
    }
    dend <- dendextend::color_branches(
        dend,
        k = k,
        col = clusterCols[seq_len(k)]
    )

    # nolint next: object_usage_linter.
    .dendrogramGgplot(
        dend,
        labelSize = label_size,
        labelColours = label_colours,
        showHeightAxis = TRUE,
        title = "Functional groups",
        subtitle = sprintf(
            "k = %d clusters · %d species",
            k,
            n_leaves
        )
    )
}

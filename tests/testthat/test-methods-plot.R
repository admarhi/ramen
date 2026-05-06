## ---- plot(CMA): heatmap (multiple) ----------------------------------------

test_that("plot(cma_multiple) returns ggplot heatmap", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    cma <- align(cms)
    p <- plot(cma)
    expect_s3_class(p, "ggplot")
})

test_that("heatmap explicit type works", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    p <- plot(cma, type = "heatmap")
    expect_s3_class(p, "ggplot")
})

test_that("heatmap has correct number of tiles", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    cma <- align(cms)
    p <- plot(cma, type = "heatmap")
    ## 3x3 = 9 tiles
    expect_equal(nrow(p$data), 9L)
})

test_that("heatmap errors for pairwise CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(
        plot(cma, type = "heatmap"),
        "multiple"
    )
})

## ---- plot(CMA): scores (both types) --------------------------------------

test_that("scores plot works for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    p <- plot(cma, type = "scores")
    expect_s3_class(p, "ggplot")
})

test_that("scores plot works for multiple", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    p <- plot(cma, type = "scores")
    expect_s3_class(p, "ggplot")
})

## ---- plot(CMA): network (pairwise) ---------------------------------------

test_that("network plot returns ggplot for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    p <- plot(cma, type = "network")
    expect_s3_class(p, "ggplot")
})

test_that("network plot has Okabe-Ito categorical legend", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    p <- plot(cma, type = "network")
    ## The categorical edge-colour scale should be present and titled.
    scale_titles <- vapply(
        p$scales$scales,
        function(s) {
            tt <- s$name
            if (
                is.null(tt) ||
                    inherits(tt, "waiver") ||
                    length(tt) != 1L ||
                    isTRUE(is.na(tt))
            ) {
                ""
            } else {
                as.character(tt)
            }
        },
        character(1L)
    )
    expect_true("Pathway type" %in% scale_titles)
})

test_that("plot(cm) returns a ggplot object", {
    cm <- synCM("test", n_species = 3, max_met = 5, seed = 1)
    p <- plot(cm, type = "Binary")
    expect_s3_class(p, "ggplot")
})

test_that("plot(cm) weighted has edge weight gradient legend", {
    cm <- synCM("test", n_species = 4, max_met = 6, seed = 1)
    p <- plot(cm, type = "nSpecies")
    expect_s3_class(p, "ggplot")
    scale_titles <- vapply(
        p$scales$scales,
        function(s) {
            tt <- s$name
            if (
                is.null(tt) ||
                    inherits(tt, "waiver") ||
                    length(tt) != 1L ||
                    isTRUE(is.na(tt))
            ) {
                ""
            } else {
                as.character(tt)
            }
        },
        character(1L)
    )
    expect_true("Weight" %in% scale_titles)
})

test_that("network errors for multiple CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    expect_error(
        plot(cma, type = "network"),
        "pairwise"
    )
})

## ---- plot(CMA): invalid type ---------------------------------------------

test_that("plot errors on invalid type", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(plot(cma, type = "invalid"))
})

## ---- B2: plotDirectedFlow preserves pre-set edge colours -------------------

test_that("plotDirectedFlow preserves pre-set E(g)$color", {
    g <- igraph::make_ring(4, directed = TRUE)
    igraph::E(g)$color <- c("red", "blue", "green", "red")
    expect_no_error(
        plotDirectedFlow(g, colourEdgesByWeight = FALSE)
    )
})

test_that("plotDirectedFlow defaults to gray when no pre-set colours", {
    g <- igraph::make_ring(4, directed = TRUE)
    expect_no_error(
        plotDirectedFlow(g, colourEdgesByWeight = FALSE)
    )
})

## ---- B8: scores plot filters non-score values ------------------------------

test_that("scores plot for multiple alignment excludes nPairs/sd", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    cma <- align(cms)
    p <- plot(cma, type = "scores")
    expect_s3_class(p, "ggplot")
    expect_false("nPairs" %in% as.character(p$data$metric))
    expect_false("sd" %in% as.character(p$data$metric))
})

## ---- Layout fidelity: 3-column structure ----------------------------------

test_that("plotDirectedFlow places sources / sinks at requested x", {
    ## Build a directed graph with one source-only node (s),
    ## two intermediates (m1, m2), and one sink-only node (t).
    edges <- data.frame(
        from = c("s", "s", "m1", "m2"),
        to = c("m1", "m2", "m2", "t")
    )
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)

    p <- plotDirectedFlow(g, sourceX = -2, mixedX = 1.5, sinkX = 5)
    expect_s3_class(p, "ggplot")

    lay <- p$data
    expect_true(all(c("name", "x", "y") %in% names(lay)))

    src_x <- lay$x[lay$name == "s"]
    sink_x <- lay$x[lay$name == "t"]
    mid_x <- lay$x[lay$name %in% c("m1", "m2")]

    expect_equal(src_x, -2)
    expect_equal(sink_x, 5)
    ## Intermediate nodes must lie strictly between source and sink x,
    ## guarding against a regression to a generic stress / kk / fr layout.
    expect_true(all(mid_x > -2 & mid_x < 5))
})

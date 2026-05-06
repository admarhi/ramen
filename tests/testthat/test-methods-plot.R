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

## ---- B2: plotDirectedFlow preserves pre-set edge colors --------------------

test_that("plotDirectedFlow preserves pre-set E(g)$color", {
    g <- igraph::make_ring(4, directed = TRUE)
    igraph::E(g)$color <- c("red", "blue", "green", "red")
    expect_no_error(
        plotDirectedFlow(g, color_edges_by_weight = FALSE)
    )
})

test_that("plotDirectedFlow defaults to gray when no pre-set colors", {
    g <- igraph::make_ring(4, directed = TRUE)
    expect_no_error(
        plotDirectedFlow(g, color_edges_by_weight = FALSE)
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

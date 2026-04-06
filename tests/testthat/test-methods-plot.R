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

test_that("network plot runs without error for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_no_error(plot(cma, type = "network"))
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

## ---- CMA class construction and validity ------------------------------------

test_that("CMA class can be instantiated with defaults", {
    cma <- ConsortiumMetabolismAlignment()
    expect_s4_class(cma, "ConsortiumMetabolismAlignment")
    expect_true(is.na(cma@Name))
    expect_true(is.na(cma@Type))
    expect_equal(cma@Metric, "FOS")
    expect_true(is.na(cma@PrimaryScore))
    expect_true(is.na(cma@Pvalue))
})

test_that("CMA validity rejects invalid PrimaryScore", {
    expect_error(
        ConsortiumMetabolismAlignment(PrimaryScore = 1.5),
        "PrimaryScore"
    )
    expect_error(
        ConsortiumMetabolismAlignment(PrimaryScore = -0.1),
        "PrimaryScore"
    )
})

test_that("CMA validity rejects invalid Type", {
    expect_error(
        ConsortiumMetabolismAlignment(Type = "invalid"),
        "Type"
    )
})

test_that("CMA validity rejects pairwise without names", {
    expect_error(
        ConsortiumMetabolismAlignment(
            Type = "pairwise",
            ReferenceName = "ref"
        ),
        "QueryName"
    )
})

test_that("CMA accepts valid pairwise configuration", {
    cma <- ConsortiumMetabolismAlignment(
        Type = "pairwise",
        QueryName = "query",
        ReferenceName = "ref",
        PrimaryScore = 0.75
    )
    expect_equal(cma@Type, "pairwise")
    expect_equal(cma@PrimaryScore, 0.75)
})

## ---- CMA accessors: pairwise ---------------------------------------------

test_that("scores() returns Scores list for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    s <- scores(cma)
    expect_true(is.list(s))
    expect_true("FOS" %in% names(s))
    expect_true("jaccard" %in% names(s))
})

test_that("pathways(type='shared') returns data.frame for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    sp <- pathways(cma, type = "shared")
    expect_true(is.data.frame(sp))
    expect_named(sp, c("consumed", "produced"))
})

test_that("pathways(type='unique') returns list for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    up <- pathways(cma, type = "unique")
    expect_true(is.list(up))
    expect_true("query" %in% names(up))
    expect_true("reference" %in% names(up))
    expect_true(is.data.frame(up$query))
})

test_that("metabolites() returns character vector for CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    m <- metabolites(cma)
    expect_true(is.character(m))
    expect_true(length(m) > 0L)
})

test_that("pathways() returns data.frame for CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    e <- pathways(cma)
    expect_true(is.data.frame(e))
    expect_named(e, c("consumed", "produced"))
    ev <- pathways(cma, verbose = TRUE)
    expect_true(is.data.frame(ev))
    expect_true(ncol(ev) >= ncol(e))
})

test_that("consortia() errors for CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(consortia(cma), "not applicable")
})

## ---- CMA accessors: type guards ------------------------------------------

test_that("pathways(type='shared') errors for multiple CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    expect_error(
        pathways(cma, type = "shared"), "pairwise"
    )
})

test_that("pathways(type='unique') errors for multiple CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    expect_error(
        pathways(cma, type = "unique"), "pairwise"
    )
})

test_that("similarityMatrix() errors for pairwise CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(similarityMatrix(cma), "multiple")
})

test_that("prevalence() errors for pairwise CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(prevalence(cma), "multiple")
})

## ---- CMA accessors: multiple ---------------------------------------------

test_that("scores() returns Scores list for multiple", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    s <- scores(cma)
    expect_true(is.list(s))
    expect_true("median" %in% names(s))
})

test_that("similarityMatrix() returns matrix for multiple", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    sm <- similarityMatrix(cma)
    expect_true(is.matrix(sm))
    expect_equal(nrow(sm), 2L)
})

test_that("prevalence() returns data.frame for multiple", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    p <- prevalence(cma)
    expect_true(is.data.frame(p))
    expect_true("nConsortia" %in% names(p))
})

test_that("pathways(type='consensus') returns data.frame for multiple", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    cp <- pathways(cma, type = "consensus")
    expect_true(is.data.frame(cp))
    expect_true("nConsortia" %in% names(cp))
    expect_true("proportion" %in% names(cp))
})

test_that("pathways(type='consensus') errors for pairwise CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(
        pathways(cma, type = "consensus"), "multiple"
    )
})

test_that("pathways(type='shared', verbose=TRUE) returns full data", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 4)
    cma <- align(cm1, cm2)
    sp_concise <- pathways(cma, type = "shared")
    sp_verbose <- pathways(
        cma, type = "shared", verbose = TRUE
    )
    expect_true(is.data.frame(sp_verbose))
    expect_true(ncol(sp_verbose) >= ncol(sp_concise))
    expect_true("querySpecies" %in% names(sp_verbose))
})

## ---- B1: species(CMA) returns only species, not flux values ----------------

test_that("species(CMA) returns only species names, not flux values", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 10)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 20)
    cma <- align(cm1, cm2)
    sp <- species(cma)
    expect_type(sp, "character")
    all_sp <- union(species(cm1), species(cm2))
    expect_true(all(sp %in% all_sp))
    expect_false(any(grepl("^[0-9]", sp)))
})

test_that("species(CMA) returns sorted unique character vector", {
    cm1 <- synCM("a", n_species = 4, max_met = 6, seed = 42)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 43)
    cma <- align(cm1, cm2)
    sp <- species(cma)
    expect_type(sp, "character")
    expect_equal(sp, sort(unique(sp)))
})

## ---- B3: similarityMatrix/prevalence on empty CMA --------------------------

test_that("similarityMatrix errors cleanly on empty CMA", {
    cma <- ConsortiumMetabolismAlignment()
    expect_error(similarityMatrix(cma), "multiple")
})

test_that("prevalence errors cleanly on empty CMA", {
    cma <- ConsortiumMetabolismAlignment()
    expect_error(prevalence(cma), "multiple")
})

## ---- B4: scores() includes MAAS composite and p-value ----------------------

test_that("scores() includes MAAS when method is MAAS", {
    cm1 <- synCM("a", n_species = 5, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 5, max_met = 8, seed = 43)
    cma <- align(cm1, cm2, method = "MAAS")
    s <- scores(cma)
    expect_true("MAAS" %in% names(s))
    expect_equal(s$MAAS, cma@PrimaryScore)
})

test_that("scores() does not include MAAS for FOS method", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2, method = "FOS")
    s <- scores(cma)
    expect_false("MAAS" %in% names(s))
})

test_that("scores() includes pvalue when computed", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 42)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 43)
    cma <- align(
        cm1, cm2,
        computePvalue = TRUE,
        nPermutations = 49L
    )
    s <- scores(cma)
    expect_true("pvalue" %in% names(s))
    expect_true(s$pvalue >= 0 && s$pvalue <= 1)
})

test_that("scores() omits pvalue when not computed", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    s <- scores(cma)
    expect_false("pvalue" %in% names(s))
})

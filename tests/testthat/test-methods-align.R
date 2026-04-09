## ---- align() dispatch tests -------------------------------------------------

test_that("align() dispatches CM,CM and returns CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cm2 <- synCM("b", n_species = 3, max_met = 5)
    cma <- align(cm1, cm2)
    expect_s4_class(cma, "ConsortiumMetabolismAlignment")
    expect_equal(cma@Type, "pairwise")
})

test_that("align() dispatches CMS,missing and returns CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cm2 <- synCM("b", n_species = 3, max_met = 5)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    expect_s4_class(cma, "ConsortiumMetabolismAlignment")
    expect_equal(cma@Type, "multiple")
    expect_equal(cma@Metric, "FOS")
})

## ---- Integration tests: align(CM, CM) --------------------------------------

test_that("align identical CMs returns score 1", {
    cm <- synCM("test", n_species = 4, max_met = 6, seed = 42)
    cma <- align(cm, cm)
    expect_s4_class(cma, "ConsortiumMetabolismAlignment")
    expect_equal(cma@PrimaryScore, 1)
    expect_equal(cma@Type, "pairwise")
    expect_equal(nrow(cma@UniqueQuery), 0)
    expect_equal(nrow(cma@UniqueReference), 0)
})

test_that("align is symmetric for FOS", {
    cm1 <- synCM("a", n_species = 4, max_met = 6, seed = 10)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 20)
    cma_ab <- align(cm1, cm2, method = "FOS")
    cma_ba <- align(cm2, cm1, method = "FOS")
    expect_equal(cma_ab@PrimaryScore, cma_ba@PrimaryScore)
})

test_that("align is symmetric for jaccard", {
    cm1 <- synCM("a", n_species = 4, max_met = 6, seed = 10)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 20)
    cma_ab <- align(cm1, cm2, method = "jaccard")
    cma_ba <- align(cm2, cm1, method = "jaccard")
    expect_equal(cma_ab@PrimaryScore, cma_ba@PrimaryScore)
})

test_that("all scores in [0, 1]", {
    cm1 <- synCM("a", n_species = 5, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 5, max_met = 8, seed = 43)
    cma <- align(cm1, cm2)
    for (nm in names(cma@Scores)) {
        if (!is.na(cma@Scores[[nm]])) {
            expect_true(
                cma@Scores[[nm]] >= 0,
                info = paste(nm, "< 0")
            )
            expect_true(
                cma@Scores[[nm]] <= 1,
                info = paste(nm, "> 1")
            )
        }
    }
})

test_that("MAAS method works", {
    cm1 <- synCM("a", n_species = 5, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 5, max_met = 8, seed = 43)
    cma <- align(cm1, cm2, method = "MAAS")
    expect_equal(cma@Metric, "MAAS")
    expect_true(!is.na(cma@PrimaryScore))
})

test_that("MAAS p-value uses composite score, not FOS", {
    cm1 <- synCM("a", n_species = 5, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 5, max_met = 8, seed = 43)
    cma <- align(
        cm1,
        cm2,
        method = "MAAS",
        computePvalue = TRUE,
        nPermutations = 49L
    )
    expect_true(!is.na(cma@Pvalue))
    expect_true(cma@Pvalue >= 0 && cma@Pvalue <= 1)
    ## p-value denominator should be nPerm + 1
    expect_true(
        cma@Pvalue %% (1 / (49L + 1L)) < .Machine$double.eps^0.5
    )
})

test_that("align with computePvalue returns valid p-value", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 42)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 43)
    cma <- align(cm1, cm2, computePvalue = TRUE, nPermutations = 99L)
    expect_true(!is.na(cma@Pvalue))
    expect_true(cma@Pvalue >= 0 && cma@Pvalue <= 1)
})

test_that("SharedPathways + UniqueQuery + UniqueRef = total pathways", {
    cm1 <- synCM("a", n_species = 4, max_met = 6, seed = 10)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 20)
    cma <- align(cm1, cm2)
    total_shared <- nrow(cma@SharedPathways)
    total_uq <- nrow(cma@UniqueQuery)
    total_ur <- nrow(cma@UniqueReference)
    expect_equal(
        total_shared + total_uq,
        sum(assays(cm1)$Binary)
    )
    expect_equal(
        total_shared + total_ur,
        sum(assays(cm2)$Binary)
    )
})

test_that("FOS from align matches CMS overlap matrix", {
    cm1 <- synCM("a", n_species = 4, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 4, max_met = 8, seed = 43)
    cma <- align(cm1, cm2, method = "FOS")
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cms_fos <- 1 - cms@OverlapMatrix[1, 2]
    expect_equal(cma@PrimaryScore, cms_fos)
})

## ---- Edge case tests -------------------------------------------------------

test_that("align errors on invalid method", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    expect_error(align(cm1, cm2, method = "invalid"))
})

test_that("align handles single-edge CMs", {
    cm1 <- ConsortiumMetabolism(
        data = tibble::tibble(
            species = c("sp1", "sp1"),
            metabolite = c("A", "B"),
            flux = c(-1, 1)
        ),
        name = "single1"
    )
    cm2 <- ConsortiumMetabolism(
        data = tibble::tibble(
            species = c("sp2", "sp2"),
            metabolite = c("A", "B"),
            flux = c(-1, 1)
        ),
        name = "single2"
    )
    cma <- align(cm1, cm2)
    expect_equal(cma@PrimaryScore, 1)
})

## ---- Multiple alignment: align(CMS, missing) ----------------------------

test_that("FOS similarity matches CMS OverlapMatrix", {
    cm1 <- synCM("a", n_species = 4, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 4, max_met = 8, seed = 43)
    cm3 <- synCM("c", n_species = 4, max_met = 8, seed = 44)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    cma <- align(cms, method = "FOS")
    om <- cms@OverlapMatrix
    sim <- cma@SimilarityMatrix
    ## Each off-diagonal entry should be 1 - distance
    for (i in seq_len(nrow(sim))) {
        for (j in seq_len(ncol(sim))) {
            if (i != j) {
                expect_equal(
                    sim[i, j],
                    1 - om[i, j],
                    info = paste("pair", i, j)
                )
            }
        }
    }
})

test_that("SimilarityMatrix is square, symmetric, diagonal=1", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    cma <- align(cms)
    sim <- cma@SimilarityMatrix
    expect_equal(nrow(sim), ncol(sim))
    expect_equal(nrow(sim), 3L)
    expect_equal(sim, t(sim))
    expect_equal(unname(diag(sim)), rep(1, 3))
    off_diag <- sim[upper.tri(sim)]
    expect_true(all(off_diag >= 0 & off_diag <= 1))
})

test_that("SimilarityMatrix has correct dimnames", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    cm_names <- names(cms@BinaryMatrices)
    expect_equal(rownames(cma@SimilarityMatrix), cm_names)
    expect_equal(colnames(cma@SimilarityMatrix), cm_names)
})

test_that("PrimaryScore is median of pairwise scores", {
    cm1 <- synCM("a", n_species = 4, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 4, max_met = 8, seed = 43)
    cm3 <- synCM("c", n_species = 4, max_met = 8, seed = 44)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    cma <- align(cms)
    vals <- cma@SimilarityMatrix[upper.tri(
        cma@SimilarityMatrix
    )]
    expect_equal(cma@PrimaryScore, stats::median(vals))
})

test_that("Scores list is complete", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    cma <- align(cms)
    expected_keys <- c(
        "mean",
        "median",
        "min",
        "max",
        "sd",
        "nPairs"
    )
    expect_true(all(expected_keys %in% names(cma@Scores)))
    expect_equal(cma@Scores$nPairs, 3L)
    expect_true(all(vapply(
        cma@Scores,
        is.finite,
        logical(1L)
    )))
})

test_that("Prevalence data.frame has correct structure", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    prev <- cma@Prevalence
    expect_true(is.data.frame(prev))
    expect_true(all(
        c("consumed", "produced", "nConsortia", "proportion") %in% names(prev)
    ))
    expect_true(all(prev$nConsortia >= 1L))
    expect_true(all(prev$nConsortia <= 2L))
    expect_true(all(
        prev$proportion > 0 &
            prev$proportion <= 1
    ))
})

test_that("ConsensusPathways matches Prevalence", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    expect_identical(cma@ConsensusPathways, cma@Prevalence)
})

test_that("Dendrogram is present and valid", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    cma <- align(cms)
    expect_equal(length(cma@Dendrogram), 1L)
    expect_true(is(cma@Dendrogram[[1L]], "dendrogram"))
})

test_that("Jaccard method produces valid result", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    cma <- align(cms, method = "jaccard")
    expect_equal(cma@Metric, "jaccard")
    sim <- cma@SimilarityMatrix
    expect_equal(sim, t(sim))
    expect_equal(unname(diag(sim)), rep(1, 3))
})

test_that("2-consortium edge case works", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    expect_equal(nrow(cma@SimilarityMatrix), 2L)
    expect_equal(cma@Scores$nPairs, 1L)
    ## PrimaryScore equals the single pairwise score
    expect_equal(
        cma@PrimaryScore,
        cma@SimilarityMatrix[1, 2]
    )
})

test_that("multiple alignment errors on invalid method", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    expect_error(align(cms, method = "invalid"))
})

test_that("BPPARAM argument is accepted", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    expect_no_error(
        align(cms, BPPARAM = BiocParallel::SerialParam())
    )
})

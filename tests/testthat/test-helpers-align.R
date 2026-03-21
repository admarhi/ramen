## ---- Unit tests: .expandMatrix() -------------------------------------------

test_that("expandMatrix maps to correct positions", {
    mat <- Matrix::sparseMatrix(
        i = 1, j = 2, x = 1,
        dims = c(2, 2),
        dimnames = list(c("A", "B"), c("A", "B"))
    )
    result <- .expandMatrix(mat, c("A", "B", "C"))
    expect_equal(dim(result), c(3, 3))
    expect_equal(result[1, 2], 1)
    expect_equal(sum(result), 1)
})

test_that("expandMatrix handles empty matrix", {
    mat <- Matrix::sparseMatrix(
        i = integer(0), j = integer(0),
        x = numeric(0), dims = c(2, 2),
        dimnames = list(c("A", "B"), c("A", "B"))
    )
    result <- .expandMatrix(mat, c("A", "B", "C"))
    expect_equal(dim(result), c(3, 3))
    expect_equal(sum(result), 0)
})

## ---- Unit tests: .harmonizeMetaboliteSpace() --------------------------------

test_that("harmonize creates correct union space", {
    cm1 <- synCM("a", n_species = 3, max_met = 4, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 4, seed = 2)
    h <- .harmonizeMetaboliteSpace(cm1, cm2)
    mets1 <- rownames(assays(cm1)$Binary)
    mets2 <- rownames(assays(cm2)$Binary)
    expect_true(all(mets1 %in% h$metabolites))
    expect_true(all(mets2 %in% h$metabolites))
    expect_equal(nrow(h$X), length(h$metabolites))
    expect_equal(nrow(h$Y), length(h$metabolites))
})

## ---- Unit tests: .functionalOverlap() --------------------------------------

test_that("FOS of identical matrices is 1", {
    m <- Matrix::sparseMatrix(
        i = c(1, 1, 2), j = c(2, 3, 3),
        x = 1, dims = c(3, 3)
    )
    expect_equal(.functionalOverlap(m, m), 1)
})

test_that("FOS of disjoint matrices is 0", {
    m1 <- Matrix::sparseMatrix(
        i = 1, j = 2, x = 1, dims = c(3, 3)
    )
    m2 <- Matrix::sparseMatrix(
        i = 2, j = 3, x = 1, dims = c(3, 3)
    )
    expect_equal(.functionalOverlap(m1, m2), 0)
})

test_that("FOS matches known value", {
    m1 <- Matrix::sparseMatrix(
        i = c(1, 1), j = c(2, 3), x = 1,
        dims = c(3, 3)
    )
    m2 <- Matrix::sparseMatrix(
        i = c(1, 2), j = c(2, 3), x = 1,
        dims = c(3, 3)
    )
    expect_equal(.functionalOverlap(m1, m2), 0.5)
})

test_that("FOS is symmetric", {
    m1 <- Matrix::sparseMatrix(
        i = c(1, 1), j = c(2, 3), x = 1,
        dims = c(3, 3)
    )
    m2 <- Matrix::sparseMatrix(
        i = c(1, 2), j = c(2, 3), x = 1,
        dims = c(3, 3)
    )
    expect_equal(
        .functionalOverlap(m1, m2),
        .functionalOverlap(m2, m1)
    )
})

test_that("FOS of empty matrices is 0", {
    m <- Matrix::sparseMatrix(
        i = integer(0), j = integer(0),
        x = numeric(0), dims = c(3, 3)
    )
    expect_equal(.functionalOverlap(m, m), 0)
})

## ---- Unit tests: .jaccardIndex() -------------------------------------------

test_that("Jaccard of identical is 1", {
    m <- Matrix::sparseMatrix(
        i = c(1, 2), j = c(2, 3), x = 1,
        dims = c(3, 3)
    )
    expect_equal(.jaccardIndex(m, m), 1)
})

test_that("Jaccard of disjoint is 0", {
    m1 <- Matrix::sparseMatrix(
        i = 1, j = 2, x = 1, dims = c(3, 3)
    )
    m2 <- Matrix::sparseMatrix(
        i = 2, j = 3, x = 1, dims = c(3, 3)
    )
    expect_equal(.jaccardIndex(m1, m2), 0)
})

test_that("Jaccard known value", {
    m1 <- Matrix::sparseMatrix(
        i = c(1, 1), j = c(2, 3), x = 1,
        dims = c(3, 3)
    )
    m2 <- Matrix::sparseMatrix(
        i = c(1, 2), j = c(2, 3), x = 1,
        dims = c(3, 3)
    )
    expect_equal(.jaccardIndex(m1, m2), 1 / 3)
})

## ---- Unit tests: .redundancyOverlap() --------------------------------------

test_that("redundancyOverlap of identical is 1", {
    m <- Matrix::sparseMatrix(
        i = c(1, 2), j = c(2, 3), x = c(3, 5),
        dims = c(3, 3)
    )
    expect_equal(.redundancyOverlap(m, m), 1)
})

test_that("redundancyOverlap of disjoint is 0", {
    m1 <- Matrix::sparseMatrix(
        i = 1, j = 2, x = 3, dims = c(3, 3)
    )
    m2 <- Matrix::sparseMatrix(
        i = 2, j = 3, x = 5, dims = c(3, 3)
    )
    expect_equal(.redundancyOverlap(m1, m2), 0)
})

test_that("redundancyOverlap known value", {
    m1 <- Matrix::sparseMatrix(
        i = c(1, 2), j = c(2, 3), x = c(3, 2),
        dims = c(3, 3)
    )
    m2 <- Matrix::sparseMatrix(
        i = c(1, 2), j = c(2, 3), x = c(5, 1),
        dims = c(3, 3)
    )
    expect_equal(.redundancyOverlap(m1, m2), 4 / 7)
})

## ---- Unit tests: .brayCurtisSimilarity() -----------------------------------

test_that("brayCurtis of identical is 1", {
    xW <- list(
        Consumption = Matrix::sparseMatrix(
            i = 1, j = 2, x = 3, dims = c(3, 3)
        ),
        Production = Matrix::sparseMatrix(
            i = 2, j = 3, x = 5, dims = c(3, 3)
        )
    )
    expect_equal(.brayCurtisSimilarity(xW, xW), 1)
})

test_that("brayCurtis of disjoint is 0", {
    xW <- list(
        Consumption = Matrix::sparseMatrix(
            i = 1, j = 2, x = 3, dims = c(3, 3)
        ),
        Production = Matrix::sparseMatrix(
            i = 1, j = 2, x = 3, dims = c(3, 3)
        )
    )
    yW <- list(
        Consumption = Matrix::sparseMatrix(
            i = 2, j = 3, x = 5, dims = c(3, 3)
        ),
        Production = Matrix::sparseMatrix(
            i = 2, j = 3, x = 5, dims = c(3, 3)
        )
    )
    result <- .brayCurtisSimilarity(xW, yW)
    expect_true(result >= 0 && result <= 1)
    expect_equal(result, 0)
})

## ---- Unit tests: .computeMAAS() --------------------------------------------

test_that("MAAS with default weights", {
    scores <- list(
        FOS = 0.8, jaccard = 0.6,
        brayCurtis = 0.7,
        redundancyOverlap = 0.5
    )
    expected <- 0.4 * 0.8 + 0.2 * 0.6 + 0.2 * 0.7 + 0.2 * 0.5
    expect_equal(.computeMAAS(scores), expected)
})

test_that("MAAS handles NA scores", {
    scores <- list(
        FOS = 0.8, jaccard = 0.6,
        brayCurtis = NA_real_,
        redundancyOverlap = NA_real_
    )
    expected <- (0.4 / 0.6) * 0.8 + (0.2 / 0.6) * 0.6
    expect_equal(.computeMAAS(scores), expected,
                 tolerance = 1e-10)
})

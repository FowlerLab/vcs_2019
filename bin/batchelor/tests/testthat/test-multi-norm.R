# This tests the multi-block normalization calculations.
# require(batchelor); require(testthat); source("test-multi-norm.R")

set.seed(20010)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- SingleCellExperiment(list(counts=dummy))
sizeFactors(X) <- runif(ncol(X))

set.seed(20011)
test_that("multiBatchNorm works properly", {
    X2 <- X
    counts(X2) <- counts(X2) * 2L
    X3 <- X
    counts(X3) <- counts(X3) * 3L

    # Checking that it works in the vanilla case.
    out <- multiBatchNorm(X, X2, X3)
    expect_equal(sizeFactors(out[[1]]), sizeFactors(out[[2]])/2)
    expect_equal(sizeFactors(out[[1]]), sizeFactors(out[[3]])/3)
    expect_equal(logcounts(out[[1]]), logcounts(out[[2]]))
    expect_equal(logcounts(out[[1]]), logcounts(out[[3]]))

    # Checking that the order does not matter.
    re.out <- multiBatchNorm(X3, X, X2)
    expect_equal(out[c(3,1,2)], re.out)

    # Single batch input just returns the same object as normalize().
    solo.out <- multiBatchNorm(X3)
    expect_equal(solo.out[[1]], normalize(X3))

    # Reverts to the library size correctly.
    Xtmp <- X
    sizeFactors(Xtmp) <- NULL
    expect_warning(out <- multiBatchNorm(Xtmp, X2), "no endogenous size factors in batch 1") 
    expect_equal(sizeFactors(out[[1]]), scater::librarySizeFactors(Xtmp))
})

set.seed(200111)
test_that("multiBatchNorm size factor centering logic is correct", {
    dummy2 <- matrix(rnbinom(ngenes*ncells, mu=means * 10, size=5), ncol=ncells, nrow=ngenes)
    rownames(dummy2) <- paste0("X", seq_len(ngenes))
    X2 <- SingleCellExperiment(list(counts=dummy2))
    sizeFactors(X2) <- runif(ncol(X2))

    # Centers of new size factors are correct.
    out <- multiBatchNorm(X, X2, min.mean=0)
    expect_equal(mean(sizeFactors(out[[1]])), 1)
    expect_equal(sizeFactors(scater::centreSizeFactors(X)), sizeFactors(out[[1]]))

    ave1 <- scater::calcAverage(X)
    ave2 <- scater::calcAverage(X2)
    M <- median(ave2/ave1)
    expect_equal(mean(sizeFactors(out[[2]])), M)
    expect_equal(sizeFactors(scater::centreSizeFactors(X2)), sizeFactors(out[[2]]) / M)

    # Renormalized values have no composition biases.
    ave1 <- scater::calcAverage(out[[1]])
    ave2 <- scater::calcAverage(out[[2]]) / M # as calcAverage automatically centers out[[2]'s SFs.
    expect_equal(1, median(ave2/ave1))
})

set.seed(20012)
test_that("multiBatchNorm behaves correctly with gene filtering", {
    dummy2 <- matrix(rnbinom(ngenes*ncells, mu=means * 10, size=5), ncol=ncells, nrow=ngenes)
    rownames(dummy2) <- paste0("X", seq_len(ngenes))
    
    X2 <- SingleCellExperiment(list(counts=dummy2))
    sizeFactors(X2) <- runif(ncol(X2))

    # Creating a reference function using calcAverage() explicitly.
    library(scater)
    REFFUN <- function(..., min.mean=1) {
	    batches <- list(...)
	    nbatches <- length(batches)
        batches <- lapply(batches, centreSizeFactors)
	    collected.ave <- lapply(batches, calcAverage, use_size_factors=TRUE)

    	collected.ratios <- rep(1, nbatches)
	    first.ave <- collected.ave[[1]]
        for (second in 2:nbatches) {
	        second.ave <- collected.ave[[second]]
            keep <- calcAverage(cbind(first.ave, second.ave)) >= min.mean
            collected.ratios[second] <- median(second.ave[keep]/first.ave[keep]) 
		}
		collected.ratios
    }

    out <- multiBatchNorm(X, X2, min.mean=1)
	refA <- REFFUN(X, X2, min.mean=1)
    expect_equal(mean(sizeFactors(out[[2]]))/mean(sizeFactors(out[[1]])), refA[2])

    out <- multiBatchNorm(X, X2, min.mean=10)
	refB <- REFFUN(X, X2, min.mean=10)
    expect_equal(mean(sizeFactors(out[[2]]))/mean(sizeFactors(out[[1]])), refB[2])

    out <- multiBatchNorm(X, X2, min.mean=100)
	refC <- REFFUN(X, X2, min.mean=100)
    expect_equal(mean(sizeFactors(out[[2]]))/mean(sizeFactors(out[[1]])), refC[2])
    
    expect_false(isTRUE(all.equal(refA, refB)))
    expect_false(isTRUE(all.equal(refA, refC)))

    # Checking that gene subsetting works correctly.
    randoms <- sample(ngenes, 500)
    sub.out <- multiBatchNorm(X, X2, subset.row=randoms)
    ref.out <- multiBatchNorm(X[randoms,], X2[randoms,])
    expect_equal(logcounts(sub.out[[1]])[randoms,], logcounts(ref.out[[1]])) 
    expect_equal(logcounts(sub.out[[2]])[randoms,], logcounts(ref.out[[2]])) 
})

set.seed(200121)
test_that("multiBatchNorm rescales spike-ins correctly", {
    isp <- rbinom(ngenes, 1, 0.1)==1
    isSpike(X, "ERCC") <- isp
    sizeFactors(X, "ERCC") <- runif(ncol(X))

    X2 <- X
    counts(X2) <- counts(X2) * 2L
    X3 <- X
    counts(X3) <- counts(X3) * 3L

    out <- multiBatchNorm(X, X2, X3)

    M <- mean(sizeFactors(X))
    expect_equal(sizeFactors(out[[1]]), sizeFactors(X)/M)
    expect_equal(sizeFactors(out[[1]]), sizeFactors(out[[2]]) / 2)
    expect_equal(sizeFactors(out[[1]]), sizeFactors(out[[3]]) / 3)

    expect_equal(sizeFactors(out[[1]], "ERCC"), sizeFactors(X, "ERCC")/M)
    expect_equal(sizeFactors(out[[1]], "ERCC"), sizeFactors(out[[2]], "ERCC")/2)
    expect_equal(sizeFactors(out[[1]], "ERCC"), sizeFactors(out[[3]], "ERCC")/3)
    
    # Missing spike-in size factors are handled correctly.
    Xtmp <- X
    sizeFactors(Xtmp, "ERCC") <- NULL
    expect_warning(mout <- multiBatchNorm(Xtmp, X2), "no ERCC size factors")
    expect_equal(sizeFactors(mout[[1]], "ERCC"), scater::librarySizeFactors(X, subset_row=isp)/M)
    expect_equal(mout[[2]], out[[2]])

    # Does nothing by default, when there are no endogenous genes.
    nout <- multiBatchNorm(X[isp,], X2[isp,], X3[isp,])

    expect_equal(sizeFactors(nout[[1]]), sizeFactors(X)/M)
    expect_equal(sizeFactors(nout[[1]]), sizeFactors(nout[[2]]))
    expect_equal(sizeFactors(nout[[1]]), sizeFactors(nout[[3]]))

    expect_equal(sizeFactors(nout[[1]], "ERCC"), sizeFactors(X, "ERCC")/M)
    expect_equal(sizeFactors(nout[[1]], "ERCC"), sizeFactors(nout[[2]], "ERCC"))
    expect_equal(sizeFactors(nout[[1]], "ERCC"), sizeFactors(nout[[3]], "ERCC"))
})

set.seed(200121)
test_that("multiBatchNorm handles separate rescaling of spike-ins correctly", {
    dummy2 <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
    rownames(dummy2) <- rownames(dummy)
    X2 <- SingleCellExperiment(list(counts=dummy2))
    sizeFactors(X2) <- runif(ncol(X2))

    Y <- X
    isp <- rbinom(ngenes, 1, 0.1)==1
    isSpike(Y, "ERCC") <- isp
    sizeFactors(Y, "ERCC") <- runif(ncol(X))

    Y2 <- X2
    isSpike(Y2, "ERCC") <- isp
    sizeFactors(Y2, "ERCC") <- runif(ncol(X2))

    ref <- multiBatchNorm(Y, Y2, separate.spikes=TRUE)

    # Non-spike-in counts are computed correctly.
    out <- multiBatchNorm(X[!isp,], X2[!isp,])
    expect_equal(logcounts(ref[[1]])[!isp,], logcounts(out[[1]]))
    expect_equal(logcounts(ref[[2]])[!isp,], logcounts(out[[2]]))

    # Spike-in counts are computed correctly.
    alt <- X[isp,]
    sizeFactors(alt) <- sizeFactors(Y, "ERCC")
    alt2 <- X2[isp,]
    sizeFactors(alt2) <- sizeFactors(Y2, "ERCC")

    sout <- multiBatchNorm(alt, alt2)
    expect_equal(sizeFactors(ref[[1]], "ERCC"), sizeFactors(sout[[1]])) 
    expect_equal(sizeFactors(ref[[2]], "ERCC"), sizeFactors(sout[[2]])) 
    expect_equal(logcounts(ref[[1]])[isp,], logcounts(sout[[1]])) 
    expect_equal(logcounts(ref[[2]])[isp,], logcounts(sout[[2]])) 

    # Missing spike-in size factors are handled correctly.
    Ytmp <- Y
    sizeFactors(Ytmp, "ERCC") <- NULL
    expect_warning(nout <- multiBatchNorm(Ytmp, Y2, separate.spikes=TRUE), "no ERCC size factors")
    expect_equal(sizeFactors(nout[[1]]), sizeFactors(out[[1]]))
    expect_equal(sizeFactors(nout[[2]]), sizeFactors(out[[2]]))

    alttmp <- X[isp,]
    sizeFactors(alttmp) <- NULL
    expect_warning(comp <- multiBatchNorm(alttmp, alt2), "no endogenous size factors")
    expect_equal(sizeFactors(comp[[1]]), sizeFactors(nout[[1]], "ERCC"))
    expect_equal(sizeFactors(comp[[2]]), sizeFactors(nout[[2]], "ERCC"))

    # Spike-in selection and subset.row interact correctly.
    randoms <- sample(ngenes, 500)
    ref <- multiBatchNorm(Y[randoms,], Y2[randoms,])
    out <- multiBatchNorm(Y, Y2, subset.row=randoms)
    expect_equal(logcounts(ref[[1]]), logcounts(out[[1]])[randoms,])
    expect_equal(logcounts(ref[[2]]), logcounts(out[[2]])[randoms,])
})

set.seed(20013)
test_that("multiBatchNorm spits the dummy correctly", {
    # Mismatching row names or dimensions.
    X2 <- X
    rownames(X2) <- sample(rownames(X2))
    expect_error(out <- multiBatchNorm(X, X2), "not the same")
    rownames(X2) <- NULL
    expect_error(out <- multiBatchNorm(X, X2), "not the same")

    X2 <- X
    expect_error(out <- multiBatchNorm(X, X2[1,]), "not the same")

    # Mismatching spike-in identities.
    Xsp <- X
    isSpike(Xsp, "ERCC") <- sample(ngenes, 200)
    expect_error(multiBatchNorm(Xsp, X2), "spike-in sets differ across batches")
    X2sp <- X
    isSpike(X2sp, "ERCC") <- sample(ngenes, 200)
    expect_error(multiBatchNorm(Xsp, X2sp), "spike-in identities differ across batches")

    Xtmp <- Xsp
    isSpike(Xtmp, "SIRV") <- sample(ngenes, 50)
    X2tmp <- X
    isSpike(X2tmp, "ERCC") <- isSpike(Xtmp, "SIRV") 
    isSpike(X2tmp, "SIRV") <- isSpike(Xtmp, "ERCC") 
    expect_error(multiBatchNorm(Xtmp, X2tmp), "ERCC spike-in identities differ across batches")

    # Input resulting in NA scaling values (add 1 to avoid NA size factors).
    X2 <- X
    counts(X2)[] <- 0
    counts(X2)[1,] <- 1
    expect_error(out <- multiBatchNorm(X, X2), "not finite")

    # Empty inputs.
    expect_error(out <- multiBatchNorm(X[0,], X2[0,]), "not finite")
    expect_error(out <- multiBatchNorm(X[,0], X2[,0]), "not finite")
    expect_error(multiBatchNorm(), "at least one SingleCellExperiment")
})


context("strain")

test_that("random from plaf", {
  
  # with 1 and 0 should be no variation
    plaf <- c(rep(1,18),rep(0,18))
    y <- test_barcode_from_PLAF(plaf, length(plaf))
    x <- as.raw(plaf)
    expect_identical(y, x)
  
  # from a differing plaf
    plaf <- dnorm(1:47,mean = 12,sd = 5)
    y <- replicate(100000,as.numeric(test_barcode_from_PLAF(plaf, length(plaf))))
    x <- rowMeans(y)
    expect_true(sum(abs(plaf - x)<mean(plaf)) == length(plaf))
  
})

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



test_that("ibd recombinants", {
  
  # with 1 and 0 should be no variation
  x <- as.raw(c(1,1,1,0,0,1,0,0))
  y <- as.raw(c(0,0,0,1,1,0,1,1))
  z <- test_recombinant_with_ibd(x,y,length(x),4,2,rep(1,4))
  expect_identical(y,z)

  })


test_that("ibd next generation", {
  
  z <- test_generate_next_ibd(bl = 96,nl = 24,ib = 2,pc = rep(1,24),id = 2)
  expect_identical(y,z)
  
})

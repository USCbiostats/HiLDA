test_that("correctly load the test data", {
  inputFile <- system.file("extdata/sample.rdata", package = "HiLDA")
  load(inputFile)

  expect_equal(as.character(class(G)), "MutationFeatureData")
  expect_equal(G@type, "independent")
  expect_equal(G@flankingBasesNum, 5)
})


test_that("running the global test and the local test", {
  inputFile <- system.file("extdata/sample.rdata", package = "HiLDA")
  load(inputFile)
  K <- 3
  Param <- pmsignature::getPMSignature(G, K = K)

  expect_equal(class(pmBarplot(G, Param, refGroup=1:4)), "list")
})

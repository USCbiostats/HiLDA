test_that("correctly load the test data", {
  inputFile <- system.file("extdata/sample.rdata", package = "HiLDA")
  load(inputFile)

  expect_equal(as.character(class(G)), "MutationFeatureData")
  expect_equal(slot(G, "type"), "independent")
  expect_equal(slot(G, "flankingBasesNum"), 5)
})


test_that("running the global test and the local test", {

  load(inputFile <- system.file("extdata/sample.rdata", package = "HiLDA"))
  Param <- pmgetSignature(G, K = 3)
  expect_equal(class(pmBarplot(G, Param, refGroup=1:4)), "list")
})

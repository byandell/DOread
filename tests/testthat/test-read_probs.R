context("Read genotype probabilities")

test_that("feather_genoprob works for reading probabilities", {
  library(DOread)
  library(qtl2feather)
  chr_id <- "10"
  datapath <- "~/Documents/Research/attie_alan/DO/data/DerivedData"
  probs <- read_probs(chr_id, start_val = 10, end_val = 50, datapath = datapath,
                              method = "calc")
  fprobs <- read_probs(chr_id, start_val = 10, end_val = 50, datapath = datapath)
  expect_equal(probs[["10"]], fprobs[["10"]])
})

context("Read mRNA expression")

test_that("feather_genoprob stored mRNA expression", {
  library(DOread)
  library(qtl2feather)
  chr_id <- "10"
  datapath <- "~/Documents/Research/attie_alan/DO/data/DerivedData"
  fprobs <- read_probs(chr_id, start_val = 10, end_val = 50, datapath = datapath)$probs
  indID <- dimnames(fprobs)$ind
  mrna <- read_mrna(indID, chr_id, start_val = 10, end_val = 50, datapath = datapath)
  expect_equal(nrow(mrna$expr), length(indID))
  expect_equal(ncol(mrna$expr), nrow(mrna$annot))
})

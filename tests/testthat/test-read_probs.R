context("Read genotype probabilities")

test_that("feather_genoprob works for reading probabilities", {
  library(DOread)
  library(qtl2feather)
  chr_id <- "10"
  datapath <- "~/Documents/Research/attie_alan/DO/data/DerivedData"

  # Whole chromosome allele=TRUE
  probs <- read_probs(chr_id, datapath = datapath,
                      method = "calc")
  fprobs <- read_probs(chr_id, datapath = datapath)
  expect_equal(probs$probs[[chr_id]], fprobs$probs[[chr_id]])
  # Subset of chromosome allele=TRUE
  probs <- read_probs(chr_id, start_val = 10, end_val = 20, datapath = datapath,
                              method = "calc")
  fprobs <- read_probs(chr_id, start_val = 10, end_val = 20, datapath = datapath)
  expect_equal(probs$probs[[chr_id]], fprobs$probs[[chr_id]])

  # Whole chromosome allele=FALSE
  probs <- read_probs(chr_id, datapath = datapath,
                      method = "calc", allele=FALSE)
  fprobs <- read_probs(chr_id, datapath = datapath, allele=FALSE)
  expect_equal(probs$probs[[chr_id]], fprobs$probs[[chr_id]])
  # Subset of chromosome allele=TRUE
  probs <- read_probs(chr_id, start_val = 10, end_val = 20, datapath = datapath,
                      method = "calc", allele=FALSE)
  fprobs <- read_probs(chr_id, start_val = 10, end_val = 20,
                       datapath = datapath, allele=FALSE)
  # This is not working. fprobs is missing some markers -- why?
  #  Error in get_dimension(mar, x$mar, type = "marker") :
  # Some markers not found: JAX00282567, UNCJPD004131, UNCJPD008279, UNCJPD004132, JAX00282796, JAX00194630, UNCJPD004135, JAX00282915, UNCJPD008017, UNCJPD008049, JAX00282986, JAX00282986q, JAX00283017, UNCJPD004137, UNCJPD004138, UNCJPD004140, UNCJPD004141, UNCJPD004142, JAX00283386, ICR268, ICR4268, JAX00283473, UNC17474037, UNCJPD004145, ICR275, ICR276, JAX00283653, JAX00283653q, JAX00015318, JAX00015318r, UNCJPD004148, JAX00283721, UNCJPD004151, JAX00283798, JAX00284009, UNCHS027600, UNCJPD004155, JAX00284145, JAX00284145q, UNCHS027617, UNCJPD004158, UNCJPD004159
  expect_equal(probs$probs[[chr_id]], fprobs$probs[[chr_id]], allele=FALSE)
})

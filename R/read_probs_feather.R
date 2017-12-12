# Read genotype probability object from file
read_probs_feather <- function(chr, datapath, allele = TRUE) {

  ## Read in feather_genotype object (small).
  probs <- readRDS(file.path(datapath, "feather",
                             ifelse(allele, "faprobs.rds", "fprobs.rds")))

  ## Modify feather directory to match current datapath.
  pr <- unclass(probs)
  oldpath <- pr$feather
  pr$feather <-
    stringr::str_replace(oldpath,
                         file.path(".*feather", basename(oldpath)),
                         file.path(datapath, "feather", basename(oldpath)))
  probs <- modify_object(probs, pr)
  subset(probs, chr = chr)
}

# Read genotype probability object from file
read_probs_feather <- function(chr, datapath, allele = TRUE) {

  ## Read in feather_genotype object (small).
  probs <- readRDS(file.path(datapath, "feather",
                             ifelse(allele, "faprobs.rds", "fprobs.rds")))
  subset(probs, chr = chr)
}

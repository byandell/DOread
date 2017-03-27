# Read genotype probability object from file
read_probs_feather <- function(chr=NULL, start_val=NULL, end_val=NULL, datapath,
                               allele = TRUE) {

  map <- readRDS(file.path(datapath, "pmap.rds"))

  ## Read in feather_genotype object (small).
  probs <- readRDS(file.path(datapath, "feather",
                             ifelse(allele, "faprobs.rds", "fprobs.rds")))

  if(!is.null(chr)) {
    map <- map[chr]
    probs <- subset(probs, chr = chr)

    if(length(chr) == 1) {
      if(!is.null(start_val) & !is.null(end_val)) {
        ## Reduce to region of interest.
        wh <- which(map[[chr]] >= start_val &
                      map[[chr]] <= end_val)
        map[[chr]] <- map[[chr]][wh]
        probs <- subset(probs, mar = wh)
      }
    }
  }
  list(probs = convert_probs(probs),
       map = map)
}

# Read genotype probability object from file
read_probs_feather <- function(chr=NULL, start_val=NULL, end_val=NULL, datapath) {

  map <- readRDS(file.path(datapath, "pmap.rds"))

  ## Read in phenotype probabilities ("probs") [very large]
  probs <- readRDS(file.path(datapath, "feather", "faprobs.rds"))

  if(!is.null(chr)) {
    map <- map[chr]
    probs <- probs[, chr]

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
# Read genotype probability object from file
read_probs36_feather <- function(chr, start_val=NULL, end_val=NULL, datapath) {

  if(missing(chr))
    stop("must specifiy one chromosome")
  if(length(chr) != 1)
    stop("must specifiy one chromosome")

  map <- readRDS(file.path(datapath, "pmap.rds"))[chr]

  # Read in feather_genoprob object for this chromosome.
  probs1 <- readRDS(file.path(datapath, "feather", paste0("fprobs_", chr, ".rds")))

  if(!is.null(start_val) & !is.null(end_val)) {
    ## Reduce to region of interest.
    wh <- which(map[[chr]] >= start_val &
                  map[[chr]] <= end_val)
    map[[chr]] <- map[[chr]][wh]
    probs1 <- subset(probs1, mar = wh)
  }

  list(probs = probs1,
       map = map)
}

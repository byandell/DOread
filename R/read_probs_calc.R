# Read genotype probability object from file
read_probs_calc <- function(chr=NULL, start_val=NULL, end_val=NULL, datapath,
                       map = readRDS(file.path(datapath, "pmap.rds"))) {
  read_all_probs <- is.null(chr)

  if(!read_all_probs) {
    ## Read all probs if requesting at least half the chromosomes.
    read_all_probs <- (2 * length(chr) >= length(map))
  }

  if(read_all_probs) {
    ## Read in phenotype probabilities ("probs") [very large]
    probs <- readRDS(file.path(datapath, "probs.rds"))
  } else {
    ## Read in probs for selected chromosomes and cbind.
    probs <- readRDS(file.path(datapath, paste0("probs_", chr[1], ".rds")))
    if(length(chr) > 1) for(chri in chr[-1]) {
      probs <- cbind(probs,
                     readRDS(file.path(datapath, paste0("probs_", chri, ".rds"))))
    } else {
      if(!is.null(start_val) & !is.null(end_val)) {
        ## Reduce to region of interest.
        wh <- which(map[[chr]] >= start_val &
                      map[[chr]] <= end_val)
        map[[chr]] <- map[[chr]][wh]
        probs1[[chr]] <- probs1[[chr]][,,wh]
      }
    }
    map <- map[chr]
  }
  list(probs = convert_probs(probs),
       map = map)
}
# Read genotype probability object from file
#' @importFrom qtl2scan interp_map
read_probs36_calc <- function(chr=NULL, start_val=NULL, end_val=NULL, datapath, map=NULL) {

  if(is.null(chr))
    stop("must supply chr")

  attieDO <- readRDS(file.path(datapath,"attieDO.rds"))  # cross object
  probs1 <- readRDS(file.path(datapath,
                              paste0("attieDO_probs", chr, ".rds")))

  if(is.null(map)) {
    map <- probs1$map
    if(is.null(map))
      stop("need physical map")
    map <- qtl2scan::interp_map(map, attieDO$gmap, attieDO$pmap)
  }

  ## Needed for transition only. Does nothing if already in new format.
  pr <- probs1 <- convert_probs(probs1)

  if(!is.null(start_val) & !is.null(end_val)) {
    ## Reduce to region of interest.
    wh <- which(map[[chr]] >= start_val &
                  map[[chr]] <= end_val)
    map[[chr]] <- map[[chr]][wh]
    probs1[[chr]] <- probs1[[chr]][,,wh]
  }

  ## Fix rownames of probs. Begin with "DO-".
  tmp <- substring(rownames(probs1[[chr]]), 4)
  rownames(probs1[[chr]]) <- tmp
  ## Sort them in increasing number order.
  probs1[[chr]] <-
    probs1[[chr]][order(as.numeric(tmp)),,]

  # Bring along attributes for calc_genoprob object.
  probs1 <- modify_probs(pr, probs1)

  list(probs = probs1,
       map = map)
}

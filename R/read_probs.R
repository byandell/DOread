# Half-way redo of read_probs to align with new R/qtl2 data format
# calc_genoprob attributes: is_x_chr, alleles, alleleprobs, and crosstype

#' Read genotype probability object from file
#'
#' Uses readRDS to read object.
#'
#' @param chr vector of chromosome identifiers
#' @param datapath name of folder with Derived Data
#' @param map genome map as list of chromosome maps
#'
#' @return list with \code{probs} = large object of class \code{\link[qtl2geno]{calc_genoprob}} and \code{map} = physical map for selected \code{chr}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{read_probs(chr, datapath)}
#'
#' @export
read_probs <- function(chr=NULL, datapath,
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
    if(length(chr) > 1) for(chri in chr[-1])
      probs <- cbind(probs,
                     readRDS(file.path(datapath, paste0("probs_", chri, ".rds"))))
    map <- map[chr]
  }
  list(probs = convert_probs(probs),
       map = map)
}
#' Read genotype probability object from file
#'
#' Uses readRDS to read object.
#'
#' @param chr vector of chromosome identifiers
#' @param start_val, end_val start and end values in Mbp
#' @param datapath name of folder with Derived Data
#'
#' @return list with \code{probs} = large object of class \code{\link[qtl2geno]{calc_genoprob}}
#'  and \code{map} = physical map for selected \code{chr} region
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{read_probs36(chr, start_val, end_val, dirpath)}
#'
#' @export
#' @rdname read_probs
#' @importFrom qtl2scan interp_map
read_probs36 <- function(chr_id, start_val=NULL, end_val=NULL, datapath, map=NULL) {
  attieDO <- readRDS(file.path(datapath,"attieDO.rds"))  # cross object
  probs1 <- readRDS(file.path(datapath,
                              paste0("attieDO_probs", chr_id, ".rds")))

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
    wh <- which(map[[chr_id]] >= start_val &
                  map[[chr_id]] <= end_val)
    map[[chr_id]] <- map[[chr_id]][wh]
    probs1[[chr_id]] <- probs1[[chr_id]][,,wh]
  }

  ## Fix rownames of probs. Begin with "DO-".
  tmp <- substring(rownames(probs1[[chr_id]]), 4)
  rownames(probs1[[chr_id]]) <- tmp
  ## Sort them in increasing number order.
  probs1[[chr_id]] <-
    probs1[[chr_id]][order(as.numeric(tmp)),,]

  # Bring along attributes for calc_genoprob object.
  probs1 <- modify_probs(pr, probs1)

  list(probs = probs1,
       map = map)
}

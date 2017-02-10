#' Read genotype probability object from file
#'
#' Uses readRDS to read object.
#'
#' @param chr vector of chromosome identifiers
#' @param datapath name of folder with Derived Data
#'
#' @return large object of class \code{\link[qtl2geno]{calc_genoprob}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{read_probs(chr, datapath)}
#'
#' @export
read_probs <- function(chr=NULL, datapath) {
  read_all_probs <- is.null(chr)

  if(!read_all_probs) {
    ## Read all probs if requesting at least half the chromosomes.
    pmap <- readRDS(file.path(datapath, "pmap.rds"))
    read_all_probs <- (2 * length(chr) >= length(pmap))
    rm(pmap)
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
  }
  probs
}
#' Read genotype probability object from file
#'
#' Uses readRDS to read object.
#'
#' @param chr vector of chromosome identifiers
#' @param start_val, end_val start and end values in Mbp
#' @param datapath name of folder with Derived Data
#'
#' @return large object of class \code{\link[qtl2geno]{calc_genoprob}}
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
read_probs36 <- function(chr_id, start_val, end_val, datapath) {
  attieDO <- readRDS(file.path(datapath,"attieDO.rds"))  # cross object
  probs1 <- readRDS(file.path(datapath,
                              paste0("attieDO_probs", chr_id, ".rds")))
  probs1$map <- qtl2scan::interp_map(probs1$map, attieDO$gmap, attieDO$pmap)
  ## Reduce to region of interest.
  wh <- which(probs1$map[[chr_id]] >= start_val &
                probs1$map[[chr_id]] <= end_val)
  probs1$map[[chr_id]] <- probs1$map[[chr_id]][wh]
  probs1$probs[[chr_id]] <- probs1$probs[[chr_id]][,,wh]
  ## Fix rownames of probs. Begin with "DO-".
  tmp <- substring(rownames(probs1$probs[[chr_id]]), 4)
  rownames(probs1$probs[[chr_id]]) <- tmp
  ## Sort them in increasing number order.
  probs1$probs[[chr_id]] <-
    probs1$probs[[chr_id]][order(as.numeric(tmp)),,]

  probs1
}

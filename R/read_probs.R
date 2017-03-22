#' Read genotype probability object from file
#'
#' Uses readRDS to read object.
#'
#' @param chr vector of chromosome identifiers
#' @param start_val, end_val start and end values in Mbp
#' @param datapath name of folder with Derived Data
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
read_probs <- function(chr=NULL, start_val=NULL, end_val=NULL, datapath) {

  map <- readRDS(file.path(datapath, "pmap.rds"))

  ## Read in phenotype probabilities ("probs") [very large]
  probs <- readRDS(file.path(datapath, "faprobs.rds"))

  if(!is.null(chr)) {
    map <- map[chr]
    probs <- probs[, chr]

    if(length(chr) == 1) {
      if(!is.null(start_val) & !is.null(end_val)) {
        ## Reduce to region of interest.
        wh <- which(map[[chr]] >= start_val &
                      map[[chr]] <= end_val)
        map[[chr]] <- map[[chr]][wh]
        probs <- subset(probs1, mar = wh)
      }
    }
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
read_probs36 <- function(chr, start_val=NULL, end_val=NULL, datapath) {

  if(missing(chr))
    stop("must specifiy one chromosome")
  if(length(chr) != 1)
    stop("must specifiy one chromosome")

  map <- readRDS(file.path(datapath, "pmap.rds"))[chr]

  # Read in feather_genoprob object for this chromosome.
  probs1 <- readRDS(file.path(datapath, paste0("fprobs_", chr, ".rds")))

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

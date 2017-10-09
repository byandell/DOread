#' Get phenotypes
#'
#' Get phenotypes using data frame of phenopypes filtered by \code{analyses_tbl}
#'
#' @param phe phenotypes in data frame
#' @param analyses_tbl table of analyses setups
#' @param transform transform phenotypes according to \code{analyses_tbl} if \code{TRUE}
#'
#' @return data frame of phenotypes
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_pheno(phe, analyses_tbl)}
#'
#' @export
#'
#' @importFrom assertthat assert_that is.number
#' @importFrom stats quantile
#'
get_pheno <- function(phe, analyses_tbl, transform = TRUE) {
  # Get phenotype names. Make sure it is character, not factor.
  phename <- as.character(unique(analyses_tbl$pheno))
  #  indID <- unlist(phe[,1]) ## IDs assumed to be first trait.
  phe <- phe[, match(phename, colnames(phe), nomatch=0), drop=FALSE]
  #  rownames(phe) <- indID

  if(transform) {
    ## Transform phenotype.
    not.id <- !((transf <- analyses_tbl$transf) %in% c("id","identity"))
    if(any(not.id)) {
      offset <- analyses_tbl$offset
      for(i in seq_along(not.id)[not.id])
        phe[,phename[i]] <- get(transf[i])(unlist(phe[,phename[i]]) + offset[i])
    }
    ## Parameter to winsorize will later be a value.
    if(any(is.logical(winz <- analyses_tbl$winsorize)))
      winz <- rep(0.02, length(phename))
    winsor <- (winz > 0)
    if(any(winsor)) {
      for(i in seq_along(winsor)[winsor])
        phe[,phename[i]] <- winsorize(unlist(phe[,phename[i]]), winz[i])
    }
  }
  phe
}

#' Read phenotypes from RDS files
#'
#' Read phenotypes previously stored as RDS filtered by \code{analyses_tbl}
#'
#' @param analyses_tbl table of analyses setups
#' @param datapath path to Derived Data
#'
#' @return data frame of phenotypes
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{read_pheno_tbl(analyses_tbl, datapath)}
#'
#' @export
#' @importFrom dplyr bind_cols
read_pheno_tbl <- function(analyses_tbl, datapath) {
  pheno_group <- unique(analyses_tbl$pheno_group)
  phe <- readRDS(file.path(datapath,
                           paste0("pheno_", pheno_group[1], ".rds")))
  indID <- rownames(phe)
  myadd <- function(df, row_names) {
    out <- as.data.frame(matrix(NA, length(row_names), ncol(df)))
    dimnames(out) <- list(row_names, names(df))
    out
  }
  if(length(pheno_group) > 1) for(phe_gp in pheno_group[-1]) {
    tmp <- readRDS(file.path(datapath, paste0("pheno_", phe_gp, ".rds")))
    tmpID <- rownames(tmp)
    m <- match(tmpID,indID)
    na.m <- is.na(m)
    if(any(na.m)) {
      ## add rows to phe if it is too short
      phe <- rbind(phe, myadd(phe, tmpID[na.m]))
      indID <- rownames(phe)
      m <- match(tmpID,indID)
    }
    na.m <- length(indID) - length(tmpID)
    ## add rows to tmp if too short
    if(na.m)
      tmp <- rbind(tmp, myadd(tmp, indID[-m]))
    phe <- dplyr::bind_cols(phe, tmp)
  }
  phe <- as.data.frame(phe)
  rownames(phe) <- indID
  phe
}

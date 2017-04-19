#' Read mrna expression object from file
#'
#' Uses feather to read mrna expression object and associated annotations.
#'
#' @param indID individual IDs to match with phenotypes
#' @param chr_id vector of chromosome identifiers
#' @param start_val, end_val start and end values in Mbp
#' @param datapath name of folder with Derived Data
#' @param allele read haplotype allele probabilities (if \code{TRUE}) or diplotype allele-pair probabilities (if \code{FALSE})
#' @param method method of genoprob storage
#'
#' @return list with \code{expr} = matrix of expression mRNA values in region and \code{annot} = data frame of annotations for mRNA.
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{read_probs(chr, datapath)}
#'
#' @export
#' @importFrom dplyr filter mutate rename
#' @importFrom feather read_feather
#'
read_mrna <- function(indID, chr_id=NULL, start_val=NULL, end_val=NULL, datapath) {

  annot.mrna <-
    dplyr::rename(
      dplyr::mutate(
        dplyr::filter(
          readRDS(file.path(datapath, "RNAseq", "annot.mrna.rds")),
          chr == chr_id,
          start >= start_val * 1e6,
          end <= end_val * 1e6),
        start = start * 1e-6,
        end = end * 1e-6,
        middle_point = middle_point * 1e-6),
      pos = middle_point)

  expr.mrna <- feather::read_feather(file.path(datapath, "RNAseq", "expr.mrna.feather"),
                                     c("Mouse.ID", annot.mrna$id))
  expr_mx <- matrix(NA, length(indID), ncol(expr.mrna) - 1)
  m <- match(expr.mrna$Mouse.ID, indID, nomatch = 0)
  if(!any(m > 0))
    stop("no individuals selected")

  expr_mx[m,] <- as.matrix(expr.mrna[m > 0, -1])
  dimnames(expr_mx) <- list(indID, colnames(expr.mrna)[-1])

  list(expr = expr_mx, annot = annot.mrna)
}

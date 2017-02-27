#' Quantile rank
#'
#' Quantile rank taken from R/qtl.
#'
#' @param x numeric vector
#' @param jitter jitter if \code{TRUE} (default \code{FALSE})
#'
#' @return transformed values
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dontrun(nqrank(x))
#'
#' @export
nqrank <- function (x, jitter = FALSE)
{
  qtl::nqrank(x, jitter)
}

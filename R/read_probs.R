#' Read genotype probability object from file
#'
#' Read object from file stored according to method.
#'
#' @param chr vector of chromosome identifiers
#' @param start_val, end_val start and end values in Mbp
#' @param datapath name of folder with Derived Data
#' @param allele read haplotype allele probabilities (if \code{TRUE}) or diplotype allele-pair probabilities (if \code{FALSE})
#' @param method method of genoprob storage
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
read_probs <- function(chr=NULL, start_val=NULL, end_val=NULL, datapath,
                       allele = TRUE, method = c("feather","calc")) {

  method <- match.arg(method)

  switch(method,
         feather = read_probs_feather(chr, start_val, end_val, datapath, allele),
         calc = {
           if(allele)
             read_probs_calc(chr, start_val, end_val, datapath)
           else
             read_probs36_calc(chr, start_val, end_val, datapath)
         })
}

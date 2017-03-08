#' Convert probs object to new format.
#'
#' @param object object of class \code{\link[qtl2geno]{calc_genoprob}} computed with old format
#'
#' @export
#'
convert_probs <- function(object) {
  # calc_genoprob attributes: is_x_chr, alleles, alleleprobs, and crosstype

  # If already an array, assume conversion already done.
  if(is.array(object[[1]]))
    return(object)

  out <- object$probs

  for(obj in c("crosstype", "is_x_chr", "alleles", "alleleprobs")) {
    if(!is.null(object[[obj]]))
    attr(out, obj) <- object[[obj]]
  }
  class(out) <- class(object)
  out
}
#' Modify probs object in new format.
#'
#' @param object object of class \code{\link[qtl2geno]{calc_genoprob}}
#' @param newprob matrix of new LOD values
#'
#' @export
#'
modify_probs <- function(object, newprob) {
  # calc_genoprob attributes: is_x_chr, alleles, alleleprobs, and crosstype

  for(obj in c("crosstype", "is_x_chr", "alleles", "alleleprobs")) {
    x_attr <- attr(object, obj)
    if(!is.null(x_attr))
      attr(newprob, obj) <- x_attr
  }
  class(newprob) <- class(object)
  newprob
}

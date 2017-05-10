#' Read mrna expression object from file
#'
#' Uses feather to read mrna expression object and associated annotations.
#'
#' @param chr_id vector of chromosome identifiers
#' @param start_val,end_val start and end values in Mbp
#' @param datapath name of folder with Derived Data
#' @param local read only mRNA values local to region if \code{TRUE};
#' otherwise include distal mRNA values that map to region
#' @param qtl read only mRNA values with QTL peak in region if \code{TRUE}
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
#' @importFrom dplyr filter group_by inner_join mutate rename summarize ungroup
#' @importFrom feather read_feather
#'
read_mrna <- function(chr_id=NULL, start_val=NULL, end_val=NULL, datapath,
                      local = TRUE, qtl = FALSE) {

  start_val6 <- 1e6 * start_val
  end_val6 <- 1e6 * end_val

  # Identify mRNA located in region or with QTL peak in region.
  peaks.mrna <- feather::read_feather(file.path(datapath, "RNAseq", "peaks.mrna.feather"))
  if(local) {
    peaks.mrna <- dplyr::filter(peaks.mrna,
                                gene_chr == chr_id,
                                pmax(gene_start, gene_end) >= start_val6,
                                pmin(gene_start, gene_end) <= end_val6)
  } else {
    peaks.mrna <- dplyr::filter(peaks.mrna,
                                ((gene_chr == chr_id &
                                    pmax(gene_start, gene_end) >= start_val6 &
                                    pmin(gene_start, gene_end) <= end_val6) |
                                   (qtl_chr == chr_id &
                                      qtl_pos >= start_val6 &
                                      qtl_pos <= end_val6)))
  }
  peaks.mrna <- dplyr::mutate(peaks.mrna,
                              gene_start = 1e-6 * gene_start,
                              gene_end = 1e-6 * gene_end,
                              qtl_pos = 1e-6 * qtl_pos)

  mrna_ids <- unique(peaks.mrna$gene_id)

  # Get annotations for unique mRNA IDs
  annot.mrna <-
    dplyr::rename(
      dplyr::mutate(
        dplyr::filter(
          readRDS(file.path(datapath, "RNAseq", "annot.mrna.rds")),
          id %in% mrna_ids),
        start = start * 1e-6,
        end = end * 1e-6,
        middle_point = middle_point * 1e-6),
      pos = middle_point)

  annot.mrna <-
    dplyr::mutate(
      dplyr::inner_join(
        annot.mrna,
        dplyr::rename(
          dplyr::ungroup(
            dplyr::summarize(
              dplyr::group_by(peaks.mrna, gene_id),
              qtl_ct = n(),
              QTL = paste0(qtl_chr, "@",
                           round(qtl_pos), ":",
                           round(lod), collapse = ","),
              qtl_pos = ifelse(any(qtl_chr == chr_id &
                                   qtl_pos >= start_val &
                                   qtl_pos <= end_val),
                               qtl_pos[qtl_chr == chr_id], NA),
              lod = ifelse(is.na(qtl_pos),
                           NA, lod))),
          id = gene_id),
        by = "id"),
      local = !is.na(qtl_pos) &
        chr == chr_id &
        pmax(start,end) >= start_val &
        pmin(start,end) <= end_val)

  if(qtl) {
    annot.mrna <- dplyr::filter(annot.mrna, !is.na(qtl_pos))
    # Reduce to mRNA having peak in region.
    expr_id <- annot.mrna$id
    peaks.mrna <- dplyr::filter(peaks.mrna, gene_id %in% expr_id)
  } else {
    expr_id <- annot.mrna$id
  }


  # Get expression data.
  expr.mrna <- as.data.frame(feather::read_feather(file.path(datapath, "RNAseq", "expr.mrna.feather"),
                                     c("Mouse.ID", expr_id)))
  rownames(expr.mrna) <- expr.mrna$Mouse.ID

  list(expr = expr.mrna[,-1], annot = annot.mrna, peaks = peaks.mrna)
}

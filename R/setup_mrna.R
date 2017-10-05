setup_peaks_mrna <- function(peaks, datapath) {
  if(dir.exists(tmp <- file.path(datapath, "RNAseq"))) {
    if(file.exists(tmp <- file.path(tmp, "peaks.mrna.feather"))) {
      chrs <- levels(peaks$chr)
      peaks_new <- dplyr::select(
        dplyr::mutate(
          dplyr::rename(
            feather::read_feather(tmp),
            chr = qtl_chr,
            pos = qtl_pos),
          chr = factor(chr, chrs),
          pos = pos * 1e-6,
          pheno = paste(gene_symbol,
                        as.integer(stringr::str_replace(gene_id, "ENSMUSG", "")),
                        sep = "."),
          longname = pheno,
          output = pheno,
          pheno_group = "Islet.mRNA",
          pheno_type = paste("Islet.mRNA", gene_chr, sep = ".")),
        pheno,longname,output,pheno_group,pheno_type,chr,lod,pos)
      peaks <- dplyr::bind_rows(peaks,
                                peaks_new[, names(peaks)])
    }
  }
  peaks
}
setup_analyses_mrna <- function(analyses_tbl, peaks, datapath) {
  peaks_mrna <- dplyr::distinct(
    dplyr::filter(peaks,
                  pheno_group == "Islet.mRNA"),
    pheno,output,pheno_group,longname,pheno_type)

  analyses_mrna <- data.frame(pheno_group = "Islet.mRNA",
                              model = "normal",
                              transf = "identity",
                              offset = 0,
                              winsorize = FALSE,
                              sex = TRUE,
                              DOwave = TRUE,
                              diet_days = TRUE,
                              stringsAsFactors = FALSE)
  analyses_mrna <- inner_join(peaks_mrna,
                              analyses_mrna, by = "pheno_group")

  analyses_tbl <- dplyr::bind_rows(analyses_tbl,
                                   analyses_mrna)

  analyses_tbl
}

setup_data_mrna <- function(pheno_data, peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "RNAseq"))) {
    if(file.exists(filename <- file.path(datapath, "peaks.mrna.feather"))) {
      peaks_new <- dplyr::mutate(
        dplyr::rename(
          dplyr::distinct(
            feather::read_feather(filename),
            gene_id, gene_symbol),
          pheno = gene_id),
        pheno_gene = paste(gene_symbol,
                           as.integer(stringr::str_replace(pheno, "ENSMUSG", "")),
                           sep = "."))
      if(file.exists(filename <- file.path(datapath, "expr.mrna.feather"))) {
        pheno_new <- feather::read_feather(filename)
        pheno_data <- setup_data_append(pheno_data, pheno_new, peaks_new)
        # Now put names we want for mrna
        m <- match(peaks_new$pheno, colnames(pheno_data))
        colnames(pheno_data)[m] <- peaks_new$pheno_gene
      }
    }
  }
  pheno_data
}





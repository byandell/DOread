setup_peaks_mrna <- function(peaks, datapath) {
  tissue <- "Islet"
  if(dir.exists(datapath <- file.path(datapath, "RNAseq"))) {
    if(file.exists(filename <- file.path(datapath, "peaks.mrna.feather"))) {
      chrs <- levels(peaks$chr)
      peaks_new <- dplyr::select(
        dplyr::mutate(
          dplyr::rename(
            feather::read_feather(filename),
            chr = qtl_chr,
            pos = qtl_pos),
          chr = factor(chr, chrs),
          pheno = pheno_tissue(tissue, symbol, gene_id),
          longname = pheno,
          output = pheno,
          pheno_group = "Islet.mRNA",
          pheno_type = paste("Islet.mRNA", gene_chr, sep = ".")),
        pheno,longname,output,pheno_group,pheno_type,chr,lod,pos)
      peaks <- dplyr::bind_rows(peaks,
                                peaks_new[, names(peaks)])
    }
    peaks <- setup_peaks_mrna_module(peaks, datapath)
    peaks <- dplyr::filter(peaks,
                           pheno != "Islet.MEgrey")
    peaks <- setup_peaks_mrna_module(peaks, datapath,
                                      "peaks_pc.csv", "Hotspot")
  }
  peaks
}
setup_peaks_mrna_module <- function(peaks, datapath,
                                    filename = "peaks_me.csv",
                                    type = "Module") {
  tissue <- "Islet"
  if(file.exists(filename <- file.path(datapath, filename))) {
    chrs <- levels(peaks$chr)
    peaks_new <- dplyr::select(
      dplyr::mutate(
        read.csv(filename),
        chr = factor(chr, chrs),
        pheno = pheno_tissue(tissue, lodcolumn),
        longname = pheno,
        output = pheno,
        pheno_group = "Islet.mRNA",
        pheno_type = paste0("Islet.mRNA.", type)),
      pheno,longname,output,pheno_group,pheno_type,chr,lod,pos)
    peaks <- dplyr::bind_rows(peaks,
                              peaks_new[, names(peaks)])
  }
  peaks
}
pheno_tissue <- function(tissue, pheno, id = NULL) {
  pheno <- paste(tissue, pheno, sep = ".")
  if(!is.null(id))
    pheno <- paste(pheno,
                   as.integer(stringr::str_replace(id, "ENSMUSG", "")),
                   sep = ".")
  pheno
}

setup_analyses_mrna <- function(analyses_tbl, peaks, datapath) {
  peaks_new <- dplyr::distinct(
    dplyr::filter(peaks,
                  pheno_group == "Islet.mRNA"),
    pheno,output,pheno_group,longname,pheno_type)

  analyses_new <- data.frame(pheno_group = "Islet.mRNA",
                              model = "normal",
                              transf = "identity",
                              offset = 0,
                              winsorize = FALSE,
                              sex = TRUE,
                              DOwave = TRUE,
                              diet_days = TRUE,
                              stringsAsFactors = FALSE)
  analyses_new <- inner_join(peaks_new,
                             analyses_new, by = "pheno_group")

  analyses_tbl <- dplyr::bind_rows(analyses_tbl,
                                   analyses_new)

  analyses_tbl
}

setup_data_mrna <- function(pheno_data, peaks, datapath) {
  tissue <- "Islet"
  if(dir.exists(datapath <- file.path(datapath, "RNAseq"))) {
    if(file.exists(filename <- file.path(datapath, "peaks.mrna.feather"))) {
      # Read peaks again but keep original pheno and new name as pheno_rename.
      peaks_new <- dplyr::mutate(
        dplyr::rename(
          dplyr::distinct(
            feather::read_feather(filename),
            gene_id, symbol),
          pheno = gene_id),
        pheno_rename = pheno_tissue(tissue, symbol, pheno))

      if(file.exists(filename <- file.path(datapath, "expr.mrna.feather"))) {
        # Read phenotypes; append; change names from pheno to pheno_rename.
        pheno_new <- data.frame(
          feather::read_feather(filename),
          check.names = FALSE)
        rownames(pheno_new) <- pheno_new$Mouse.ID
        pheno_new$Mouse.ID <- NULL
        pheno_data <- setup_data_append(pheno_data, pheno_new, peaks_new)
        # Now put names we want for mrna
        m <- match(peaks_new$pheno, colnames(pheno_data))
        colnames(pheno_data)[m] <- peaks_new$pheno_rename
      }
    }
    pheno_data <- setup_data_mrna_module(pheno_data, peaks, datapath)
    pheno_data <- setup_data_mrna_module(pheno_data, peaks, datapath,
                                         "peaks_pc.csv", "pheno_pc.csv")
  }
  pheno_data
}
setup_data_mrna_module <- function(pheno_data, peaks, datapath,
                                   peakname = "peaks_me.csv",
                                   dataname = "ME.mrna.csv") {
  tissue <- "Islet"
  if(file.exists(filename <- file.path(datapath, peakname))) {
    peaks_new <- dplyr::mutate(
      dplyr::rename(
        dplyr::distinct(
          read.csv(filename),
          lodcolumn),
        pheno = lodcolumn),
      pheno_rename = pheno_tissue(tissue, pheno))

    if(file.exists(filename <- file.path(datapath, dataname))) {
      pheno_new <- read.csv(filename, row.names = 1)
      rownames(pheno_new) <- as.integer(stringr::str_replace(rownames(pheno_new), "DO-", ""))
      pheno_data <- setup_data_append(pheno_data, pheno_new, peaks_new)
      # Now put names we want for mrna
      m <- match(peaks_new$pheno, colnames(pheno_data))
      colnames(pheno_data)[m] <- peaks_new$pheno_rename
    }
  }
  pheno_data
}





setup_peaks_bileacid <- function(peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "otu"))) {
    ## Replace old bile acid with new BileAcid
    if(file.exists(filename <- file.path(datapath, "peaks_BileAcid.rds"))) {
      peaks <- dplyr::filter(peaks,
                             !(pheno_type == "bile acid"))
      peaks_new <- readRDS(filename)
      peaks <- dplyr::bind_rows(peaks,
                                peaks_new[, names(peaks)])
    }
  }
  peaks
}

setup_analyses_bileacid <- function(analyses_tbl, peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "otu"))) {
    ## Replace bile acid with new BileAcid
    if(file.exists(filename <- file.path(datapath, "analyses_BileAcid.rds"))) {
      analyses_tbl <- dplyr::filter(analyses_tbl,
                                    !(pheno_type == "bile acid"))
      analyses_new <- dplyr::filter(readRDS(filename),
                                    output %in% peaks$output)
      analyses_new$model <- "normal"
      analyses_tbl <- dplyr::bind_rows(analyses_tbl,
                                       analyses_new[, names(analyses_tbl)])
    }
  }
  analyses_tbl
}

setup_data_bileacid <- function(pheno_data, peaks, datapath) {
  ## Remove old bile acid
  peaks_old <- dplyr::distinct(
    dplyr::filter(
      readRDS(file.path(datapath, "peaks.rds")),
      pheno_type == "bile acid"),
    pheno)
  m <- match(peaks_old$pheno, colnames(pheno_data), nomatch = 0)
  pheno_data <- pheno_data[, -m[m>0]]

  if(dir.exists(datapath <- file.path(datapath, "otu"))) {
    ## Replace old bile acid with new BileAcid

    ## Find new bile acid phenotype data.
    peaks_new <- dplyr::filter(
      dplyr::distinct(peaks, pheno_group, pheno),
      pheno_group == "BileAcid")
    pheno_new <- readRDS(file.path(datapath, "pheno_BileAcid.rds"))
    # Add total.ba.
    pheno_new$total.ba <- rowSums(
      pheno_new[, seq_len(grep("primary", colnames(pheno_new)) - 1)],
      na.rm = TRUE)
    # Remove any columns from pheno_data that still match.
    m <- match(colnames(pheno_new), colnames(pheno_data), nomatch = 0)
    pheno_data <- pheno_data[, -m[m>0]]

    # Append new data.
    pheno_data <- setup_data_append(pheno_data, pheno_new, peaks_new)
  }
  pheno_data
}


setup_peaks_otucr <- function(peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "otu"))) {
    ## Add Closed Reference OTUs, OTU Modules, and Bile Acids
    peaks <- read_otu(peaks, datapath, "peaks_OTU_CR.rds")
    peaks <- read_otu(peaks, datapath, "peaks_OTU_Module.rds")
  }
  peaks
}
read_otu <- function(peaks, datapath, filename) {
  if(file.exists(filename <- file.path(datapath, filename))) {
    peaks_new <- readRDS(filename)
    peaks <- dplyr::bind_rows(peaks,
                              peaks_new[, names(peaks)])
  }
  peaks
}

setup_analyses_otucr <- function(analyses_tbl, peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "otu"))) {
    ## Add Closed Reference OTUs, OTU Modules, and Bile Acids
    analyses_tbl <- read_otua(analyses_tbl, peaks,
                              datapath, "analyses_OTU_CR.rds")
    analyses_tbl <- read_otua(analyses_tbl, peaks,
                              datapath, "analyses_OTU_Module.rds")

    # Add model column to analyses_tbl
    analyses_tbl <-
      dplyr::mutate(
        analyses_tbl,
        model = ifelse(pheno_type == "OTU_Bin",
                       "binary",
                       "normal"))
    analyses_tbl <-
      dplyr::select(
        analyses_tbl,
        pheno:pheno_type, model, transf:ncol(analyses_tbl))
  }
  analyses_tbl
}
read_otua <- function(analyses_tbl, peaks, datapath, filename) {
  if(file.exists(filename <- file.path(datapath, filename))) {
    analyses_otu <-
      dplyr::filter(
        readRDS(filename),
        output %in% peaks$output)
    analyses_tbl <-
      dplyr::bind_rows(analyses_tbl,
                       analyses_otu[, names(analyses_tbl)])
  }
  analyses_tbl
}

setup_data_otucr <- function(pheno_data, peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "otu"))) {
    peaks_new <- dplyr::filter(
      dplyr::distinct(peaks, pheno_group, pheno),
      pheno_group %in% c("OTU_Closed_Ref", "OTU_Module"))
    pheno_new <- readRDS(file.path(datapath, "pheno_OTU_CR.rds"))
    pheno_data <- setup_data_append(pheno_data, pheno_new, peaks_new)
    pheno_new <- readRDS(file.path(datapath, "pheno_OTU_Module.rds"))
    pheno_data <- setup_data_append(pheno_data, pheno_new, peaks_new)
  }
  pheno_data
}


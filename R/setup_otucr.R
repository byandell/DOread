setup_peaks_otucr <- function(peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "otu"))) {
    ## Add Closed Reference OTUs, OTU Modules, and Bile Acids
    if(file.exists(filename <- file.path(datapath, "peaks_OTU_CR.rds"))) {
      peaks_new <- readRDS(filename)
      peaks <- dplyr::bind_rows(peaks,
                                peaks_new[, names(peaks)])
    }
    peaks <- setup_peaks_otucr_module(peaks, datapath)
  }
  peaks
}
setup_peaks_otucr_module <- function(peaks, datapath) {
  if(file.exists(filename <- file.path(datapath, "peaks_OTU_Module.rds"))) {
    peaks_new <- dplyr::mutate(
      readRDS(filename),
      pheno = pheno_tissue("OTU", pheno), # prefix MEname by OTU.
      longname = pheno,
      output = pheno,
      pheno_group = "OTU_Closed_Ref") # Change group from OTU_Module
    # Drop grey module.
    peaks_new <- dplyr::filter(
      peaks_new,
      pheno != "OTU.MEgrey")

    peaks <- dplyr::bind_rows(peaks,
                              peaks_new[, names(peaks)])
  }
  peaks
}

setup_analyses_otucr <- function(analyses_tbl, peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "otu"))) {
    ## Add Closed Reference OTUs, OTU Modules, and Bile Acids
    if(file.exists(filename <- file.path(datapath, "analyses_OTU_CR.rds"))) {
      analyses_new <-
        dplyr::filter(
          readRDS(filename),
          output %in% peaks$output)
      analyses_tbl <-
        dplyr::bind_rows(analyses_tbl,
                         analyses_new[, names(analyses_tbl)])
    }

    analyses_tbl <- setup_analyses_otucr_module(analyses_tbl, peaks, datapath)

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
setup_analyses_otucr_module <- function(analyses_tbl, peaks, datapath) {
  if(file.exists(filename <- file.path(datapath, "analyses_OTU_Module.rds"))) {
    analyses_new <-
      dplyr::filter(
        dplyr::mutate(
          readRDS(filename),
          pheno = pheno_tissue("OTU", pheno),
          longname = pheno,
          output = pheno,
          pheno_group = "OTU_Closed_Ref"),
        output %in% peaks$output)
    analyses_tbl <-
      dplyr::bind_rows(analyses_tbl,
                       analyses_new[, names(analyses_tbl)])
  }
  analyses_tbl
}

setup_data_otucr <- function(pheno_data, peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "otu"))) {
    if(file.exists(filename <- file.path(datapath, "pheno_OTU_CR.rds"))) {
      peaks_new <- dplyr::distinct(
        dplyr::filter(
          peaks,
          pheno_group == "OTU_Closed_Ref",
          pheno_type != "OTU_Module"),
        pheno)

      pheno_new <- readRDS(filename)
      pheno_data <- setup_data_append(pheno_data, pheno_new, peaks_new)
    }

    pheno_data <- setup_data_otucr_module(pheno_data, peaks, datapath)
  }
  pheno_data
}
setup_data_otucr_module <- function(pheno_data, peaks, datapath) {
  if(file.exists(filename <- file.path(datapath, "peaks_OTU_Module.rds"))) {
    peaks_new <- dplyr::mutate(
      dplyr::distinct(
        dplyr::filter(
          readRDS(filename),
          pheno_type == "OTU_Module"),
        pheno),
      pheno_rename = pheno_tissue("OTU", pheno))

    pheno_new <- readRDS(file.path(datapath, "pheno_OTU_Module.rds"))
    pheno_data <- setup_data_append(pheno_data, pheno_new, peaks_new)
    # Now put names we want for mrna
    m <- match(peaks_new$pheno, colnames(pheno_data))
    colnames(pheno_data)[m] <- peaks_new$pheno_rename
  }
  pheno_data
}


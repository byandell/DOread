setup_peaks_molecules <- function(peaks, datapath) {
  if(dir.exists(tmp <- file.path(datapath, "molecule"))) {
    if(file.exists(tmp <- file.path(tmp, "peaks_small.rds"))) {
      peaks_new <- readRDS(tmp)
      peaks <- dplyr::bind_rows(peaks,
                                peaks_new[, names(peaks)])
    }
  }
  peaks
}

setup_analyses_molecules <- function(analyses_tbl, peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "molecule"))) {
    if(file.exists(filename <- file.path(datapath, "analyses_small.rds"))) {
      ## Add rows for small molecules here.
      analyses_tbl <- dplyr::bind_rows(analyses_tbl,
                                       readRDS(filename))
    }
  }
  analyses_tbl
}

setup_data_molecules <- function(pheno_data, peaks, datapath) {
  if(dir.exists(datapath <- file.path(datapath, "molecule"))) {
    if(file.exists(filename <- file.path(datapath, "pheno_small.rds"))) {
      peaks_new <- dplyr::filter(
        dplyr::distinct(peaks, pheno_group, pheno),
        pheno_group == "Molecule")
      pheno_new <- readRDS(filename)
      pheno_data <- setup_data_append(pheno_data, pheno_new, peaks_new)
    }
  }
  pheno_data
}



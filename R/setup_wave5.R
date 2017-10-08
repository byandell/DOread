# See inst/DerivedData/convertRDS.Rmd for conversion from RData to RDS files.

setup_peaks_wave5 <- function(datapath) {

  peaks <- readRDS(file.path(datapath, "peaks.rds"))

  ## Filter peaks and analyses_tbl to "best" analysis (first encountered in table).
  peak_info <- dplyr::ungroup(
    dplyr::summarize(
      dplyr::arrange(
        dplyr::distinct(
          dplyr::group_by(peaks, pheno),
          output),
        dplyr::desc(output)),
      output = output[1]))$output

  peaks <- dplyr::filter(peaks, output %in% peak_info)
  dplyr::mutate(peaks, chr = factor(chr, c(1:19,"X")))
}

setup_analyses_wave5 <- function(peaks, datapath) {
  dplyr::filter(readRDS(file.path(datapath, "analyses.rds")),
                output %in% peaks$output)
}

#' @export
setup_data_wave5 <- function(analyses_tbl, peaks, datapath) {
  ## Want to filter to "best" analysis based -- anal1 or anal2
  analyses_std <- dplyr::filter(analyses_tbl,
                                pheno_group %in% c("clin","otu","otufam","gutMB"))
  pheno_data <- read_pheno_tbl(analyses_std, datapath)
  dplyr::select(pheno_data,
                which(names(pheno_data) %in% peaks$pheno))
}

setup_covar_wave5 <- function(datapath) {
  # Read covar matrix prepared by Karl Broman.
  readRDS(file.path(datapath, "covar.rds"))
}

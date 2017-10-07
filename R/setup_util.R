# Want to use file.path(datapath, "source_pheno.csv").
# Need also files for peaks and analyses.
#
#' @export
#' @importFrom dplyr arrange bind_rows desc distinct filter funs group_by
#'                   mutate mutate_if select summarize ungroup
setup_peaks <- function(datapath) {

  ## Filter peaks and analyses_tbl to "best" analysis (first encountered in table).
  peaks <- readRDS(file.path(datapath, "peaks.rds"))

  peak_info <- dplyr::ungroup(
    dplyr::summarize(
      dplyr::arrange(
        dplyr::distinct(
          dplyr::group_by(peaks, pheno),
          output),
        dplyr::desc(output)),
      output = output[1]))$output

  peaks <- dplyr::filter(peaks, output %in% peak_info)
  peaks <- dplyr::mutate(peaks, chr = factor(chr, c(1:19,"X")))

  ## Add in newer peaks.
  peaks <- setup_peaks_otucr(peaks, datapath)
  peaks <- setup_peaks_bileacid(peaks, datapath)
  peaks <- setup_peaks_molecules(peaks, datapath)
  peaks <- setup_peaks_mrna(peaks, datapath)

  peaks
}

#' @export
setup_analyses <- function(peaks, datapath) {
  analyses_tbl <- dplyr::filter(readRDS(file.path(datapath, "analyses.rds")),
                                output %in% peaks$output)

  ## Add OTU_CR and BileAcid analyses. These have DOwaven columns.
  analyses_tbl <- setup_analyses_otucr(analyses_tbl, peaks, datapath)
  analyses_tbl <- setup_analyses_bileacid(analyses_tbl, peaks, datapath)

  ## Change DOwave columns to one character column; add small molecule batches
  analyses_tbl <- dplyr::select(
    dplyr::mutate(analyses_tbl,
                  DOwave = DOwave2 | DOwave3 | DOwave4 | DOwave5),
    -(DOwave2:DOwave5))

  analyses_tbl <- setup_analyses_molecules(analyses_tbl, peaks, datapath)
  analyses_tbl <- setup_analyses_mrna(analyses_tbl, peaks, datapath)

  ## Change NA to FALSE.
  analyses_tbl <- dplyr::mutate_if(analyses_tbl,
                                   is.logical,
                                   dplyr::funs(ifelse(is.na(.), FALSE, .)))

  analyses_tbl
}

#' @export
setup_data <- function(analyses_tbl, peaks, datapath) {
  ## Want to filter to "best" analysis based -- anal1 or anal2
  analyses_std <- dplyr::filter(analyses_tbl,
                                pheno_group %in% c("clin","otu","otufam","gutMB"))
  pheno_data <- read_pheno_tbl(analyses_std, datapath)
  pheno_data <- dplyr::select(pheno_data,
                              which(names(pheno_data) %in% peaks$pheno))

  pheno_data <- setup_data_otucr(pheno_data, peaks, datapath)
  pheno_data <- setup_data_bileacid(pheno_data, peaks, datapath)
  pheno_data <- setup_data_molecules(pheno_data, peaks, datapath)
  pheno_data <- setup_data_mrna(pheno_data, peaks, datapath)

  pheno_data
}
setup_data_append <- function(pheno_data, pheno_new, peaks_new) {
  pheno_new <- dplyr::select(
    as.data.frame(pheno_new),
    which(colnames(pheno_new) %in% peaks_new$pheno))
  # join, but need to reestablish rownames
  pheno_new$id <- rownames(pheno_new)
  pheno_data$id <- rownames(pheno_data)
  out <- dplyr::full_join(pheno_data, pheno_new, by = "id")
  rownames(out) <- out$id
  out$id <- NULL
  out
}

#' @export
setup_covar <- function(datapath) {
  # Read covar matrix prepared by Karl Broman.
  covar <- readRDS(file.path(datapath, "covar.rds"))

  # Convert to data frame.
  covar <- as.data.frame(covar)

  # Condense DOwave into one factor.
  if(!("DOwave" %in% colnames(covar))) {
    covar <- cbind(covar, DOwave = factor(1 + covar[,"DOwave2"] +
                                            2 * covar[,"DOwave3"] +
                                            3 * covar[,"DOwave4"] +
                                            4 * covar[,"DOwave5"]))
    covar <- covar[, -grep("DOwave[1-5]", colnames(covar))]
  }

  # Convert sex to F=0,M=1
  covar$sex = c("F","M")[1 + covar$sex]

  # Get small molecule variables
  covar_small <- readRDS(file.path(datapath, "molecule", "covar_small.rds"))

  # Combine these together.
  covar$id <- rownames(covar)
  covar_small$id <- rownames(covar_small)
  covar <- dplyr::left_join(covar, covar_small, by = "id")
  rownames(covar) <- covar$id
  covar$id <- NULL

  covar
}

#' @export
setup_type <- function(analyses_tbl) {
  c("all", sort(unique(analyses_tbl$pheno_type)))
}

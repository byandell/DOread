#' @export
#' @importFrom dplyr arrange bind_rows desc distinct filter group_by
#'                   mutate select summarize ungroup
setup_peaks <- function(datapath) {

  ## Filter peaks and analyses_tbl to "best" analysis.
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

  if(dir.exists(file.path(datapath, "otu"))) {
    ## Add Closed Reference OTUs, OTU Modules, and Bile Acids
    read_otu <- function(peaks, datapath, filename) {
      if(file.exists(file.path(datapath, "otu", filename))) {
        peaks_otu <- readRDS(file.path(datapath, "otu", filename))
        peaks <- dplyr::bind_rows(peaks,
                                  peaks_otu[, names(peaks)])
      }
      peaks
    }
    peaks <- read_otu(peaks, datapath, "peaks_OTU_CR.rds")
    peaks <- read_otu(peaks, datapath, "peaks_OTU_Module.rds")
    ## Replace bile acid with new BileAcid
    if(file.exists(file.path(datapath, "otu", "peaks_BileAcid.rds")))
      peaks <- dplyr::filter(peaks, !(pheno_type == "bile acid"))
    peaks <- read_otu(peaks, datapath, "peaks_BileAcid.rds")
  }
  if(dir.exists(file.path(datapath, "molecule"))) {
    if(file.exists(file.path(datapath, "molecule", "peaks_small.rds"))) {
      peaks_small <- readRDS(file.path(datapath, "molecule", "peaks_small.rds"))
      peaks <- dplyr::bind_rows(peaks,
                                peaks_small[, names(peaks)])
    }
  }
  peaks
}

#' @export
setup_analyses <- function(peaks, datapath) {
  analyses_tbl <- dplyr::filter(readRDS(file.path(datapath, "analyses.rds")),
                                output %in% peaks$output)
  if(dir.exists(file.path(datapath, "otu"))) {
    ## Add Closed Reference OTUs, OTU Modules, and Bile Acids
    read_otua <- function(analyses_tbl, peaks, datapath, filename) {
      if(file.exists(file.path(datapath, "otu", filename))) {
        analyses_otu <-
          dplyr::filter(
            readRDS(file.path(datapath, "otu", filename)),
            output %in% peaks$output)
        analyses_tbl <-
          dplyr::bind_rows(analyses_tbl,
                           analyses_otu[, names(analyses_tbl)])
      }
      analyses_tbl
    }
    analyses_tbl <- read_otua(analyses_tbl, peaks, datapath, "analyses_OTU_CR.rds")
    analyses_tbl <- read_otua(analyses_tbl, peaks, datapath, "analyses_OTU_Module.rds")
    ## Replace bile acid with new BileAcid
    analyses_tbl <- dplyr::filter(analyses_tbl, !(pheno_type == "bile acid"))
    analyses_tbl <- read_otua(analyses_tbl, peaks, datapath, "analyses_BileAcid.rds")

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

  ## Change DOwave columns to one character column; add small molecule batches
  analyses_tbl <- dplyr::select(
    dplyr::mutate(analyses_tbl,
                  DOwave = DOwave2 | DOwave3 | DOwave4 | DOwave5),
    -(DOwave2:DOwave5))

  if(dir.exists(file.path(datapath, "molecule"))) {
    if(file.exists(file.path(datapath, "molecule", "analyses_small.rds"))) {
      ## Add rows for small molecules here.
      analyses_molecule <- readRDS(file.path(datapath, "molecule", "analyses_small.rds"))
      analyses_tbl <- dplyr::bind_rows(analyses_tbl, analyses_molecule)
    }
  }

  analyses_tbl[is.na(analyses_tbl)] <- FALSE

  analyses_tbl
}

#' @export
setup_data <- function(analyses_tbl, peaks, datapath) {
  ## Want to filter to "best" analysis based -- anal1 or anal2
  analyses_std <- dplyr::filter(analyses_tbl,
                                pheno_group %in% c("clin","otu","otufam"))
  pheno_data <- read_pheno_tbl(analyses_std, datapath)
  pheno_data <- dplyr::select(pheno_data,
                              which(names(pheno_data) %in% peaks$pheno))

  tmpfn <- function(pheno_data, pheno_otu, peaks_otu) {
    pheno_otu <- dplyr::select(
      as.data.frame(pheno_otu),
      which(colnames(pheno_otu) %in% peaks_otu$pheno))
    # join, but need to reestablish rownames
    pheno_otu$id <- rownames(pheno_otu)
    pheno_data$id <- rownames(pheno_data)
    out <- dplyr::full_join(pheno_data, pheno_otu, by = "id")
    rownames(out) <- out$id
    out$id <- NULL
    out
  }
  ## Add Closed Reference OTUs
  if(dir.exists(file.path(datapath, "otu"))) {
    peaks_otu <- dplyr::filter(
      dplyr::distinct(peaks, pheno_group, pheno),
      !(pheno_group %in% c("clin","gutMB","otu","otufam")))
    pheno_otu <- readRDS(file.path(datapath, "otu", "pheno_OTU_CR.rds"))
    pheno_data <- tmpfn(pheno_data, pheno_otu, peaks_otu)
  }
  if(dir.exists(file.path(datapath, "molecule"))) {
    if(file.exists(file.path(datapath, "molecule", "pheno_small.rds"))) {
      peaks_molecule <- dplyr::filter(
        dplyr::distinct(peaks, pheno_group, pheno),
        pheno_group == "Molecule")
      pheno_molecule <- readRDS(file.path(datapath, "molecule", "pheno_small.rds"))
      pheno_data <- tmpfn(pheno_data, pheno_molecule, peaks_molecule)
    }
  }

  pheno_data
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

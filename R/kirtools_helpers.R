#' @include kir_allele_list.R

#' @keywords internal
get_kir_path <- function() {
  system.file(package = "kirtools")
}

#' @keywords internal
get_feature_mat <- function(x) {
  exmat.cols <- paste0(c("utr", rep(c("Exon ", "Intron "), 8), "Exon ", "utr"),
                       c(5, rep(1:8, each = 2), 9, 3))
  ms         <- hlatools::sequences(x)
  all.ranges <- hlatools::ranges(hlatools::features(x))
  intex      <- t(sapply(all.ranges, function(x) exmat.cols %in% names(x)))
  intex
}

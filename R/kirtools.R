#' @include kirtools_helpers.R
#'

#' @export
KIRGeneList <- function(kirGene = "all") {
  pdb              <- PythonDB(gene = kirGene)
  ka               <- getAllelesFromPython(pdb)
  ka
}

#' @export
KIRGeneNames <- function(kirGene = "all") {
  pdb              <- PythonDB(gene = "all")
  ka               <- getAllelesFromPython(pdb)
  na               <- getKIRGeneNames(ka)
  nm               <- getKIRAlleleCounts(ka)
  names(nm) <- na
  nm
}

#' @export
KIRGeneAlleles <- function(kirGene = "KIR2DL1") {
  gb <- KIRGeneList(kirGene = kirGene)
  getKIRGene(gb)
}


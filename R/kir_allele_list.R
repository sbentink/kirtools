#' @include python_tools.R

#' @title Class \code{KIRAlleleList}
#'
#' @description An object that stores KIRAlleles for multiple Genes.
#' @slot alleles  A \code{list} with one \code{KIRAllele} object per KIR gene.
#' @slot data_source A \code{character} with raw data information. Currently defaults to 'genebank'
#' @slot data_name A \code{character} with the raw data file name. Current default is the name of a genebank file.
setClass("KIRAlleleList",
         slots = list(alleles     = "list",
                      data_source = "character",
                      data_name   = "character"))

setMethod(f = "show", signature = "KIRAlleleList", definition = function(object) {
  kil <- object
  le  <- length(kil@alleles)
  lea <- sapply(kil@alleles,length)
  cat("KIRAlleleList that contains", le, "KIR genes\n")
  cat("Alleles per gene:\n")
  print(lea)}
)

setMethod(f = "getKIRGene", signature = "KIRAlleleList", definition = function(kil, kirgene = NA) {
  if (is.na(kirgene)) return(kil@alleles[[1]])
  geneFound <- any(names(kil@alleles) == kirgene)
  if (!geneFound) stop(paste("Gene",kirgene,"not found"))
  return(kil@alleles[[kirgene]])
})

setMethod(f = "getKIRGeneNames", signature = "KIRAlleleList", definition = function(kil) {
  return(names(kil@alleles))})
setMethod(f = "getKIRAlleleCounts", signature = "KIRAlleleList", definition = function(kil) {
  return(sapply(kil@alleles,length))})

KIRAllele <- function(my_idx, my_sequences, my_features, my_metaData) {
  ma           <- hlatools::HLAAllele()
  ma@sequence  <- Biostrings::DNAStringSet(my_sequences)[my_idx]
  ma@features  <- hlatools::HLARangesList(lapply(my_features,function(x) x$feat))[my_idx]
  ma@metadata  <- S4Vectors::DataFrame(my_metaData)[my_idx,]
  ma
}

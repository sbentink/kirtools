% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kir_allele_list.R
\docType{class}
\name{KIRAlleleList-class}
\alias{KIRAlleleList-class}
\title{Class \code{KIRAlleleList}}
\description{
An object that stores KIRAlleles for multiple Genes.
}
\section{Slots}{

\describe{
\item{\code{alleles}}{A \code{list} with one \code{KIRAllele} object per KIR gene.}

\item{\code{data_source}}{A \code{character} with raw data information. Currently defaults to 'genebank'}

\item{\code{data_name}}{A \code{character} with the raw data file name. Current default is the name of a genebank file.}
}}


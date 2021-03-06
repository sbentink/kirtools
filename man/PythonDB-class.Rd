% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/python_db.R, R/python_tools.R
\docType{class}
\name{PythonDB-class}
\alias{PythonDB-class}
\alias{getAllelesFromPython,PythonDB-method}
\title{Class \code{PythonDB}}
\usage{
\S4method{getAllelesFromPython}{PythonDB}(db, ncores = parallel::detectCores()
  - 4)
}
\arguments{
\item{db}{A \linkS4class{PythonDB} object from which we want to extract the alleles}

\item{nccores}{The number of CPU cores that should be assigned this job. It is not huge computations carried out by this function,
but parallel makes it a little bit faster anyway.}
}
\value{
A \linkS4class{KIRAlleleList}
}
\description{
S4 representing a KIR data base parsed from a Genebank file using a Python script.
}
\section{Methods (by generic)}{
\itemize{
\item \code{getAllelesFromPython}: Get list of \linkS4class{KIRAllele}s from a \linkS4class{PythonDB} object.
}}
\section{Slots}{

\describe{
\item{\code{pydb}}{A \code{list_Python} object that hold the raw data obtained from a Genebank file.}

\item{\code{locusNames}}{A \code{character} vector with the KIR gene name of each allele in \code{pydb}}

\item{\code{dbPath}}{A \code{character} with a reference to the Genebank data file used to initialize the object}

\item{\code{locusType}}{\code{character} that defaults to 'KIR'. The python script could also handle locuses other than KIR.
However we do not use this feature in this package.}
}}
\keyword{internal}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/python_db.R
\name{PythonDB}
\alias{PythonDB}
\title{Constructor for \linkS4class{PythonDB}}
\usage{
PythonDB(gene, my_py_db = new("PythonDB"))
}
\arguments{
\item{gene}{Giving a String with an actual gene name returns a \linkS4class{PythonDB} with a single entry that holds alleles for one gene. Setting this to 'all'
will cause the parser to return a \linkS4class{PythonDB} that contains the entire data base}

\item{Full}{path to KIR gene bank file}

\item{locus}{type. The python script can handle more than just KIR data. However this package does not implement anything other than KIR alles.
Hence this parameter defaults to 'KIR'. Don't change it.}

\item{pythonPath}{the path to the Python script (the core parser)}
}
\value{
A \linkS4class{PythonDB} object
}
\description{
Constructor for \linkS4class{PythonDB}
}
\examples{
\dontrun{
dbPath=file.path("extdata","2016_10_05_KIR.dat")
pythonPath="./exec"
my_db_raw  <- PythonDB(dbPath,"KIR","all",pythonPath)
my_db_list <- getAllelesFromPython(my_db_raw)
my_db_list
getKIRGene(kil=my_db_list,kirgene="KIR3DL2")
getKIRGene(kil=my_db_list,kirgene=NA)
}
}
\keyword{internal}


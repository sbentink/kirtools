#' @include ctrl_env.R
#'

setClass("PythonParams",
         slots =
           list(dbUrl      = "character",
                dbPath     = "character",
                dbName     = "character",
                dbDate     = "character",
                locusType  = "character"))

setMethod("initialize", "PythonParams",
          function(.Object,
                   dbUrl      = getOption("kt_dbUrl"),
                   dbDate     = getOption("kt_dbDate"),
                   dbPath     = getOption("kt_dbPath"),
                   dbName     = getOption("kt_dbName"),
                   ##The python scripts support more than one type of hla gene, that explains this option,
                   ##which is not really an option. Here it is all about kir.
                   locusType  = getOption("kt_locusType")) {
            KIRSetPars(kt_dbDate    = dbDate,
                       kt_dbPath    = dbPath,
                       kt_dbName    = dbName)
            .Object@dbUrl      = getOption("kt_dbUrl")
            .Object@dbPath     = getOption("kt_dbPath")
            .Object@dbName     = getOption("kt_dbName")
            .Object@dbDate     = getOption("kt_dbDate")
            .Object@locusType  = getOption("kt_locusType")
            .Object
})

#' @title Class \code{PythonDB}
#' @description S4 representing a KIR data base parsed from a Genebank file using a Python script.
#' @slot pydb A \code{list_Python} object that hold the raw data obtained from a Genebank file.
#' @slot locusNames A \code{character} vector with the KIR gene name of each allele in \code{pydb}
#' @slot dbPath A \code{character} with a reference to the Genebank data file used to initialize the object
#' @slot locusType \code{character} that defaults to 'KIR'. The python script could also handle locuses other than KIR.
#' However we do not use this feature in this package.
#' @importClassesFrom XRPython list_Python
setClass("PythonDB",
         slots =
           list(pydb       = "list_Python",
                locusNames = "character"),
         contains = "PythonParams"
)
#' Constructor for \linkS4class{PythonDB}
#' @keywords internal
#' @param Full path to KIR gene bank file
#' @param locus type. The python script can handle more than just KIR data. However this package does not implement anything other than KIR alles.
#' Hence this parameter defaults to 'KIR'. Don't change it.
#' @param gene Giving a String with an actual gene name returns a \linkS4class{PythonDB} with a single entry that holds alleles for one gene. Setting this to 'all'
#' will cause the parser to return a \linkS4class{PythonDB} that contains the entire data base
#' @param pythonPath the path to the Python script (the core parser)
#' @return A \linkS4class{PythonDB} object
#' @examples
#' \dontrun{
#' dbPath=file.path("extdata","2016_10_05_KIR.dat")
#' pythonPath="./exec"
#' my_db_raw  <- PythonDB(dbPath,"KIR","all",pythonPath)
#' my_db_list <- getAllelesFromPython(my_db_raw)
#' my_db_list
#' getKIRGene(kil=my_db_list,kirgene="KIR3DL2")
#' getKIRGene(kil=my_db_list,kirgene=NA)
#' }
PythonDB <- function(gene, my_py_db = new("PythonDB")) {
  dbPath     <- my_py_db@dbPath
  dbName     <- my_py_db@dbName
  locusType  <- my_py_db@locusType
  ev <- XRPython::RPython()
  XRPython::pythonAddToPath(evaluator = ev)
  ev$Import("hla_embl_parser", "read_dat_file_simple_locus", "read_dat_file_simple")
  myParseAll   <- XRPython::PythonFunction("read_dat_file_simple")
  myParseGene  <- XRPython::PythonFunction("read_dat_file_simple_locus")
  if (gene == "all") {
    pydb    <- myParseAll(file.path(dbPath, dbName), locusType)
  }
  else {
    pydb    <- myParseGene(file.path(dbPath, dbName), locusType, gene)
  }
  ##next we could use pop however pop currently removes the last element
  ##but it does not reduce the size.
  ##I explicity access the last element where I know that is contains the list with gene
  ##names
  py_size    <- pydb$size() - 1
  py_el      <- pydb$el(py_size)
  py_raw     <- XRPython::pythonGet(py_el)
  locusNames <- unlist(py_raw)
  if (is.null(locusNames)) stop(paste("Gene",gene,"not found"))
  my_py_db <- new("PythonDB")
  my_py_db@pydb        <- pydb
  my_py_db@locusNames  <- locusNames
  my_py_db
}

.onLoad <- function(libname, pkgname) {
  kt_dbPath     <- file.path(system.file(package = pkgname), "extdata")
  kt_dbName     <- "2016_10_05_KIR.dat"
  kt_locusType  <- "KIR"
  kt_pythonPath <- "exec"
  options(kt_dbPath     = kt_dbPath,
          kt_dbName     = kt_dbName,
          kt_locusType  = kt_locusType,
          kt_pythonPath = kt_pythonPath)
}

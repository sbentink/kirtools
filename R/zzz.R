.onLoad <- function(libname, pkgname) {
  kt_dbUrl      <- NA
  kt_dbDate     <- NA
  kt_dbPath     <- NA
  kt_dbName     <- NA
  kt_locusType  <- NA
  options(kt_dbUrl      = kt_dbUrl,
          kt_dbDate     = kt_dbDate,
          kt_dbPath     = kt_dbPath,
          kt_dbName     = kt_dbName,
          kt_locusType  = kt_locusType)
}

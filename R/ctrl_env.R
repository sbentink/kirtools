#' @include generics.R

#' @export
KIRGetPars     <- function() {
  curr_opts       <- c(getOption("kt_dbUrl"),
                       getOption("kt_dbDate"),
                       getOption("kt_dbPath"),
                       getOption("kt_dbName"),
                       getOption("kt_locusType"))
  localCopyPresent <- file.exists(file.path(curr_opts[3],curr_opts[4]))
  curr_opts        <- c(curr_opts,localCopyPresent)
  names(curr_opts) <- c("kt_dbUrl", "kt_dbDate",
                        "kt_dbPath", "kt_dbName",
                        "kt_locusType","localCopyPresent")
  curr_opts
}

#' @export
KIRSetPars    <- function(kt_dbUrl     = "ftp://ftp.ebi.ac.uk/pub/databases/ipd/kir/KIR.dat",
                          kt_dbDate    = NA,
                          kt_dbPath    = NA,
                          kt_dbName    = NA,
                          kt_locusType = "KIR") {
  valid_file <- file.exists(file.path(kt_dbPath,kt_dbName))
  valid_date <- !is.na(kt_dbDate)
  if (!(valid_file & valid_date)) {
    KIRCleanEnv()
    curr_date  <- Sys.time()
    kt_dbDate     <- as.character(curr_date)
    kt_dbName     <- paste(format(curr_date, "%Y-%m-%d|%H:%M"), "KIR.dat", sep = "_")
    kt_dbPath     <- tempdir()
    cat("Pulling new local copy of KIR data base from", kt_dbUrl, "\n")
    kt_dbObj      <- RCurl::getURL(kt_dbUrl)
    cat("Saving file", kt_dbName, "to directory", kt_dbPath, "\n")
    writeLines(kt_dbObj,file.path(kt_dbPath,kt_dbName))
  }
  options(kt_dbUrl      = kt_dbUrl,
          kt_dbDate     = kt_dbDate,
          kt_dbPath     = kt_dbPath,
          kt_dbName     = kt_dbName,
          kt_locusType  = kt_locusType)
}

#' @export
KIREditDb <- function(my_editor = "emacs",
                      kt_dbPath = getOption("kt_dbPath"),
                      kt_dbName = getOption("kt_dbName")) {
    valid_file <- file.exists(file.path(kt_dbPath,kt_dbName))
    if (!valid_file) {
      cat("No valid data base loaded")
    } else {
      my_command <- paste(my_editor, file.path(kt_dbPath,
                                     kt_dbName), "&")
      cat("Calling", my_command, "\n")
      system(my_command)
    }
}

#' @export
KIRCleanEnv <- function() {
  options(kt_dbUrl      = NA,
          kt_dbDate     = NA,
          kt_dbPath     = NA,
          kt_dbName     = NA,
          kt_locusType  = NA)
  my_temp     <- tempdir()
  my_kirfiles <- dir(my_temp, pattern = "*KIR\\.dat$")
  cat("Setting existing kirtools options to NA\n")
  cat("Found", length(my_kirfiles), "existing file(s) in directory", my_temp,"\n")
  if (length(my_kirfiles) > 0) {
    frm <- unlink(file.path(my_temp, my_kirfiles)) + 1
    if (sum(frm) > 0)  cat(paste("File",my_kirfiles[frm],"deleted\n"))
    if (sum(!frm) > 0) cat(paste("File",my_kirfiles[!frm],"not deleted\n"))
  }
}


#' @include python_db.R

#' @keywords internal
#' @importClassesFrom XRPython dict_Python
#' @importClassesFrom hlatools HLARanges
setMethod(f = "getAlleleList", signature = "dict_Python",
          definition = function(dictPy) {
            rawObj      <- XRPython::pythonGet(dictPy)
            DString     <- Biostrings::DNAString(rawObj[4])
            exonIndex      <- grep("Exon", rawObj)
            intronIndex    <- grep("Intron", rawObj)
            utrIndex       <- grep("utr", rawObj)
            featIndex      <- grep("Exon|Intron|utr", rawObj)
            featsplit      <- strsplit(rawObj[featIndex], " ")
            feat           <- sapply(featsplit, function(x) as.numeric(x[2:3]))
            featraw        <- sapply(featsplit, function(x) x[1])
            colnames(feat) <- featraw
            fid <- sapply(featsplit, function(x)
              if (x[1] %in% c("Exon", "Intron")) paste(x[1], x[4]) else x[1]
            )
            r        <- IRanges::IRanges(feat[1,], feat[2,])
            names(r) <- fid
            o        <- order(order(r))
            snames   <- rep(rawObj[1], ncol(feat))
            ty       <- gsub("utr3|utr5", "UTR", colnames(feat))
            sid     <- paste0(snames,"_", fid)
            seqnames_raw  <- rep(rawObj[1], ncol(feat))
            seqnames_s4   <- S4Vectors::Rle(seqnames_raw)
            outranges_raw <-
              hlatools::HLARanges(
                seqnames = seqnames_s4,
                ranges   = r,
                id       = sid,
                order    = o,
                type     = ty,
                status   = rep(rawObj[7],length(r)),
                frame    = integer(length(r)))
            outranges  = sort(outranges_raw)
            returnList  <- list(
              ID      = rawObj[1],
              locus   = rawObj[2],
              name    = rawObj[3],
              seq     = DString,
              len     = rawObj[5],
              isRef   = as.logical(rawObj[6]),
              fullSeq = rawObj[7],
              feat    = outranges)
            return(returnList)
          }
)

#' @keywords internal
#' @describeIn PythonDB Get list of \linkS4class{KIRAllele}s from a \linkS4class{PythonDB} object.
#' @param db A \linkS4class{PythonDB} object from which we want to extract the alleles
#' @param nccores The number of CPU cores that should be assigned this job. It is not huge computations carried out by this function,
#' but parallel makes it a little bit faster anyway.
#' @return A \linkS4class{KIRAlleleList}
#' @importClassesFrom XRPython dict_Python
setMethod(f = "getAllelesFromPython",signature = "PythonDB",
          definition = function(db, ncores = parallel::detectCores() - 4) {

            ##Do not process last element because it is the vector with all gene names
            n.chunk    <- ncores
            db_len     <- db@pydb$size() - 1
            proc.chunk <- round(seq(1, n.chunk,len = db_len))
            proc.list  <- tapply(1:db_len, proc.chunk, function(i) i)
            result     <- unlist(parallel::mclapply(proc.list,
                                                    function(n) {l = list(); for (i in 1:length(n))
                                                      l[[i]] = getAlleleList(db@pydb$el(n[i] - 1));l},
                                                    mc.cores = ncores), recursive = FALSE)
            my_sequences_raw <- sapply(result, function(x) x$seq)
            my_sequences <- Biostrings::DNAStringSet(my_sequences_raw)

            meta.raw <- t(sapply(result,function(x) x[c("ID","name","fullSeq","locus")]))
            names(my_sequences) <- unlist(meta.raw[,1])

            my_metaData <-  S4Vectors::DataFrame(allele_name   = unlist(meta.raw[,2]),
                                                 allel_id      = unlist(meta.raw[,1]),
                                                 date_assigned = rep(NA,nrow(meta.raw)),
                                                 cwd_status    = rep(NA,nrow(meta.raw)),
                                                 complete      = as.logical(unlist(meta.raw[,3])),
                                                 pmid          = rep(NA,nrow(meta.raw)),
                                                 ethnicity     = rep(NA,nrow(meta.raw)),
                                                 locus         = unlist(meta.raw[,4]))

            my.alleles.ind <- tapply(1:nrow(my_metaData), my_metaData$locus, function(i) i)
            my.alleles     <- parallel::mclapply(my.alleles.ind,KIRAllele, my_sequences,result,
                                                 my_metaData, mc.cores = ncores)
            ka <- new("KIRAlleleList")
            ka@alleles     <- my.alleles
            ka@data_source <- "genebank"
            ka@data_name   <- sapply(strsplit(db@dbPath, .Platform$file.sep), function(x) x[length(x)])
            ka
          })

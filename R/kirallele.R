#Generics:
setGeneric(name="getAlleleList",def=function(dictPy)  {standardGeneric("getAlleleList")})
setGeneric(name="getAllelesFromPython",def=function(db,...) {standardGeneric("getAllelesFromPython")})
setGeneric(name="getKIRGene",def=function(kil,...) {standardGeneric("getKIRGene")})
setGeneric(name="getKIRGeneNames",def=function(kil,...) {standardGeneric("getKIRGeneNames")})


##Class definitions

#' Class \code{PythonDB}
#'
#' @description {S4 representing a KIR data base parsed from a Genebank file using a Python script.}
#' @slot {pydb A \code{list_Python} object that hold the raw data obtained from a Genebank file.}
#' @slot {locusNames A \code{character} vector with the KIR gene name of each allele in \code{pydb}}
#' @slot {dbPath A \code{character} with a reference to the Genebank data file used to initialize the object}
#' @slot {locusType \code{character} that defaults to 'KIR'. The python script could also handle locuses other than KIR.
#' However we do not use this feature in this package.}
#' @importClassesFrom XRPython list_Python
setClass("PythonDB",
         slots=
           list(pydb="list_Python",
                locusNames="character",
                dbPath="character",
                locusType="character",
                pythonPath="character")
         )

#' Class \code{KIRAllele}
#'
#' @description {A wrapper extending \code{HLAAllele-class} imported from \code{hlatools}.
#' No methods or attributes in addition to the ones inherrited from \code{HLAAllele} are implemented in this implementation.
#' This class should allows to handle KIR alleles in a separate object.}
#' @importClassesFrom hlatools HLAAllele
setClass("KIRAllele",
         contains="HLAAllele")

#' Class \code{KIRAlleleList}
#'
#' @description {An object that stores KIRAlleles for multiple Genes.}
#' @slot alleles  A \code{list} with one \code{KIRAllele} object per KIR gene.
#' @slot data_source A \code{character} with raw data information. Currently defaults to 'genebank'
#' @slot data_name A \code{character} with the raw data file name. Current default is the name of a genebank file.
setClass("KIRAlleleList",
         slots=list(alleles="list",
                    data_source="character",
                    data_name="character"))

PythonDB <- function(dbPath,locusType,gene,pythonPath) {
            ev <- XRPython::RPython()
            XRPython::pythonAddToPath(pythonPath,evaluator = ev)
            ev$Import("hla_embl_parser","read_dat_file_simple_locus","read_dat_file_simple")
            if (gene=="all") {
              myParse <- XRPython::PythonFunction("read_dat_file_simple")
              pydb=myParse(dbPath,locusType)
            }
            else {
              myParse <- XRPython::PythonFunction("read_dat_file_simple_locus")
              pydb=myParse(dbPath,locusType,gene)
            }
            ##next we could use pop however pop currently removes the last element
            ##but it does not reduce the size.
            ##I explicity access the last element where I know that is contains the list with gene
            ##names
            locusNames <- unlist(XRPython::pythonGet(pydb$el(pydb$size()-1)))
            if (is.null(locusNames)) stop(paste("Gene",gene,"not found"))
            new("PythonDB",
                  pydb=pydb,
                  locusNames=locusNames,
                  dbPath=dbPath,
                  locusType=locusType,
                  pythonPath=pythonPath)}

#' @importClassesFrom XRPython dict_Python
setMethod(f="getAlleleList",signature="dict_Python",
          definition=function(dictPy) {
            rawObj      <- XRPython::pythonGet(dictPy)
            DString     <- Biostrings::DNAString(rawObj[4])
            exonIndex   <- grep("Exon",rawObj)
            intronIndex <- grep("Intron",rawObj)
            utrIndex    <- grep("utr",rawObj)
            featIndex   <- grep("Exon|Intron|utr",rawObj)
            featsplit   <- strsplit(rawObj[featIndex]," ")
            feat        <- sapply(featsplit,function(x) as.numeric(x[2:3]))
            featraw     <- sapply(featsplit,function(x) x[1])
            colnames(feat) <- featraw
            fid <- sapply(featsplit,function(x) if (x[1] %in% c("Exon","Intron")) paste(x[1],x[4]) else {x[1]})
            r <- IRanges::IRanges(feat[1,],feat[2,])
            names(r) <- fid
            o <- order(order(r))
            snames <- rep(rawObj[1],ncol(feat))
            ty <- gsub("utr3|utr5","UTR",colnames(feat))
            ##sid    <- paste0(snames,"_",colnames(feat),"_",sapply(paste0("00",o),function(x) substr(x,nchar(x)-2,nchar(x))))
            sid <- paste0(snames,"_",fid)
            outranges=sort(
              hlatools::HLARanges(seqnames=S4Vectors::Rle(rep(rawObj[1],ncol(feat))),
              ranges=r,
              id=sid,
              order=o,
              type=ty,
              status=rep(rawObj[7],length(r)),
              frame=integer(length(r))))
              returnList  <- list(
                ID=rawObj[1],
                locus=rawObj[2],
                name=rawObj[3],
                seq=DString,
                len=rawObj[5],
                isRef=as.logical(rawObj[6]),
                fullSeq=rawObj[7],
                feat=outranges)
              return(returnList)
}
)

#' Get list of \code{KIRAllele}s from a \code{PythonDB} object.
#' @param {db A \code{PythonDB object from which we want to extract the allels}}
#' @param {nccores The number of CPU cores that should be assigned this job. It is not huge computations carried out by this function,
#' but parallel makes it a little bit faster anyway.}
#' @return {A \code{list} of \code{KIRAllele} objects}
setMethod(f="getAllelesFromPython",signature="PythonDB",
          definition=function(db,ncores=parallel::detectCores()-4) {

            ##Do not process last element because it is the vector with all gene names
            n.chunk <- ncores
            proc.chunk <- round(seq(1,n.chunk,len=db@pydb$size()-1))
            proc.list  <- tapply(1:(db@pydb$size()-1),proc.chunk,function(i) i)
            result     <- unlist(parallel::mclapply(proc.list,
                                        function(n) {l=list(); for (i in 1:length(n))
                                          l[[i]] <- getAlleleList(db@pydb$el(n[i]-1));l},
                                          mc.cores=ncores),recursive=FALSE)
            sequences <- Biostrings::DNAStringSet(sapply(result,function(x) x$seq))
            ##names(sequences) <- sapply(result,function(x) x$ID)

            meta.raw <- t(sapply(result,function(x) x[c("ID","name","fullSeq","locus")]))
            names(sequences) <- unlist(meta.raw[,1])

            metaData <-  S4Vectors::DataFrame(allele_name=unlist(meta.raw[,2]),
                                              allel_id=unlist(meta.raw[,1]),
                                              date_assigned=rep(NA,nrow(meta.raw)),
                                              cwd_status=rep(NA,nrow(meta.raw)),
                                              complete=as.logical(unlist(meta.raw[,3])),
                                              pmid=rep(NA,nrow(meta.raw)),
                                              ethnicity=rep(NA,nrow(meta.raw)),
                                              locus=unlist(meta.raw[,4]))
            my.alleles.ind <- tapply(1:nrow(metaData),metaData$locus,function(i) i)
            my.alleles     <- parallel::mclapply(my.alleles.ind,KIRAllele,sequences,result,metaData,mc.cores=ncores)
            my.alleles})

setMethod(f="show", signature="KIRAlleleList",definition=function(object) {
  kil=object
  le  <- length(kil@alleles)
  lea <- sapply(kil@alleles,length)

  cat("KIRAlleleList that contains", le, "KIR alleles\n")
  cat("Alleles per gene:\n")
  print(lea)}
  )

setMethod(f="getKIRGene",signature="KIRAlleleList",definition=function(kil,kirgene=NA) {
  if (is.na(kirgene)) return(kil@alleles[[1]])
  geneFound <- any(names(kil@alleles)==kirgene)
  if (!geneFound) stop(paste("Gene",kirgene,"not found"))
  return(kil@alleles[[kirgene]])
})

setMethod(f="getKIRGeneNames",signature="KIRAlleleList",definition=function(kil) {
  return(names(kil@alleles))})

KIRAllele <- function(idx,sequences,features,metaData) {
  ma           <- new("KIRAllele")
  ma@sequence  <- sequences[idx]
  ma@features  <- hlatools::HLARangesList(lapply(features,function(x) x$feat))[idx]
  ma@metadata  <- metaData[idx,]
  ma
}

KIRAlleleList_gb <- function(
  kirGene="all",
  dbPath=file.path("extdata","2016_10_05_KIR.dat"),
  pythonPath="./exec",
  locusType="KIR") {
  ka      <- new("KIRAlleleList")
  my_db   <- PythonDB(dbPath,locusType,kirGene,pythonPath)
  ka@alleles     <- getAllelesFromPython(my_db)
  ka@data_source <- "genebank"
  ka@data_name   <- sapply(strsplit(dbPath,.Platform$file.sep),function(x) x[length(x)])
  ka
}

#' @export
KIRAlleleFromDB <- function(genename="KIR2DL1") {
  getKIRGene(KIRAlleleList_gb(kirGene = genename))
}

calc_common_exon_distance <- function(x, selex = c("Exon 3","Exon 4","Exon 5"), verbose = TRUE) {
  stopifnot(requireNamespace("DECIPHER", quietly = TRUE))
  ms <- hlatools::sequences(x)
  all.ranges <- hlatools::ranges(hlatools::features(x))
  intex <- lapply(all.ranges,function(x) x[names(x) %in% selex])
  exseq <- Biostrings::DNAStringSet(sapply(1:length(ms),
                                           function(i) unlist(GenomicFeatures::extractTranscriptSeqs(ms[[i]],IRanges::IRangesList(intex[[i]])))))
  aln <- DECIPHER::AlignSeqs(exseq, iterations = 0, refinements = 0,
                             restrict = -500, verbose = verbose)
  DECIPHER::DistanceMatrix(aln, includeTerminalGaps = TRUE, verbose = verbose)
}


##calc_common_exon_distance(k@alleles[[1]])
##myDBsingle   <- PythonDB(dbPath,locusType,"KIR2DL1",pythonPath)
##myDBall      <- PythonDB(dbPath,locusType,"all",pythonPath)
##singleAllele <- getAllelesFromPython(myDBsingle)
##allAlleles   <- getRangesFromPython(myDBall)

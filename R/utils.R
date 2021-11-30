#' @importFrom stats complete.cases
NULL

N_MAX_MISMATCHES <- 3

utils::globalVariables(c('.', "SpCas9", "AsCas12a"))


#' @importFrom methods is
.validateCrisprNuclease <- function(crisprNuclease){
    if (is.null(crisprNuclease)){
        crisprNuclease <- .getDefaultCrisprNuclease()
    } else {
        if (!is(crisprNuclease, "CrisprNuclease")){
            stop("Provided nuclease must be a 'CrisprNuclease' object. ")
        }
    }
    return(crisprNuclease)
}


.getDefaultCrisprNuclease <- function(type=c("Cas9", "Cas12a")){
    type <- match.arg(type)
    if (type=="Cas9"){
        data("SpCas9",
             package="crisprBase",
             envir=environment())
        nuc <- SpCas9
    } else {
        data("AsCas12a",
             package="crisprBase",
             envir=environment())
        nuc <- AsCas12a
    }
    return(nuc)
}


#' @importFrom methods is
.validateBSGenome <- function(genome){
    if (is.null(genome)){
        bsgenome <- .getDefaultBSGenome()
    } else {
        if (!is(bsgenome, "BSgenome")){
            stop("Provided genome must be a 'BSgenome' object. ")
        }
    }
    return(bsgenome)
}

.getDefaultBSGenome <- function(){
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38")){
        x <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38 
    } else {
        x <- NULL
    }
    return(x)
}



# Takes a character vector of sequences
# and writes to disk the sequences 
# in a fastq format. If temporary=TRUE,
# the fastq file is written in a 
# temporary folder. 
#' @importFrom utils write.table
.fastqfy <- function(sequences,
                     temporary=TRUE,
                     file=NULL
){
    lines <- list()
    lines[[1]] <- paste0("@", sequences)
    lines[[2]] <- sequences
    lines[[3]] <- paste0("+", sequences)
    lines[[4]] <- vapply(nchar(sequences), function(x){
                      paste0(rep('~', x), collapse='')
                  }, FUN.VALUE="a")    
    temp <- split(do.call(cbind, lines),
                  f=sequences)
    temp <- matrix(unlist(temp), ncol=1)
    if (temporary){
        file <- tempfile()
    } else {
        if (is.null(file)){
            stop("If temporary=FALSE, 'file' must be provided.")
        }
    }
    write.table(temp, 
                file=file,
                quote=FALSE,
                row.names=FALSE,
                col.names=FALSE)
    return(file)
}





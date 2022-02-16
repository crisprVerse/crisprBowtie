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

#' @importFrom utils data
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
.validateBSGenome <- function(bsgenome){
    if (is.null(bsgenome)){
        stop("Provided bsgenome argument cannot be NULL.")
    } else {
        if (!is(bsgenome, "BSgenome")){
            stop("Provided bsgenome argument must be a 'BSgenome' object. ")
        }
    }
    return(bsgenome)
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




# Takes a character vector of sequences
# and writes to disk the sequences 
# in a fasta format. If temporary=TRUE,
# the fasta file is written in a 
# temporary folder. 
#' @importFrom utils write.table
.fastafy <- function(sequences,
                     temporary=TRUE,
                     file=NULL
){
    lines <- list()
    lines[[1]] <- paste0(">", sequences)
    lines[[2]] <- sequences
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





.validateBowtieIndex <- function(bowtie_index){
    suffixes_fwd <- paste0(".", seq_len(4), ".ebwt")
    suffixes_rev <- paste0(".rev.", seq_len(2), ".ebwt")
    suffixes <- c(suffixes_fwd, suffixes_rev)
    files   <- paste0(bowtie_index, suffixes)
    missing <- files[!file.exists(files)]
    if (length(missing)>0){
        missing <- paste0(missing, collapse="\n") 
        stop("The following files needed for bowtie",
             " are missing: \n", missing)
    }
    return(bowtie_index)
}




#' @importFrom methods is
.checkBSGenomeOrNull <- function(bsgenome){
    if (!is.null(bsgenome)){
        if (!is(bsgenome, "BSgenome")){
            stop("Provided bsgenome argument must be a 'BSgenome' object or NULL. ")
        }
    }
    invisible(NULL)
}




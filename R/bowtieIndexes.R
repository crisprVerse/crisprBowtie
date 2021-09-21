#' @title List available bowtie indexes
#' @description List all available bowtie indexes found in
#'     the directory returned by \code{getBowtieIndexesDir}.
#' 
#' @return Character vector ofall available bowtie indexes found
#' in the directory returned by \code{getBowtieIndexesDir}.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' listBowtieIndexes()
#' 
#' @export
listBowtieIndexes <- function(){
    dir <- getBowtieIndexesDir()
    if (is.null(dir)){
        message("A Bowtie indexes path does not seem to be set.",
                "Consider setting setBowtieIndexesDir.")
        return(NULL)
    } else {
        out <- list.files(dir,
                          pattern="ebwt",
                          recursive=TRUE)
        out <- gsub("\\.rev","", out)
        out <- gsub("\\.[0-4]\\.ebwt|","", out)
        out <- unique(out)
    }
    if (length(out)==0){
        message("A Bowtie indexes path is found, but there does ",
                "not seem to be any valid indexes in it.")
        out <- NULL
    }
    return(out)
}




#' @title Get bowtie indexes directory  
#' @description Return bowtie indexes directory name if set up by user. 
#' 
#' @return If the environment path `BOWTIE_INDEXES` is found,
#'     \code{getBowtieDir} will return the directory name.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' getBowtieIndexesDir()
#' 
#' @export
getBowtieIndexesDir <- function(){
    dir <- Sys.getenv("BOWTIE_INDEXES")
    if (dir==""){
        dir <- NULL
    } else if (!file.exists(dir)){
        warning("A BOWTIE_INDEXES variable was found (`", dir, "`), 
            but the directory does not seem to exist.", 
            "See setBowtieIndexesDir() to set up the bowtie indexes",
            "directory manually.")
        dir <- NULL
    }
    return(dir)
}


#' @title Set bowtie indexes directory 
#' @description Set bowtie indexes directory in an R session. 
#' 
#' @param dir String specifying directory where bowtie indexes are stored. 
#' @return No object is returned. Instead, the BOWTIE_INDEXES
#'     environment variable is set to 'dir' for the remaining of the session.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' setBowtieIndexesDir("./")
#' 
#' @export
setBowtieIndexesDir <- function(dir){
    if (!dir.exists(dir)){
        stop("The provided path does not seem to exist on disk.")
    } else {
        Sys.setenv(BOWTIE_INDEXES=dir)
        message("The Bowtie indexes path is set up successfully.")
    }
}










.validateBowtieIndex <- function(bowtie_index,
                                 verbose=TRUE){
    if (is.null(bowtie_index)){
        bowtie_index <- .getDefaultBowtieIndex(verbose=verbose)
    }
    bowtie_index <- .appendBowtieDir(bowtie_index)
    bowtie_index <- .checkBowtieIndexFiles(bowtie_index)
    if (is.null(bowtie_index)){
        stop("Default bowtie_index could not be found and",
             " an explicit bowtie_index must be therefore provided.",
             " See getBowtieIndexesDir()",
             " and documentation to set up a default",
             " bowtie indexes directory.")
    } 
    return(bowtie_index)
}


.appendBowtieDir <- function(bowtie_index){
    if (!is.null(bowtie_index) & !dir.exists(dirname(bowtie_index))){
        bowtie_index <- file.path(getBowtieIndexesDir(),
                                  bowtie_index)
    }
    return(bowtie_index)
}



.checkBowtieIndexFiles <- function(bowtie_index){

    if (!is.null(bowtie_index)){
        if (!file.exists(dirname(bowtie_index))){
            stop("Provided bowtie_index could not be found.")
        } 
        suffixes_fwd <- paste0(".", seq_len(4), ".ebwt")
        suffixes_rev <- paste0(".rev.", seq_len(2), ".ebwt")
        suffixes <- c(suffixes_fwd, suffixes_rev)
        files   <- paste0(bowtie_index, suffixes)
        missing <- files[!file.exists(files)]
        if (length(missing)>0){
            missing <- paste0(missing, collapse="\n") 
            stop("bowtie_index found, but the following files",
                 " are missing in the directory: \n", missing)
        }
        
    }
    return(bowtie_index)
}

.getDefaultBowtieIndex <- function(first=FALSE,
                                   verbose=TRUE){
    choices <- listBowtieIndexes()
    if (is.null(choices)){
        return(NULL)
    }
    if (first){
        return(choices[[1]])
    }
    if ("hg38/hg38" %in% choices){
        index="hg38/hg38"
    } else if ("GRCh38/GRCh38" %in% choices){
        index="GRCh38/GRCh38"
    } else if (grepl("human", choices, ignore.case=TRUE)){
        wh <- which(grepl("human", choices, ignore.case=TRUE))
        index=choices[wh]
    } else if (grepl("sapiens", choices, ignore.case=TRUE)){
        wh <- which(grepl("sapiens", choices, ignore.case=TRUE))
        index=choices[wh]
    } else {
        index <- NULL
    }
    if (verbose){
        cat(paste0("[runBowtie] Using index ", index, "\n"))
    }
    return(index)
}



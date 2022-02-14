#' @title Perform short sequence alignment with bowtie 
#' @description Perform short sequence alignment with bowtie. 
#' 
#' @param sequences Character vector of DNA sequences.
#' @param bowtie_index String specifying path to a bowtie index.
#' @param n_mismatches Integer between 0 and 3 specifying maximum number
#'     of mismatches allowed between query sequences and target DNA.
#'     0 by default. 
#' @param all_alignments Should all possible alignments be returned?
#'     TRUE by default. 
#' @param n_max_alignments Maximum number of alignments to return if
#'     \code{all_alignments} is FALSE. 1000 by default. 
#' @param verbose Should messages be printed to the console?
#'     TRUE by default.
#' 
#' @return A data.frame of the alignments with the following columns:
#'    \itemize{
#'        \item \code{query} — string specifying query DNA sequence
#'        \item \code{target} — string specifying target DNA sequence
#'        \item \code{chr} - string specifying chromosome name
#'        \item \code{pos} - string specifying genomic coordinate of the start
#'              of the target DNA sequence
#'        \item \code{strand} - string specifying strand ("+" or "-") 
#'        \item \code{n_mismatches} - integer specifying number of mismatches
#'              between query and target sequences
#'        \item \code{mm1} - position of the first mismatch. NA if none.
#'        \item \code{mm2} - position of the second mismatch. NA if none.
#'        \item \code{mm3} - position of the third mismatch. NA if none.
#'        \item \code{mmnuc1} - nucleotide in target space of the first mismatch.
#'        \item \code{mmnuc2} - nucleotide in target space of the second mismatch.
#'        \item \code{mmnuc3} - nucleotide in target space of the third mismatch.
#'    }
#'     
#' 
#' @details \code{runBowtie} can be used to map short DNA sequences 
#'     to a reference genome. To search for sequences while imposing
#'     constraints on PAM sequences (such as gRNA spacer sequences), see
#'     \code{runCrisprBowtie} instead.  
#' 
#' @examples
#' fasta <- system.file(package="crisprBowtie", "example/chr1.fa")
#' outdir <- tempdir()
#' Rbowtie::bowtie_build(fasta,outdir=outdir, force=TRUE, prefix="tempIndex")
#' index <- file.path(outdir, "tempIndex")
#' seqs <- c("GGAAGT",
#'           "GTGGAC",
#'           "GTGTGC") 
#' results <- runBowtie(seqs, bowtie_index=index, n_mismatches=2)
#' 
#' 
#' @seealso \code{\link{runCrisprBowtie}} to map gRNA spacer sequences.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @importFrom Biostrings DNAStringSet DNAString replaceLetterAt
#' @importFrom Rbowtie bowtie
#' @importFrom stringr str_split
#' @export
runBowtie <- function(sequences, 
                      bowtie_index,
                      n_mismatches=0,
                      all_alignments=TRUE,
                      n_max_alignments=1000,
                      verbose=TRUE
){    
    .checkNMismatches(n_mismatches)
    # Note: bowtie is based on python,
    # so index is 0-based instead of 1-based
    sequences <- unique(sequences)
    if (is.null(bowtie_index)){
        stop("bowtie_index must be provided.")
    }
    bowtie_index <- .validateBowtieIndex(bowtie_index)
                                         
    
    # Generating alignments:
    input <- .fastqfy(sequences,
                      temporary=TRUE)
    results <- bowtie(sequences=input,
                      type="single",
                      index=bowtie_index,
                      v=n_mismatches,
                      a=all_alignments,
                      k=crisprBase:::.makeLongInteger(n_max_alignments), 
                      force=TRUE)

    if (file.exists(input)){
        file.remove(input)
    }
    # In case of no alignments:
    if (length(results)==0){ 
        return(NULL) 
    }  

    # Reformatting data:
    results <- str_split(results, pattern="\t")
    results <- do.call("rbind", results)
    results <- as.data.frame(results, stringsAsFactors=FALSE)
    colnames(results) <- c("query", "strand", "chr", "pos",
                           "target", "qc", "idk", "mismatches")
    cols    <- c("query", "chr", "pos", "strand", "mismatches")
    results <- results[,cols]
    results$query  <- as.character(results$query)
    results$chr    <- as.character(results$chr)
    results$strand <- as.character(results$strand)
    results$pos    <- as.numeric(results$pos)
    results$pos    <- results$pos+1 #since bowtie is 0-based

    # Processing mismatch info:
    results$mismatches <- as.character(results$mismatches)
    results$mismatches[results$mismatches==""]  <- "none"
    results$mismatches[is.na(results$mismatches)] <- "none"
    mismatches <- .processBowtieMismatches(results$mismatches)
    results    <- data.frame(results,
                             mismatches,
                             stringsAsFactors=FALSE)
    results$mismatches   <- NULL
    mm_cols <- paste0("mm", seq_len(N_MAX_MISMATCHES))
    results$n_mismatches <- rowSums(!is.na(results[, mm_cols]))



    # Complementing mismatches for negative strand:
    wh.neg <- which(results$strand=="-")
    if (length(wh.neg)>0){
        for (i in seq_len(N_MAX_MISMATCHES)){
            col <- paste0("mmnuc", i)
            .complement <- function(x){
                chartr("ATGC","TACG",x)
            }
            results[[col]][wh.neg] <- .complement(results[[col]][wh.neg])  
        }
    }

    # Getting target sequences in the absence of names:
    if (is.numeric(results$query)){
        results$query  <- sequences[as.integer(results$query)+1]
    }
    results$target <- results$query

    # Getting target genomic DNA sequence:
    sameLength <- length(unique(nchar(results$target)))==1
    if (sameLength){
        # Fast replacement with rectangular array:
        gRNA <- as.matrix(DNAStringSet(results$target))
        for (k in seq_len(N_MAX_MISMATCHES)){
            i <- seq_len(nrow(gRNA))
            j <- as.numeric(results[, paste0("mm",k)])
            x <- as.character(results[, paste0("mmnuc",k)])
            ind <- cbind(i,j)
            wh.complete <- complete.cases(ind)
            ind <- ind[wh.complete,,drop=FALSE]
            x <- x[wh.complete]
            if (length(x)>0){
                gRNA[ind] <- x
            }
        }
        DNA <- apply(gRNA, 1, paste, collapse="")
    } else {
        # Slow replacement with query sequences of different lengths:
        gRNAs <- lapply(results$target, DNAString)
        for (k in seq_len(N_MAX_MISMATCHES)){
            pos    <- as.numeric(results[, paste0("mm",k)])
            letter <- as.numeric(results[, paste0("mmnuc",k)])
            good <- which(!is.na(letter))
            if (length(good)>0){
                for (wh in good){
                    gRNAs[[wh]] <- replaceLetterAt(gRNAs[[wh]],
                                                   at=pos[wh],
                                                   letter=letter[wh])
                }
            }
        }
        DNA  <- as.character(DNAStringSet(gRNAs))
    }
    results$target <- DNA

    # Cleaning up:
    cols <- c("query", "target", "chr", "pos", "strand", "n_mismatches",
              paste0("mm", seq_len(N_MAX_MISMATCHES)),
              paste0("mmnuc", seq_len(N_MAX_MISMATCHES)))
    results <- results[, cols, drop=FALSE]
    return(results)
}


.checkNMismatches <- function(n_mismatches){
    if (n_mismatches>3){
        stop("Number of mismatches must be either 0, 1, 2 or 3.")
    }
}

.processBowtieMismatches <- function(mismatches){
    out <- data.frame(.processMismatchPositions(mismatches),
                      .processMismatchLetters(mismatches),
                      stringsAsFactors=FALSE,
                      check.names=FALSE)
    return(out)
}


#' @importFrom stringr str_extract
.processMismatchPositions <- function(mismatches){
    mismatches <- strsplit(mismatches, split=",") 
    
    # Getting mismatch positions:
    mismatches.pos <- lapply(mismatches, function(x){
        x <- as.numeric(str_extract(x, "[0-9]+"))
        x <- x + 1 #since bowtie are 0-based
        for (i in seq_len(N_MAX_MISMATCHES-1)){
            if (length(x)==i){
                x <- c(x, rep(NA, N_MAX_MISMATCHES-i))
                break
            }
        }
        return(x)
    })
    mismatches.pos     <- do.call("rbind", mismatches.pos)
    for (i in seq_len(N_MAX_MISMATCHES)){
        mismatches.pos[,i] <- as.numeric(mismatches.pos[,i])  
    }
    colnames(mismatches.pos) <- paste0("mm",
                                       seq_len(N_MAX_MISMATCHES))
    return(mismatches.pos)
}


#' @importFrom stringr str_extract
.processMismatchLetters <- function(mismatches){
    mismatches <- strsplit(mismatches, split=",") 
    
    # Getting mismatch positions:
    mismatches.letters <- lapply(mismatches, function(x){
        x <- str_extract(x, "\\:[A,C,G,T]") 
        x <- gsub("\\:","",x)
        x <- as.character(x)
        for (i in seq_len(N_MAX_MISMATCHES-1)){
            if (length(x)==i){
                x <- c(x, rep(NA, N_MAX_MISMATCHES-i))
                break
            }
        }
        return(x)
    })
    mismatches.letters <- do.call("rbind",
                                  mismatches.letters)
    colnames(mismatches.letters) <- paste0("mmnuc",
                                           seq_len(N_MAX_MISMATCHES))
    return(mismatches.letters)
}







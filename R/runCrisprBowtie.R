#' @title Find gRNA spacer alignments with bowtie
#' @description Return bowtie alignments for a list of gRNA spacer sequences.
#' 
#' @param spacers Character vector of DNA sequences corresponding
#'     to gRNA spacer sequences. Must all be of equal length. 
#' @param mode String specifying which alignment mode should be used:
#'     \code{protospacer} or \code{spacer}. In spacer mode, sequences are
#'     aligned to the genome without PAM (spacer only). In protospacer mode,
#'     sequences are aligned with all valid PAM sequences appended
#'     (spacer + PAM). Valid PAMs depend on the \code{canonical} and 
#'     \code{ignore_pam} user inputs. Spacer mode is recommended
#'     for alignments with a large number of potental PAM sequences
#'     such as non-canonical Cas9 PAM sequences, and especially
#'     non-canonical Cas12a PAM sequences. 
#' @param bsgenome BSgenome object to be used in spacer mode. 
#' @param bowtie_index Path to the bowtie index to be used for alignment.
#' @param crisprNuclease \code{CrisprNuclease} object. 
#' @param canonical Should only canonical PAM sequences be considered?
#'     TRUE by default.
#' @param ignore_pam If TRUE, will return all matches regardless of
#'     PAM sequence. FALSE by default.
#' @param n_mismatches Integer between 0 and 3 specifying maximum
#'     number of mismatches allowed between spacer and protospacer sequences.
#' @param all_alignments Should all possible alignments be returned?
#'     TRUE by default. 
#' @param n_max_alignments Maximum number of alignments to return if
#'     \code{all_alignments} is FALSE. 1000 by default. 
#' @param force_spacer_length Should the spacer length be overwritten in the
#'     crisprNuclease object? FALSE by default. 
#' @param verbose Should messages be printed to the consolde? TRUE by default.
#' @return \strong{runBowtie} returns spacer alignment data, including genomic 
#'     coordinates and sequence, and position of mismatches relative
#'     to \code{pam_site}.
#' 
#' @details \code{runCrisprBowtie} is similar to \code{runBowtie}, with the 
#'     addition of imposing constraints on PAM sequences such that the query
#'     sequences are valid protospacer sequences in the searched genome. 
#' 
#' @examples
#' fasta <- system.file(package="crisprBowtie", "example/chr1.fa")
#' outdir <- tempdir()
#' Rbowtie::bowtie_build(fasta,outdir=outdir, force=TRUE, prefix="tempIndex")
#' index <- file.path(outdir, "tempIndex")
#' seqs <- c("GGAAATTCCCCCAGTGGCGC",
#'           "ACACAGCTGCGGACAGGGCC")
#' data(SpCas9, package="crisprBase")
#' results <- runCrisprBowtie(seqs,
#'                            bowtie_index=index,
#'                            n_mismatches=2,
#'                            crisprNuclease=SpCas9)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
#' @importFrom BSgenome getSeq
#' @importFrom crisprBase extractPamFromProtospacer
#' @importFrom crisprBase extractSpacerFromProtospacer
#' @importFrom crisprBase nucleaseName pams pamLength pamIndices
#' @importFrom crisprBase spacerLength spacerLength<- spacerSide
runCrisprBowtie <- function(spacers, 
                            mode=c("protospacer", "spacer"),
                            bowtie_index=NULL,
                            bsgenome=NULL,
                            crisprNuclease=NULL,
                            canonical=TRUE,
                            ignore_pam=FALSE,
                            n_mismatches=0, 
                            all_alignments=TRUE,
                            n_max_alignments=1000,
                            force_spacer_length=FALSE,
                            verbose=TRUE
){ 
    #Checking inputs:
    mode     <- match.arg(mode)
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    if (mode=="spacer"){
        bsgenome <- .validateBSGenome(bsgenome)
    }
    

    # Checking mode arguments:
    if (mode=="spacer" & is.null(bsgenome)){
        stop("In 'spacer' mode, a BSgenome object must be provided using",
             " the argument 'bsgenome'.")
    }
    if (mode=="spacer" & verbose){
        cat(paste0("[runCrisprBowtie] Using ", bsgenome@pkgname, " \n"))
    }
    if (verbose){
        cat(paste0("[runCrisprBowtie] Searching for ",
                   nucleaseName(crisprNuclease), " protospacers \n"))
    }

    # Getting nuclease info:
    pam.len     <- pamLength(crisprNuclease)
    spacer.side <- spacerSide(crisprNuclease) 
    spacer.len  <- unique(nchar(spacers))

    if (length(spacer.len)>1){
        stop("Spacers must all have the same length.")
    }
    spacer.len.from.nuclease <- spacerLength(crisprNuclease)
    if (!force_spacer_length){
        if (spacer.len.from.nuclease != spacer.len){
            stop("Length of provided spacers is ",
                 spacer.len,
                 ", but spacer length for the provided nuclease is ",
                 spacer.len.from.nuclease, 
                 ". Consider force_spacer_length=TRUE to overwrite the",
                 " nuclease spacer length. ")
        }
    } else {
        if (spacer.len.from.nuclease != spacer.len){
            message("Setting spacer length to be ", spacer.len)
            spacerLength(crisprNuclease) <- spacer.len
        }
    }
    
    possiblePams      <- pams(crisprNuclease, 
                              primary=canonical,
                              ignore_pam=ignore_pam,
                              as.character=TRUE)
    pams.canonical    <- pams(crisprNuclease,
                              primary=TRUE,
                              as.character=TRUE)
    pams.noncanonical <- pams(crisprNuclease,
                              primary=FALSE,
                              as.character=TRUE)
    pam.indices <- pamIndices(crisprNuclease)

     

    # Getting input sequences:
    if (mode=="spacer"){
        sequences <- spacers
    } else {
        if (spacer.side=="5prime"){
            protospacers <- lapply(possiblePams, function(x){
                paste0(spacers, x)
            })
        } else {
            protospacers <- lapply(possiblePams, function(x){
                paste0(x, spacers)
            })
        }
        sequences <- do.call("c", protospacers)
    }

    # Performing alignment:
    aln <- runBowtie(sequences=sequences, 
                     bowtie_index=bowtie_index, 
                     n_mismatches=n_mismatches,
                     n_max_alignments=n_max_alignments,
                     all_alignments=all_alignments,
                     verbose=verbose)
    if (is.null(aln)){
        return(.emptyAlignments())
    }

    aln$pam_site <- .getPamSiteFromBowtieOutput(pos=aln$pos, 
                                                strand=aln$strand, 
                                                spacer.len=spacer.len,
                                                crisprNuclease=crisprNuclease,
                                                mode=mode)

    if (mode=="spacer"){
        if (nrow(aln)>0){
            pam.gr <- .getBowtiePamRanges(chr=aln$chr,
                                          pam_site=aln$pam_site,
                                          strand=aln$strand,
                                          crisprNuclease=crisprNuclease)
            aln$pam <- as.character(getSeq(bsgenome, pam.gr))
            aln$canonical <- aln$pam %in% pams.canonical
            if (canonical & !ignore_pam){
                aln <- aln[aln$canonical,,drop=FALSE]
            } else if (!canonical & !ignore_pam){
                aln  <- aln[aln$pam %in% pams.noncanonical,,drop=FALSE]
            }
        } 
    } else {
        bad <- which(aln$mm1 %in% pam.indices|
                     aln$mm2 %in% pam.indices| 
                     aln$mm3 %in% pam.indices)
        if (length(bad)>0){
            aln <- aln[-bad,,drop=FALSE]
        }
        if (nrow(aln)>0){
            aln$pam       <- extractPamFromProtospacer(aln$target,crisprNuclease)
            aln$target    <- extractSpacerFromProtospacer(aln$target, crisprNuclease)
            aln$query     <- extractSpacerFromProtospacer(aln$query, crisprNuclease)
            aln$canonical <- aln$pam %in% pams.canonical
        } 
    }


    if (nrow(aln)==0){
        aln <- .emptyAlignments()
    } else {
        #Changing relative position of mismatch
        if (spacer.side=="3prime" & mode=="protospacer"){
            wh <- paste0("mm", seq_len(N_MAX_MISMATCHES))
            aln[, wh] <- aln[, wh]-pam.len
        }
        aln <- aln[order(aln$query, aln$n_mismatches),,drop=FALSE]
        aln$pos <- NULL
        colnames(aln)[colnames(aln)=="query"]  <- "spacer"
        colnames(aln)[colnames(aln)=="target"] <- "protospacer"
        aln <- aln[ ,.outputColumns(), drop=FALSE]
        rownames(aln) <- NULL
    }
    return(aln)
}



.outputColumns <- function(){
     cols <- c("spacer",
               "protospacer",
               "pam",
               "chr",
               "pam_site",
               "strand", 
               "n_mismatches", 
               "mm1", "mm2", "mm3", 
               "mmnuc1", "mmnuc2", "mmnuc3", 
               "canonical")
    return(cols)
}


.emptyAlignments <- function(){
    cols <- .outputColumns()
    out <- data.frame(matrix(0, ncol=length(cols)))
    colnames(out) <- cols
    out <- out[-1,,drop=FALSE]
    return(out)
}


#' @importFrom crisprBase pamLength spacerSide
.getPamSiteFromBowtieOutput <- function(pos, 
                                        strand,
                                        spacer.len,
                                        crisprNuclease=NULL,
                                        mode=c("spacer", "protospacer")
){
    mode <- match.arg(mode)
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    pam.len  <- pamLength(crisprNuclease)
    spacer.side <- spacerSide(crisprNuclease)
    wh <- which(strand=="-")
    if (spacer.side=="5prime"){
        pam_site <- pos + spacer.len
        if (length(wh)>0){
            pam_site[wh] <- pos[wh] - 1
            if (mode=="protospacer"){
                pam_site[wh] <- pam_site[wh] + pam.len
            }
        }
    } else {
        pam_site <- pos - pam.len
        if (mode=="protospacer"){
            pam_site <- pos
        }
        if (length(wh)>0){
            pam_site[wh] <- pos[wh] +spacer.len + pam.len - 1
        }
    }
    return(pam_site)
}




#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges promoters
#' @importFrom crisprBase pamLength spacerSide
.getBowtiePamRanges <- function(chr,
                                pam_site,
                                strand,
                                crisprNuclease=NULL
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    pam.len  <- pamLength(crisprNuclease)
    spacer.side <- spacerSide(crisprNuclease)
    gr <- GRanges(chr, IRanges(pam_site, width=1), strand=strand)
    if (spacer.side=="5prime"){
        gr <- promoters(gr, downstream=pam.len, upstream=0)
    } else {
        gr <- promoters(gr, downstream=pam.len, upstream=0)
    }
    return(gr)
}




#' @importFrom S4Vectors split mcols<- mcols
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
.converToGRanges <- function(aln, collapse=TRUE){
    out <- GRanges(aln$chr,
                   IRanges(start=aln$pam_site,
                           width=1),
                   strand=aln$strand)
    mcols(out) <- aln
    if (!collapse){
        out <- split(out, f=mcols(out)[["spacer"]])
    }
    return(out)
}   






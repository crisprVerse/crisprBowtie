#' @title Perform CRISPR gRNA spacer alignment with bowtie
#' @description Perform CRISPR gRNA spacer alignment with bowtie.
#' 
#' @param spacers Character vector specifying gRNA spacer sequences.
#'     Sequences must all be of equal length. 
#' @param mode String specifying which alignment mode should be used:
#'     \code{protospacer} or \code{spacer}. 
#'     For RNA-targeting nucleases such as CasRx, only the 
#'     protospacer mode can be used. 
#' @param bowtie_index String specifying path to a bowtie index.
#' @param bsgenome A \linkS4class{BSgenome} object.
#'     Must be provided if \code{mode} is "spacer".
#'     Ignore
#' @param crisprNuclease A \linkS4class{CrisprNuclease} object. 
#' @param canonical Should only canonical PAM sequences be considered?
#'     TRUE by default.
#' @param ignore_pam Should PAM sequences be ignore?
#'     If TRUE, all alignments are returned regardless of PAM tolerance.
#'     FALSE by default.
#' @param n_mismatches Integer between 0 and 3 specifying maximum number
#'     of mismatches allowed between spacer sequences and target DNA.
#'     0 by default. 
#' @param all_alignments Should all possible alignments be returned?
#'     TRUE by default. 
#' @param n_max_alignments Maximum number of alignments to return if
#'     \code{all_alignments} is FALSE. 1000 by default. 
#' @param force_spacer_length Should the spacer length be overwritten in the
#'     \code{crisprNuclease} object? FALSE by default. 
#' @param rna_strict_directionality Should only protospacers found in the 
#'     original direction of the RNA be considered for RNA-targeting
#'     nucleases? TRUE by default.
#' @param verbose Should messages be printed to the console? TRUE by default.
#' 
#' @return A data.frame of the spacer alignments with the following columns:
#'    \itemize{
#'        \item \code{spacer} — string specifying gRNA spacer sequence
#'        \item \code{protospacer} — string specifying target protospacer sequence
#'        \item \code{pam} — string specifying target PAM sequence
#'        \item \code{chr} - string specifying chromosome name
#'        \item \code{pam_site} - string specifying genomic coordinate of the
#'              first nucleotide of the PAM sequence.
#'        \item \code{strand} - string specifying strand ("+" or "-") 
#'        \item \code{n_mismatches} - integer specifying number of mismatches
#'              between spacer and protospacer sequences
#'        \item \code{canonical} - logical indicating whether or not PAM sequence
#'              is canonical. 
#'    }
#'     
#' 
#' @details When \code{mode} is "spacer", spacer sequences are aligned to the
#'     reference index without appending PAM sequences first. This requires the
#'     specification of a \linkS4class{BSgenome} object through the argument
#'     \code{bsgenome} to validate that the aligned spacer sequences are
#'     adjacent to valid PAM sequences. 
#'     
#'     When \code{mode} is "protospacer", sequences are aligned with all
#'     valid PAM sequences appended (spacer + PAM). The set of valid PAM
#'     sequences depend on the inputs  \code{canonical} and \code{ignore_pam}.
#'     This is faster than the "spacer" mode if the number of possible
#'     PAM sequences is small (e.g. SpCas9).
#' 
#'     For RNA-targeting nucleases, such as RfxCas13d (CasRx), the bowtie
#'     index should be built on a transcriptome. For such applications,
#'     only the "protospacer" mode can be used as there is no
#'     corresponding bsgenome package. The protospacer sequences
#'     searched in the reference index will be the reverse complement
#'     of the input spacer sequences.
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
#' @seealso \code{\link{runBowtie}} to map general DNA sequences.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom BSgenome getSeq
#' @importFrom crisprBase extractPamFromTarget
#' @importFrom crisprBase extractProtospacerFromTarget
#' @importFrom crisprBase nucleaseName pams pamLength pamIndices
#' @importFrom crisprBase spacerLength spacerLength<- pamSide
#' @importFrom crisprBase hasSpacerGap isRnase
#' @importFrom crisprBase getTargetRanges
#' @importFrom Biostrings reverseComplement DNAStringSet
#' @importFrom BiocGenerics start end
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
                            rna_strict_directionality=TRUE,
                            verbose=TRUE
){ 
    #Checking inputs:
    mode     <- match.arg(mode)
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    if (hasSpacerGap(crisprNuclease)){
        stop("CRISPR nucleases with spacer gaps are not ",
             "supported at the moment.")
    }
    if (isRnase(crisprNuclease) & mode=="spacer"){
        stop("For RNA-targeting nucleases, the alignment mode must ",
             "be set to 'protospacer'.")
    }
    if (mode=="spacer"){
        bsgenome <- .validateBSGenome(bsgenome)
    }
    

    if (is.null(bowtie_index)){
        stop("bowtie_index must be provided.")
    }
    bowtie_index <- .validateBowtieIndex(bowtie_index)
    

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
    pam.side    <- pamSide(crisprNuclease) 
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
        if (isRnase(crisprNuclease)){
            sequences <- reverseComplement(DNAStringSet(spacers))
            sequences <- as.character(sequences)
        } else {
            sequences <- spacers
        }

        if (pam.side=="3prime"){
            protospacers <- lapply(possiblePams, function(x){
                paste0(sequences, x)
            })
        } else {
            protospacers <- lapply(possiblePams, function(x){
                paste0(x, sequences)
            })
        }
        sequences <- do.call("c", protospacers)
    }

    # Performing alignment:
    aln <- runBowtie(sequences=sequences, 
                     bowtie_index=bowtie_index, 
                     bsgenome=bsgenome,
                     n_mismatches=n_mismatches,
                     n_max_alignments=n_max_alignments,
                     all_alignments=all_alignments,
                     verbose=verbose)
    if (is.null(aln)){
        return(.emptyAlignments())
    }

    #cat("Getting pam site \n")
    aln$pam_site <- .getPamSiteFromBowtieOutput(pos=aln$pos, 
                                                strand=aln$strand, 
                                                spacer.len=spacer.len,
                                                crisprNuclease=crisprNuclease,
                                                mode=mode)
    aln <- aln[aln$pam_site>0,,drop=FALSE]

    #cat("Getting PAM sequences \n")
    if (mode=="spacer"){
        if (nrow(aln)>0){

            # Filtering out PAMs falling outside of chrs
            protoRanges <- getTargetRanges(seqnames=aln$chr,
                                           pam_site=aln$pam_site,
                                           strand=aln$strand,
                                           nuclease=crisprNuclease)
            chr_lens <- seqlengths(bsgenome)[as.character(seqnames(protoRanges))]
            valid <- BiocGenerics::start(protoRanges)>0 & BiocGenerics::end(protoRanges) <= chr_lens
            aln <- aln[valid,,drop=FALSE]
        }
        if (nrow(aln)>0){
            pam.gr <- .getBowtiePamRanges(chr=aln$chr,
                                          pam_site=aln$pam_site,
                                          strand=aln$strand,
                                          crisprNuclease=crisprNuclease)
            aln$pam <- as.character(getSeq(bsgenome, pam.gr))
            aln$canonical <- aln$pam %in% pams.canonical
        } 
    } else {
        if (nrow(aln)>0){
            aln$pam    <- extractPamFromTarget(aln$target,
                                                    crisprNuclease)
            aln$target <- extractProtospacerFromTarget(aln$target,
                                                       crisprNuclease)
            aln$query  <- extractProtospacerFromTarget(aln$query,
                                                       crisprNuclease)
            aln$canonical <- aln$pam %in% pams.canonical
        } 
        # To remove mismatches occuring in pam sequences:
        # Ordering by number of mismatches is important
        # as higher number of mismatches for identical 
        # spacer and protospacer indicates that the 
        # last-ordered row has a mismatch in the PAM
        aln <- aln[order(aln$query,
                         aln$target,
                         aln$chr,
                         aln$pam_site,
                         aln$n_mismatches),,drop=FALSE]

        ids <- paste0(aln$query, aln$target, aln$chr, aln$pam_site)
        bad <- which(duplicated(ids))
        #bad <- which(aln$mm1 %in% pam.indices|
        #             aln$mm2 %in% pam.indices| 
        #             aln$mm3 %in% pam.indices)
        if (length(bad)>0){
            aln <- aln[-bad,,drop=FALSE]
        }
    }
    if (canonical & !ignore_pam){
        aln <- aln[aln$canonical,,drop=FALSE]
    } else if (!canonical & !ignore_pam){
        aln  <- aln[aln$pam %in% pams.noncanonical,,drop=FALSE]
    }


    if (nrow(aln)==0){
        aln <- .emptyAlignments()
    } else {
        aln <- aln[order(aln$query, aln$n_mismatches),,drop=FALSE]
        aln$pos <- NULL
        colnames(aln)[colnames(aln)=="query"]  <- "spacer"
        colnames(aln)[colnames(aln)=="target"] <- "protospacer"
        aln <- aln[ ,.outputColumns(), drop=FALSE]
        rownames(aln) <- NULL
    }

    # RNAse considerations:
    if (isRnase(crisprNuclease)){
        spacers <- reverseComplement(DNAStringSet(aln$spacer))
        aln$spacer <- as.character(spacers)
        if (rna_strict_directionality){
            aln <- aln[aln$strand=="+",,drop=FALSE]
        }
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


#' @importFrom crisprBase pamLength pamSide
.getPamSiteFromBowtieOutput <- function(pos, 
                                        strand,
                                        spacer.len,
                                        crisprNuclease=NULL,
                                        mode=c("spacer", "protospacer")
){
    mode <- match.arg(mode)
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    pam.len  <- pamLength(crisprNuclease)
    pam.side <- pamSide(crisprNuclease)
    wh <- which(strand=="-")
    if (pam.side=="3prime"){
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
#' @importFrom crisprBase pamLength
.getBowtiePamRanges <- function(chr,
                                pam_site,
                                strand,
                                crisprNuclease=NULL
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    pam.len  <- pamLength(crisprNuclease)
    gr <- GRanges(chr, IRanges(pam_site, width=1), strand=strand)
    gr <- promoters(gr, downstream=pam.len, upstream=0)
    return(gr)
}









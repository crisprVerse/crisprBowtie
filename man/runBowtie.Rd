% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runBowtie.R
\name{runBowtie}
\alias{runBowtie}
\title{Perform short sequence alignment with bowtie}
\usage{
runBowtie(
  sequences,
  bowtie_index,
  bsgenome = NULL,
  n_mismatches = 0,
  all_alignments = TRUE,
  n_max_alignments = 1000,
  verbose = TRUE
)
}
\arguments{
\item{sequences}{Character vector of DNA sequences.}

\item{bowtie_index}{String specifying path to a bowtie index.}

\item{bsgenome}{BSgenome object.}

\item{n_mismatches}{Integer between 0 and 3 specifying maximum number
of mismatches allowed between query sequences and target DNA.
0 by default.}

\item{all_alignments}{Should all possible alignments be returned?
TRUE by default.}

\item{n_max_alignments}{Maximum number of alignments to return if
\code{all_alignments} is FALSE. 1000 by default.}

\item{verbose}{Should messages be printed to the console?
TRUE by default.}
}
\value{
A data.frame of the alignments with the following columns:
   \itemize{
       \item \code{query} — string specifying query DNA sequence
       \item \code{target} — string specifying target DNA sequence
       \item \code{chr} - string specifying chromosome name
       \item \code{pos} - string specifying genomic coordinate of the start
             of the target DNA sequence
       \item \code{strand} - string specifying strand ("+" or "-") 
       \item \code{n_mismatches} - integer specifying number of mismatches
             between query and target sequences
   }
}
\description{
Perform short sequence alignment with bowtie.
}
\details{
fasta <- system.file(package="crisprBowtie", "example/chr1.fa")
outdir <- tempdir()
Rbowtie::bowtie_build(fasta,outdir=outdir, force=TRUE, prefix="tempIndex")

\code{runBowtie} can be used to map short DNA sequences 
    to a reference genome. To search for sequences while imposing
    constraints on PAM sequences (such as gRNA spacer sequences), see
    \code{runCrisprBowtie} instead.
}
\examples{
fasta <- system.file(package="crisprBowtie", "example/chr1.fa")
outdir <- tempdir()
Rbowtie::bowtie_build(fasta,outdir=outdir, force=TRUE, prefix="tempIndex")
index <- file.path(outdir, "tempIndex")
seqs <- c("GGAAGT",
          "GTGGAC",
          "GTGTGC") 
results <- runBowtie(seqs, bowtie_index=index, n_mismatches=2)


}
\seealso{
\code{\link{runCrisprBowtie}} to map gRNA spacer sequences.
}
\author{
Jean-Philippe Fortin
}

library(crisprBowtie)
library(crisprBase)
library(Biostrings)
data(CasRx, package="crisprBase")
crisprNuclease <- CasRx

spacers <- c("GCCTTGAGAGCACCTCCTACACG",
             "GGAGAATTCCCCTACTCCTGCAT",
             "CTCCTTCTTCCTTTCCAGGCTTT",
             "CTGTCGTGCAGTCTCTTCCACAT")

bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/ensembl_human_104/ensembl_human_104"
aln <- runCrisprBowtie(spacers,
                       crisprNuclease=crisprNuclease,
                       bowtie_index=bowtie_index,
                       n_mismatches=2,
                       canonical=FALSE)

aln

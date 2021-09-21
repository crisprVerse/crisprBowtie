library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
bsgenome=BSgenome.Hsapiens.UCSC.hg38
seq <- DNAStringSet(getSeq(bsgenome, "chr1", 1000000,1001000-1))
names(seq) <- "chr1"
writeXStringSet(seq, "chr1.fa")

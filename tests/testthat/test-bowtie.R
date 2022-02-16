context('Testing runBowtie function')

# Creating bowtie index:
fasta  <- system.file(package="crisprBowtie","example/chr1.fa")
outdir <- tempdir()
Rbowtie::bowtie_build(fasta,
                      outdir=outdir,
                      force=TRUE,
                      prefix="tempIndex")
index <- file.path(outdir, "tempIndex")


test_that('Testing runBowtie', {
    seqs <- c("GGAAGTA", "GTGGACC", "GTGTGCG") 
    long_seq <- "AGGGAGGAGAGAGAGAGAGGGAGGTGGGTTAGGATTTAGGATTAGTTA"
    results_bowtie_0 <- runBowtie(seqs, bowtie_index=index, n_mismatches=0)
    results_bowtie_1 <- runBowtie(seqs, bowtie_index=index, n_mismatches=1)
    results_bowtie_2 <- runBowtie(seqs, bowtie_index=index, n_mismatches=2)
    results_bowtie_3 <- runBowtie(seqs, bowtie_index=index, n_mismatches=3)
    results_null <- runBowtie(long_seq, bowtie_index=index, n_mismatches=0)
    update <- FALSE

    expect_equal_to_reference(results_bowtie_0,
                              file=file.path("objects/results_bowtie_0.rds"),
                              update=update)
    expect_equal_to_reference(results_bowtie_1,
                              file=file.path("objects/results_bowtie_1.rds"),
                              update=update)
    expect_equal_to_reference(results_bowtie_2,
                              file=file.path("objects/results_bowtie_2.rds"),
                              update=update)
    expect_equal_to_reference(results_bowtie_3,
                              file=file.path("objects/results_bowtie_3.rds"),
                              update=update)
    expect_true(is.null(results_null))
})




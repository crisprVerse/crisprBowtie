context('Testing runCrisprBowtie functions')

library(crisprBase)
library(Rbowtie)

# Creating bowtie index:
fasta  <- system.file(package="crisprBowtie", "example/chr1.fa")
outdir <- tempdir()
Rbowtie::bowtie_build(fasta,
                      outdir=outdir,
                      force=TRUE,
                      prefix="tempIndex")
index <- file.path(outdir, "tempIndex")

Cas9 <- SpCas9
Cas12a <- AsCas12a


spacers_cas9 <- c("GGAAATTCCCCCAGTGGCGC",
                  "GGAACTTCCCCCCGTGGCGC",
                  "GGAGCACGGAGGGCTGCAGA",
                  "GCGGGCGAGCGGCGAGCGCG",
                  "GCGGGCGAGCAGCGAGGGCG",
                  "CGCCCTCCCGCCTGCCGACC",
                  "GGGCCTCACGCCTGCCGACC")

spacers_cas12a <- c("AGGCCTCGCCGCCGCGGGGTGTG",
                    "CCACACTCGCCCCAGCCAATCGA",
                    "TGCAGCCCTCCGTGCTCCAGGCC",
                    "TGCAGCCCTCGGTGCTCCAGGCC",
                    "TGCAGCCCTCGGTACTCCAGGCC",
                    "TGCAGCCCTCGGTACTGCAGGCC")



test_that('Testing runCrisprBowtie with Cas9 canonical', {
    results_cas9_mm0 <- runCrisprBowtie(spacers_cas9,
                                        bowtie_index=index,
                                        n_mismatches=0,
                                        crisprNuclease=Cas9,
                                        canonical=TRUE)
    results_cas9_mm1 <- runCrisprBowtie(spacers_cas9,
                                        bowtie_index=index,
                                        n_mismatches=1,
                                        crisprNuclease=Cas9,
                                        canonical=TRUE)
    results_cas9_mm2 <- runCrisprBowtie(spacers_cas9,
                                        bowtie_index=index,
                                        n_mismatches=2,
                                        crisprNuclease=Cas9,
                                        canonical=TRUE)
    results_cas9_mm3 <- runCrisprBowtie(spacers_cas9,
                                        bowtie_index=index,
                                        n_mismatches=3,
                                        crisprNuclease=Cas9,
                                        canonical=TRUE)
    expect_equal_to_reference(results_cas9_mm0,
                              file=file.path("objects/results_cas9_mm0.rds"))
    expect_equal_to_reference(results_cas9_mm1,
                              file=file.path("objects/results_cas9_mm1.rds"))
    expect_equal_to_reference(results_cas9_mm2,
                              file=file.path("objects/results_cas9_mm2.rds"))
    expect_equal_to_reference(results_cas9_mm3,
                              file=file.path("objects/results_cas9_mm3.rds"))
})


test_that('Testing runCrisprBowtie with Cas9 non-canonical', {
    results_cas9_mm0_nc <- runCrisprBowtie(spacers_cas9,
                                        bowtie_index=index,
                                        n_mismatches=0,
                                        crisprNuclease=Cas9,
                                        canonical=FALSE)
    results_cas9_mm1_nc <- runCrisprBowtie(spacers_cas9,
                                        bowtie_index=index,
                                        n_mismatches=1,
                                        crisprNuclease=Cas9,
                                        canonical=FALSE)
    results_cas9_mm2_nc <- runCrisprBowtie(spacers_cas9,
                                        bowtie_index=index,
                                        n_mismatches=2,
                                        crisprNuclease=Cas9,
                                        canonical=FALSE)
    results_cas9_mm3_nc <- runCrisprBowtie(spacers_cas9,
                                        bowtie_index=index,
                                        n_mismatches=3,
                                        crisprNuclease=Cas9,
                                        canonical=FALSE)
    expect_equal_to_reference(results_cas9_mm0_nc,
                              file=file.path("objects/results_cas9_mm0_nc.rds"))
    expect_equal_to_reference(results_cas9_mm1_nc,
                              file=file.path("objects/results_cas9_mm1_nc.rds"))
    expect_equal_to_reference(results_cas9_mm2_nc,
                              file=file.path("objects/results_cas9_mm2_nc.rds"))
    expect_equal_to_reference(results_cas9_mm3_nc,
                              file=file.path("objects/results_cas9_mm3_nc.rds"))
})


spacers_cas9_short <- c("AATTCCCCCAGTGGCGC",
                        "ACTTCCCCCCGTGGCGC",
                        "GCACGGAGGGCTGCAGA",
                        "GGCGAGCGGCGAGCGCG",
                        "GGCGAGCAGCGAGGGCG",
                        "CCTCCCGCCTGCCGACC",
                        "CCTCACGCCTGCCGACC")


test_that('Testing runCrisprBowtie with Cas9 (short spacers)', {
    expect_error(runCrisprBowtie(spacers_cas9_short,
                                 bowtie_index=index,
                                 n_mismatches=0,
                                 crisprNuclease=Cas9,
                                 canonical=FALSE))
    results_short_cas9 <- runCrisprBowtie(spacers_cas9_short,
                                          bowtie_index=index,
                                          n_mismatches=3,
                                          crisprNuclease=Cas9,
                                          canonical=FALSE,
                                          force_spacer_length=TRUE)
    expect_equal_to_reference(results_short_cas9,
                              file=file.path("objects/results_short_cas9.rds"))
})



test_that('Testing runCrisprBowtie with Cas12a', {
                      
    results_cas12a_mm0 <- runCrisprBowtie(spacers_cas12a,
                                          bowtie_index=index,
                                          n_mismatches=0,
                                          crisprNuclease=Cas12a,
                                          canonical=TRUE)
    results_cas12a_mm1 <- runCrisprBowtie(spacers_cas12a,
                                          bowtie_index=index,
                                          n_mismatches=1,
                                          crisprNuclease=Cas12a,
                                          canonical=TRUE)
    results_cas12a_mm2 <- runCrisprBowtie(spacers_cas12a,
                                          bowtie_index=index,
                                          n_mismatches=2,
                                          crisprNuclease=Cas12a,
                                          canonical=TRUE)
    results_cas12a_mm3 <- runCrisprBowtie(spacers_cas12a,
                                          bowtie_index=index,
                                          n_mismatches=3,
                                          crisprNuclease=Cas12a,
                                          canonical=TRUE)
    expect_equal_to_reference(results_cas12a_mm0,
                              file=file.path("objects/results_cas12a_mm0.rds"))
    expect_equal_to_reference(results_cas12a_mm1,
                              file=file.path("objects/results_cas12a_mm1.rds"))
    expect_equal_to_reference(results_cas12a_mm2,
                              file=file.path("objects/results_cas12a_mm2.rds"))
    expect_equal_to_reference(results_cas12a_mm3,
                              file=file.path("objects/results_cas12a_mm3.rds"))
})


spacers_cas12a_short <- c("AGGCCTCGCCGCCGCGGGGT",
                          "CCACACTCGCCCCAGCCAAT",
                          "TGCAGCCCTCCGTGCTCCAG",
                          "TGCAGCCCTCGGTGCTCCAG",
                          "TGCAGCCCTCGGTACTCCAG",
                          "TGCAGCCCTCGGTACTGCAG")


test_that('Testing runCrisprBowtie with Cas12a (short spacers)', {
    expect_error(runCrisprBowtie(spacers_cas12a_short,
                                 bowtie_index=index,
                                 n_mismatches=0,
                                 crisprNuclease=Cas12a,
                                 canonical=FALSE))
    results_short_cas12a <- runCrisprBowtie(spacers_cas12a_short,
                                            bowtie_index=index,
                                            n_mismatches=3,
                                            crisprNuclease=Cas12a,
                                            canonical=FALSE,
                                            force_spacer_length=TRUE)
    expect_equal_to_reference(results_short_cas12a,
                              file=file.path("objects/results_short_cas12a.rds"))
})




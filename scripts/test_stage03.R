source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/pathway_sets.R"))
source(path("R/pathway_scoring.R"))
source(path("R/pathway_qc.R"))
run_stage("stage03_tests", {
  rejected <- function(code, decision) {
    tryCatch({force(code); FALSE}, stage03_error = function(e) identical(e$decision, decision))
  }
  bad <- setNames(data.frame(x = 1), "Early.pCR")
  stopifnot(rejected(assert_blinded(bad), "STOP_REVIEW"))
  bad <- setNames(data.frame(x = 1), "residual_disease")
  stopifnot(rejected(assert_blinded(bad), "STOP_REVIEW"))
  stopifnot(length(forbidden_fields(data.frame(patient_id = "test", arm = "test"))) == 0L)
  authority <- blinded_identity()
  scores <- read_tsv(path(primary_files[2]))
  for (a in names(axes)) scores[[a]] <- as.numeric(scores[[a]])
  bad <- scores
  bad$patient_id[match(authority$pairs$baseline_sample_id[1], bad$sample_id)] <- "wrong_patient"
  stopifnot(rejected(pair_pathway_scores(bad, authority$pairs), "STOP_REVIEW"))
  constant <- scores
  constant$apoptosis <- 0
  stopifnot(!score_qc(constant)$passed[6])
  nonfinite <- scores
  nonfinite$apoptosis[1] <- Inf
  stopifnot(!score_qc(nonfinite)$passed[6])
  expected <- pair_pathway_scores(scores, authority$pairs)
  bad_pairs <- authority$pairs
  bad_pairs$baseline_sample_id[1] <- "wrong_sample"
  stopifnot(rejected(link_pathway_composition(expected, bad_pairs), "MODIFY_REVIEW"))
  before <- sha256(path("metadata/pathways/confirmatory_gene_sets_v1.tsv"))
  freeze_panel()
  stopifnot(identical(before, sha256(path("metadata/pathways/confirmatory_gene_sets_v1.tsv"))))
  verify_stage02()
  tests <- c("case_insensitive_leakage", "residual_disease_leakage", "arm_allowed",
    "patient_mismatch_rejected", "constant_scores_rejected", "nonfinite_scores_rejected",
    "composition_sample_mismatch_rejected", "frozen_membership_reused", "stage01_02_preserved")
  write_tsv(data.frame(check = tests, passed = TRUE), path("results/stage03_pathways/stage03_regression_tests.tsv"))
})
cat("Checks passed\n")

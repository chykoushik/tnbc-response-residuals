source("D:/TNBC-ResponseResiduals/R/paths.R")
source(path("R/composition_checks.R"))
status <- list(decision = "PROCEED", paired_complete = 150L, mapping_failures = 0L, metadata_conflicts = 0L)
stopifnot(stage01_status_ok(status))
for (field in names(status)) {
  invalid <- status
  invalid[[field]] <- if (field == "decision") "MODIFY" else 1L
  stopifnot(!stage01_status_ok(invalid))
}
identity <- data.frame(sample_id = c("A_B", "A_D"), gsm = c("GSM1", "GSM2"),
  patient_id = "A", timepoint = c("Baseline", "D1C2"))
fractions <- matrix(c(0.4, 0.6, 0.3, 0.7), nrow = 2, byrow = TRUE, dimnames = list(identity$sample_id, c("X", "Y")))
stopifnot(all(composition_qc(fractions, identity, identity$sample_id)$passed))
stopifnot(!all(composition_qc(fractions[2:1, ], identity, identity$sample_id)$passed))
for (value in c(NA_real_, Inf, -0.1, 1.1)) {
  invalid <- fractions
  invalid[1, 1] <- value
  stopifnot(!all(composition_qc(invalid, identity, identity$sample_id)$passed))
}
invalid <- fractions
rownames(invalid) <- rep("A_B", 2)
stopifnot(!all(composition_qc(invalid, identity, identity$sample_id)$passed))
composition <- cbind(identity, cancer_epithelial = c(60, 30), stromal_fraction = c(10, 20), immune_fraction = c(20, 40))
pairs <- data.frame(patient_id = "A", baseline_sample_id = "A_B", d1c2_sample_id = "A_D")
changes <- pair_composition(composition, pairs)
stopifnot(changes$delta_cancer_epithelial == -30, changes$delta_stromal_fraction == 10,
  changes$delta_immune_fraction == 20, !any(grepl("pcr|arm|tnbctype", names(changes), ignore.case = TRUE)))
stopifnot(inherits(try(pair_composition(composition[1, ], pairs), silent = TRUE), "try-error"))
cat("Checks complete\n")

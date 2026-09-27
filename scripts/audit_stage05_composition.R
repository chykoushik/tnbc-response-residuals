source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/outcome_manifest.R"))
verify_stage04_05()
spec <- read_spec05()
check <- read_numeric05("results/stage05_primary/composition_constraint_check.tsv", c("patient_id", "passed"))
full <- read_numeric05("data/processed/composition/neotrip_instaprism_full.tsv.gz",
  c("sample_id", "gsm", "patient_id", "timepoint"))
pairs <- read_tsv(path("data/processed/composition/neotrip_pair_composition_change.tsv"))
i <- match(check$patient_id, pairs$patient_id)
b <- match(pairs$baseline_sample_id[i], full$sample_id)
d <- match(pairs$d1c2_sample_id[i], full$sample_id)
if (anyNA(c(i, b, d)) || anyDuplicated(pairs$patient_id) ||
    any(full$patient_id[b] != check$patient_id | full$patient_id[d] != check$patient_id) ||
    any(full$timepoint[b] != "Baseline" | full$timepoint[d] != "D1C2"))
  stop05("Composition audit identity failed")
normal_delta <- full$Normal_Epithelial[d] - full$Normal_Epithelial[b]
report <- data.frame(patient_id = check$patient_id, delta_three_components = check$delta_sum,
  delta_normal_epithelial = normal_delta, four_component_sum = check$delta_sum + normal_delta,
  tolerance = spec$composition_constraint_tolerance,
  four_component_balance = abs(check$delta_sum + normal_delta) <= spec$composition_constraint_tolerance)
write_tsv(report, path("results/stage05_primary/composition_definition_audit.tsv"))
jsonlite::write_json(list(stage = 5L, decision = "MODIFY_REVIEW", patients = nrow(report),
  three_component_failures = sum(abs(report$delta_three_components) > spec$composition_constraint_tolerance),
  fourth_component = "Normal_Epithelial",
  four_component_balance_failures = sum(!report$four_component_balance),
  conclusion = "The three requested aggregates exclude Normal_Epithelial. Their paired changes do not sum to zero; adding the omitted normal-epithelial change verifies the original complete composition. This is a constraint diagnosis only; no reparameterization or renormalization was applied.",
  source_units = "Stage 2 pair changes: percentage points; full cell-type fractions: 0-1",
  model_action = "No logistic models fitted; explicit methodology review required"),
  path("metadata/outcome/stage05_composition_review.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Composition definition audited\n")

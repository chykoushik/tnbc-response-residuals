source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/stage06_common.R"))
run06("18_freeze_stage06_spec", {
  stage05_gate06()
  spec <- list(version = "v1", freeze_timestamp = utc_now(), primary_stage5_version = "v2",
    primary_residual_version = "v1", master_seed = 319641L, permutation_seed = 319642L,
    random_gene_seed = 319643L, permutations = 10000L, random_gene_sets_per_axis = 100L,
    primary_robustness_variants = variants06,
    additional_analyses = c("arm_specific_residuals", "outcome_arm_interaction", "leave_one_out_influence",
      "permutation_test", "random_matched_gene_sets", "complete_case_check"),
    axes = axes06, ridge_alpha = 0, standardize = TRUE, intercept = TRUE,
    nlambda = 100L, lambda_min_ratio = 1e-4, outer_folds = 10L, inner_folds = 10L,
    arm_specific_inner_folds = 5L,
    arm_inner_seed_rule = "319641 + 100 * arm_index + outer_fold; CT=1, CT/A=2",
    random_matching = "20 type-7 quantile bins on eligible mean expression; outer limits infinite; cut right=TRUE; confirmatory genes excluded; exact within-bin counts without replacement",
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    null_baseline_adjustment = "Corresponding null baseline z, consistent with replacing only the target pathway",
    ssgsea = list(alpha = 0.25, normalize = TRUE, minSize = 1, maxSize = "Inf"),
    ilr_component_order = c("cancer_epithelial", "normal_epithelial", "immune", "stromal"),
    ilr_basis = "stats::contr.helmert(4), columns divided by their Euclidean norms",
    stability = "Same beta sign as primary AND beta within primary ordinary log-odds 95% CI",
    permutation_quantile_type = 7L, permutation_statistic = "Frozen directional residual beta",
    random_percentile = "100 * mean(abs(null beta) <= abs(primary beta))",
    influence_zero_tolerance = 1e-8, influence_gate = 0.25,
    GSVA_version = as.character(packageVersion("GSVA")),
    glmnet_version = as.character(packageVersion("glmnet")),
    scoring_documentation = "https://bioconductor.posit.co/packages/3.22/bioc/manuals/GSVA/man/GSVA.pdf",
    external_data_access = "PROHIBITED", exploratory_pathways = "PROHIBITED")
  freeze_yaml06(spec, spec06_file)
})
cat("Stage 6 frozen\n")

novelty_audit <- function() {
  output <- data.frame(analysis_component = character(), neotrip_2026_status = character(),
    paper3_status = character(), overlap_risk = character(), evidence_source = character(),
    evidence_location = character(), decision = character())
  base <- paste0(repo_url, "/blob/", repo_commit, "/")
  add <- function(component, location, interpretation, done = TRUE, intended = "Supporting analysis only; not a novelty claim") {
    file <- sub(":.*$", "", location)
    line <- sub("^[^:]+:", "", location)
    output <<- rbind(output, data.frame(analysis_component = component,
      neotrip_2026_status = if (done) paste0("ALREADY_PERFORMED: ", interpretation) else paste0("NOT_IDENTIFIED_IN_REVIEWED_SOURCES: ", interpretation),
      paper3_status = intended, overlap_risk = if (done) "HIGH for this component alone" else "No direct residual-construct overlap found; limited by inaccessible full text/supplements",
      evidence_source = paste0(base, file), evidence_location = paste0(file, ":", line),
      decision = if (done) "DO_NOT_CLAIM_AS_NOVEL" else "RETAIN_FOR_PROTOCOL_REVIEW"))
  }
  add("baseline pathway scoring", "scripts/01-data_download_and_singscore.R:34-97", "All-timepoint TPM matrix is scored with multiScore; baseline analyses consume these scores.")
  add("D1C2 pathway scoring", "scripts/01-data_download_and_singscore.R:34-97", "Same all-timepoint scoring pipeline; D1C2 analysis in scripts/05-Figure4.R:24-35.")
  add("baseline vs D1C2 pathway change", "scripts/10-SuppFigure3.R:21-44", "Paired patient limma models test on-treatment versus baseline contrasts.")
  add("pCR association of baseline signatures", "scripts/03-Figure2.R:25-36", "Baseline signature logistic models use pCR.num as response.")
  add("pCR association of D1C2 signatures", "scripts/05-Figure4.R:24-35", "On-treatment signature logistic models use pCR.num as response.")
  add("pCR association of delta signatures", "scripts/06-Figure5.R:27-56", "Observed D1C2-minus-baseline scores enter pCR logistic models; no expected-change subtraction.")
  add("Hallmark signatures", "data/geneSets_Supplementary_Table1.txt:28-77", "HM gene-set definitions are loaded and scored by scripts/01-data_download_and_singscore.R:61-97; article abstract Methods corroborates Hallmark analysis.")
  add("immune signatures", "data/geneSets_Supplementary_Table1.txt:5-22", "cTME definitions enter the shared scoring pipeline; article abstract Methods corroborates immune signatures.")
  add("ferroptosis signatures", "data/geneSets_Supplementary_Table1.txt:2-4;23-27;78-88", "Explicit ferroptosis annotations and iron/redox definitions enter the shared scoring pipeline.")
  add("immune deconvolution", "scripts/12-SuppFigure5.R:21-39;85-100", "BRCA InstaPrism fractions are generated and analysed by timepoint; article abstract Results corroborates deconvolution.")
  add("tumor epithelial fraction", "scripts/13-SuppFigure6.R:22-30", "Cancer epithelial theta is exposed as Cancer_Epithelial and Tumor_fraction; scripts/12-SuppFigure5.R:435-439 analyses the epithelial fraction.")
  add("tumor-purity/composition adjustment", "scripts/13-SuppFigure6.R:58-84;300-318;374-386", "Stromal/immune fractions adjust signature-outcome and longitudinal models; this is transcriptomic composition adjustment, not independent pathology purity.")
  add("treatment-arm interaction analyses", "scripts/00a-helper_functions.R:67-106", "Explicit Gene*Covariate_interaction logistic formula; callers specify Arm in scripts/03-Figure2.R:33-36 and scripts/06-Figure5.R:48-51.")
  add("XGBoost prediction", "scripts/07-Figure6.R:44-105;127-173", "Repeated cross-validation optimizes binary pCR classification using baseline/on-treatment signatures; outcome is label=pCR.num.")
  add("prediction of expected treatment-induced change from pretreatment state without pCR",
    "scripts/07-Figure6.R:22-44;76-105", "The prediction target is pCR, not outcome-blind molecular change.", FALSE,
    "Intended central construct; not implemented in Stage 2")
  add("cross-fitted expected molecular response", "scripts/07-Figure6.R:76-105;127-173",
    "Cross-validation evaluates pCR classifiers, not held-out expected molecular response models.", FALSE, "Intended central construct; not implemented")
  add("patient-specific expected-versus-observed residual", "scripts/06-Figure5.R:27-50",
    "deltaMat is an observed timepoint difference with no expected-change prediction. scripts/02-Figure1.R:64-92 uses variance-partition residual variance, a different object.", FALSE, "Intended observed-minus-expected molecular residual; not implemented")
  add("out-of-fold residual phenotype", "scripts/07-Figure6.R:127-181",
    "CV produces AUC, accuracy and classifier importance; no per-patient out-of-fold molecular residual phenotype.", FALSE, "Intended patient-level cross-fitted residual; not implemented")
  add("pCR analysis after residual freezing", "scripts/07-Figure6.R:76-105",
    "pCR labels enter model optimization directly; no frozen outcome-blind residual artifact in reviewed pipeline.", FALSE, "Intended later analysis after residual freezing; prohibited in Stage 2")
  add("composition-adjusted residual-pCR association", "scripts/13-SuppFigure6.R:58-84",
    "Adjustment is applied to observed signature associations; no expected-change residual is constructed.", FALSE, "Intended later residual-specific association; prohibited in Stage 2")
  add("residual-specific negative controls", "scripts/07-Figure6.R:127-181",
    "Repeated supervised classifier CV is not a residual null/permutation analysis; complete tracked script review found no residual-specific control pipeline.", FALSE, "Intended residual null/permutation controls; not implemented")
  add("external residual transportability", "scripts/01-data_download_and_singscore.R:14-18;34-97",
    "Pipeline uses NeoTRIP GSE319641 and saves within-cohort scores; no frozen residual model or external-cohort residual transport pipeline identified.", FALSE, "Intended transportability of frozen residual construct; not implemented")
  write_tsv(output, path("metadata/novelty/neotrip_already_done_vs_paper3.tsv"))
  scripts <- git_output("external/neotrip_2026_repo", c("ls-files", "scripts"))
  coverage <- data.frame(source = scripts,
    sha256 = vapply(path("external/neotrip_2026_repo", scripts), sha256, ""),
    line_count = vapply(path("external/neotrip_2026_repo", scripts), function(x) length(readLines(x, warn = FALSE)), integer(1)),
    review = "Static read-only review of inputs, transformations, model targets, predictions and exported objects; not executed.")
  write_tsv(coverage, path("metadata/novelty/review_coverage.tsv"))
  write_tsv(data.frame(item = c("scope", "semantic_review", "access_limitation", "novelty_limit", "residual_method_overlap_found"),
    value = c("Pinned complete public repository; public abstract and bibliography; restricted dataset descriptions.",
      "Distinguished observed delta, variance-partition residual variance and pCR-classifier cross-validation from outcome-blind expected-change residuals.",
      "No publisher full article or supplement body retrieved; no claim that uninspected sources lack the method.",
      "Bounded source audit, not a proof of universal novelty or a manuscript claim.", "FALSE")),
    path("metadata/novelty/review_scope.tsv"))
  output
}

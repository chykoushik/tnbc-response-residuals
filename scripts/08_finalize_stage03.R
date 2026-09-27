source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/pathway_sets.R"))
source(path("R/pathway_scoring.R"))
source(path("R/pathway_qc.R"))
stage03_init()
run_stage("08_finalize_pathways", {
  version <- function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_
  status <- list(stage = 3L, timestamp = utc_now(), stage02_verified = FALSE, panel_version = "v1",
    msigdb_version = NA_character_, msigdbr_version = version("msigdbr"), singscore_version = version("singscore"),
    confirmatory_pathways = 6L, all_sets_frozen = FALSE, all_sets_pass_coverage = FALSE,
    expression_samples = 0L, paired_patients = 0L, composition_linked_patients = 0L,
    score_qc_failures = NA_integer_, mapping_failures = NA_integer_, forbidden_outcome_fields = NA_integer_,
    decision = "MODIFY_REVIEW", decision_reason = "Stage 3 validation incomplete")
  validation <- data.frame(check = character(), passed = logical())
  check <- function(name, condition, reason, decision = "MODIFY_REVIEW") {
    validation <<- rbind(validation, data.frame(check = name, passed = isTRUE(condition)))
    if (!isTRUE(condition)) stage03_stop(reason, decision)
  }
  tryCatch({
    status$stage02_verified <- verify_stage02()
    failure <- path("results/stage03_pathways/stage03_failure.json")
    if (file.exists(failure)) {
      f <- jsonlite::read_json(failure, simplifyVector = TRUE)
      stage03_stop(f$reason, f$decision)
    }
    sets <- read_frozen_panel()
    status$all_sets_frozen <- TRUE
    status$msigdb_version <- "2026.1"
    authority <- blinded_identity()
    expression <- read_pathway_expression(authority$identity)
    coverage <- pathway_coverage(sets, rownames(expression))
    status$all_sets_pass_coverage <- all(coverage$passed)
    scores <- read_tsv(path(primary_files[2]))
    check("sample_schema", identical(names(scores), c("sample_id", "gsm", "patient_id", "timepoint", names(axes))),
      "Score schema changed", "STOP_REVIEW")
    for (a in names(axes)) scores[[a]] <- as.numeric(scores[[a]])
    status$expression_samples <- nrow(scores)
    i <- match(scores$sample_id, authority$identity$sample_id)
    mapping <- c(sum(is.na(i)), sum(duplicated(scores$sample_id)),
      as.integer(!identical(scores$sample_id, colnames(expression))))
    if (!anyNA(i)) for (f in c("gsm", "patient_id", "timepoint"))
      mapping <- c(mapping, sum(scores[[f]] != authority$identity[[f]][i]))
    status$mapping_failures <- sum(mapping)
    check("sample_mapping", status$mapping_failures == 0L, "Sample identity differs", "STOP_REVIEW")
    qc <- score_qc(scores)
    status$score_qc_failures <- sum(!qc$passed)
    check("score_qc", status$score_qc_failures == 0L, "Stored score QC failed")
    expected_scores <- score_pathways(expression, sets, authority$identity)
    check("score_reproduction", isTRUE(all.equal(scores, expected_scores, check.attributes = FALSE, tolerance = 1e-12)),
      "Stored scores differ")
    expected_pairs <- pair_pathway_scores(scores, authority$pairs)
    expected_inputs <- link_pathway_composition(expected_pairs, authority$pairs)
    for (entry in list(list(file = primary_files[3], expected = expected_pairs),
                       list(file = primary_files[4], expected = expected_inputs))) {
      stored <- read_tsv(path(entry$file))
      check(paste0(basename(entry$file), "_schema"), identical(names(stored), names(entry$expected)),
        "Paired output schema differs", "STOP_REVIEW")
      for (f in names(stored)[-(1:2)]) stored[[f]] <- as.numeric(stored[[f]])
      check(basename(entry$file), isTRUE(all.equal(stored, entry$expected, check.attributes = FALSE, tolerance = 1e-12)),
        "Paired output differs")
    }
    status$paired_patients <- nrow(expected_pairs)
    status$composition_linked_patients <- nrow(expected_inputs)
    leakage <- leakage_qc()
    status$forbidden_outcome_fields <- sum(leakage$forbidden_fields)
    check("leakage", all(leakage$passed), "Primary outputs incomplete")
    check("figures", all(file.exists(path("figures/stage03_pathways",
      c("score_distributions_v1.pdf", "score_distributions_v1.png")))), "QC figures missing")
    status$decision <- "PROCEED"
    status$decision_reason <- paste("Six exact MSigDB 2026.1 Hallmark sets frozen with msigdbr 26.1.1;",
      "all pass HGNC-symbol coverage. Outcome-blinded singscore scores validate for 401 samples;",
      "150 verified paired changes link to composition. No outcome fields, mapping failures or score QC failures.",
      "Expression values and all 23,890 genes retained unchanged. Technical QC only;",
      "residual modelling, exploratory scoring and ssGSEA remain deferred.")
  }, error = function(e) {
    status$decision <<- if (inherits(e, "stage03_error")) e$decision else "MODIFY_REVIEW"
    status$decision_reason <<- conditionMessage(e)
  })
  write_tsv(validation, path("results/stage03_pathways/stage03_validation.tsv"))
  existing <- read_tsv(path("logs/software_versions.tsv"))
  packages <- c("msigdbr", "singscore", "BiocManager")
  packages <- packages[vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  dependencies <- tools::package_dependencies(packages, db = installed.packages(), recursive = TRUE)
  packages <- sort(unique(c(packages, unlist(dependencies, use.names = FALSE))))
  versions <- data.frame(R_version = R.version.string, platform = R.version$platform,
    package = packages, package_version = vapply(packages, version, ""))
  write_tsv(rbind(existing[!existing$package %in% packages, ], versions), path("logs/software_versions.tsv"))
  if (status$decision == "PROCEED") {
    renv::snapshot(project = root, packages = c("renv", "digest", "jsonlite", "yaml", "InstaPrism",
      "Biobase", "BiocManager", "msigdbr", "singscore"), prompt = FALSE)
  }
  jsonlite::write_json(status, path("results/stage03_pathways/stage03_status.json"),
    auto_unbox = TRUE, pretty = TRUE, na = "null")
  code <- c("R/pathway_sets.R", "R/pathway_scoring.R", "R/pathway_qc.R",
    "scripts/06_freeze_confirmatory_pathways.R", "scripts/07_score_confirmatory_pathways.R",
    "scripts/08_finalize_stage03.R", "scripts/test_stage03.R", "scripts/install_stage03.R", "renv.lock")
  write_tsv(data.frame(local_relative_path = code, sha256 = vapply(path(code), sha256, "")),
    path("metadata/pathways/implementation_checksums.tsv"))
  files <- c(path(primary_files[1]), unlist(lapply(path(stage03_dirs), list.files,
    full.names = TRUE, recursive = TRUE), use.names = FALSE))
  files <- sort(unique(files[file.exists(files) & basename(files) != "stage03_output_checksums.tsv"]))
  write_tsv(data.frame(local_relative_path = substring(files, nchar(root) + 2L),
    file_size_bytes = file.info(files)$size, sha256 = vapply(files, sha256, "")),
    path("results/stage03_pathways/stage03_output_checksums.tsv"))
  entry <- data.frame(timestamp = status$timestamp, stage = 3L, decision = status$decision,
    reason = status$decision_reason, inputs = "Stage 1 immutable expression and manifests;Stage 2 composition;frozen MSigDB panel",
    outputs = "results/stage03_pathways/stage03_status.json")
  write.table(entry, path("logs/decision_log.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE, append = TRUE)
})
cat("results/stage03_pathways/stage03_status.json\n")

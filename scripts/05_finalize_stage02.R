source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/repository_audit.R"))
source(path("R/composition_checks.R"))
source(path("R/novelty_audit.R"))
run_stage("05_finalize", {
  authority <- verify_stage01()
  pin_repository("external/neotrip_2026_repo", repo_url, repo_commit)
  novelty <- novelty_audit()
  overlap <- any(novelty$decision == "STOP_REVIEW")
  feasibility <- jsonlite::read_json(path("results/stage02_composition/instaprism_feasibility.json"), simplifyVector = TRUE)
  pathology <- read_tsv(path("metadata/composition/pathology_variable_audit.tsv"))
  pathology_available <- any(pathology$decision == "AVAILABLE_PUBLICLY" & pathology$linkable == "TRUE")
  qc <- read_tsv(path("results/stage02_composition/composition_qc.tsv"))
  failures <- sum(as.integer(qc$failures), na.rm = TRUE)
  qc_pass <- nrow(qc) > 0L && all(qc$passed == "TRUE") && !any(missing_value(qc$failures))
  composition_file <- path("data/processed/composition/neotrip_composition.tsv")
  pair_file <- path("data/processed/composition/neotrip_pair_composition_change.tsv")
  reproduced <- isTRUE(feasibility$reproducible_from_public_inputs) && qc_pass &&
    file.exists(composition_file) && file.exists(pair_file)
  sample_count <- pair_count <- 0L
  if (reproduced) {
    composition <- read_tsv(composition_file)
    changes <- read_tsv(pair_file)
    forbidden <- c("pCR", "pCR.num", "Early.pCR", "Arm", "TNBCtype", "pcr", "arm")
    if (any(names(composition) %in% forbidden) || any(names(changes) %in% forbidden)) stop("Forbidden composition fields")
    index <- match(composition$sample_id, authority$identity$sample_id)
    if (anyNA(index) || anyDuplicated(composition$sample_id) ||
        !setequal(composition$sample_id, authority$identity$sample_id) ||
        any(composition$patient_id != authority$identity$patient_id[index]) ||
        any(composition$timepoint != authority$identity$timepoint[index]) ||
        any(composition$gsm != authority$identity$gsm[index])) stop("Composition identity changed")
    for (field in c("cancer_epithelial", "stromal_fraction", "immune_fraction")) {
      composition[[field]] <- as.numeric(composition[[field]])
      if (any(!is.finite(composition[[field]])) || any(composition[[field]] < 0 | composition[[field]] > 100)) stop("Composition values changed")
    }
    expected <- pair_composition(composition, authority$pairs)
    if (!identical(changes$patient_id, expected$patient_id) ||
        !identical(changes$baseline_sample_id, expected$baseline_sample_id) ||
        !identical(changes$d1c2_sample_id, expected$d1c2_sample_id)) stop("Composition pair order changed")
    for (field in names(expected)[-(1:3)]) {
      if (any(abs(as.numeric(changes[[field]]) - expected[[field]]) > 1e-10)) stop("Composition delta changed")
    }
    full_connection <- gzfile(path("data/processed/composition/neotrip_instaprism_full.tsv.gz"), "rt")
    full <- read.delim(full_connection, check.names = FALSE, stringsAsFactors = FALSE)
    close(full_connection)
    fractions <- as.matrix(full[, -(1:4)])
    rownames(fractions) <- full$sample_id
    expression_connection <- gzfile(path("data/raw/geo/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz"), "rt")
    expected_samples <- strsplit(readLines(expression_connection, n = 1L), "\t", fixed = TRUE)[[1L]][-(1:3)]
    close(expression_connection)
    recheck <- composition_qc(fractions, authority$identity, expected_samples)
    if (any(!recheck$passed)) stop("Stored composition QC failed")
    if (!identical(full$sample_id, composition$sample_id) ||
        any(abs(composition$cancer_epithelial - fractions[, "Cancer_Epithelial"] * 100) > 1e-10) ||
        any(abs(composition$stromal_fraction - rowSums(fractions[, 1:3]) * 100) > 1e-10) ||
        any(abs(composition$immune_fraction - rowSums(fractions[, 4:6]) * 100) > 1e-10)) stop("Composition aggregation changed")
    sample_count <- nrow(composition)
    pair_count <- nrow(changes)
  }
  if (overlap || (isTRUE(feasibility$reproducible_from_public_inputs) && !reproduced)) {
    decision <- "STOP_REVIEW"
    reason <- if (overlap) "Direct evidence of substantial prior residual-method implementation requires review." else "Composition reproduction or sample-level QC failed; no downstream analysis is permitted."
  } else if (!isTRUE(feasibility$reproducible_from_public_inputs)) {
    decision <- if (isTRUE(feasibility$reference_public) && isTRUE(feasibility$package_identified)) "MODIFY_REVIEW" else "STOP_REVIEW"
    reason <- feasibility$reason
  } else {
    decision <- if (pathology_available) "PROCEED" else "PROCEED_WITH_MODIFICATION"
    reason <- paste("Published InstaPrism method reproduced from pinned public inputs;", sample_count, "samples and", pair_count,
      "Stage 1 pairs pass QC. No qualifying public pathology tumour-content measurement is linkable in inspected sources; transcriptomic composition is the available measurement.",
      "Negative ComBat input values were retained exactly as in the published code; no undocumented clipping was applied.",
      "No substantial residual construct was identified in reviewed code and abstract. Publisher full text/supplements could not be retrieved; pathology availability and novelty conclusions remain limited to inspected sources. No pCR analysis performed.")
  }
  status <- list(stage = 2L, timestamp = utc_now(), stage01_verified = TRUE,
    pathology_variable_found = pathology_available, pathology_variable_linkable = pathology_available,
    instaprism_exact_method_identified = all(unlist(feasibility[c("package_identified", "package_version_identified",
      "reference_identified", "input_scale_identified", "generation_code_identified", "cell_types_identified", "aggregation_rules_identified")])),
    instaprism_public_inputs = isTRUE(feasibility$reference_public), instaprism_reproduced = reproduced,
    composition_samples = sample_count, paired_composition_patients = pair_count,
    composition_qc_failures = failures, residual_method_overlap_found = overlap, decision = decision, decision_reason = reason)
  packages <- c("renv", "digest", "jsonlite", "yaml", "InstaPrism", "Biobase", "BiocVersion")
  packages <- packages[vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  dependencies <- tools::package_dependencies(packages, db = installed.packages(), recursive = TRUE)
  packages <- sort(unique(c(packages, unlist(dependencies, use.names = FALSE))))
  versions <- data.frame(R_version = R.version.string, platform = R.version$platform, package = packages,
    package_version = vapply(packages, function(x) as.character(utils::packageVersion(x)), ""))
  if (file.exists(path("metadata/composition/auxiliary_software.tsv"))) {
    versions <- rbind(versions, read_tsv(path("metadata/composition/auxiliary_software.tsv")))
  }
  write_tsv(versions, path("logs/software_versions.tsv"))
  renv::snapshot(project = root, packages = c("renv", "digest", "jsonlite", "yaml", "InstaPrism", "Biobase"), prompt = FALSE)
  jsonlite::write_json(status, path("results/stage02_composition/stage02_status.json"), auto_unbox = TRUE, pretty = TRUE)
  outputs <- c(list.files(path("metadata/composition"), full.names = TRUE),
    list.files(path("metadata/novelty"), full.names = TRUE),
    list.files(path("results/stage02_composition"), full.names = TRUE),
    list.files(path("data/processed/composition"), full.names = TRUE))
  outputs <- outputs[!file.info(outputs)$isdir & basename(outputs) != "stage02_output_checksums.tsv"]
  write_tsv(data.frame(local_relative_path = substring(outputs, nchar(root) + 2L),
    file_size_bytes = file.info(outputs)$size, sha256 = vapply(outputs, sha256, "")),
    path("results/stage02_composition/stage02_output_checksums.tsv"))
  log <- data.frame(timestamp = status$timestamp, stage = 2L, decision = decision, reason = reason,
    inputs = "metadata/composition/stage01_input_checksums.tsv;metadata/data_inventory.tsv;metadata/composition/instaprism_method_audit.tsv",
    outputs = "results/stage02_composition/stage02_status.json;metadata/composition;metadata/novelty;data/processed/composition")
  write.table(log, path("logs/decision_log.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
})
cat(path("results/stage02_composition/stage02_status.json"), "\n", sep = "")

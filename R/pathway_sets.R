axes <- c(immune_interferon = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  proliferation = "HALLMARK_E2F_TARGETS",
  stromal_emt = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  dna_repair = "HALLMARK_DNA_REPAIR", hypoxia_stress = "HALLMARK_HYPOXIA",
  apoptosis = "HALLMARK_APOPTOSIS")
stage03_dirs <- c("metadata/pathways", "data/processed/pathways",
  "results/stage03_pathways", "figures/stage03_pathways")
stage03_init <- function() {
  for (x in stage03_dirs) dir.create(path(x), recursive = TRUE, showWarnings = FALSE)
}
stage03_stop <- function(reason, decision = "MODIFY_REVIEW") {
  stop(structure(list(message = reason, call = NULL, decision = decision),
    class = c("stage03_error", "error", "condition")))
}
verify_stage02 <- function() {
  s <- jsonlite::read_json(path("results/stage02_composition/stage02_status.json"), simplifyVector = TRUE)
  if (!identical(s$decision, "PROCEED_WITH_MODIFICATION") ||
      !isTRUE(s$stage01_verified) || !isTRUE(s$instaprism_reproduced) ||
      !isTRUE(s$composition_samples == 401) || !isTRUE(s$paired_composition_patients == 150) ||
      !isTRUE(s$composition_qc_failures == 0) || !identical(s$residual_method_overlap_found, FALSE))
    stage03_stop("Stage 2 prerequisite failed", "STOP_REVIEW")
  for (f in c("metadata/composition/stage01_input_checksums.tsv",
              "results/stage02_composition/stage02_output_checksums.tsv")) {
    x <- read_tsv(path(f))
    files <- if ("path" %in% names(x)) x$path else x$local_relative_path
    if (any(!file.exists(path(files))) || !identical(unname(vapply(path(files), sha256, "")), x$sha256))
      stage03_stop("Stage 1/2 provenance failed", "STOP_REVIEW")
  }
  tryCatch(verify_inventory(), error = function(e) stage03_stop("Source provenance failed", "STOP_REVIEW"))
  TRUE
}
read_restricted <- function(file, fields) {
  header <- names(read_tsv_header(file))
  if (!all(fields %in% header)) stage03_stop("Required identity fields missing", "STOP_REVIEW")
  classes <- ifelse(header %in% fields, "character", "NULL")
  x <- read.delim(file, sep = "\t", quote = "", comment.char = "", check.names = FALSE,
    colClasses = classes, na.strings = NULL)
  x[, fields, drop = FALSE]
}
read_tsv_header <- function(file) {
  read.delim(file, nrows = 0L, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
}
forbidden_fields <- function(x) {
  names(x)[grepl("pcr|early[.]pcr|early_pcr|response|outcome|residual_disease", names(x), ignore.case = TRUE)]
}
assert_blinded <- function(x) {
  if (length(forbidden_fields(x))) stage03_stop("Forbidden outcome fields", "STOP_REVIEW")
  invisible(x)
}
blinded_identity <- function() {
  fields <- c("patient_id", "baseline_gsm", "baseline_sample_id", "d1c2_gsm", "d1c2_sample_id",
    "arm", "baseline_present", "d1c2_present", "paired", "metadata_complete")
  manifest <- read_restricted(path("metadata/neotrip_pair_manifest_v1.tsv"), fields)
  assert_blinded(manifest)
  write_tsv(manifest, path("metadata/neotrip_analysis_manifest_blinded_v1.tsv"))
  identity <- read_restricted(path("metadata/neotrip_mapping_audit.tsv"),
    c("geo_library_name", "gsm", "geo_patient_id", "geo_timepoint", "mapping_status"))
  names(identity) <- c("sample_id", "gsm", "patient_id", "timepoint", "mapping_status")
  if (nrow(identity) != 401L || anyDuplicated(identity$sample_id) || anyDuplicated(identity$gsm) ||
      anyDuplicated(paste(identity$patient_id, identity$timepoint)) ||
      any(identity$mapping_status != "exact") || any(missing_value(as.matrix(identity))))
    stage03_stop("Sample identity unreliable", "STOP_REVIEW")
  pairs <- manifest[manifest$paired == "TRUE", , drop = FALSE]
  if (nrow(pairs) != 150L || anyDuplicated(pairs$patient_id) ||
      any(pairs$baseline_present != "TRUE" | pairs$d1c2_present != "TRUE" |
          pairs$metadata_complete != "TRUE") || any(missing_value(pairs$arm)))
    stage03_stop("Paired identity unreliable", "STOP_REVIEW")
  for (tp in c("baseline", "d1c2")) {
    ids <- pairs[[paste0(tp, "_sample_id")]]
    i <- match(ids, identity$sample_id)
    if (anyNA(i) || anyDuplicated(ids) || any(identity$patient_id[i] != pairs$patient_id) ||
        any(identity$gsm[i] != pairs[[paste0(tp, "_gsm")]]) ||
        any(identity$timepoint[i] != if (tp == "baseline") "Baseline" else "D1C2"))
      stage03_stop("Paired sample identity unreliable", "STOP_REVIEW")
  }
  list(identity = identity[, 1:4], pairs = pairs)
}
read_frozen_panel <- function() {
  file <- path("metadata/pathways/confirmatory_gene_sets_v1.tsv")
  panel <- yaml::read_yaml(path("metadata/pathways/confirmatory_panel_v1.yml"))
  checksum <- readLines(path("metadata/pathways/confirmatory_gene_sets_v1.sha256"))
  if (!identical(sha256(file), panel$sha256) || !identical(checksum, panel$sha256))
    stage03_stop("Frozen membership checksum failed", "STOP_REVIEW")
  x <- read_tsv(file)
  if (!setequal(unique(x$axis), names(axes)) || any(x$gs_name != unname(axes[x$axis])) ||
      any(!x$db_version %in% "2026.1.Hs") || any(missing_value(x$gene_symbol)) ||
      anyDuplicated(paste(x$axis, x$gene_symbol)) ||
      !identical(panel$msigdbr_version, "26.1.1") || !identical(panel$msigdb_version, "2026.1"))
    stage03_stop("Frozen membership invalid", "STOP_REVIEW")
  if (!identical(order(x$axis, x$gene_symbol, method = "radix"), seq_len(nrow(x))) ||
      !identical(panel$species, "Homo sapiens") || !identical(panel$collection, "H") ||
      !identical(panel$primary_scoring, "singscore") || !identical(panel$sensitivity_scoring, "ssGSEA") ||
      !identical(panel$sensitivity_status, "DEFERRED") || length(panel$pathways) != 6L)
    stage03_stop("Frozen panel specification changed", "STOP_REVIEW")
  for (p in panel$pathways) {
    if (!p$axis %in% names(axes) || !identical(p$gs_name, axes[[p$axis]]) ||
        p$gene_count != sum(x$axis == p$axis)) stage03_stop("Frozen gene count changed", "STOP_REVIEW")
  }
  x
}
freeze_panel <- function() {
  file <- path("metadata/pathways/confirmatory_gene_sets_v1.tsv")
  if (file.exists(file) || file.exists(path("metadata/pathways/confirmatory_panel_v1.yml"))) {
    read_frozen_panel()
    return(invisible(NULL))
  }
  if (!requireNamespace("msigdbr", quietly = TRUE) ||
      as.character(utils::packageVersion("msigdbr")) != "26.1.1")
    stage03_stop("Exact msigdbr release unavailable")
  x <- msigdbr::msigdbr(db_species = "HS", species = "Homo sapiens", collection = "H")
  if (!identical(unique(x$db_version), "2026.1.Hs"))
    stage03_stop("Exact MSigDB release unavailable")
  x <- x[x$gs_name %in% unname(axes), , drop = FALSE]
  x$axis <- names(axes)[match(x$gs_name, axes)]
  x$entrez_gene <- x$ncbi_gene
  x$source <- "msigdbr 26.1.1; MSigDB human 2026.1; collection H"
  x <- x[, c("axis", "gs_name", "gene_symbol", "entrez_gene", "ensembl_gene", "db_version", "source")]
  x <- unique(as.data.frame(x))
  x <- x[order(x$axis, x$gene_symbol, method = "radix"), , drop = FALSE]
  if (!setequal(unique(x$axis), names(axes)) || anyDuplicated(paste(x$axis, x$gene_symbol)) ||
      any(missing_value(x$gene_symbol))) stage03_stop("Gene membership ambiguous")
  write_tsv(x, file)
  hash <- sha256(file)
  writeLines(hash, path("metadata/pathways/confirmatory_gene_sets_v1.sha256"))
  panel <- list(panel_version = "v1", freeze_timestamp = utc_now(), msigdbr_version = "26.1.1",
    msigdb_version = "2026.1", species = "Homo sapiens", collection = "H",
    primary_scoring_method = "singscore", primary_scoring = "singscore",
    sensitivity_scoring = "ssGSEA", sensitivity_status = "DEFERRED",
    direction = "positive", sha256 = hash, decision_status = "FROZEN",
    source_documentation = "https://igordot.github.io/msigdbr/news/index.html",
    pathways = lapply(names(axes), function(a) list(axis = a, gs_name = axes[[a]],
      gene_count = sum(x$axis == a))))
  yaml::write_yaml(panel, path("metadata/pathways/confirmatory_panel_v1.yml"))
  read_frozen_panel()
  invisible(NULL)
}
stage03_run <- function(step, code) {
  stage03_init()
  run_stage(step, {
    tryCatch({
      force(code)
      failure <- path("results/stage03_pathways/stage03_failure.json")
      if (file.exists(failure) && identical(jsonlite::read_json(failure)$step, step)) unlink(failure)
    }, error = function(e) {
      jsonlite::write_json(list(step = step, timestamp = utc_now(),
        decision = if (inherits(e, "stage03_error")) e$decision else "MODIFY_REVIEW",
        reason = conditionMessage(e)), path("results/stage03_pathways/stage03_failure.json"),
        auto_unbox = TRUE, pretty = TRUE)
      stop(e)
    })
  })
}

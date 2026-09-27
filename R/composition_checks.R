stage01_status_ok <- function(status) {
  identical(status$decision, "PROCEED") && isTRUE(status$paired_complete == 150) &&
    isTRUE(status$mapping_failures == 0) && isTRUE(status$metadata_conflicts == 0)
}
verify_stage01 <- function() {
  required <- c("metadata/neotrip_pair_manifest_v1.tsv", "metadata/neotrip_mapping_audit.tsv",
    "metadata/data_inventory.tsv", "results/stage01_manifest/expression_audit.tsv",
    "results/stage01_manifest/stage01_status.json",
    "data/raw/geo/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz",
    "data/raw/author_repo/NeoTRIP_baseline_D1C2_pheno_and_clinical.txt")
  if (!all(file.exists(path(required)))) stop("Stage 1 inputs missing")
  status <- jsonlite::read_json(path("results/stage01_manifest/stage01_status.json"), simplifyVector = TRUE)
  if (!stage01_status_ok(status)) stop("Stage 1 gate failed")
  verify_inventory()
  audit <- read_tsv(path("metadata/neotrip_mapping_audit.tsv"))
  identity <- audit[, c("gsm", "geo_library_name", "geo_patient_id", "geo_timepoint", "mapping_status")]
  names(identity) <- c("gsm", "sample_id", "patient_id", "timepoint", "mapping_status")
  if (nrow(identity) != status$expression_samples || any(identity$mapping_status != "exact") ||
      anyDuplicated(identity$sample_id) || anyDuplicated(identity$gsm) || any(missing_value(identity$patient_id)) ||
      anyDuplicated(paste(identity$patient_id, identity$timepoint)) ||
      any(!identity$timepoint %in% c("Baseline", "D1C2")) ||
      any(audit$geo_library_name != audit$author_sample_id) || any(audit$geo_patient_id != audit$author_patient_id) ||
      any(audit$geo_timepoint != audit$author_timepoint)) stop("Stage 1 identity failed")
  pairs <- read_tsv(path("metadata/neotrip_pair_manifest_v1.tsv"))
  pairs <- pairs[, c("patient_id", "baseline_sample_id", "d1c2_sample_id", "baseline_gsm", "d1c2_gsm", "paired")]
  pairs <- pairs[pairs$paired == "TRUE", , drop = FALSE]
  if (nrow(pairs) != 150L || anyDuplicated(pairs$patient_id)) stop("Stage 1 pairing failed")
  for (tp in c("baseline", "d1c2")) {
    index <- match(pairs[[paste0(tp, "_sample_id")]], identity$sample_id)
    expected <- if (tp == "baseline") "Baseline" else "D1C2"
    if (anyNA(index) || any(identity$patient_id[index] != pairs$patient_id) ||
        any(identity$timepoint[index] != expected) || any(identity$gsm[index] != pairs[[paste0(tp, "_gsm")]])) stop("Stage 1 pair mismatch")
  }
  frozen <- c(required[!grepl("data_inventory", required)], "config/project.yml")
  hashes <- data.frame(path = frozen, sha256 = vapply(path(frozen), sha256, ""))
  file <- path("metadata/composition/stage01_input_checksums.tsv")
  if (file.exists(file)) {
    previous <- read_tsv(file)
    if (!identical(previous$path, hashes$path) || !identical(previous$sha256, hashes$sha256)) stop("Stage 1 inputs changed")
  } else write_tsv(hashes, file)
  list(identity = identity[, 1:4], pairs = pairs[, 1:5])
}
pathology_audit <- function() {
  table <- data.frame(variable = c("pathology_tumour_cell_presence_or_percentage", "Early.pCR", "pCR", "TNBCtype", "Arm", "Cancer_Epithelial", "IMC epithelial masks", "publisher supplementary pathology"),
    source = c("GEO SOFT and public author data schema", rep("public author phenotype", 4), "12-SuppFigure5.R:36;13-SuppFigure6.R:24-30", "Zenodo 7990870 metadata.description and metadata.access_right", "Publisher article/supplement access audit"),
    level = c("sample/patient", rep("sample", 5), "image/cell", "unknown"),
    link_key = c("none for a qualifying variable", rep("Sample_ID", 5), "unverified", "unverified"),
    public = c("FALSE", rep("TRUE", 4), "generation_code", "metadata_only", "UNKNOWN_ACCESS"),
    linkable = c(FALSE, rep(TRUE, 5), FALSE, FALSE),
    coverage = c("0/401 samples; 0/251 patients", rep("401/401 sample identifiers; values not used", 4), "401 samples eligible for reproduction", "unverified: restricted files", "unverified: full text/supplements not retrieved"),
    definition = c("Independent pathology-derived tumour presence/content measurement", "Author combined early response category; prohibited substitute", "Surgical outcome; prohibited substitute", "Molecular subtype; prohibited substitute", "Treatment allocation; prohibited substitute", "Transcriptomic epithelial fraction; not pathology", "Dataset description lists epithelial masks; measurements and keys not publicly inspected", "No accessible keyed measurement identified"),
    decision = c("UNAVAILABLE_PUBLICLY", rep("EXCLUDED_BY_PROTOCOL", 5), "UNAVAILABLE_PUBLICLY", "UNAVAILABLE_PUBLICLY"),
    reason = c("No qualifying measurement in inspected public data; restricted and inaccessible sources cannot establish public linkage.",
      rep("Not accepted as pathology tumour content; no recoding or proxy construction.", 5),
      "Public record only; data access requires request. No request made.",
      "Access limitation retained; no inference of biological absence or unpublished-variable absence."))
  if (file.exists(path("metadata/composition/imc_key_audit.tsv"))) {
    table <- rbind(table, data.frame(variable = "isEpithelial;cellAnnotation;cell_cluster",
      source = "Zenodo 14290163 nt_preprocessed.h5ad /obs", level = "cell/image; baseline only",
      link_key = "Patient_ID;Tissue_ID;fov;_index", public = "TRUE", linkable = FALSE,
      coverage = "119 IMC baseline samples described; zero exact Stage 1 key matches",
      definition = "Image-derived epithelial/invasive cell annotations; not an independently verified histopathologic tumour percentage",
      decision = "UNLINKABLE_EXACT_KEYS",
      reason = "All candidate public keys tested against Stage 1 verbatim; zero matches. No inferred NT-to-TNBC conversion, no cell-fraction reconstruction, no outcome values read."))
  }
  write_tsv(table, path("metadata/composition/pathology_variable_audit.tsv"))
  table
}
composition_qc <- function(fractions, identity, expected_samples) {
  checks <- c(sample_count = nrow(fractions) != length(expected_samples),
    sample_order = !identical(rownames(fractions), expected_samples),
    duplicate_sample_ids = anyDuplicated(rownames(fractions)) > 0L,
    missing_fractions = sum(is.na(fractions)), non_finite_values = sum(!is.finite(fractions)),
    negative_fractions = sum(fractions < 0, na.rm = TRUE), fractions_above_one = sum(fractions > 1, na.rm = TRUE),
    component_sum = sum(abs(rowSums(fractions) - 1) > 1e-8, na.rm = TRUE),
    unmatched_samples = sum(!rownames(fractions) %in% identity$sample_id),
    patient_mapping = sum(missing_value(identity$patient_id[match(rownames(fractions), identity$sample_id)])),
    timepoint_mapping = sum(!identity$timepoint[match(rownames(fractions), identity$sample_id)] %in% c("Baseline", "D1C2")))
  data.frame(check = names(checks), failures = as.integer(checks), passed = checks == 0,
    detail = "Full cell-type theta fractions in [0,1]; sum tolerance 1e-8; Stage 1 exact keys.")
}
pair_composition <- function(composition, pairs) {
  baseline <- match(pairs$baseline_sample_id, composition$sample_id)
  d1c2 <- match(pairs$d1c2_sample_id, composition$sample_id)
  if (anyNA(baseline) || anyNA(d1c2)) stop("Composition pairing failed")
  output <- pairs[, c("patient_id", "baseline_sample_id", "d1c2_sample_id")]
  for (field in c("cancer_epithelial", "stromal_fraction", "immune_fraction")) {
    output[[paste0("baseline_", field)]] <- composition[[field]][baseline]
    output[[paste0("d1c2_", field)]] <- composition[[field]][d1c2]
    output[[paste0("delta_", field)]] <- composition[[field]][d1c2] - composition[[field]][baseline]
  }
  output
}

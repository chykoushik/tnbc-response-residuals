audit_instaprism <- function() {
  directory <- "external/instaprism_source"
  package_url <- "https://github.com/humengying0907/InstaPrism"
  pin_repository(directory, package_url, instaprism_commit)
  for (file in git_output(directory, "ls-files")) {
    inventory_source(paste0(directory, "/", file), paste0(package_url, "/blob/", instaprism_commit, "/", file),
      "InstaPrism source", package_url, instaprism_commit)
  }
  reference_url <- paste0("https://raw.githubusercontent.com/humengying0907/InstaPrismSourceCode/", reference_commit, "/refPhi/BRCA_refPhi.RDS")
  public <- fetch_public(reference_url, "external/instaprism_reference/BRCA_refPhi.RDS", "InstaPrism BRCA reference",
    "https://github.com/humengying0907/InstaPrismSourceCode", reference_commit)
  fetch_public("https://api.github.com/repos/humengying0907/InstaPrismSourceCode/commits?path=refPhi/BRCA_refPhi.RDS&per_page=100",
    "external/instaprism_reference/reference_history.json", "BRCA reference history", "InstaPrismSourceCode")
  if (public) {
    ref <- readRDS(path("external/instaprism_reference/BRCA_refPhi.RDS"))
    cells <- attr(ref, "map")
    cell_types <- names(cells)
    write_tsv(data.frame(cell_type = rep(names(cells), lengths(cells)), cell_state = unlist(cells, use.names = FALSE)),
      path("metadata/composition/instaprism_reference_cells.tsv"))
  } else cell_types <- character()
  source <- "external/neotrip_2026_repo/"
  lines <- readLines(path(source, "scripts/12-SuppFigure5.R"), warn = FALSE)
  aggregation <- readLines(path(source, "scripts/13-SuppFigure6.R"), warn = FALSE)
  lock <- jsonlite::read_json(path(source, "renv.lock"), simplifyVector = TRUE)$Packages$InstaPrism
  command_ok <- identical(lines[29L], "deconv_tpm = InstaPrism(bulk_Expr = expData, refPhi_cs = ref)")
  reference_ok <- identical(lines[28L], "ref = InstaPrism_reference('BRCA')")
  pheno_columns <- names(read_tsv(path("data/raw/author_repo/NeoTRIP_baseline_D1C2_pheno_and_clinical.txt")))
  expected_pheno <- c("Sample_ID", "Patient_ID", "Timepoint", "Arm", "pCR", "pCR.num", "Early.pCR", "TNBCtype", "Batch_correlation")
  columns_ok <- identical(pheno_columns, expected_pheno)
  aggregation_ok <- columns_ok && length(cell_types) == 8L &&
    identical(aggregation[29L], "scoreSS$Stromal_fraction <- rowSums(pData(scoreSS)[, c(10:12)])") &&
    identical(aggregation[30L], "scoreSS$Immune_fraction <- rowSums(pData(scoreSS)[, c(13:15)])") &&
    identical(gsub(" ", "_", cell_types[8L]), "Cancer_Epithelial")
  method <- data.frame(item = character(), value = character(), source_file = character(), source_line = character(), verified = logical(), notes = character())
  add <- function(item, value, file, line, verified = TRUE, notes = "") {
    method <<- rbind(method, data.frame(item = item, value = value, source_file = file, source_line = line, verified = verified, notes = notes))
  }
  add("package", "InstaPrism", paste0(source, "renv.lock"), "884-933")
  add("package_version", lock$Version, paste0(source, "renv.lock"), "886", lock$Version == "0.1.6")
  add("package_commit", lock$RemoteSha, paste0(source, "renv.lock"), "932", lock$RemoteSha == instaprism_commit)
  add("generation_command", lines[29], paste0(source, "scripts/12-SuppFigure5.R"), "21-30", command_ok)
  add("reference", "InstaPrism_reference('BRCA')", paste0(source, "scripts/12-SuppFigure5.R"), "28", reference_ok)
  add("reference_public", reference_url, paste0(directory, "/R/InstaPrism.R"), "951-968", public,
    "Original function resolves mutable main; archived file revision pinned from public history, last changed 2024-05-15. Original run reference checksum not supplied.")
  add("input_scale", "Linear ComBat-adjusted TPM", paste0(source, "scripts/12-SuppFigure5.R"), "21-29", TRUE,
    "Stage 1 official expression source; no outcome-dependent filtering or log transform.")
  add("preprocessing", "HGNC symbols become row names; remove three annotation columns; no external gene filter", paste0(source, "scripts/12-SuppFigure5.R"), "21-27")
  add("source_file_resolution", "Official GSE319641 supplementary TPM file retained with GEO accession prefix", paste0(source, "scripts/01-data_download_and_singscore.R"), "14-18;34-39", TRUE,
    "Authors' local read path omits the GEO prefix; reproduction reads the official gzip directly without renaming, unzipping to disk or editing source bytes.")
  add("negative_input_handling", "Published code supplies linear ComBat values unchanged; neither caller nor default bpPrepare clips negative bulk values",
    paste0(source, "scripts/12-SuppFigure5.R"), "21-29;external/instaprism_source/R/InstaPrism.R:86-161", TRUE,
    "Source input QC records negative entries. Retained for literal reproduction; this is a methodological limitation, not proof of biological negative abundance.")
  add("internal_preprocessing", "Default bpPrepare: shared genes, bulk outlier filtering, reference normalization; n.iter=max(100,2*cell_states)",
    paste0(directory, "/R/InstaPrism.R"), "208-285;bpPrepare 84-181")
  add("cell_types", paste(cell_types, collapse = ";"), "external/instaprism_reference/BRCA_refPhi.RDS", "map slot; ordered names", length(cell_types) == 8L)
  add("Cancer_Epithelial", if (length(cell_types) == 8L) cell_types[8L] else "unverified", paste0(source, "scripts/13-SuppFigure6.R"), "24-28", aggregation_ok,
    "Cell type 8 in public BRCA reference; phenodata has nine columns before binding theta. Name normalization only; not pathology.")
  add("Stromal_fraction", paste(head(cell_types, 3L), collapse = "+"), paste0(source, "scripts/13-SuppFigure6.R"), "26;29", aggregation_ok,
    "Phenodata columns 10:12 correspond to first three reference cell types; exact positional rule.")
  add("Immune_fraction", if (length(cell_types) == 8L) paste(cell_types[4:6], collapse = "+") else "unverified", paste0(source, "scripts/13-SuppFigure6.R"), "26;30", aggregation_ok,
    "Phenodata columns 13:15 correspond to reference cell types 4:6; normal epithelial component is not added.")
  add("fraction_units", "Full theta: 0-1; summary and paired changes: percent/percentage points", paste0(source, "scripts/13-SuppFigure6.R"), "24-30", TRUE,
    "Published aggregation multiplies theta by 100. The full eight-type vector sums to 1; the three summaries exclude normal epithelium and need not sum to 100.")
  add("normalization_evidence", "ppguess divided by sum(ppguess); cell-state theta summed by map", paste0(directory, "/src/cppFuncs.cpp"), "55;R/InstaPrism.R:302")
  write_tsv(method, path("metadata/composition/instaprism_method_audit.tsv"))
  feasibility <- list(package_identified = TRUE, package_version_identified = lock$Version == "0.1.6" && lock$RemoteSha == instaprism_commit,
    reference_identified = reference_ok, reference_public = public, input_scale_identified = TRUE,
    generation_code_identified = command_ok, cell_types_identified = length(cell_types) == 8L,
    aggregation_rules_identified = aggregation_ok)
  feasibility$reproducible_from_public_inputs <- all(unlist(feasibility))
  feasibility$reason <- if (feasibility$reproducible_from_public_inputs) {
    "Exact package commit, published command, BRCA reference, linear TPM input and positional aggregation identified. Public reference revision pinned; historical run checksum not supplied. No alternative reference or method substituted."
  } else "Published method has an unresolved component; no deconvolution or substitute implementation is permitted."
  jsonlite::write_json(feasibility, path("results/stage02_composition/instaprism_feasibility.json"), auto_unbox = TRUE, pretty = TRUE)
  feasibility
}
run_instaprism <- function(authority) {
  feasibility <- jsonlite::read_json(path("results/stage02_composition/instaprism_feasibility.json"), simplifyVector = TRUE)
  if (!isTRUE(feasibility$reproducible_from_public_inputs)) return(invisible(FALSE))
  if (!requireNamespace("InstaPrism", quietly = TRUE)) stop("Pinned InstaPrism missing")
  description <- utils::packageDescription("InstaPrism")
  if (description$Version != "0.1.6" || description$RemoteSha != instaprism_commit) stop("InstaPrism version mismatch")
  suppressPackageStartupMessages(library(InstaPrism))
  dir.create(path("data/processed/composition"), recursive = TRUE, showWarnings = FALSE)
  expression_file <- path("data/raw/geo/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz")
  connection <- gzfile(expression_file, "rt")
  expression <- read.delim(connection, check.names = FALSE, stringsAsFactors = FALSE)
  close(connection)
  genes <- expression$HGNC_approved_symbol
  if (anyDuplicated(genes) || any(missing_value(genes))) stop("Expression genes invalid")
  expression <- as.matrix(expression[, -(1:3)])
  rownames(expression) <- genes
  samples <- colnames(expression)
  if (anyDuplicated(samples) || !setequal(samples, authority$identity$sample_id)) stop("Expression identity mismatch")
  input_qc <- data.frame(check = c("nonfinite_values", "negative_values", "samples_with_negative_values", "minimum_value", "nonpositive_sample_totals"),
    value = c(sum(!is.finite(expression)), sum(expression < 0), sum(colSums(expression < 0) > 0), min(expression), sum(colSums(expression) <= 0)),
    interpretation = c("Must be zero", "Retained unchanged as supplied in published linear ComBat TPM input",
      "No samples excluded; published code does not clip negative values", "Observed source minimum", "Must be zero"))
  write_tsv(input_qc, path("results/stage02_composition/expression_input_qc.tsv"))
  if (any(!is.finite(expression)) || any(colSums(expression) <= 0)) stop("Expression values invalid")
  reference <- readRDS(path("external/instaprism_reference/BRCA_refPhi.RDS"))
  expected_types <- c("Endothelial", "CAFs", "PVL", "B-cells", "T-cells", "Myeloid", "Normal Epithelial", "Cancer Epithelial")
  if (!identical(names(reference@map), expected_types)) stop("Reference cell order changed")
  warnings <- character()
  started <- utc_now()
  result <- withCallingHandlers(
    InstaPrism::InstaPrism(bulk_Expr = expression, refPhi_cs = reference),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  saveRDS(result, path("data/processed/composition/neotrip_instaprism_result.rds"))
  fractions <- t(result@Post.ini.ct@theta)
  colnames(fractions) <- gsub(" ", "_", colnames(fractions))
  qc <- composition_qc(fractions, authority$identity, samples)
  qc <- rbind(qc, data.frame(check = "cell_type_order", failures = as.integer(!identical(colnames(fractions), gsub(" ", "_", expected_types))),
    passed = identical(colnames(fractions), gsub(" ", "_", expected_types)), detail = "Pinned BRCA map order."))
  qc <- rbind(qc, data.frame(check = "method_warnings", failures = length(warnings), passed = !length(warnings),
    detail = if (length(warnings)) paste(unique(warnings), collapse = ";") else "No method warnings."))
  write_tsv(qc, path("results/stage02_composition/composition_qc.tsv"))
  write_tsv(data.frame(warning = unique(warnings)), path("results/stage02_composition/instaprism_warnings.tsv"))
  identity <- authority$identity[match(samples, authority$identity$sample_id), c("sample_id", "gsm", "patient_id", "timepoint")]
  write_gzip <- function(table, name) {
    connection <- gzfile(path("data/processed/composition", name), "wt")
    on.exit(close(connection))
    write_tsv(table, connection)
  }
  write_gzip(cbind(identity, as.data.frame(fractions, check.names = FALSE)), "neotrip_instaprism_full.tsv.gz")
  states <- t(result@Post.ini.cs@theta)
  if (!identical(rownames(states), samples)) stop("Cell state order mismatch")
  write_gzip(cbind(identity, as.data.frame(states, check.names = FALSE)), "neotrip_instaprism_cell_states.tsv.gz")
  provenance <- list(started = started, completed = utc_now(), package_version = description$Version,
    package_commit = description$RemoteSha, repository_commit = repo_commit,
    expression_sha256 = sha256(expression_file),
    reference_sha256 = sha256(path("external/instaprism_reference/BRCA_refPhi.RDS")),
    input_genes = nrow(expression), input_samples = ncol(expression),
    used_genes = nrow(result@initial.reference@phi.cs), cell_types = colnames(fractions), cell_states = ncol(states),
    command = "InstaPrism(bulk_Expr = expression, refPhi_cs = reference)",
    preprocessing = "HGNC row names; discard annotation columns; unchanged linear TPM including negative ComBat values. Published code and package default bpPrepare contain no clipping step. Package default bulk outlier filtering and shared-gene reference normalization retained; negative-input limitation documented, no correction invented.",
    defaults = list(n.iter = max(100L, 2L * ncol(reference@phi.cs)), filter = TRUE, outlier.cut = 0.01,
      outlier.fraction = 0.1, pseudo.min = 1e-8, n.core = 1L, input_type = "refPhi_cs"),
    units = list(full = "0-1", summary = "percent", change = "percentage_points"),
    identities = "Stage 1 exact mapping; all expression samples; no clinical covariates or outcomes used.",
    environment = "Pinned InstaPrism 0.1.6 commit; project R and dependency versions recorded separately; no claim of bitwise reproduction of original authors' results.")
  jsonlite::write_json(provenance, path("metadata/composition/reproduction_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
  if (any(!qc$passed)) stop("Composition QC failed")
  composition <- cbind(identity, data.frame(cancer_epithelial = fractions[, "Cancer_Epithelial"] * 100,
    stromal_fraction = rowSums(fractions[, 1:3, drop = FALSE]) * 100,
    immune_fraction = rowSums(fractions[, 4:6, drop = FALSE]) * 100))
  write_tsv(composition, path("data/processed/composition/neotrip_composition.tsv"))
  changes <- pair_composition(composition, authority$pairs)
  write_tsv(changes, path("data/processed/composition/neotrip_pair_composition_change.tsv"))
  invisible(TRUE)
}

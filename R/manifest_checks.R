empty_issues <- function() {
  data.frame(issue = character(), severity = character(), patient_id = character(),
             sample_id = character(), gsm = character(), detail = character())
}
duplicate_values <- function(x) {
  unique(x[!missing_value(x) & (duplicated(x) | duplicated(x, fromLast = TRUE))])
}
expression_check <- function(file, phenotype_ids) {
  connection <- gzfile(file, "rt")
  on.exit(close(connection))
  split_line <- function(x) {
    fields <- strsplit(paste0(x, "\t.END"), "\t", fixed = TRUE)[[1L]]
    head(fields, -1L)
  }
  header <- split_line(readLines(connection, n = 1L, warn = FALSE))
  annotations <- c("HGNC_approved_symbol", "HGNC_id", "entrez_gene_id")
  if (!identical(header[1:3], annotations)) stop("Expression schema changed")
  samples <- header[-(1:3)]
  genes <- lapply(annotations, function(x) character())
  malformed <- 0L
  repeat {
    lines <- readLines(connection, n = 1000L, warn = FALSE)
    if (!length(lines)) break
    fields <- strsplit(paste0(lines, "\t.END"), "\t", fixed = TRUE)
    widths <- lengths(fields) - 1L
    malformed <- malformed + sum(widths != length(header))
    for (j in 1:3) genes[[j]] <- c(genes[[j]], vapply(fields, function(x) x[j], ""))
  }
  audit <- data.frame(check = character(), count = integer(), detail = character())
  add <- function(check, count, detail = "") {
    audit <<- rbind(audit, data.frame(check = check, count = count, detail = detail))
  }
  add("expression_sample_columns", length(samples))
  add("gene_rows", length(genes[[1L]]))
  add("duplicate_sample_names", sum(duplicated(samples)), paste(duplicate_values(samples), collapse = ";"))
  add("missing_sample_names", sum(missing_value(samples)))
  add("malformed_expression_rows", malformed)
  for (j in 1:3) {
    values <- genes[[j]][!missing_value(genes[[j]])]
    add(paste0("duplicate_", annotations[j]), sum(duplicated(values)), paste(duplicate_values(values), collapse = ";"))
    add(paste0("missing_", annotations[j]), sum(missing_value(genes[[j]])))
  }
  unmatched_expression <- setdiff(samples, phenotype_ids)
  unmatched_phenotype <- setdiff(phenotype_ids, samples)
  add("unmatched_expression_samples", length(unmatched_expression), paste(unmatched_expression, collapse = ";"))
  add("unmatched_phenotype_samples", length(unmatched_phenotype), paste(unmatched_phenotype, collapse = ";"))
  add("phenotype_expression_expectation", length(phenotype_ids), "All supplied phenotype rows are expected to have expression.")
  list(samples = samples, audit = audit, unmatched_expression = unmatched_expression,
    unmatched_phenotype = unmatched_phenotype,
    fatal = anyDuplicated(samples) > 0L || any(missing_value(samples)) || malformed > 0L || !length(genes[[1L]]))
}
build_manifest <- function(geo, author, expression, config) {
  issues <- empty_issues()
  add_issue <- function(issue, severity = "error", patient = "", sample = "", gsm = "", detail = "") {
    issues <<- rbind(issues, data.frame(issue = issue, severity = severity,
      patient_id = patient, sample_id = sample, gsm = gsm, detail = detail))
  }
  report_rows <- function(index, issue, severity = "error") {
    for (i in which(index)) add_issue(issue, severity, author$Patient_ID[i], author$Sample_ID[i])
  }
  duplicated_gsm <- duplicate_values(geo$gsm)
  duplicated_library <- duplicate_values(geo$geo_library_name)
  duplicated_sample <- duplicate_values(author$Sample_ID)
  for (i in which(missing_value(geo$gsm))) add_issue("missing_gsm", "fatal", sample = geo$geo_library_name[i])
  geo_patient_time <- paste(geo$geo_patient_id, geo$geo_timepoint, sep = "\r")
  geo_duplicate_time <- geo_patient_time %in% duplicate_values(geo_patient_time) &
    !missing_value(geo$geo_patient_id) & !missing_value(geo$geo_timepoint)
  for (i in which(geo_duplicate_time)) {
    add_issue("duplicate_geo_patient_timepoint", "fatal", geo$geo_patient_id[i], geo$geo_library_name[i], geo$gsm[i])
  }
  for (gsm in duplicated_gsm) add_issue("duplicate_gsm", "fatal", gsm = gsm)
  for (id in duplicated_library) add_issue("duplicate_geo_library", "fatal", sample = id)
  report_rows(author$Sample_ID %in% duplicated_sample, "duplicate_sample_id", "fatal")
  report_rows(missing_value(author$Sample_ID), "missing_sample_id", "fatal")
  report_rows(missing_value(author$Patient_ID), "missing_patient_id")
  report_rows(missing_value(author$Arm), "missing_treatment_arm")
  report_rows(!missing_value(author$Arm) & !author$Arm %in% config$arms, "invalid_treatment_arm")
  report_rows(missing_value(author$pCR), "missing_pcr")
  report_rows(!missing_value(author$pCR) & !author$pCR %in% config$pcr_values, "invalid_pcr")
  report_rows(missing_value(author$pCR.num), "missing_pcr_num")
  report_rows(!missing_value(author$pCR.num) & !author$pCR.num %in% config$pcr_numeric, "invalid_pcr_num")
  report_rows(!author$Timepoint %in% config$timepoints, "invalid_timepoint")
  pcr_agrees <- (author$pCR == "pCR" & author$pCR.num == "1") | (author$pCR == "RD" & author$pCR.num == "0")
  report_rows(!missing_value(author$pCR) & !missing_value(author$pCR.num) & !pcr_agrees, "pcr_encoding_conflict")
  patient_time <- paste(author$Patient_ID, author$Timepoint, sep = "\r")
  duplicated_patient_time <- patient_time %in% duplicate_values(patient_time)
  report_rows(duplicated_patient_time, "duplicate_patient_timepoint", "fatal")
  report_rows(duplicated_patient_time & author$Timepoint == "Baseline", "duplicate_baseline", "fatal")
  report_rows(duplicated_patient_time & author$Timepoint == "D1C2", "duplicate_d1c2", "fatal")
  for (i in which(!geo$geo_timepoint %in% config$timepoints)) {
    add_issue("invalid_geo_timepoint", patient = geo$geo_patient_id[i], gsm = geo$gsm[i])
  }
  mapping <- data.frame(gsm = character(), geo_library_name = character(), geo_patient_id = character(),
    geo_timepoint = character(), geo_batch = character(), author_sample_id = character(),
    author_patient_id = character(), author_timepoint = character(), author_batch = character(),
    arm = character(), pcr = character(), pcr_num = character(), early_pcr = character(),
    tnbc_type = character(), mapping_status = character())
  append_mapping <- function(g, a, status) {
    gv <- function(field) if (length(g)) geo[[field]][g] else ""
    av <- function(field) if (length(a)) author[[field]][a] else ""
    mapping <<- rbind(mapping, data.frame(gsm = gv("gsm"), geo_library_name = gv("geo_library_name"),
      geo_patient_id = gv("geo_patient_id"), geo_timepoint = gv("geo_timepoint"), geo_batch = gv("geo_batch"),
      author_sample_id = av("Sample_ID"), author_patient_id = av("Patient_ID"), author_timepoint = av("Timepoint"),
      author_batch = av("Batch_correlation"), arm = av("Arm"), pcr = av("pCR"), pcr_num = av("pCR.num"),
      early_pcr = av("Early.pCR"), tnbc_type = av("TNBCtype"), mapping_status = status))
  }
  used <- integer()
  for (g in seq_len(nrow(geo))) {
    matches <- which(!missing_value(geo$geo_library_name[g]) & author$Sample_ID == geo$geo_library_name[g])
    if (!length(matches)) {
      append_mapping(g, integer(), "missing_author")
      add_issue("missing_author", patient = geo$geo_patient_id[g], sample = geo$geo_library_name[g], gsm = geo$gsm[g])
      next
    }
    for (a in matches) {
      used <- c(used, a)
      status <- character()
      if (geo$gsm[g] %in% duplicated_gsm || geo$geo_library_name[g] %in% duplicated_library ||
          author$Sample_ID[a] %in% duplicated_sample || duplicated_patient_time[a]) status <- c(status, "duplicate")
      if (any(missing_value(c(geo$geo_patient_id[g], geo$geo_timepoint[g], geo$geo_batch[g])))) status <- c(status, "missing_geo")
      if (any(missing_value(c(author$Patient_ID[a], author$Timepoint[a], author$Batch_correlation[a])))) status <- c(status, "missing_author")
      if (geo$geo_patient_id[g] != author$Patient_ID[a]) status <- c(status, "patient_conflict")
      if (geo$geo_timepoint[g] != author$Timepoint[a]) status <- c(status, "timepoint_conflict")
      if (geo$geo_batch[g] != author$Batch_correlation[a]) status <- c(status, "batch_conflict")
      for (problem in status) add_issue(problem, patient = author$Patient_ID[a], sample = author$Sample_ID[a], gsm = geo$gsm[g])
      append_mapping(g, a, if (length(status)) status[1L] else "exact")
    }
  }
  for (a in setdiff(seq_len(nrow(author)), used)) {
    append_mapping(integer(), a, "missing_geo")
    add_issue("missing_geo", patient = author$Patient_ID[a], sample = author$Sample_ID[a])
  }
  write_tsv(mapping, path("metadata/neotrip_mapping_audit.tsv"))
  mapping$patient_id <- ifelse(mapping$mapping_status == "exact", mapping$geo_patient_id, "")
  mapping$timepoint <- ifelse(mapping$mapping_status == "exact", mapping$geo_timepoint, "")
  mapping$expression_present <- mapping$author_sample_id %in% expression$samples
  write_tsv(mapping, path("metadata/neotrip_sample_manifest_v1.tsv"))
  for (id in expression$unmatched_expression) add_issue("unmatched_expression_sample", sample = id)
  for (id in expression$unmatched_phenotype) add_issue("unmatched_phenotype_sample", sample = id)
  for (i in which(expression$audit$count > 0L & grepl("^(duplicate_|missing_|malformed_)", expression$audit$check))) {
    check <- expression$audit$check[i]
    severity <- if (check %in% c("duplicate_sample_names", "missing_sample_names", "malformed_expression_rows")) "fatal" else "warning"
    add_issue(check, severity, detail = paste("Count:", expression$audit$count[i]))
  }
  conflicts <- data.frame(patient_id = character(), field = character(), values = character(), sample_ids = character())
  for (patient in unique(author$Patient_ID[!missing_value(author$Patient_ID)])) {
    rows <- author[author$Patient_ID == patient, , drop = FALSE]
    for (field in c("Arm", "pCR", "pCR.num")) {
      values <- unique(rows[[field]][!missing_value(rows[[field]])])
      if (length(values) > 1L) {
        conflicts <- rbind(conflicts, data.frame(patient_id = patient, field = field,
          values = paste(values, collapse = ";"), sample_ids = paste(rows$Sample_ID, collapse = ";")))
        add_issue(paste0("patient_", tolower(field), "_conflict"), patient = patient, detail = paste(values, collapse = ";"))
      }
    }
  }
  write_tsv(conflicts, path("results/stage01_manifest/patient_metadata_conflicts.tsv"))
  canonical <- mapping[mapping$mapping_status == "exact", , drop = FALSE]
  pairs <- data.frame(patient_id = character(), baseline_gsm = character(), baseline_sample_id = character(),
    d1c2_gsm = character(), d1c2_sample_id = character(), arm = character(), pcr = character(), pcr_num = character(),
    baseline_present = logical(), d1c2_present = logical(), paired = logical(), metadata_complete = logical(), pair_status = character())
  for (patient in unique(canonical$patient_id)) {
    rows <- canonical[canonical$patient_id == patient, , drop = FALSE]
    baseline <- rows[rows$timepoint == "Baseline" & rows$expression_present, , drop = FALSE]
    d1c2 <- rows[rows$timepoint == "D1C2" & rows$expression_present, , drop = FALSE]
    clinical <- author[author$Patient_ID == patient, , drop = FALSE]
    one <- function(x) {
      values <- unique(x[!missing_value(x)])
      if (length(values) == 1L) values else ""
    }
    arm <- one(clinical$Arm)
    pcr <- one(clinical$pCR)
    pcr_num <- one(clinical$pCR.num)
    problems <- issues$issue[issues$patient_id == patient & issues$severity != "warning"]
    valid <- !length(problems)
    paired <- nrow(baseline) == 1L && nrow(d1c2) == 1L
    complete <- valid && arm %in% config$arms && pcr %in% config$pcr_values && pcr_num %in% config$pcr_numeric
    status <- if (!valid) "invalid_metadata" else if (!complete) "incomplete_metadata" else if (paired) "paired" else if (!nrow(baseline)) "missing_baseline" else "missing_d1c2"
    extract <- function(x, field) if (nrow(x) == 1L) x[[field]] else ""
    pairs <- rbind(pairs, data.frame(patient_id = patient,
      baseline_gsm = extract(baseline, "gsm"), baseline_sample_id = extract(baseline, "author_sample_id"),
      d1c2_gsm = extract(d1c2, "gsm"), d1c2_sample_id = extract(d1c2, "author_sample_id"),
      arm = arm, pcr = pcr, pcr_num = pcr_num, baseline_present = nrow(baseline) == 1L,
      d1c2_present = nrow(d1c2) == 1L, paired = paired, metadata_complete = complete, pair_status = status))
  }
  write_tsv(pairs, path("metadata/neotrip_pair_manifest_v1.tsv"))
  write_tsv(issues, path("results/stage01_manifest/metadata_issues.tsv"))
  usable <- pairs$paired & pairs$metadata_complete
  counts <- c(geo_samples = nrow(geo), expression_samples = length(expression$samples), phenotype_rows = nrow(author),
    unique_patients = length(unique(author$Patient_ID[!missing_value(author$Patient_ID)])),
    baseline_samples = sum(canonical$timepoint == "Baseline" & canonical$expression_present),
    d1c2_samples = sum(canonical$timepoint == "D1C2" & canonical$expression_present),
    patients_with_baseline = sum(pairs$baseline_present), patients_with_d1c2 = sum(pairs$d1c2_present),
    matched_baseline_d1c2_pairs = sum(pairs$paired),
    paired_with_arm = sum(pairs$paired & pairs$arm %in% config$arms & !pairs$patient_id %in% conflicts$patient_id),
    paired_with_pcr = sum(pairs$paired & pairs$pcr %in% config$pcr_values & !pairs$patient_id %in% conflicts$patient_id),
    paired_complete = sum(usable), ct_patients = sum(usable & pairs$arm == "CT"),
    ct_a_patients = sum(usable & pairs$arm == "CT/A"), pcr_patients = sum(usable & pairs$pcr == "pCR"),
    rd_patients = sum(usable & pairs$pcr == "RD"))
  scopes <- c(rep("source", 4), rep("verified_identity_with_expression", 5),
              rep("verified_pairs", 3), rep("usable_complete_pairs", 4))
  write_tsv(data.frame(metric = names(counts), count = unname(counts), scope = scopes),
            path("results/stage01_manifest/cohort_counts.tsv"))
  mapping_failures <- sum(mapping$mapping_status != "exact")
  fatal <- expression$fatal || any(issues$severity == "fatal")
  if (fatal || !any(usable)) {
    decision <- "STOP"
    reason <- if (fatal) "Core sample identity or matrix structure checks failed; inspect metadata_issues.tsv." else "No complete, conflict-free, explicitly mapped baseline-D1C2 pairs can be reconstructed."
  } else if (any(issues$severity != "warning")) {
    decision <- "MODIFY"
    reason <- paste("Documented exceptions are isolated in the issue reports;", sum(usable),
      "complete pairs remain reconstructable by exact metadata agreement. Exceptions require methodology review; no sample-size adequacy claim is made.")
  } else {
    decision <- "PROCEED"
    reason <- paste("All sample mappings agree exactly;", sum(usable),
      "complete pairs have consistent treatment and pCR. No unresolved systematic mismatch was detected. Gene annotation warnings, if present, do not affect sample identity. Observed counts require methodology review; sample-size adequacy is not inferred.")
  }
  status <- list(stage = 1L, timestamp = utc_now(), geo_accession = config$geo_accession,
    geo_samples = unname(counts["geo_samples"]), expression_samples = unname(counts["expression_samples"]),
    phenotype_rows = unname(counts["phenotype_rows"]), unique_patients = unname(counts["unique_patients"]),
    paired_patients = sum(pairs$paired), paired_with_arm = unname(counts["paired_with_arm"]),
    paired_with_pcr = unname(counts["paired_with_pcr"]), paired_complete = sum(usable),
    mapping_failures = mapping_failures,
    metadata_conflicts = sum(grepl("conflict", issues$issue)), decision = decision, decision_reason = reason)
  list(status = status, issues = issues, pairs = pairs, mapping = mapping)
}
save_status <- function(status) {
  jsonlite::write_json(status, path("results/stage01_manifest/stage01_status.json"),
                       auto_unbox = TRUE, pretty = TRUE, na = "null")
  log <- data.frame(timestamp = status$timestamp, stage = status$stage,
    decision = status$decision, reason = status$decision_reason,
    inputs = "metadata/data_inventory.tsv;data/raw/geo;data/raw/author_repo",
    outputs = "metadata/neotrip_mapping_audit.tsv;metadata/neotrip_pair_manifest_v1.tsv;results/stage01_manifest")
  file <- path("logs/decision_log.tsv")
  write.table(log, file, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = !file.exists(file), append = file.exists(file), na = "")
}

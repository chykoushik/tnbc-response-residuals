read_pathway_expression <- function(identity) {
  file <- path("data/raw/geo/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz")
  con <- gzfile(file, "rt")
  on.exit(close(con))
  x <- read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, quote = "",
    comment.char = "", colClasses = c(rep("character", 3), rep("numeric", 401)))
  symbols <- x$HGNC_approved_symbol
  m <- as.matrix(x[, -(1:3), drop = FALSE])
  rownames(m) <- symbols
  expected <- read_restricted(path("metadata/neotrip_sample_manifest_v1.tsv"), "geo_library_name")
  checks <- c(expression_genes_unique = anyDuplicated(x$HGNC_id) > 0,
    hgnc_symbols_unique = anyDuplicated(symbols) > 0,
    missing_symbols = sum(missing_value(symbols)), duplicated_symbols = sum(duplicated(symbols)),
    expression_samples = ncol(m) != 401L, expression_genes = nrow(m) != 23890L,
    sample_identity = !setequal(colnames(m), identity$sample_id),
    duplicate_samples = anyDuplicated(colnames(m)) > 0,
    finite_expression = sum(!is.finite(m)),
    sample_order = !identical(colnames(m), expected$geo_library_name),
    missing_entrez_retained = sum(missing_value(x$entrez_gene_id)) != 617L)
  qc <- data.frame(check = names(checks), failures = as.integer(checks), passed = checks == 0)
  qc$detail <- c("HGNC IDs", "HGNC approved symbols", "No empty symbol", "No duplicate symbol",
    as.character(ncol(m)), as.character(nrow(m)), "Exact Stage 1 sample IDs", "Unique sample IDs",
    "All values finite", "Stage 1 expression order",
    paste(sum(missing_value(x$entrez_gene_id)), "genes retained"))
  qc <- rbind(qc, data.frame(check = "negative_values_unchanged", failures = 0L, passed = TRUE,
    detail = paste(sum(m < 0), "negative values retained; no transformation")))
  write_tsv(qc, path("results/stage03_pathways/pathway_input_qc.tsv"))
  if (any(!qc$passed)) stage03_stop("Expression input validation failed")
  m
}
pathway_coverage <- function(sets, symbols) {
  coverage <- do.call(rbind, lapply(names(axes), function(a) {
    genes <- sets$gene_symbol[sets$axis == a]
    matched <- sum(genes %in% symbols)
    data.frame(axis = a, gs_name = axes[[a]], expected_genes = length(genes),
      matched_genes = matched, missing_genes = length(genes) - matched,
      coverage_fraction = matched / length(genes),
      passed = matched >= 20L && matched / length(genes) >= 0.80)
  }))
  missing <- sets[!sets$gene_symbol %in% symbols, , drop = FALSE]
  write_tsv(coverage, path("results/stage03_pathways/pathway_coverage.tsv"))
  write_tsv(missing, path("results/stage03_pathways/pathway_missing_genes.tsv"))
  if (any(!coverage$passed)) stage03_stop("Pathway coverage failed")
  coverage
}
score_pathways <- function(m, sets, identity) {
  if (!requireNamespace("singscore", quietly = TRUE) ||
      as.character(BiocManager::version()) != "3.22") stage03_stop("Compatible singscore unavailable")
  ranks <- singscore::rankGenes(m, tiesMethod = "min")
  if (!identical(dimnames(ranks), dimnames(m))) stage03_stop("Rank identity changed")
  output <- identity[match(colnames(m), identity$sample_id), , drop = FALSE]
  for (a in names(axes)) {
    genes <- intersect(sets$gene_symbol[sets$axis == a], rownames(m))
    score <- singscore::simpleScore(ranks, upSet = genes, centerScore = TRUE, knownDirection = TRUE)
    if (!identical(rownames(score), colnames(m))) stage03_stop("Score identity changed")
    output[[a]] <- score$TotalScore
  }
  assert_blinded(output)
  provenance <- list(timestamp = utc_now(), singscore_version = as.character(packageVersion("singscore")),
    bioconductor_version = as.character(BiocManager::version()), R_version = R.version.string,
    expression_sha256 = sha256(path("data/raw/geo/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz")),
    frozen_gene_sets_sha256 = sha256(path("metadata/pathways/confirmatory_gene_sets_v1.tsv")),
    gene_universe = nrow(m), samples = ncol(m), matching_key = "HGNC_approved_symbol",
    rankGenes = list(tiesMethod = "min", stable_genes = "not supplied"),
    simpleScore = list(upSet = "frozen genes matched by HGNC symbol", downSet = "not supplied",
      centerScore = TRUE, knownDirection = TRUE, score_column = "TotalScore"),
    expression_transformation = "none", patient_score_standardization = "none",
    missing_entrez_genes_retained = 617L,
    documentation = "https://bioconductor.org/packages/3.22/bioc/vignettes/singscore/inst/doc/singscore.html")
  jsonlite::write_json(provenance, path("metadata/pathways/scoring_provenance.json"),
    auto_unbox = TRUE, pretty = TRUE)
  output
}
pair_pathway_scores <- function(scores, pairs) {
  b <- match(pairs$baseline_sample_id, scores$sample_id)
  d <- match(pairs$d1c2_sample_id, scores$sample_id)
  if (anyNA(b) || anyNA(d) || anyDuplicated(b) || anyDuplicated(d) ||
      any(scores$patient_id[b] != pairs$patient_id | scores$patient_id[d] != pairs$patient_id) ||
      any(scores$timepoint[b] != "Baseline" | scores$timepoint[d] != "D1C2"))
    stage03_stop("Paired score identity failed", "STOP_REVIEW")
  output <- pairs[, c("patient_id", "arm")]
  for (a in names(axes)) {
    output[[paste0("baseline_", a)]] <- scores[[a]][b]
    output[[paste0("d1c2_", a)]] <- scores[[a]][d]
    output[[paste0("delta_", a)]] <- scores[[a]][d] - scores[[a]][b]
  }
  if (nrow(output) != 150L || anyDuplicated(output$patient_id) ||
      any(!is.finite(as.matrix(output[, -(1:2)]))) || !identical(output$arm, pairs$arm))
    stage03_stop("Paired score validation failed")
  assert_blinded(output)
  output
}
link_pathway_composition <- function(pair_scores, pairs) {
  composition <- read_tsv(path("data/processed/composition/neotrip_pair_composition_change.tsv"))
  i <- match(pairs$patient_id, composition$patient_id)
  if (nrow(composition) != 150L || anyDuplicated(composition$patient_id) || anyNA(i) ||
      any(composition$baseline_sample_id[i] != pairs$baseline_sample_id) ||
      any(composition$d1c2_sample_id[i] != pairs$d1c2_sample_id))
    stage03_stop("Composition linkage incomplete")
  output <- pair_scores[, c("patient_id", "arm", paste0("baseline_", names(axes)), paste0("delta_", names(axes)))]
  for (a in c("cancer_epithelial", "stromal_fraction", "immune_fraction")) {
    for (tp in c("baseline", "delta")) {
      field <- paste0(tp, "_", a)
      output[[field]] <- as.numeric(composition[[field]][i])
    }
  }
  if (nrow(output) != 150L || anyDuplicated(output$patient_id) ||
      any(!is.finite(as.matrix(output[, -(1:2)])))) stage03_stop("Composition values incomplete")
  assert_blinded(output)
  output
}

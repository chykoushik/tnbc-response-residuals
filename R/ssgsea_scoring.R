ssgsea06 <- function(m, context) {
  sets <- read06("metadata/pathways/confirmatory_gene_sets_v1.tsv", blinded = TRUE)
  genes <- setNames(lapply(axes06, function(a) sets$gene_symbol[sets$axis == a]), axes06)
  param <- GSVA::ssgseaParam(m, genes, alpha = 0.25, normalize = TRUE, minSize = 1, maxSize = Inf)
  scored <- GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
  assert06(identical(colnames(scored), colnames(m)) && setequal(rownames(scored), axes06) &&
    all(is.finite(scored)), "ssGSEA score identity invalid")
  identity <- context$identity[match(colnames(m), context$identity$sample_id), ]
  assert06(identical(identity$sample_id, colnames(m)), "ssGSEA sample mapping failed")
  output <- cbind(identity, as.data.frame(t(scored[axes06, , drop = FALSE])))
  write06(output, "data/processed/sensitivity/ssgsea_confirmatory_scores_v1.tsv.gz")
  jsonlite::write_json(list(GSVA_version = as.character(packageVersion("GSVA")),
    alpha = 0.25, normalize = TRUE, minSize = 1, maxSize = "Inf", gene_matching = "HGNC symbols",
    expression_transformation = "none", genes = nrow(m), samples = ncol(m),
    negative_values_retained = sum(m < 0),
    frozen_membership_sha256 = sha256(path("metadata/pathways/confirmatory_gene_sets_v1.tsv"))),
    path("metadata/sensitivity/ssgsea_provenance_v1.json"), auto_unbox = TRUE, pretty = TRUE)
  output
}

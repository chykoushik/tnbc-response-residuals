counts09 <- function(gsm) {
  folder <- path("data/raw/single_cell/GSE263995")
  find <- function(pattern) {
    f <- list.files(folder, pattern = paste0("^", gsm, ".*", pattern), full.names = TRUE)
    assert09(length(f) == 1, "Sample matrix file ambiguous"); f
  }
  features <- read.delim(gzfile(find("(features|genes)[.]tsv.gz$")), header = FALSE)
  barcodes <- read.delim(gzfile(find("barcodes[.]tsv.gz$")), header = FALSE)[[1]]
  con <- gzfile(find("matrix[.]mtx.gz$")); on.exit(try(close(con), silent = TRUE))
  m <- as(Matrix::readMM(con), "CsparseMatrix")
  assert09(nrow(m) == nrow(features) && ncol(m) == length(barcodes) &&
    all(m@x >= 0 & m@x == floor(m@x)), "Counts dimensions or values invalid")
  genes <- as.character(features[[2]])
  if (anyDuplicated(genes)) {
    unique_genes <- unique(genes)
    collapse <- Matrix::sparseMatrix(i = match(genes, unique_genes), j = seq_along(genes), x = 1,
      dims = c(length(unique_genes),length(genes)))
    m <- collapse %*% m; genes <- unique_genes
  }
  rownames(m) <- genes; colnames(m) <- barcodes
  m
}
process09 <- function() {
  eligible <- read09("results/stage09_final/single_cell_eligibility.tsv")
  if (!eligible$eligible[eligible$accession == "GSE263995"]) return(invisible(NULL))
  manifest <- read09("metadata/single_cell/GSE263995_tumor_manifest.tsv")
  specification <- list(seed = 319660L, cell_qc = "scuttle quickPerCellQC independently by sample; default 3 MAD",
    mitochondria_regex = "^MT-", doublets = "scDblFinder per sample with default parameters and SerialParam",
    normalization = "Seurat LogNormalize scale.factor 10000; no integration or regression",
    annotation = "SingleR HumanPrimaryCellAtlasData label.main; pruned NA mapped to other",
    minimum_compartment_cells = 50L, scoring = "UCell maxRank 1500; average ties; missing genes imputed zero",
    attribution = "Highest median of evaluable patient-compartment medians at both timepoints; require at least two evaluable patients per compartment and concordant top compartment. Different top compartments are MULTICOMPARTMENT; absent sufficient support or substantial scoring discordance is INCONCLUSIVE.",
    substantial_discordance = "More than half of evaluable paired patient-compartment directions disagree between UCell and pseudobulk singscore",
    pseudobulk = "sum raw counts; edgeR TMM across profiles; logCPM prior.count 2; singscore positive-set scores",
    no_clinical_response = TRUE, no_significance_tests = TRUE)
  spec_file <- path("metadata/single_cell/stage09_sc_spec_v1.yml")
  if (!file.exists(spec_file)) {
    yaml::write_yaml(specification, spec_file); writeLines(sha256(spec_file), paste0(spec_file, ".sha256"))
  } else assert09(identical(yaml::read_yaml(spec_file), specification) &&
    identical(readLines(paste0(spec_file, ".sha256")), sha256(spec_file)), "Single-cell specification changed")
  qc_rows <- dbl_rows <- list()
  for (i in seq_len(nrow(manifest))) {
    info <- manifest[i, ]; gsm <- info$GSM
    cache <- path("models/stage09", paste0(gsm, "_qc.rds"))
    if (file.exists(cache)) {
      cached <- readRDS(cache); qc_rows[[gsm]] <- cached$qc; dbl_rows[[gsm]] <- cached$doublet
      next
    }
    m <- counts09(gsm)
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = m))
    qc <- scuttle::perCellQCMetrics(sce, subsets = list(Mito = grepl("^MT-", rownames(sce))))
    filter <- scuttle::quickPerCellQC(qc, sub.fields = "subsets_Mito_percent")
    row <- data.frame(GSM = gsm, patient_id = info$patient_id, timepoint = info$timepoint,
      cells_raw = ncol(sce), cells_after_qc = sum(!filter$discard), cells_removed = sum(filter$discard),
      median_features = median(qc$detected), median_counts = median(qc$sum), median_mito = median(qc$subsets_Mito_percent))
    write09(data.frame(cell_barcode = colnames(sce), as.data.frame(qc), as.data.frame(filter)),
      paste0("metadata/single_cell/", gsm, "_cell_qc.tsv"))
    sce <- sce[, !filter$discard]
    set.seed(319660)
    sce <- scDblFinder::scDblFinder(sce, BPPARAM = BiocParallel::SerialParam(RNGseed = 319660))
    doublet <- data.frame(GSM = gsm, cells_after_qc = ncol(sce),
      predicted_doublets = sum(sce$scDblFinder.class == "doublet"),
      retained_singlets = sum(sce$scDblFinder.class == "singlet"), seed = 319660)
    write09(data.frame(cell_barcode = colnames(sce), doublet_class = sce$scDblFinder.class,
      doublet_score = sce$scDblFinder.score), paste0("metadata/single_cell/",gsm,"_doublet_calls.tsv"))
    sce <- sce[, sce$scDblFinder.class == "singlet"]
    seurat <- Seurat::CreateSeuratObject(counts = SummarizedExperiment::assay(sce, "counts"),
      min.cells = 0, min.features = 0)
    seurat <- Seurat::NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    normalized <- SeuratObject::LayerData(seurat, layer = "data")
    assert09(identical(rownames(normalized), rownames(sce)) && identical(colnames(normalized),colnames(sce)),
      "Normalization identity mismatch")
    SummarizedExperiment::assay(sce, "logcounts") <- normalized
    saveRDS(list(sce = sce, qc = row, doublet = doublet), cache)
    qc_rows[[gsm]] <- row; dbl_rows[[gsm]] <- doublet
    cat("QC complete", gsm, ncol(sce), "singlets\n")
    rm(m, sce, seurat, normalized); gc()
  }
  write09(do.call(rbind,qc_rows), "results/stage09_final/single_cell_qc_by_sample.tsv")
  write09(do.call(rbind,dbl_rows), "results/stage09_final/single_cell_doublet_qc.tsv")
}

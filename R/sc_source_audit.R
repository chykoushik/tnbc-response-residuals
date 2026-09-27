dirs09 <- c("data/raw/single_cell", "data/processed/single_cell", "metadata/single_cell",
  "results/stage09_final", "figures/final", "tables/final", "models/stage09")
read09 <- function(f) read.delim(path(f), check.names = FALSE, stringsAsFactors = FALSE)
write09 <- function(x, f) write_tsv(x, path(f))
json09 <- function(f) jsonlite::read_json(path(f), simplifyVector = TRUE)
assert09 <- function(x, reason) if (!isTRUE(x)) stop(reason, call. = FALSE)
init09 <- function() for (d in dirs09) dir.create(path(d), recursive = TRUE, showWarnings = FALSE)
prerequisite09 <- function() {
  verify08()
  m <- read_tsv(path("results/stage08_final_gate/stage08_output_checksums.tsv"))
  assert09(identical(unname(vapply(path(m$local_relative_path), sha256, "")), m$sha256), "Stage 8 checksum failure")
  s <- json09("results/stage08_final_gate/stage08_status.json")
  assert09(s$decision == "PROCEED_NEOTRIP_ONLY" && s$patients == 150 &&
    s$confirmatory_pathways == 6 && s$multiplicity_confirmed_pathways == 0 &&
    s$full_transportability_eligible_pathways == 0 && s$GSE122630_status == "SKIP" &&
    !s$consistent_incremental_prediction_gain && !s$fundamental_stop_condition &&
    file.exists(path("metadata/final_gate/analysis_scope_freeze_v1.yml")) &&
    file.exists(path("metadata/final_gate/claim_boundaries_v1.tsv")), "Stage 8 prerequisite failed")
}
soft09 <- function(acc) {
  lines <- readLines(path("metadata/single_cell", paste0(acc, "_family.soft")), warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE =", lines)
  blocks <- lapply(seq_along(starts), function(i) lines[starts[i]:(if (i < length(starts)) starts[i+1] - 1 else length(lines))])
  value <- function(b, key) sub(paste0("^", key, " = "), "", grep(paste0("^", key, " = "), b, value = TRUE))
  samples <- do.call(rbind, lapply(blocks, function(b) data.frame(GSM = sub("^\\^SAMPLE = ", "", b[1]),
    title = value(b, "!Sample_title"), tissue = value(b, "!Sample_source_name_ch1"),
    timepoint_evidence = paste(value(b, "!Sample_characteristics_ch1")[grepl("treatment:", value(b, "!Sample_characteristics_ch1"))], collapse = ";"))))
  list(lines = lines, samples = samples, value = value)
}
audit_sources09 <- function() {
  prerequisite09()
  download <- read09("metadata/single_cell/download_manifest.tsv")
  assert09(all(file.exists(path(download$local_relative_path))) &&
    identical(unname(vapply(path(download$local_relative_path), sha256, "")), download$sha256), "Single-cell download integrity failed")
  a <- soft09("GSE263995"); b <- soft09("GSE289825")
  expected <- data.frame(GSM = c("GSM8208510","GSM8208511","GSM8391107","GSM8391110","GSM8391113","GSM8391116"),
    patient_id = rep(c("P0102","P02","P03"), each = 2), timepoint = rep(c("Pre","Post"), 3))
  observed <- a$samples[match(expected$GSM, a$samples$GSM), ]
  assert09(nrow(a$samples) == 14 && !anyNA(observed) &&
    all(observed$tissue == "triple-negative breast cancer") &&
    all(mapply(grepl, expected$patient_id, observed$title, MoreArgs = list(fixed = TRUE))) &&
    all(mapply(grepl, tolower(expected$timepoint), observed$timepoint_evidence, MoreArgs = list(fixed = TRUE))),
    "MODIFY_REVIEW: official tumour sample mapping differs")
  manifest <- cbind(expected, observed[, -1])
  manifest$repository_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", manifest$GSM)
  write09(manifest, "metadata/single_cell/GSE263995_tumor_manifest.tsv")
  write09(a$samples, "metadata/single_cell/GSE263995_all_sample_metadata.tsv")
  write09(b$samples, "metadata/single_cell/GSE289825_sample_metadata.tsv")
  sources <- data.frame(accession = c("GSE263995","GSE289825"), public_status = c(
    a$value(a$lines,"!Series_status"), b$value(b$lines,"!Series_status")),
    repository_URL = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",c("GSE263995","GSE289825")),
    platform = c("GPL24676; Illumina NovaSeq 6000", "GPL30173; NextSeq 2000"),
    whole_transcriptome_or_targeted = c("whole_transcriptome","BD Rhapsody Human Immune Response Targeted Panel"),
    patients = c(3,2), sample_identifiers = c(paste(a$samples$GSM,collapse=";"),paste(b$samples$GSM,collapse=";")),
    tissue = c("tumour; TdLN; PBMC (tumour only analysed)", "PBMC and tumour"),
    timepoint = c("Pre; Post", "before NAC; day 3; day 21"),
    processed_matrices_available = TRUE, raw_sequencing_available = TRUE, public_without_request = TRUE,
    download_decision = c("SIX_TUMOUR_SAMPLES_ONLY", "METADATA_AND_NINE_PROCESSED_MATRICES_FOR_COVERAGE"),
    analysis_role = "Optional orthogonal cell-type interpretation; no clinical validation")
  write09(sources, "metadata/single_cell/single_cell_source_audit.tsv")
  inv <- read_tsv(path("metadata/data_inventory.tsv"))
  add <- data.frame(source = "NCBI GEO", url = download$url, retrieval_timestamp = download$download_date,
    local_relative_path = download$local_relative_path, file_size_bytes = download$file_size_bytes,
    sha256 = download$sha256, accession_or_repository = ifelse(nzchar(download$GSM),download$GSM,download$accession),
    public_access_status = "public_unauthenticated", repository_commit = "", acquisition_status = "downloaded",
    checksum_timestamp = download$download_date)
  write_tsv(rbind(inv, add[!add$local_relative_path %in% inv$local_relative_path, ]), path("metadata/data_inventory.tsv"))
  sets <- read09("metadata/pathways/confirmatory_gene_sets_v1.tsv")
  features <- list()
  for (gsm in manifest$GSM) {
    f <- list.files(path("data/raw/single_cell/GSE263995"), pattern = paste0("^",gsm,".*(features|genes).*gz$"), full.names = TRUE)
    tab <- read.delim(gzfile(f), header = FALSE)
    features[[gsm]] <- unique(as.character(tab[[2]]))
  }
  zips <- list.files(path("data/raw/single_cell/GSE289825"), pattern = "[.]zip$", full.names = TRUE)
  targeted <- lapply(zips, function(f) {
    members <- unzip(f, list = TRUE)$Name
    member <- members[grepl("_DBEC_MolsPerCell[.]csv$", members)]
    assert09(length(member) == 1, "Targeted archive feature schema ambiguous")
    con <- unz(f, member); on.exit(close(con), add = TRUE)
    header <- readLines(con, n = 8)
    assert09(startsWith(header[8], "Cell_Index,"), "Targeted gene header differs")
    strsplit(header[8], ",", fixed = TRUE)[[1]][-1]
  })
  # Common measured features are required across samples; no expression-based filtering.
  available <- list(GSE263995 = Reduce(intersect, features), GSE289825 = Reduce(intersect, targeted))
  coverage <- do.call(rbind, lapply(names(available), function(acc) do.call(rbind, lapply(axes07, function(axis) {
    genes <- unique(sets$gene_symbol[sets$axis == axis]); n <- sum(genes %in% available[[acc]])
    data.frame(accession = acc, axis = axis, expected_genes = length(genes), available_genes = n,
      coverage_fraction = n / length(genes))
  }))))
  write09(coverage, "results/stage09_final/single_cell_gene_coverage.tsv")
  eligibility <- data.frame(accession = names(available), eligible = vapply(names(available),
    function(a) all(coverage$coverage_fraction[coverage$accession == a] >= .8), logical(1)))
  eligibility$dataset_status <- ifelse(eligibility$eligible, "ELIGIBLE", "OMIT_TECHNICAL")
  eligibility$reason <- ifelse(eligibility$eligible, "All six frozen pathways meet 80% coverage in every sample",
    "At least one frozen pathway has less than 80% common measured-gene coverage; no subset rescue")
  write09(eligibility, "results/stage09_final/single_cell_eligibility.tsv")
  write09(data.frame(claim = c("patient-level pCR validation","external residual replication",
    "clinical prediction validation","mechanistic proof","cell-type interpretation"),
    status = c(rep("NOT_PERFORMED",3),"NOT_SUPPORTED","DESCRIPTIVE_ORTHOGONAL_SUPPORT")),
    "metadata/single_cell/single_cell_claim_boundary.tsv")
}

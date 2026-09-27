source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("renv/activate.R"))
source(path("R/checksums.R"))
source(path("R/downloads.R"))
run_stage("01_download", {
  options(timeout = 1200)
  config <- yaml::read_yaml(path("config/project.yml"))
  force <- "--force" %in% commandArgs(trailingOnly = TRUE)
  commit <- repository_sha(config$repository)
  ref <- if (nzchar(commit)) commit else "main"
  base <- paste0("https://raw.githubusercontent.com/mdugo/Early_GEP_NeoTRIP/", ref)
  files <- c("GSE319641_family.soft.gz", config$expression_file,
    "NeoTRIP_baseline_D1C2_pheno_and_clinical.txt", "geneSets_Supplementary_Table1.txt")
  urls <- c(paste0(config$geo_base, "/soft/", files[1]),
    paste0(config$geo_base, "/suppl/", files[2]), paste0(base, "/data/", files[3:4]))
  relative <- c(paste0("data/raw/geo/", files[1:2]), paste0("data/raw/author_repo/", files[3:4]))
  for (i in seq_along(files)) {
    download_source(if (i <= 2L) "NCBI GEO" else "Author repository", urls[i], relative[i],
      if (i <= 2L) config$geo_accession else config$repository,
      if (i <= 2L) "" else commit, force)
  }
  notes <- data.frame(timestamp = utc_now(), source = c(config$repository, config$repository,
    config$repository, config$geo_accession),
    topic = c("clinical_availability_discrepancy", "repository_revision", "gene_sets", "expression_scale"),
    note = c(paste("The README states clinical data are not publicly available due to privacy restrictions",
      "and may be obtained upon reasonable request; the public repository contains",
      "data/NeoTRIP_baseline_D1C2_pheno_and_clinical.txt. Discrepancy recorded without interpretation."),
      if (nzchar(commit)) paste("git ls-remote HEAD:", commit) else "Commit unavailable; main URL used for new acquisitions.",
      "Downloaded for provenance and later work only; not used in Stage 1.",
      "Official processed ComBat-adjusted TPM; no additional transformation or filtering in Stage 1."))
  write_tsv(notes, path("metadata/source_notes.tsv"))
  verify_inventory()
})
cat("Download complete\n")

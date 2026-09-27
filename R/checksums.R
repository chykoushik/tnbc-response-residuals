sha256 <- function(file) digest::digest(file = file, algo = "sha256", serialize = FALSE)
verify_inventory <- function() {
  inventory <- read_tsv(path("metadata", "data_inventory.tsv"))
  required <- c("data/raw/geo/GSE319641_family.soft.gz",
    "data/raw/geo/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz",
    "data/raw/author_repo/NeoTRIP_baseline_D1C2_pheno_and_clinical.txt",
    "data/raw/author_repo/geneSets_Supplementary_Table1.txt")
  if (!all(required %in% inventory$local_relative_path) || anyDuplicated(inventory$local_relative_path)) {
    stop("Invalid source inventory")
  }
  for (i in seq_len(nrow(inventory))) {
    file <- path(inventory$local_relative_path[i])
    if (!file.exists(file) || sha256(file) != inventory$sha256[i]) {
      stop("Source checksum mismatch")
    }
  }
  invisible(inventory)
}

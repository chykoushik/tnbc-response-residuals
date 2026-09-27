source("D:/TNBC-ResponseResiduals/R/paths.R")
source(path("R/manifest_checks.R"))
source(path("R/geo_soft.R"))
source(path("R/checksums.R"))
source(path("R/downloads.R"))
config <- yaml::read_yaml(path("config/project.yml"))
test_root <- tempfile("stage01-tests-", tmpdir = path("external/tmp"))
dir.create(test_root)
root <- test_root
invisible(lapply(c("metadata", "results/stage01_manifest", "external/tmp"),
  function(x) dir.create(path(x), recursive = TRUE)))
author <- data.frame(Sample_ID = c("A_B", "A_D", "B_B", "B_D"),
  Patient_ID = rep(c("A", "B"), each = 2), Timepoint = rep(c("Baseline", "D1C2"), 2),
  Arm = rep(c("CT", "CT/A"), each = 2), pCR = rep(c("pCR", "RD"), each = 2),
  pCR.num = rep(c("1", "0"), each = 2), Early.pCR = "RD", TNBCtype = "M", Batch_correlation = "Batch 1")
geo <- data.frame(gsm = paste0("GSM", 1:4), sample_title = "test", geo_library_name = author$Sample_ID,
  geo_patient_id = author$Patient_ID, geo_timepoint = author$Timepoint,
  geo_batch = author$Batch_correlation, organism = "Homo sapiens", platform = "GPL1")
expression <- list(samples = author$Sample_ID, audit = data.frame(check = character(), count = integer()),
                   unmatched_expression = character(), unmatched_phenotype = character(), fatal = FALSE)
run <- function(g = geo, a = author, e = expression) build_manifest(g, a, e, config)
result <- run()
stopifnot(result$status$decision == "PROCEED", result$status$paired_complete == 2L)
a <- author
a$Arm[2] <- "CT/A"
result <- run(a = a)
stopifnot(result$status$decision == "MODIFY", result$status$paired_complete == 1L,
  result$pairs$pair_status[1] == "invalid_metadata", result$pairs$arm[1] == "")
a <- author
a$pCR[2] <- "RD"
a$pCR.num[2] <- "0"
result <- run(a = a)
stopifnot(result$status$paired_complete == 1L, result$pairs$pcr[1] == "")
a <- author
a$Sample_ID[2] <- a$Sample_ID[1]
stopifnot(run(a = a)$status$decision == "STOP")
g <- geo
g$geo_patient_id[1] <- "C"
result <- run(g = g)
stopifnot(result$status$mapping_failures == 1L, result$status$paired_complete == 1L,
  result$mapping$patient_id[1] == "", "patient_conflict" %in% result$issues$issue)
g <- geo
g$geo_timepoint[1] <- "baseline"
result <- run(g = g)
stopifnot("timepoint_conflict" %in% result$issues$issue, result$mapping$patient_id[1] == "")
g <- geo
g$geo_batch[1] <- "Batch 2"
stopifnot("batch_conflict" %in% run(g = g)$issues$issue)
g <- geo
g$geo_library_name[1] <- "similar_A_B"
result <- run(g = g)
stopifnot(all(c("missing_geo", "missing_author") %in% result$mapping$mapping_status))
a <- author
a$pCR <- "NA"
stopifnot(run(a = a)$status$decision == "STOP")
a <- author
a$pCR[1] <- "unexpected"
stopifnot("invalid_pcr" %in% run(a = a)$issues$issue)
g <- geo
g$gsm[2] <- g$gsm[1]
stopifnot(run(g = g)$status$decision == "STOP")
a <- author
a$Timepoint[2] <- "Baseline"
stopifnot(run(a = a)$status$decision == "STOP")
fixture <- path("expression.txt.gz")
connection <- gzfile(fixture, "wt")
writeLines(c(paste(c("HGNC_approved_symbol", "HGNC_id", "entrez_gene_id", "A_B", "A_B", "", "OTHER"), collapse = "\t"),
  "X\tHGNC:1\t1\t1\t2\t3\t4", "X\tHGNC:1\t1\t2\t3\t4\t5"), connection)
close(connection)
checked <- expression_check(fixture, author$Sample_ID)
stopifnot(checked$fatal, "OTHER" %in% checked$unmatched_expression,
  "A_D" %in% checked$unmatched_phenotype,
  checked$audit$count[checked$audit$check == "duplicate_HGNC_approved_symbol"] == 1L)
fixture <- path("sample.soft.gz")
connection <- gzfile(fixture, "wt")
writeLines(c("^SAMPLE = GSM1", "!Sample_title = patient A at Baseline", "!Sample_description = Library name: A_B",
  "!Sample_characteristics_ch1 = patient id: A", "!Sample_characteristics_ch1 = timepoint: Baseline",
  "!Sample_characteristics_ch1 = batch: Batch 1", "^SAMPLE = GSM2", "!Sample_title = patient B at D1C2",
  "!Sample_description = Library name: B_D"), connection)
close(connection)
parsed <- parse_geo_soft(fixture)
stopifnot(parsed$geo_patient_id[2] == "", parsed$geo_timepoint[2] == "", parsed$geo_library_name[2] == "B_D")
writeLines("unchanged", path("source.txt"))
before <- sha256(path("source.txt"))
invisible(download_source("test", "https://invalid.invalid/source", "source.txt", "fixture"))
inventory <- read_tsv(path("metadata/data_inventory.tsv"))
stopifnot(inventory$retrieval_timestamp == "", inventory$sha256 == before)
invisible(download_source("test", "https://invalid.invalid/source", "source.txt", "fixture"))
stopifnot(sha256(path("source.txt")) == before)
cat("Checks complete\n")

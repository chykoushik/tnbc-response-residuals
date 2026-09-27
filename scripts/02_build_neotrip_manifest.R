source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("renv/activate.R"))
source(path("R/checksums.R"))
source(path("R/geo_soft.R"))
source(path("R/manifest_checks.R"))
run_stage("02_manifest", {
  config <- yaml::read_yaml(path("config/project.yml"))
  result <- tryCatch({
    verify_inventory()
    geo <- parse_geo_soft(path("data/raw/geo/GSE319641_family.soft.gz"))
    write_tsv(geo, path("metadata/neotrip_geo_samples.tsv"))
    author <- read_tsv(path("data/raw/author_repo/NeoTRIP_baseline_D1C2_pheno_and_clinical.txt"))
    required <- c("Sample_ID", "Patient_ID", "Timepoint", "Arm", "pCR", "pCR.num", "Early.pCR", "TNBCtype", "Batch_correlation")
    if (!all(required %in% names(author)) || anyDuplicated(names(author))) stop("Phenotype schema changed")
    expression <- expression_check(path("data/raw/geo", config$expression_file), author$Sample_ID)
    write_tsv(expression$audit, path("results/stage01_manifest/expression_audit.tsv"))
    build_manifest(geo, author, expression, config)
  }, error = function(e) {
    issues <- empty_issues()
    issues[1L, ] <- c("core_input_failure", "fatal", "", "", "", conditionMessage(e))
    write_tsv(issues, path("results/stage01_manifest/metadata_issues.tsv"))
    status <- list(stage = 1L, timestamp = utc_now(), geo_accession = config$geo_accession,
      geo_samples = NA_integer_, expression_samples = NA_integer_, phenotype_rows = NA_integer_,
      unique_patients = NA_integer_, paired_patients = NA_integer_, paired_with_arm = NA_integer_,
      paired_with_pcr = NA_integer_, paired_complete = NA_integer_, mapping_failures = NA_integer_,
      metadata_conflicts = NA_integer_, decision = "STOP", decision_reason = conditionMessage(e))
    list(status = status)
  })
  save_status(result$status)
  software_versions()
  if (result$status$decision == "STOP") stop("Core identity checks failed")
  renv::snapshot(project = root, packages = c("renv", "digest", "jsonlite", "yaml"), prompt = FALSE)
})
cat(path("results/stage01_manifest/stage01_status.json"), "\n", sep = "")

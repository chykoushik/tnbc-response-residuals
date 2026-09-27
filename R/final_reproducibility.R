stage_manifests09 <- c("metadata/composition/stage01_input_checksums.tsv",
  "results/stage02_composition/stage02_output_checksums.tsv","results/stage03_pathways/stage03_output_checksums.tsv",
  "results/stage04_residuals/stage04_output_checksums.tsv","results/stage05_primary/stage05_output_checksums_v2.tsv",
  "results/stage06_robustness/stage06_output_checksums.tsv","results/stage07_prediction/stage07_output_checksums.tsv",
  "results/stage08_final_gate/stage08_output_checksums.tsv")
stage_status09 <- c("results/stage01_manifest/stage01_status.json","results/stage02_composition/stage02_status.json",
  "results/stage03_pathways/stage03_status.json","results/stage04_residuals/stage04_status.json",
  "results/stage05_primary/stage05_status_v2.json","results/stage06_robustness/stage06_status.json",
  "results/stage07_prediction/stage07_status.json","results/stage08_final_gate/stage08_status.json")
manifest_check09 <- function(file) {
  m<-read_tsv(path(file)); field<-intersect(c("local_relative_path","path"),names(m))[1]
  assert09(!is.na(field) && "sha256"%in%names(m),"Checksum manifest schema unrecognized")
  actual_path<-m[[field]]
  replacement<-actual_path=="renv.lock" &
    file.exists(path("metadata/single_cell/pre_stage09_renv.lock"))
  actual_path[replacement]<-"metadata/single_cell/pre_stage09_renv.lock"
  exists<-file.exists(path(actual_path))
  current<-vapply(seq_along(actual_path),function(i)if(exists[i])sha256(path(actual_path[i])) else "MISSING","")
  data.frame(manifest=file,original_path=m[[field]],verified_path=actual_path,expected_sha256=m$sha256,
    current_sha256=current,passed=exists & current==m$sha256,
    applicability=ifelse(replacement,"Prior environment preserved byte-for-byte; root lock superseded by explicitly requested Stage 9 snapshot","CURRENT_FROZEN_FILE"))
}
environment09 <- function() {
  provenance <- path("metadata/single_cell/environment_snapshot_provenance.tsv")
  if(file.exists(provenance)) {
    record<-read_tsv(provenance)
    assert09(sha256(path("renv.lock"))==record$final_sha256 &&
      sha256(path(record$historical_copy))==record$historical_sha256,"Frozen environment snapshot changed")
    versions<-read09("results/stage09_final/final_software_versions.tsv")
    r_versions<-versions[!versions$package%in%c("Python","PyMuPDF"),]
    assert09(all(vapply(r_versions$package,function(p)as.character(packageVersion(p)),"")==r_versions$version),
      "Package version changed after final snapshot")
    return(invisible(NULL))
  }
  Sys.setenv(RENV_DOWNLOAD_METHOD="libcurl")
  options(timeout=60)
  previous<-path("metadata/single_cell/pre_stage09_renv.lock")
  if(!file.exists(previous)) {
    old<-read_tsv(path("metadata/pathways/implementation_checksums.tsv"))
    assert09(sha256(path("renv.lock"))==old$sha256[old$local_relative_path=="renv.lock"],
      "Historical lock checksum failed before final snapshot")
    file.copy(path("renv.lock"),previous)
  }
  # Include the historical analysis package set and all packages used by Stage 9.
  old<-renv::lockfile_read(path("metadata/prediction/stage07_environment_v1.lock"))
  packages<-unique(c(names(old$Packages),"Seurat","SeuratObject","SingleCellExperiment","scuttle","scDblFinder",
    "UCell","SingleR","celldex","edgeR","singscore","Matrix","BiocParallel"))
  renv::snapshot(project=root,lockfile=path("renv.lock"),packages=packages,prompt=FALSE)
  lock<-renv::lockfile_read(path("renv.lock"))
  packages<-names(lock$Packages)
  versions<-data.frame(package=packages,version=vapply(packages,function(p)as.character(packageVersion(p)),""),
    R_version=R.version.string,OS=Sys.info()[["sysname"]],OS_release=Sys.info()[["release"]],platform=R.version$platform,
    role="Actual analysis package or runtime dependency; union of stage environment records")
  versions<-rbind(versions,read09("metadata/single_cell/auxiliary_software.tsv"))
  write09(versions,"results/stage09_final/final_software_versions.tsv")
  write09(data.frame(path="renv.lock",historical_copy="metadata/single_cell/pre_stage09_renv.lock",
    historical_sha256=sha256(previous),final_sha256=sha256(path("renv.lock")),
    authorization="Stage 9 explicitly requires final renv snapshot; earlier checksum remains applicable to preserved historical bytes"),
    "metadata/single_cell/environment_snapshot_provenance.tsv")
}
claims09 <- function() {
  boundaries<-read09("metadata/final_gate/claim_boundaries_v1.tsv")
  assert09(nrow(boundaries)==13 && sum(boundaries$status=="NOT_SUPPORTED")==7,"Stage 8 claim boundary mismatch")
  complete<-read09("results/stage09_final/single_cell_completion.tsv")$status=="COMPLETED"
  claims<-data.frame(claim=c("Cross-fitted residual response technically measurable",
    "Residual associations after composition adjustment","Multiplicity-confirmed residual pathway",
    "Independent replication","GSE122630 transportability","Incremental predictive improvement","Clinical utility",
    "Orthogonal cell-type interpretation","Mechanistic proof","Tumour-intrinsic residual biology",
    "Residuals contain structured biological variation","Composition adjustment materially affects residual associations",
    "Findings are single-cohort and require independent validation"),
    status=c("SUPPORTED","HYPOTHESIS_GENERATING_ONLY","NOT_SUPPORTED","NOT_SUPPORTED","NOT_PERFORMED",
      "NOT_SUPPORTED","NOT_SUPPORTED",if(complete)"SUPPORTED_DESCRIPTIVELY" else "NOT_PERFORMED_TECHNICAL",
      "NOT_SUPPORTED","NOT_SUPPORTED_WITHOUT_QUALIFICATION","HYPOTHESIS_GENERATING_ONLY",
      "SUPPORTED_DESCRIPTIVELY","SUPPORTED"),
    boundary_check="PASS",qualification=c("Technical measurement only","Single cohort; no multiplicity-confirmed pathway",
      "Zero of six pass frozen multiplicity rule","No independent validation cohort","Skipped by eligibility gate",
      "No consistent mean improvement","DCA not run","Three patients; no clinical outcome analysis or cell-level inference",
      "Descriptive associations only","Bulk expression and composition dependence",
      "Structured residual variation is descriptive; cell-type scores are not residual replication or mechanistic proof",
      "Frozen adjusted versus base-model coefficients differ; no new hypothesis test",
      "Only NeoTRIP supplies clinical residual inference; single-cell support is orthogonal and descriptive"))
  write09(claims,"results/stage09_final/final_claim_audit.tsv")
}
audit09 <- function() {
  checks<-do.call(rbind,lapply(unique(c(stage_manifests09,
    "results/stage05_primary/stage05_output_checksums.tsv",
    "metadata/modeling/stage04_input_checksums.tsv",
    "metadata/single_cell/source_metadata_checksums.tsv",
    "metadata/pathways/implementation_checksums.tsv")),manifest_check09))
  sidecars<-list.files(path("metadata"),pattern="[.]sha256$",recursive=TRUE,full.names=TRUE)
  side_results<-lapply(sidecars,function(f) {
    stem<-sub("[.]sha256$","",f); candidates<-c(stem,paste0(stem,c(".yml",".tsv",".tsv.gz")))
    target<-candidates[file.exists(candidates)][1];expected<-readLines(f,warn=FALSE)[1]
    data.frame(manifest=substring(f,nchar(root)+2),original_path=substring(target,nchar(root)+2),
      verified_path=substring(target,nchar(root)+2),expected_sha256=expected,
      current_sha256=if(!is.na(target))sha256(target) else "MISSING",
      passed=!is.na(target) && sha256(target)==expected,applicability="CURRENT_FROZEN_FILE")
  })
  checks<-rbind(checks,do.call(rbind,side_results))
  write09(checks,"results/stage09_final/frozen_checksum_validation.tsv")
  inv<-read_tsv(path("metadata/data_inventory.tsv"))
  raw<-inv[nzchar(inv$sha256),]
  raw$current_sha256<-vapply(path(raw$local_relative_path),function(f)if(file.exists(f))sha256(f) else "MISSING","")
  raw$passed<-raw$current_sha256==raw$sha256
  write09(raw,"results/stage09_final/raw_data_integrity.tsv")
  checks_existing <- c(
    stage04_regression = all(read09("results/stage04_residuals/stage04_regression_tests.tsv")$passed),
    stage05_saved_models = all(read09("results/stage05_primary/model_reproduction_check_v2.tsv")$passed),
    stage06_crossfit_reproduction = all(read09("results/stage06_robustness/crossfit_reproduction_validation.tsv")$passed),
    stage07_saved_models = all(read09("results/stage07_prediction/model_reproduction_check.tsv")$passed),
    stage07_no_leakage = all(read09("results/stage07_prediction/leakage_audit.tsv")$violations == 0))
  exports<-read09("results/stage09_final/figure_export_validation.tsv")
  assert09(nrow(exports)==6 && all(toupper(as.character(exports$passed))=="TRUE") &&
    all(exports$png_dpi>=300) && all(exports$raster_images==0),"Final figure export validation failed")
  write09(data.frame(check=names(checks_existing),passed=unname(checks_existing),
    method="Read and verify frozen validation results; no model fitting or statistic recomputation"),
    "results/stage09_final/existing_reproduction_checks.tsv")
  assert09(all(checks_existing),"Existing reproduction or leakage verification failed")
  if(read09("results/stage09_final/single_cell_completion.tsv")$status=="COMPLETED") {
    cells<-read09("data/processed/single_cell/GSE263995_ucell_scores.tsv.gz")
    qc<-read09("results/stage09_final/single_cell_doublet_qc.tsv")
    annotation<-do.call(rbind,lapply(qc$GSM,function(g)read09(paste0("metadata/single_cell/",g,"_annotation.tsv"))))
    per_sample<-table(cells$GSM)
    sc_checks<-c(expected_six_samples=nrow(qc)==6,
      singlet_cell_coverage=all(as.integer(per_sample[qc$GSM])==qc$retained_singlets),
      unique_cell_identity=!anyDuplicated(cells[,c("GSM","cell_barcode")]),
      all_six_scores=all(axes07%in%names(cells)),
      finite_bounded_scores=all(is.finite(as.matrix(cells[,axes07])) & as.matrix(cells[,axes07])>=0 & as.matrix(cells[,axes07])<=1),
      no_clinical_fields=!any(grepl("pcr|response|responder",names(cells),ignore.case=TRUE)),
      annotation_coverage=nrow(annotation)==nrow(cells))
    write09(data.frame(check=names(sc_checks),passed=unname(sc_checks)),"results/stage09_final/single_cell_validation.tsv")
    assert09(all(sc_checks),"Single-cell identity or scoring QC failed")
  }
  expected<-c("PROCEED","PROCEED_WITH_MODIFICATION",rep("PROCEED",5),"PROCEED_NEOTRIP_ONLY")
  rows<-lapply(1:8,function(i) {
    s<-json09(stage_status09[i]);c<-checks[checks$manifest==stage_manifests09[i],]
    passed<-s$decision==expected[i] && all(c$passed)
    data.frame(stage=as.character(i),status=if(passed)"PASS" else "FAIL",status_file=stage_status09[i],
      checksum_file=stage_manifests09[i],decision=s$decision,required_outputs_present=all(c$current_sha256!="MISSING"),
      checksum_validation=all(c$passed),notes="Frozen outputs verified; existing model reproduction and leakage audits retained; no primary models rerun")
  })
  stage_files<-c(list.files(path("figures/final"),full.names=TRUE),
    list.files(path("tables/final"),full.names=TRUE),
    list.files(path("data/processed/single_cell"),full.names=TRUE),
    list.files(path("metadata/single_cell"),full.names=TRUE),
    list.files(path("results/stage09_final"),full.names=TRUE))
  stage_files<-stage_files[!basename(stage_files)%in%c("stage09_status.json","stage09_output_checksums.tsv",
    "stage09_analysis_checksums.tsv","final_reproducibility_audit.tsv","final_output_registry.tsv")]
  stage_files<-stage_files[!file.info(stage_files)$isdir]
  write09(data.frame(local_relative_path=substring(stage_files,nchar(root)+2),
    file_size_bytes=file.info(stage_files)$size,sha256=vapply(stage_files,sha256,"")),
    "results/stage09_final/stage09_analysis_checksums.tsv")
  complete<-read09("results/stage09_final/single_cell_completion.tsv")$status
  stage9<-data.frame(stage="9",status=if(complete%in%c("COMPLETED","OMITTED"))"PASS" else "FAIL",
    status_file="results/stage09_final/stage09_status.json",
    checksum_file="results/stage09_final/stage09_analysis_checksums.tsv",decision="AUDIT_COMPLETE",
    required_outputs_present=length(list.files(path("figures/final"),pattern="[.](pdf|png)$"))==12 &&
      length(list.files(path("tables/final"),pattern="^Table.*tsv$"))==6,
    checksum_validation=all(manifest_check09("results/stage09_final/stage09_analysis_checksums.tsv")$passed),
    notes="Stage 9 analysis artifact manifest verified; final status and final manifest issued after this audit to avoid circular hashing")
  audit<-rbind(do.call(rbind,rows),stage9)
  write09(audit,"results/stage09_final/final_reproducibility_audit.tsv")
  assert09(all(checks$passed) && all(raw$passed) && all(audit$status=="PASS") &&
    all(audit$required_outputs_present),"STOP_REVIEW: reproducibility or integrity failure")
}
registry09 <- function() {
  folders<-c("figures/final","tables/final","results/stage05_primary","results/stage06_robustness",
    "results/stage07_prediction","results/stage08_final_gate","results/stage09_final",
    "metadata/final_gate","metadata/single_cell","metadata/pathways","metadata/modeling","metadata/outcome",
    "metadata/sensitivity","metadata/prediction","data/processed/single_cell")
  files<-unique(unlist(lapply(path(folders),list.files,recursive=TRUE,full.names=TRUE)))
  files<-files[!basename(files)%in%c("final_output_registry.tsv","stage09_output_checksums.tsv")]
  rel<-substring(files,nchar(root)+2)
  role<-ifelse(grepl("^figures/",rel),"final figure",ifelse(grepl("^tables/",rel),"final table",
    ifelse(grepl("^metadata/",rel),"metadata freeze",ifelse(grepl("stage05",rel),"primary result",
      ifelse(grepl("stage06",rel),"robustness result",ifelse(grepl("stage07",rel),"prediction result","single-cell or final audit result"))))))
  stage<-ifelse(grepl("stage05|outcome",rel),"5",ifelse(grepl("stage06|sensitivity",rel),"6",
    ifelse(grepl("stage07|prediction",rel),"7",ifelse(grepl("stage08",rel),"8",
      ifelse(grepl("metadata/modeling",rel),"4",ifelse(grepl("metadata/pathways",rel),"3","9"))))))
  write09(data.frame(artifact=basename(files),path=rel,stage=stage,role=role,
    primary_or_secondary=ifelse(role=="primary result","PRIMARY",ifelse(role%in%c("prediction result",
      "single-cell or final audit result"),"SECONDARY_OR_AUDIT","SUPPORTING")),
    checksum=vapply(files,sha256,"")),"results/stage09_final/final_output_registry.tsv")
}

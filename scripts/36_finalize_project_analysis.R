source("R/paths.R")
source("R/checksums.R")
source("R/prediction_folds.R")
source("R/sc_source_audit.R")
source("R/final_reproducibility.R")
run_stage("36_finalize_project_analysis", {
  timestamp<-utc_now()
  s<-list(stage=9L,timestamp=timestamp,stage08_verified=FALSE,GSE263995_audited=FALSE,
    GSE263995_eligible=FALSE,GSE263995_patients=3L,GSE263995_tumor_pairs=3L,
    GSE289825_audited=FALSE,GSE289825_eligible=FALSE,single_cell_status="NOT_COMPLETED",
    final_figures=0L,final_tables=0L,reproducibility_failures=0L,checksum_failures=0L,
    raw_integrity_failures=0L,claim_audit_failures=0L,external_transportability="SKIPPED",
    multiplicity_confirmed_pathways=0L,incremental_prediction_gain=FALSE,
    decision="STOP_REVIEW",decision_reason="Final audit not completed")
  tryCatch({
    audit<-read09("results/stage09_final/final_reproducibility_audit.tsv")
    checks<-read09("results/stage09_final/frozen_checksum_validation.tsv")
    raw<-read09("results/stage09_final/raw_data_integrity.tsv")
    claims<-read09("results/stage09_final/final_claim_audit.tsv")
    eligibility<-read09("results/stage09_final/single_cell_eligibility.tsv")
    sources<-read09("metadata/single_cell/single_cell_source_audit.tsv")
    completion<-read09("results/stage09_final/single_cell_completion.tsv")
    s$stage08_verified<-audit$status[audit$stage==8]=="PASS"
    s$GSE263995_audited<-"GSE263995"%in%sources$accession
    s$GSE289825_audited<-"GSE289825"%in%sources$accession
    s$GSE263995_eligible<-eligibility$eligible[eligibility$accession=="GSE263995"]
    s$GSE289825_eligible<-eligibility$eligible[eligibility$accession=="GSE289825"]
    s$single_cell_status<-completion$status
    s$final_figures<-length(list.files(path("figures/final"),pattern="^Figure[1-6][.]pdf$"))
    s$final_tables<-length(list.files(path("tables/final"),pattern="^Table[1-6]_.*[.]tsv$"))
    s$reproducibility_failures<-sum(audit$status!="PASS" | !audit$required_outputs_present)
    s$checksum_failures<-sum(!checks$passed)+sum(!manifest_check09("results/stage09_final/stage09_analysis_checksums.tsv")$passed)
    s$raw_integrity_failures<-sum(!raw$passed)
    s$claim_audit_failures<-sum(claims$boundary_check!="PASS")
    assert09(s$stage08_verified && s$GSE263995_audited && s$GSE289825_audited &&
      s$final_figures==6 && s$final_tables==6 && s$single_cell_status%in%c("COMPLETED","OMITTED") &&
      sum(unlist(s[c("reproducibility_failures","checksum_failures","raw_integrity_failures","claim_audit_failures")]))==0,
      "Final prerequisite or integrity check failed")
    s$decision<-"PROCEED_TO_MANUSCRIPT"
    s$decision_reason<-paste("Frozen Stage 1-8 outputs and raw-data hashes verified; six final figures and six tables",
      "are frozen. Single-cell support is descriptive orthogonal cell-type interpretation only;",
      "the targeted panel is omitted for inadequate gene coverage. Claims remain NeoTRIP-only,",
      "hypothesis-generating, without confirmatory multiplicity success, independent replication,",
      "incremental predictive improvement, clinical utility or mechanistic proof. Manuscript writing not started.")
  },error=function(e) {s$decision<<-"STOP_REVIEW";s$decision_reason<<-conditionMessage(e)})
  freeze<-list(primary_cohort="GSE319641",analytic_patients=150L,confirmatory_pathways=6L,
    residual_version="Stage4_v1",primary_outcome_version="Stage5_v2",robustness_version="Stage6_v1",
    prediction_version="Stage7_v1",external_transportability="SKIPPED",single_cell_status=s$single_cell_status,
    single_cell_dataset_used=if(s$single_cell_status=="COMPLETED")"GSE263995" else "NONE",
    multiplicity_confirmed_pathways=0L,clinical_prediction_gain=FALSE,
    manuscript_ready=s$decision=="PROCEED_TO_MANUSCRIPT",freeze_timestamp=timestamp)
  file<-path("metadata/final_gate/final_analysis_freeze_v1.yml")
  yaml::write_yaml(freeze,file);writeLines(sha256(file),paste0(file,".sha256"))
  jsonlite::write_json(s,path("results/stage09_final/stage09_status.json"),auto_unbox=TRUE,pretty=TRUE)
  if(exists("audit")) {
    audit$decision[audit$stage==9]<-s$decision
    audit$status[audit$stage==9]<-if(s$decision=="PROCEED_TO_MANUSCRIPT")"PASS" else "FAIL"
    write09(audit,"results/stage09_final/final_reproducibility_audit.tsv")
  }
  write.table(data.frame(timestamp=timestamp,stage=9L,decision=s$decision,reason=s$decision_reason,
    inputs="Frozen Stages 1-8; official GEO single-cell matrices",outputs="results/stage09_final/stage09_status.json"),
    path("logs/decision_log.tsv"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
  registry09()
})
code<-c("R/sc_source_audit.R","R/sc_qc.R","R/sc_annotation.R","R/sc_pathway_scoring.R","R/final_figures.R",
  "R/final_tables.R","R/final_reproducibility.R",paste0("scripts/",c("31_audit_single_cell_sources.R",
  "32_process_single_cell.R","33_score_single_cell_programs.R","34_build_final_outputs.R",
  "35_run_final_reproducibility_audit.R","36_finalize_project_analysis.R")),
  "scripts/download_stage09.py","scripts/install_stage09.R","scripts/prepare_stage09_reference.R",
  "scripts/archive_stage09.py","scripts/validate_stage09_figures.py","renv.lock","metadata/data_inventory.tsv",
  "metadata/final_gate/final_analysis_freeze_v1.yml","metadata/final_gate/final_analysis_freeze_v1.yml.sha256",
  "metadata/final_gate/final_figure_plan_v1.tsv")
files<-unique(c(path(code),unlist(lapply(path(c("metadata/single_cell","results/stage09_final",
  "data/processed/single_cell","figures/final","tables/final","models/stage09")),list.files,recursive=TRUE,full.names=TRUE)),
  list.files(path("logs"),pattern="^(31_|32_|33_|34_|35_|36_).*\\.log$",full.names=TRUE)))
files<-sort(files[file.exists(files) & basename(files)!="stage09_output_checksums.tsv"])
write09(data.frame(local_relative_path=substring(files,nchar(root)+2),
  file_size_bytes=file.info(files)$size,sha256=vapply(files,sha256,"")),"results/stage09_final/stage09_output_checksums.tsv")
result<-system2("python",shQuote(path("scripts/archive_stage09.py")))
assert09(result==0,"Final archive construction failed")
cat("results/stage09_final/stage09_status.json\n")

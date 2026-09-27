tables09 <- function() {
  eligibility <- read09("results/stage09_final/single_cell_eligibility.tsv")
  complete <- read09("results/stage09_final/single_cell_completion.tsv")$status == "COMPLETED"
  feasibility <- data.frame(accession=c("GSE319641","GSE122630","GSE263995","GSE289825"),
    role=c("primary","transportability planned but skipped by eligibility gate",
      "orthogonal descriptive cell-type interpretation","targeted-panel coverage audit"),
    status=c("COMPLETED","SKIPPED",if(complete)"COMPLETED" else "OMITTED",
      eligibility$dataset_status[eligibility$accession=="GSE289825"]),
    clinical_validation=c("NOT_EXTERNAL_VALIDATION","NOT_PERFORMED","NOT_PERFORMED","NOT_PERFORMED"),
    evidence=c("results/stage01_manifest/stage01_status.json",
      "results/stage08_final_gate/transportability_decision.json",
      "results/stage09_final/single_cell_completion.tsv","results/stage09_final/single_cell_gene_coverage.tsv"))
  write09(feasibility,"tables/final/Table1_dataset_feasibility.tsv")
  write09(read09("results/stage01_manifest/cohort_counts.tsv"),"tables/final/Table2_neotrip_cohort.tsv")
  sets <- read09("metadata/pathways/confirmatory_gene_sets_v1.tsv")
  primary <- read09("results/stage05_primary/primary_residual_associations_v2.tsv")
  pathway <- do.call(rbind,lapply(axes07,function(a) {
    z<-sets[sets$axis==a,]
    data.frame(axis=a,gene_set=unique(z$gs_name),frozen_genes=length(unique(z$gene_symbol)),
      prespecified_direction=primary$prespecified_direction[primary$axis==a],
      version=unique(z$db_version),source=unique(z$source))
  }))
  write09(pathway,"tables/final/Table3_confirmatory_pathways.tsv")
  write09(primary[match(axes07,primary$axis),],"tables/final/Table4_primary_residual_results.tsv")
  summary <- read_tsv(path("results/stage07_prediction/prediction_summary.tsv"))
  paired <- read_tsv(path("results/stage07_prediction/paired_metric_summary.tsv"))
  summary$result_type <- "repeat_metric_summary"; paired$result_type <- "paired_repeat_difference"
  columns <- union(names(summary),names(paired))
  fill <- function(x) {for(n in setdiff(columns,names(x)))x[[n]]<-"";x[,columns]}
  write09(rbind(fill(summary),fill(paired)),"tables/final/Table5_prediction_results.tsv")
  if(complete) {
    support<-read09("results/stage09_final/pathway_compartment_interpretation.tsv")
    support<-support[order(match(support$axis,axes07),support$broad_compartment),]
  } else support<-read09("results/stage09_final/single_cell_completion.tsv")
  write09(support,"tables/final/Table6_orthogonal_support.tsv")
  write09(data.frame(note="No additional protocol dataset roles were present in locally frozen config, metadata or scope; no roles inferred or replacement datasets searched.",
    evidence="config/project.yml;metadata/final_gate/analysis_scope_freeze_v1.yml"),
    "metadata/single_cell/additional_dataset_role_audit.tsv")
}

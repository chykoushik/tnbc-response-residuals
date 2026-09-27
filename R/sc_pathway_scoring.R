score09 <- function() {
  stopifnot(requireNamespace("SingleCellExperiment", quietly = TRUE),
    requireNamespace("SummarizedExperiment", quietly = TRUE))
  eligible <- read09("results/stage09_final/single_cell_eligibility.tsv")
  if (!eligible$eligible[eligible$accession == "GSE263995"]) {
    write09(data.frame(status = "OMITTED", reason = "Six-pathway coverage gate failed"),
      "results/stage09_final/single_cell_completion.tsv")
    return(invisible(NULL))
  }
  reference <- readRDS(path("models/stage09/HPCA_reference.rds"))
  mapping <- mapping09(reference)
  manifest <- read09("metadata/single_cell/GSE263995_tumor_manifest.tsv")
  genes <- read09("metadata/pathways/confirmatory_gene_sets_v1.tsv")
  sets <- lapply(axes07, function(a) unique(genes$gene_symbol[genes$axis == a])); names(sets) <- axes07
  cells <- summaries <- marker_rows <- pseudobulk <- list(); profile_info <- list()
  for (i in seq_len(nrow(manifest))) {
    info <- manifest[i, ]; gsm <- info$GSM
    checkpoint <- path("models/stage09",paste0(gsm,"_annotated.rds"))
    if (file.exists(checkpoint)) sce <- readRDS(checkpoint) else {
      sce <- readRDS(path("models/stage09", paste0(gsm,"_qc.rds")))$sce
      sce <- annotate09(sce,gsm,reference,mapping); saveRDS(sce,checkpoint)
    }
    score <- UCell::ScoreSignatures_UCell(SummarizedExperiment::assay(sce,"logcounts"),
      features = sets, maxRank = 1500, ties.method = "average", missing_genes = "impute",
      BPPARAM = BiocParallel::SerialParam(), ncores = 1)
    assert09(identical(rownames(score), colnames(sce)), "UCell barcode mismatch")
    score <- as.data.frame(score)
    names(score) <- sub("_UCell$", "", names(score))
    cells[[gsm]] <- data.frame(GSM = gsm, patient_id = info$patient_id, timepoint = info$timepoint,
      cell_barcode = colnames(sce), broad_compartment = sce$broad_compartment, score[,axes07])
    marker_rows[[gsm]] <- markers09(sce,info)
    for (compartment in sort(unique(sce$broad_compartment))) {
      keep <- sce$broad_compartment == compartment; n <- sum(keep)
      for (a in axes07) {
        v <- score[[a]][keep]
        summaries[[length(summaries)+1L]] <- data.frame(patient_id = info$patient_id,
          timepoint = info$timepoint, broad_compartment = compartment, axis = a, cell_count = n,
          median_UCell = median(v), mean_UCell = mean(v), IQR = IQR(v),
          status = if (n >= 50) "EVALUABLE" else "LOW_CELL_COUNT")
      }
      if (n >= 50) {
        id <- paste(info$patient_id,info$timepoint,compartment,sep="|")
        pseudobulk[[id]] <- Matrix::rowSums(SummarizedExperiment::assay(sce,"counts")[,keep,drop=FALSE])
        profile_info[[id]] <- data.frame(profile_id = id, patient_id = info$patient_id,
          timepoint = info$timepoint,broad_compartment = compartment,cell_count = n)
      }
    }
    cat("Annotated and scored",gsm,ncol(sce),"cells\n")
    rm(sce); gc()
  }
  con <- gzfile(path("data/processed/single_cell/GSE263995_ucell_scores.tsv.gz"),"wt")
  tryCatch(write_tsv(do.call(rbind,cells),con),finally=close(con))
  summary <- do.call(rbind,summaries)
  write09(summary,"data/processed/single_cell/GSE263995_compartment_summary.tsv")
  write09(do.call(rbind,marker_rows),"results/stage09_final/cell_annotation_qc.tsv")
  common <- Reduce(intersect,lapply(pseudobulk,names))
  counts <- do.call(cbind,lapply(pseudobulk,function(v) v[common]))
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts),method="TMM")
  logcpm <- edgeR::cpm(dge,log=TRUE,prior.count=2)
  ranks <- singscore::rankGenes(logcpm)
  pb <- do.call(rbind, lapply(axes07,function(a) {
    scores <- singscore::simpleScore(ranks,upSet=sets[[a]],centerScore=TRUE,knownDirection=TRUE)
    data.frame(do.call(rbind,profile_info),axis=a,singscore=scores$TotalScore)
  }))
  write09(pb,"data/processed/single_cell/GSE263995_pseudobulk_scores.tsv")
  saveRDS(list(counts=counts,normalization_factors=dge$samples$norm.factors,logCPM=logcpm),
    path("models/stage09/pseudobulk_profiles.rds"))
  pre <- summary[summary$timepoint=="Pre",]; post <- summary[summary$timepoint=="Post",]
  keys <- c("patient_id","broad_compartment","axis")
  paired <- merge(pre,post,by=keys,suffixes=c("_Pre","_Post"),all=TRUE)
  paired$evaluable <- !is.na(paired$cell_count_Pre) & !is.na(paired$cell_count_Post) &
    paired$cell_count_Pre>=50 & paired$cell_count_Post>=50
  paired$delta_median_UCell <- ifelse(paired$evaluable,paired$median_UCell_Post-paired$median_UCell_Pre,NA_real_)
  paired$status <- ifelse(paired$evaluable,"DESCRIPTIVE_ONLY","LOW_CELL_COUNT_OR_MISSING")
  write09(paired,"results/stage09_final/single_cell_paired_program_change.tsv")
  ppre <- pb[pb$timepoint=="Pre",c(keys,"singscore")]
  ppost <- pb[pb$timepoint=="Post",c(keys,"singscore")]
  pc <- merge(ppre,ppost,by=keys,suffixes=c("_Pre","_Post"))
  pc$delta_singscore <- pc$singscore_Post-pc$singscore_Pre
  concord <- merge(paired[,c(keys,"evaluable","delta_median_UCell")],pc,by=keys,all.x=TRUE)
  concord$direction_agrees <- ifelse(concord$evaluable & !is.na(concord$delta_singscore),
    sign(concord$delta_median_UCell)==sign(concord$delta_singscore),NA)
  write09(concord,"results/stage09_final/single_cell_scoring_concordance.tsv")
  interpretations <- list()
  for (a in axes07) {
    s <- summary[summary$axis==a & summary$status=="EVALUABLE",]
    distributions <- aggregate(median_UCell~broad_compartment+timepoint,s,median)
    nn <- aggregate(patient_id~broad_compartment+timepoint,s,length)
    distributions$n_patients <- nn$patient_id
    top <- function(time) {
      d <- distributions[distributions$timepoint==time & distributions$n_patients>=2,]
      if (!nrow(d)) return(NA_character_)
      d$broad_compartment[d$median_UCell==max(d$median_UCell)]
    }
    before <- top("Pre"); after <- top("Post")
    c <- concord[concord$axis==a & concord$evaluable & !is.na(concord$direction_agrees),]
    discordant <- if(nrow(c)) mean(!c$direction_agrees)>.5 else TRUE
    label <- if(anyNA(c(before,after)) || discordant) "INCONCLUSIVE" else if(
      length(before)!=1 || length(after)!=1 || before!=after) "MULTICOMPARTMENT" else
      switch(before,immune="IMMUNE_ENRICHED",stromal="STROMAL_ENRICHED",
        epithelial="EPITHELIAL_ENRICHED",other="INCONCLUSIVE")
    for (comp in c("immune","stromal","epithelial","other")) {
      d <- distributions[distributions$broad_compartment==comp,]
      pp <- paired[paired$axis==a & paired$broad_compartment==comp,]
      value <- function(time) {v<-d$median_UCell[d$timepoint==time];if(length(v))v else NA_real_}
      interpretations[[length(interpretations)+1L]] <- data.frame(axis=a,broad_compartment=comp,
        baseline_median_across_patients=value("Pre"),post_median_across_patients=value("Post"),
        paired_patient_directions=paste(paste(pp$patient_id,ifelse(pp$evaluable,
          sign(pp$delta_median_UCell),"NOT_EVALUABLE"),sep=":"),collapse=";"),
        evaluable_paired_patients=sum(pp$evaluable),interpretation=label,
        scoring_discordance_fraction=if(nrow(c))mean(!c$direction_agrees) else NA_real_,
        role="orthogonal cell-type interpretation; descriptive only")
    }
  }
  write09(do.call(rbind,interpretations),"results/stage09_final/pathway_compartment_interpretation.tsv")
  write09(data.frame(status="COMPLETED",reason="Three-patient orthogonal descriptive cell-type interpretation; no clinical outcomes or hypothesis tests"),
    "results/stage09_final/single_cell_completion.tsv")
}

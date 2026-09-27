mapping09 <- function(reference) {
  groups <- list(
    immune = c("B_cell","BM","BM & Prog.","CMP","DC","Erythroblast","GMP","HSC_-G-CSF","HSC_CD34+",
      "Macrophage","MEP","Monocyte","Myelocyte","Neutrophils","NK_cell","Platelets",
      "Pre-B_cell_CD34-","Pro-B_cell_CD34+","Pro-Myelocyte","T_cells"),
    stromal = c("Chondrocytes","Endothelial_cells","Fibroblasts","MSC","Osteoblasts","Smooth_muscle_cells"),
    epithelial = c("Epithelial_cells","Keratinocytes"),
    other = c("Astrocyte","Embryonic_stem_cells","Gametocytes","Hepatocytes","iPS_cells",
      "Neuroepithelial_cell","Neurons","Tissue_stem_cells","UNASSIGNED"))
  mapping <- do.call(rbind, lapply(names(groups), function(g)
    data.frame(reference_label = groups[[g]], broad_compartment = g)))
  assert09(all(unique(reference$label.main) %in% mapping$reference_label), "MODIFY_REVIEW: unmapped reference label")
  file <- path("metadata/single_cell/singler_compartment_mapping_v1.tsv")
  if (!file.exists(file)) {
    write_tsv(mapping, file); writeLines(sha256(file), paste0(file, ".sha256"))
  } else assert09(identical(readLines(paste0(file, ".sha256")),sha256(file)) &&
    identical(read_tsv(file)$reference_label, mapping$reference_label) &&
    identical(read_tsv(file)$broad_compartment, mapping$broad_compartment), "Annotation mapping changed")
  mapping
}
annotate09 <- function(sce, gsm, reference, mapping) {
  labels <- SingleR::SingleR(test = sce, ref = reference, labels = reference$label.main,
    assay.type.test = "logcounts", BPPARAM = BiocParallel::SerialParam())
  label <- labels$pruned.labels; label[is.na(label)] <- "UNASSIGNED"
  assert09(all(label %in% mapping$reference_label), "MODIFY_REVIEW: observed label missing explicit mapping")
  sce$SingleR_label <- labels$labels; sce$SingleR_pruned_label <- label
  sce$broad_compartment <- mapping$broad_compartment[match(label, mapping$reference_label)]
  write09(data.frame(GSM = gsm, cell_barcode = colnames(sce), SingleR_label = labels$labels,
    SingleR_pruned_label = label, broad_compartment = sce$broad_compartment),
    paste0("metadata/single_cell/",gsm,"_annotation.tsv"))
  sce
}
markers09 <- function(sce, info) {
  markers <- c("PTPRC","CD3D","CD3E","CD79A","MS4A1","LYZ","LST1","EPCAM","KRT8","KRT18","KRT19",
    "COL1A1","COL1A2","DCN","COL3A1","PECAM1","VWF")
  m <- SummarizedExperiment::assay(sce,"logcounts")
  do.call(rbind, lapply(sort(unique(sce$broad_compartment)), function(compartment) {
    keep <- sce$broad_compartment == compartment
    do.call(rbind,lapply(markers,function(gene) {
      present <- gene %in% rownames(m); v <- if (present) as.numeric(m[gene,keep]) else NA_real_
      data.frame(GSM = info$GSM, patient_id = info$patient_id, timepoint = info$timepoint,
        broad_compartment = compartment, marker = gene, present = present, cells = sum(keep),
        mean_log_normalized = if (present) mean(v) else NA_real_,
        fraction_detected = if (present) mean(v > 0) else NA_real_)
    }))
  }))
}

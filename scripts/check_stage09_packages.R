p <- c('Seurat','SingleCellExperiment','scuttle','scDblFinder','UCell','SingleR','celldex','edgeR')
print(sapply(p,requireNamespace,quietly=TRUE))

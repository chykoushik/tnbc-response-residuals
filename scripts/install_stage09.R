source('R/paths.R')
options(repos=BiocManager::repositories(version='3.22'),timeout=1200)
renv::install(c('Seurat','bioc::scuttle','bioc::scDblFinder','bioc::UCell','bioc::SingleR','bioc::celldex'),prompt=FALSE)

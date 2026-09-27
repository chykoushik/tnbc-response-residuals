source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
options(repos = BiocManager::repositories(version = "3.22"), timeout = 1200)
renv::install("bioc::GSVA", prompt = FALSE)

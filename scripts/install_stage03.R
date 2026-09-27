source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
options(repos = BiocManager::repositories(version = "3.22"), timeout = 1200)
renv::install(c("msigdbr@26.1.1", "bioc::singscore"), prompt = FALSE)

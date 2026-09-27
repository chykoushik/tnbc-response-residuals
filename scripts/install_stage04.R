source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1200)
renv::install("glmnet", prompt = FALSE)

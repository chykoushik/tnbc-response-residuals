source("D:/TNBC-ResponseResiduals/R/paths.R")
directories <- c("config", "data/raw/geo", "data/raw/author_repo", "data/processed",
  "external", "external/tmp", "external/R-library", "figures", "logs", "metadata",
  "models", "R", "results/stage01_manifest", "scripts", "tables")
invisible(lapply(directories, function(x) dir.create(path(x), recursive = TRUE, showWarnings = FALSE)))
setwd(root)
run_stage("00_setup", {
  options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1200)
  .libPaths(c(path("external/R-library"), .libPaths()))
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv", lib = path("external/R-library"), quiet = TRUE)
  }
  if (!file.exists(path("renv/activate.R"))) renv::init(project = root, bare = TRUE, restart = FALSE)
  renv::load(project = root, quiet = TRUE)
  packages <- c("digest", "jsonlite", "yaml")
  needed <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(needed)) renv::install(needed, project = root)
  config <- yaml::read_yaml(path("config/project.yml"))
  if (!identical(config$root, root)) stop("Project root mismatch")
  software_versions()
})
cat("Setup complete\n")

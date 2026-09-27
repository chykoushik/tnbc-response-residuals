root <- "D:/TNBC-ResponseResiduals"
path <- function(...) file.path(root, ...)
utc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
read_tsv <- function(file) {
  read.delim(file, sep = "\t", header = TRUE, quote = "", comment.char = "",
             colClasses = "character", check.names = FALSE, na.strings = NULL)
}
write_tsv <- function(x, file) {
  write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
missing_value <- function(x) is.na(x) | trimws(x) %in% c("", "NA", "N/A")
run_stage <- function(name, code) {
  previous_options <- options(warn = 1)
  on.exit(options(previous_options), add = TRUE)
  log <- file(path("logs", paste0(name, ".log")), "at")
  sink(log)
  sink(log, type = "message")
  error <- tryCatch({ force(code); NULL }, error = function(e) e)
  if (!is.null(error)) cat(utc_now(), conditionMessage(error), "\n")
  sink(type = "message")
  sink()
  close(log)
  if (!is.null(error)) stop("Checks failed", call. = FALSE)
}
software_versions <- function() {
  packages <- c("renv", "yaml", "digest", "jsonlite")
  write_tsv(data.frame(R_version = R.version.string, platform = R.version$platform,
    package = c("base", packages), package_version = c(as.character(getRversion()),
    vapply(packages, function(x) as.character(utils::packageVersion(x)), ""))),
    path("logs", "software_versions.tsv"))
}

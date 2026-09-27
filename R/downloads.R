download_source <- function(source, url, relative, accession, commit = "", force = FALSE) {
  destination <- path(relative)
  inventory_file <- path("metadata", "data_inventory.tsv")
  old <- if (file.exists(inventory_file)) read_tsv(inventory_file) else NULL
  previous <- if (!is.null(old)) old[old$local_relative_path == relative, , drop = FALSE] else NULL
  retrieved <- ""
  provenance <- "preexisting_retrieval_unknown"
  if (file.exists(destination)) {
    hash <- sha256(destination)
    if (!force && !is.null(previous) && nrow(previous) == 1L) {
      if (hash != previous$sha256) stop("Existing checksum mismatch")
      return(previous)
    }
  }
  if (!file.exists(destination) || force) {
    temporary <- tempfile("download-", tmpdir = path("external", "tmp"))
    on.exit(unlink(temporary), add = TRUE)
    status <- utils::download.file(url, temporary, mode = "wb", quiet = TRUE,
                                  method = "libcurl")
    if (status != 0L || file.info(temporary)$size == 0) stop("Source download failed")
    hash <- sha256(temporary)
    if (!file.copy(temporary, destination, overwrite = force)) stop("Source write failed")
    retrieved <- utc_now()
    provenance <- "downloaded"
  }
  row <- data.frame(source = source, url = url, retrieval_timestamp = retrieved,
    local_relative_path = relative, file_size_bytes = file.info(destination)$size,
    sha256 = hash, accession_or_repository = accession,
    public_access_status = "public_unauthenticated", repository_commit = commit,
    acquisition_status = provenance, checksum_timestamp = utc_now())
  if (!is.null(old)) old <- old[old$local_relative_path != relative, , drop = FALSE]
  write_tsv(rbind(old, row), inventory_file)
  row
}
repository_sha <- function(repository) {
  tryCatch({
    result <- suppressWarnings(system2("git", c("ls-remote", repository, "HEAD"),
                                      stdout = TRUE, stderr = TRUE))
    line <- grep("^[0-9a-f]{40}[[:space:]]+HEAD$", result, value = TRUE)
    if (length(line) != 1L) "" else sub("[[:space:]].*$", "", line)
  }, error = function(e) "")
}

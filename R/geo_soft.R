parse_geo_soft <- function(file) {
  connection <- gzfile(file, "rt")
  on.exit(close(connection))
  lines <- readLines(connection, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE = ", lines)
  if (!length(starts)) stop("GEO samples missing")
  entities <- grep("^\\^", lines)
  rows <- lapply(starts, function(start) {
    following <- entities[entities > start]
    end <- if (length(following)) following[1L] - 1L else length(lines)
    block <- lines[start:end]
    get <- function(prefix) {
      values <- unique(substring(block[startsWith(block, prefix)], nchar(prefix) + 1L))
      if (length(values) > 1L) stop("Repeated GEO field conflict")
      if (length(values)) values else ""
    }
    data.frame(gsm = get("^SAMPLE = "), sample_title = get("!Sample_title = "),
      geo_library_name = get("!Sample_description = Library name: "),
      geo_patient_id = get("!Sample_characteristics_ch1 = patient id: "),
      geo_timepoint = get("!Sample_characteristics_ch1 = timepoint: "),
      geo_batch = get("!Sample_characteristics_ch1 = batch: "),
      organism = get("!Sample_organism_ch1 = "), platform = get("!Sample_platform_id = "))
  })
  do.call(rbind, rows)
}

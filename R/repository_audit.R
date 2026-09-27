repo_url <- "https://github.com/mdugo/Early_GEP_NeoTRIP"
repo_commit <- "26b765cc816fcb04794797f22d78d944e834a25b"
instaprism_commit <- "7d3b57cde9342c5eac7282f4c2d8bff09f21e693"
reference_commit <- "33d5e0cb76fbe5d06f4375e7b69de27532d4f075"
git_output <- function(directory, arguments) {
  value <- system2("git", c("-C", shQuote(path(directory)), arguments), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(value, "status"))) stop("Repository command failed")
  value
}
pin_repository <- function(directory, url, commit) {
  if (!dir.exists(path(directory, ".git"))) {
    status <- system2("git", c("clone", url, shQuote(path(directory))), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(status, "status"))) stop("Repository clone failed")
  }
  if (length(git_output(directory, c("status", "--porcelain", "--untracked-files=no")))) stop("Source repository modified")
  if (!identical(git_output(directory, c("rev-parse", "HEAD")), commit)) {
    git_output(directory, c("checkout", "--detach", commit))
  }
  if (!identical(git_output(directory, c("rev-parse", "HEAD")), commit)) stop("Repository revision mismatch")
}
inventory_source <- function(relative, url, source, accession, commit = "", retrieved = utc_now()) {
  inventory_file <- path("metadata/data_inventory.tsv")
  inventory <- read_tsv(inventory_file)
  previous <- inventory[inventory$local_relative_path == relative, , drop = FALSE]
  hash <- sha256(path(relative))
  if (nrow(previous)) {
    if (nrow(previous) != 1L || previous$sha256 != hash) stop("Source checksum mismatch")
    return(invisible(previous))
  }
  row <- data.frame(source = source, url = url, retrieval_timestamp = retrieved,
    local_relative_path = relative, file_size_bytes = file.info(path(relative))$size,
    sha256 = hash, accession_or_repository = accession, public_access_status = "public_unauthenticated",
    repository_commit = commit, acquisition_status = "downloaded", checksum_timestamp = utc_now())
  write_tsv(rbind(inventory, row), inventory_file)
  invisible(row)
}
fetch_public <- function(url, relative, source = "Publication metadata", accession = "10.1016/j.annonc.2026.05.694", commit = "") {
  access_file <- path("metadata/composition/public_access_attempts.tsv")
  attempts <- if (file.exists(access_file)) read_tsv(access_file) else NULL
  result <- "cached"
  if (!file.exists(path(relative))) {
    temporary <- tempfile("public-", tmpdir = path("external/tmp"))
    on.exit(unlink(temporary), add = TRUE)
    result <- tryCatch({
      code <- suppressWarnings(utils::download.file(url, temporary, quiet = TRUE, mode = "wb", method = "libcurl"))
      if (code != 0L || file.info(temporary)$size == 0) stop("Empty public response")
      dir.create(dirname(path(relative)), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(temporary, path(relative), overwrite = FALSE)) stop("Public source write failed")
      "downloaded"
    }, error = function(e) gsub("[\r\n\t]", " ", conditionMessage(e)))
  }
  ok <- file.exists(path(relative))
  if (ok) inventory_source(relative, url, source, accession, commit)
  row <- data.frame(timestamp = utc_now(), url = url, local_path = if (ok) relative else "",
    downloaded = ok, access_result = result)
  write_tsv(rbind(attempts, row), access_file)
  ok
}
audit_repository <- function() {
  directory <- "external/neotrip_2026_repo"
  pin_repository(directory, repo_url, repo_commit)
  files <- git_output(directory, c("ls-files"))
  terms <- c("InstaPrism", "instaprism", "deconvolution", "Tumor_fraction", "Cancer_Epithelial",
    "Stromal_fraction", "Immune_fraction", "tumor cells", "tumour cells", "tumor content", "tumour content",
    "tumor purity", "tumour purity", "pathology", "pathologic", "early pCR", "Early.pCR", "pCR", "sTIL", "PD-L1")
  hits <- data.frame(term = character(), file = character(), line_number = integer(), matched_text = character())
  manifest <- list()
  for (file in files) {
    relative <- paste0(directory, "/", file)
    inventory_source(relative, paste0(repo_url, "/blob/", repo_commit, "/", file), "NeoTRIP source repository", repo_url, repo_commit)
    lines <- readLines(path(relative), warn = FALSE, encoding = "UTF-8")
    for (term in terms) {
      index <- grep(term, lines, fixed = TRUE)
      if (length(index)) hits <- rbind(hits, data.frame(term = term, file = file, line_number = index,
        matched_text = substr(gsub("[\t\r\n]", " ", trimws(lines[index])), 1L, 180L)))
    }
    has <- function(pattern) any(grepl(pattern, lines, ignore.case = TRUE))
    manifest[[file]] <- data.frame(relative_path = file, size_bytes = file.info(path(relative))$size,
      sha256 = sha256(path(relative)), file_type = tools::file_ext(file),
      contains_instaprism = has("instaprism"), contains_tumor_fraction = has("Tumor_fraction|Cancer_Epithelial"),
      contains_pathology = has("patholog|tumou?r cells|tumou?r content|tumou?r purity"),
      contains_pcr = has("pCR"), notes = "Text occurrence flags only; variable availability requires content review.")
  }
  write_tsv(do.call(rbind, manifest), path("metadata/neotrip_public_repo_files.tsv"))
  write_tsv(hits, path("results/stage02_composition/repository_term_audit.tsv"))
  inventory <- read_tsv(path("metadata/data_inventory.tsv"))
  retrieved <- inventory$retrieval_timestamp[match(paste0(directory, "/", files[1L]), inventory$local_relative_path)]
  write_tsv(data.frame(repository_url = repo_url, commit_sha = repo_commit, retrieval_timestamp = retrieved),
    path("metadata/composition/repository_pin.tsv"))
}
audit_public_sources <- function() {
  paper <- "external/neotrip_2026_paper/"
  endpoints <- c(
    crossref = "https://api.crossref.org/works/10.1016/j.annonc.2026.05.694",
    europepmc = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:42162664%20AND%20SRC:MED&format=json&resultType=core",
    elsevier = "https://api.elsevier.com/content/article/pii/S0923753426008781?httpAccept=text/xml",
    zenodo_rna = "https://zenodo.org/api/records/17871576",
    zenodo_imc = "https://zenodo.org/api/records/7990870",
    publisher = "https://www.annalsofoncology.org/article/S0923-7534(26)00878-1/fulltext",
    sciencedirect = "https://www.sciencedirect.com/science/article/pii/S0923753426008781")
  filenames <- c("crossref.json", "europepmc.json", "elsevier.xml", "zenodo_rna.json", "zenodo_imc.json", "publisher.html", "sciencedirect.html")
  for (i in seq_along(endpoints)) fetch_public(endpoints[i], paste0(paper, filenames[i]))
  if (!file.exists(path(paper, "crossref.json")) || !file.exists(path(paper, "europepmc.json"))) stop("Bibliographic evidence missing")
  crossref <- jsonlite::fromJSON(path(paper, "crossref.json"))$message
  europe <- jsonlite::fromJSON(path(paper, "europepmc.json"), simplifyVector = FALSE)$resultList$result[[1L]]
  write_tsv(data.frame(pmid = "42162664", doi = crossref$DOI, title = crossref$title[1L],
    journal = crossref$`container-title`[1L], year = "2026", volume = crossref$volume,
    issue = crossref$issue, pages = crossref$page), path(paper, "bibliography.tsv"))
  writeLines(gsub("</?[^>]+>", "", europe$abstractText), path(paper, "abstract.txt"))
  sources <- data.frame(source_type = character(), source_name = character(), public = character(), downloaded = logical(),
    requires_request = character(), contains_sample_level_data = character(), contains_composition_data = character(),
    contains_pathology_data = character(), sample_link_possible = character(), evidence = character(), local_path = character(), sha256 = character())
  add <- function(type, name, public, downloaded, request, sample, composition, pathology, link, evidence, local = "") {
    sources <<- rbind(sources, data.frame(source_type = type, source_name = name, public = public, downloaded = downloaded,
      requires_request = request, contains_sample_level_data = sample, contains_composition_data = composition,
      contains_pathology_data = pathology, sample_link_possible = link, evidence = evidence,
      local_path = local, sha256 = if (nzchar(local) && file.exists(path(local))) sha256(path(local)) else ""))
  }
  add("repository", repo_url, "TRUE", TRUE, "FALSE", "TRUE", "generation_code_only", "FALSE", "TRUE",
    "Pinned tracked-file audit; public phenotype has nine fields, no independent pathology tumour-content field.", "external/neotrip_2026_repo/data/NeoTRIP_baseline_D1C2_pheno_and_clinical.txt")
  add("GEO", "GSE319641", "TRUE", TRUE, "FALSE", "TRUE", "bulk_expression_only", "FALSE", "TRUE",
    "SOFT SAMPLE characteristics: tissue, patient id, timepoint, batch; no pathology tumour-content measurement.", "data/raw/geo/GSE319641_family.soft.gz")
  add("publication_abstract", "PMID 42162664", "TRUE", TRUE, "FALSE", "FALSE", "aggregate_description", "aggregate_description", "FALSE",
    "Abstract Results reports on-treatment tumour-cell absence and deconvolution; no keyed sample-level measurements in abstract.", paste0(paper, "europepmc.json"))
  add("publication_metadata", "DOI 10.1016/j.annonc.2026.05.694", "TRUE", TRUE, "FALSE", "FALSE", "FALSE", "FALSE", "FALSE",
    "Crossref identifies CC BY-NC-ND open access; full article body was not downloaded through publisher endpoints.", paste0(paper, "crossref.json"))
  add("publication_fulltext_and_supplements", endpoints["publisher"], "UNKNOWN_ACCESS", file.exists(path(paper, "publisher.html")), "UNKNOWN", "UNKNOWN", "UNKNOWN", "UNKNOWN", "UNVERIFIED",
    "Publisher and ScienceDirect automated access attempts logged; no accessible supplement file URL was identified. Unretrieved material is not evidence of absence.")
  for (id in c("rna", "imc")) {
    local <- paste0(paper, "zenodo_", id, ".json")
    if (!file.exists(path(local))) next
    record <- jsonlite::fromJSON(path(local), simplifyVector = FALSE)
    if (!identical(record$metadata$access_right, "restricted") || length(record$files)) stop("Public source changed")
    add("dataset_record", record$doi, "metadata_only", TRUE, "TRUE", "restricted_unverified", "restricted_unverified",
      "restricted_unverified", "UNVERIFIED", paste("Public API metadata.access_right = restricted; files list empty.",
        if (id == "rna") "Description directs requests for additional clinical metadata." else "IMC description lists epithelial masks and clinical data; files require access request."), local)
  }
  add("repository_gene_sets", "Supplementary Table 1 gene sets", "TRUE", TRUE, "FALSE", "FALSE", "FALSE", "FALSE", "FALSE",
    "Public gene set definitions; no patient-level pathology measurements; used only for evidence audit.", "data/raw/author_repo/geneSets_Supplementary_Table1.txt")
  write_tsv(sources, path("metadata/neotrip_public_source_audit.tsv"))
}
audit_source_schemas <- function() {
  connection <- gzfile(path("data/raw/geo/GSE319641_family.soft.gz"), "rt")
  lines <- readLines(connection, warn = FALSE)
  close(connection)
  characteristics <- unique(sub(":.*$", "", sub("^!Sample_characteristics_ch1 = ", "",
    lines[startsWith(lines, "!Sample_characteristics_ch1 = ")])))
  descriptions <- unique(sub(":.*$", "", sub("^!Sample_description = ", "",
    lines[startsWith(lines, "!Sample_description = ")])))
  author <- names(read_tsv(path("data/raw/author_repo/NeoTRIP_baseline_D1C2_pheno_and_clinical.txt")))
  schema <- data.frame(source = c(rep("GEO characteristics", length(characteristics)),
    rep("GEO description", length(descriptions)), rep("Author phenotype", length(author))),
    field = c(characteristics, descriptions, author),
    review = "Public field schema inspected directly; no unlisted pathology measurement inferred.")
  write_tsv(schema, path("metadata/composition/public_field_schema.tsv"))
  article_url <- "https://iris.unisr.it/retrieve/7983debc-661f-45e8-8d0d-11427faad894/PIIS0923753426008781.pdf"
  fetch_public("https://iris.unisr.it/handle/20.500.11768/206316", "external/neotrip_2026_paper/institutional_record.html")
  article_ok <- fetch_public(article_url, "external/neotrip_2026_paper/PIIS0923753426008781.pdf", "Open access article")
  access <- read_tsv(path("metadata/neotrip_public_source_audit.tsv"))
  local <- if (article_ok) "external/neotrip_2026_paper/PIIS0923753426008781.pdf" else ""
  row <- data.frame(source_type = "open_access_article_copy", source_name = article_url, public = "listed_open_access",
    downloaded = article_ok, requires_request = "FALSE", contains_sample_level_data = "UNVERIFIED",
    contains_composition_data = "aggregate_description", contains_pathology_data = "aggregate_description",
    sample_link_possible = "UNVERIFIED",
    evidence = "Institutional record lists published PDF as accesso aperto, Creative commons. Download endpoint access result is logged; no credentials used.",
    local_path = local, sha256 = if (article_ok) sha256(path(local)) else "")
  access <- access[access$source_type != "open_access_article_copy", , drop = FALSE]
  write_tsv(rbind(access, row), path("metadata/neotrip_public_source_audit.tsv"))
  write_tsv(data.frame(source = c("Publisher supplement landing page", "Repository Supplementary Table 1"),
    url = c("https://doi.org/10.1016/j.annonc.2026.05.694", paste0(repo_url, "/blob/", repo_commit, "/data/geneSets_Supplementary_Table1.txt")),
    access_result = c("Publisher fulltext access failed; no direct supplement file URL identified in accessible metadata",
      "Public repository gene-set table downloaded unchanged"), downloaded = c(FALSE, TRUE)),
    path("metadata/composition/supplement_access_audit.tsv"))
}
audit_extended_sources <- function() {
  fetch_public("https://zenodo.org/api/records/14290163", "external/neotrip_2026_paper/quiche_zenodo.json",
    "Derived NeoTRIP IMC record", "10.5281/zenodo.14290163")
  downloaded <- fetch_public("https://zenodo.org/api/records/14290163/files/nt_preprocessed.h5ad/content",
    "external/neotrip_2026_paper/nt_preprocessed.h5ad", "Derived NeoTRIP IMC data", "10.5281/zenodo.14290163")
  fetch_public("https://zenodo.org/api/records/19226009", "external/neotrip_2026_paper/trajectory_zenodo.json",
    "Related IMC record", "10.5281/zenodo.19226009")
  if (!downloaded) stop("IMC source audit incomplete")
  status <- system2("python", shQuote(path("scripts/audit_imc_keys.py")), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(status, "status"))) stop("IMC key audit failed")
  keys <- read_tsv(path("metadata/composition/imc_key_audit.tsv"))
  if (any(as.integer(keys$stage01_exact_key_matches) != 0L)) stop("New pathology linkage found")
  sources <- read_tsv(path("metadata/neotrip_public_source_audit.tsv"))
  sources <- sources[!sources$source_type %in% c("derived_IMC_dataset", "related_trajectory_dataset"), , drop = FALSE]
  local <- "external/neotrip_2026_paper/nt_preprocessed.h5ad"
  row <- data.frame(source_type = "derived_IMC_dataset", source_name = "https://doi.org/10.5281/zenodo.14290163",
    public = "TRUE", downloaded = TRUE, requires_request = "FALSE", contains_sample_level_data = "TRUE",
    contains_composition_data = "image_cell_annotations", contains_pathology_data = "epithelial_image_candidate_not_histopathologic_percentage",
    sample_link_possible = "FALSE",
    evidence = "QUICHE Data availability and Zenodo description identify 119 baseline NeoTRIP IMC samples; HDF5 obs contains isEpithelial and cell annotations. Exact Patient_ID/Tissue_ID/fov/index overlap with Stage 1 is zero; no crosswalk supplied. Arm and pCR values not read.",
    local_path = local, sha256 = sha256(path(local)))
  extra <- data.frame(source_type = "related_trajectory_dataset", source_name = "https://doi.org/10.5281/zenodo.19226009",
    public = "UNKNOWN_ACCESS", downloaded = FALSE, requires_request = "UNKNOWN", contains_sample_level_data = "UNVERIFIED",
    contains_composition_data = "UNVERIFIED", contains_pathology_data = "UNVERIFIED", sample_link_possible = "UNVERIFIED",
    evidence = "PMC13060115 Data Availability cites this processed-data record; public API returned unavailable. No keyed data inspected.",
    local_path = "", sha256 = "")
  write_tsv(rbind(sources, row, extra), path("metadata/neotrip_public_source_audit.tsv"))
}

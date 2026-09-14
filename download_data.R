#-------------------------------------------------------------------------------
# download_data.R
#
# Fetches the large processed data objects for this analysis from Zenodo.
#
# The repository holds metadata, phylogenies, reference files and analysis
# products. The bulk processed data objects (~270 MB) are deposited on Zenodo
# instead, and are downloaded into data/ by this script.
#
# Usage, from the repository root:
#
#   Rscript download_data.R                 # fetch everything that is missing
#   Rscript download_data.R --list          # show what is needed and what is present
#   Rscript download_data.R --only NW,CML   # fetch only files matching these patterns
#   Rscript download_data.R --force         # re-download even if present
#   Rscript download_data.R --check         # verify checksums of what is already present
#
# Raw sequencing data is not here and not on Zenodo: it is deposited in the
# EGA. See the Data Availability statement in the manuscript.
#-------------------------------------------------------------------------------

# Zenodo record for the processed data: DOI 10.5281/zenodo.22754723
# Override at runtime with ZENODO_RECORD_ID=<id> if pointing at a new version.
ZENODO_RECORD_ID <- Sys.getenv("ZENODO_RECORD_ID", unset = "22754723")

args <- commandArgs(trailingOnly = TRUE)
opt_list  <- "--list"  %in% args
opt_force <- "--force" %in% args
opt_check <- "--check" %in% args
only_pat  <- if (any(grepl("^--only", args))) {
  v <- sub("^--only=?", "", args[grep("^--only", args)][1])
  if (!nchar(v)) v <- args[which(args == "--only") + 1]
  strsplit(v, ",")[[1]]
} else NULL

manifest_file <- "data/zenodo_manifest.csv"
if (!file.exists(manifest_file)) stop("Run this from the repository root (", manifest_file, " not found)")
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)

if (!is.null(only_pat)) {
  keep <- Reduce(`|`, lapply(only_pat, function(p) grepl(p, manifest$filename, fixed = TRUE)))
  manifest <- manifest[keep, ]
  if (!nrow(manifest)) stop("No files in the manifest match: ", paste(only_pat, collapse = ", "))
}

md5_ok <- function(path, want) file.exists(path) && as.character(tools::md5sum(path)) == want
status <- vapply(seq_len(nrow(manifest)), function(i) {
  if (!file.exists(manifest$path[i])) "missing"
  else if (md5_ok(manifest$path[i], manifest$md5[i])) "ok" else "CORRUPT"
}, character(1))

if (opt_list || opt_check) {
  cat(sprintf("%-56s %9s  %s\n", "file", "MB", "status"))
  for (i in seq_len(nrow(manifest)))
    cat(sprintf("%-56s %9.1f  %s\n", manifest$filename[i], manifest$bytes[i]/1e6, status[i]))
  cat("\n", sum(status == "ok"), " of ", nrow(manifest), " present and verified\n", sep = "")
  if (any(status == "CORRUPT")) cat("Re-download corrupt files with: Rscript download_data.R --force\n")
  quit(save = "no")
}

todo <- if (opt_force) seq_len(nrow(manifest)) else which(status != "ok")
if (!length(todo)) { cat("All ", nrow(manifest), " files already present and verified.\n", sep = ""); quit(save = "no") }

if (!nchar(ZENODO_RECORD_ID)) {
  stop("ZENODO_RECORD_ID is not set.\n",
       "  Set it at the top of this script, or pass it in the environment:\n",
       "    ZENODO_RECORD_ID=12345678 Rscript download_data.R\n",
       "  The record id is the number in the Zenodo URL and in the DOI 10.5281/zenodo.<id>.")
}

base_url <- paste0("https://zenodo.org/records/", ZENODO_RECORD_ID, "/files/")
cat("Downloading ", length(todo), " file(s) from Zenodo record ", ZENODO_RECORD_ID, "\n", sep = "")

# Large files over a slow link: no timeout, and resume is not assumed.
old_timeout <- getOption("timeout"); options(timeout = 60 * 60); on.exit(options(timeout = old_timeout))

failed <- character(0)
for (i in todo) {
  dest <- manifest$path[i]
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  url  <- paste0(base_url, utils::URLencode(manifest$filename[i]), "?download=1")
  cat(sprintf("  [%d/%d] %s (%.1f MB)\n", match(i, todo), length(todo), manifest$filename[i], manifest$bytes[i]/1e6))
  tmp <- paste0(dest, ".part")
  ok <- tryCatch({ utils::download.file(url, tmp, mode = "wb", quiet = TRUE); TRUE },
                 error = function(e) { cat("        download failed: ", conditionMessage(e), "\n", sep = ""); FALSE })
  if (!ok) { unlink(tmp); failed <- c(failed, manifest$filename[i]); next }
  if (!md5_ok(tmp, manifest$md5[i])) {
    cat("        checksum mismatch - file not kept\n"); unlink(tmp)
    failed <- c(failed, manifest$filename[i]); next
  }
  file.rename(tmp, dest)
}

if (length(failed)) {
  cat("\n", length(failed), " file(s) failed:\n  ", paste(failed, collapse = "\n  "), "\n", sep = "")
  cat("Re-run to retry; files that succeeded are not downloaded again.\n")
  quit(save = "no", status = 1)
}
cat("\nDone. Verify at any time with: Rscript download_data.R --check\n")

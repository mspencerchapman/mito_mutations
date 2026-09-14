#-------------------------------------------------------------------------------
# make_zenodo_manifest.R
#
# Regenerates data/zenodo_manifest.csv - the list of large data objects that are
# hosted on Zenodo rather than in this git repository.
#
# Run from the repository root whenever the deposited data changes, then
# re-upload the staged files (see prepare_zenodo_upload.sh) and publish a new
# Zenodo version.
#
# The split is by role, not merely by size:
#   on Zenodo - processed per-cohort data objects and simulation outputs
#   in repo   - metadata, references, phylogenies, and analysis products
#               (including the ABC posterior tables)
#-------------------------------------------------------------------------------

root_dir <- normalizePath(".")

zenodo_files <- c(
  "data/mito_data.Rds",
  "data/mito_data_foetal.Rds",
  "data/whole_genome_coverage_pileup_and_bedtools_annotated.csv",
  list.files("data/nonblood", pattern = "^mito_mutation_data_.*\\.RDS$", full.names = TRUE),
  list.files("data/colony_drift_simulations", pattern = "^model[0-9]+_res\\.Rds$", full.names = TRUE)
)

descriptions <- c(
  mito_data.Rds = "Processed blood/normal haematopoiesis cohort: VAF matrices, phylogenies (tree, tree.ultra), shearwater calls, germline assignments",
  mito_data_foetal.Rds = "Processed foetal cohort",
  whole_genome_coverage_pileup_and_bedtools_annotated.csv = "Per-sample mtDNA copy number from whole-genome coverage",
  model1_res.Rds = "Colony growth drift simulation, model 1",
  model2_res.Rds = "Colony growth drift simulation, model 2",
  model3_res.Rds = "Colony growth drift simulation, model 3",
  model4_res.Rds = "Colony growth drift simulation, model 4"
)
cohort_labels <- c(blood="blood", CML="CML", HL="HL (colon)", KY="KY (lung organoid)",
                   LM="LM (endometrium)", lymph="lymphoid", NW="NW (blood MPN)",
                   PR="PR (MUTYH mutant)", SO="SO (colon IBD)")

describe <- function(f) {
  b <- basename(f)
  if (b %in% names(descriptions)) return(unname(descriptions[b]))
  coh <- sub("^mito_mutation_data_(.*)\\.RDS$", "\\1", b)
  if (coh %in% names(cohort_labels)) return(paste0("Processed mtDNA mutation data, ", cohort_labels[[coh]], " cohort"))
  ""
}

stopifnot(all(file.exists(zenodo_files)))
if (any(duplicated(basename(zenodo_files)))) stop("Zenodo stores files by name only - basenames must be unique")

cat("Computing checksums for", length(zenodo_files), "files...\n")
manifest <- data.frame(
  filename    = basename(zenodo_files),                    # name as stored on Zenodo
  path        = zenodo_files,                              # destination, relative to repo root
  bytes       = file.size(zenodo_files),
  md5         = vapply(zenodo_files, function(f) as.character(tools::md5sum(f)), character(1)),
  description = vapply(zenodo_files, describe, character(1)),
  row.names   = NULL, stringsAsFactors = FALSE
)
manifest <- manifest[order(-manifest$bytes), ]

write.csv(manifest, "data/zenodo_manifest.csv", row.names = FALSE)
cat("Wrote data/zenodo_manifest.csv:", nrow(manifest), "files,",
    round(sum(manifest$bytes)/1e6), "MB\n")

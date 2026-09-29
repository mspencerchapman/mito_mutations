#-----------------------------------------------------------------------------------#
# Generate_SuppTable_donor_metadata.R
#
# Donor-level metadata for every individual in the study, blood and comparator
# tissues together, written to tables/donor_level_metadata.csv.
#
# Supplementary Table 1 summarises the studies. This is the per-donor companion:
# the covariate-relevant characteristics (age, sex, disease status) that the
# reporting summary asks for, alongside the unit of study for each donor - the
# number of clonal samples sequenced and the number of mtDNA mutations called.
#
# Sources:
#   data/non_blood_metadata.xlsx   age, sex, disease, for all cohorts including blood
#   data/mito_data.Rds             the blood phylogenies and matrices
#   data/nonblood/mito_mutation_data_<cohort>.RDS   the comparator phylogenies
#-----------------------------------------------------------------------------------#

suppressMessages({library(dplyr); library(tidyr); library(readxl); library(readr)})
options(stringsAsFactors = FALSE)
source(here::here("config.R"))   #root_dir

#Blacklisted recurrent artefacts, as elsewhere in the analysis
exclude_muts <- c("MT_302_A_C","MT_311_C_T","MT_456_C_T","MT_567_A_C","MT_574_A_C",
                  "MT_8270_C_T","MT_16170_A_C","MT_16181_A_C","MT_16182_A_C",
                  "MT_16183_A_C","MT_16189_T_C")
mutCN_cutoff <- 25

#Cohort codes to the dataset names used in Supplementary Table 1
cohort_names <- c(
  blood = "Normal aging hematopoiesis / fetal hemopoiesis",
  lymph = "Lymphoid",
  NW    = "Myeloproliferative neoplasms",
  LM    = "Endometrium",
  KY    = "Bronchial epithelium",
  HL    = "Normal colon",
  SO    = "IBD-affected colon",
  PR    = "MUTYH-mutated colon")

tissue_names <- c(
  Adult_Cord = "Blood (adult and cord)",
  Foetal     = "Blood (fetal liver / bone marrow)",
  MPN        = "Blood",
  Endometrium= "Endometrial glands",
  Bronchial  = "Bronchial epithelium",
  Colon_normal = "Colonic crypts",
  Colon_IBD    = "Colonic crypts (IBD-affected)",
  Colon_MAP    = "Colonic crypts (MUTYH-associated)")

#-----------------------------------------------------------------------------------#
### Per-donor counts from the data itself
#-----------------------------------------------------------------------------------#

#Distinct mtDNA mutations called in a donor, after the standard filters
count_muts <- function(list) {
  m <- list$matrices
  keep <- m$implied_mutCN > mutCN_cutoff |
    (matrix(!rownames(m$vaf) %in% list$CN_correlating_muts, ncol = 1) %*%
     matrix(rep(1, ncol(m$vaf)), nrow = 1))
  filt <- m$vaf * m$SW * keep
  if ("global" %in% colnames(filt)) filt <- filt[, setdiff(colnames(filt), "global"), drop = FALSE]
  muts <- rownames(filt)[rowSums(filt > 0, na.rm = TRUE) > 0]
  length(setdiff(muts[!grepl("DEL|INS", muts)], exclude_muts))
}

counts_from <- function(dat, cohort) {
  bind_rows(Map(list = dat, id = names(dat), function(list, id) {
    if (is.null(list$tree)) return(NULL)
    data.frame(Cohort = cohort, ID = id,
               n_samples = length(list$tree$tip.label),
               n_mtDNA_mutations = count_muts(list))
  }))
}

cat("Counting samples and mutations per donor\n")
blood_counts <- counts_from(readRDS(paste0(root_dir, "/data/mito_data.Rds")), "blood")

nonblood_counts <- bind_rows(lapply(c("lymph","NW","LM","KY","HL","SO","PR"), function(ds) {
  f <- paste0(root_dir, "/data/nonblood/mito_mutation_data_", ds, ".RDS")
  if (!file.exists(f)) { cat("  no data for", ds, "- skipped\n"); return(NULL) }
  counts_from(readRDS(f), ds)
}))

all_counts <- bind_rows(blood_counts, nonblood_counts)

#-----------------------------------------------------------------------------------#
### Join to the recorded metadata
#-----------------------------------------------------------------------------------#

#The blood donors are "8 pcw"/"18 pcw" in the metadata and "8pcw"/"18pcw" in the
#data objects; match on a whitespace-stripped key
meta <- read_excel(paste0(root_dir, "/data/non_blood_metadata.xlsx")) %>%
  dplyr::select(Cohort, Tissue_type, ID, Sex, Age, Disease) %>%
  mutate(key = gsub("\\s", "", ID))

#Join on cohort as well as ID: the same donor identifier occurs in more than
#one cohort, so matching on ID alone duplicates those donors
donor_metadata <- all_counts %>%
  mutate(key = gsub("\\s", "", ID)) %>%
  left_join(meta %>% dplyr::select(-ID), by = c("Cohort", "key")) %>%
  transmute(
    Dataset   = unname(cohort_names[Cohort]),
    Tissue    = ifelse(is.na(Tissue_type), NA_character_,
                       ifelse(Tissue_type %in% names(tissue_names),
                              unname(tissue_names[Tissue_type]),
                              gsub("_", " ", Tissue_type))),
    Donor_ID  = ID,
    Sex       = recode(Sex, F = "Female", M = "Male", .missing = NA_character_),
    #Fetal donors carry a negative age, in years relative to birth. Their
    #gestational age is taken from the donor identifier rather than converted
    #from that figure, which is approximate and gives 8.8 and 18.6 pcw.
    Age_years = ifelse(Age < 0, NA_real_, round(Age, 1)),
    Gestation_pcw = ifelse(grepl("^[0-9]+ ?pcw$", ID),
                           as.numeric(sub("^([0-9]+) ?pcw$", "\\1", ID)), NA_real_),
    Disease_or_status = Disease,
    N_clonal_samples  = n_samples,
    N_mtDNA_mutations = n_mtDNA_mutations) %>%
  arrange(Dataset, Age_years, Donor_ID)

#-----------------------------------------------------------------------------------#
### Report and write
#-----------------------------------------------------------------------------------#

cat("\nDonors:", nrow(donor_metadata),
    "| clonal samples:", sum(donor_metadata$N_clonal_samples),
    "| mtDNA mutations:", sum(donor_metadata$N_mtDNA_mutations), "\n\n")

print(donor_metadata %>%
        group_by(Dataset) %>%
        summarise(donors = n(),
                  samples = sum(N_clonal_samples),
                  mutations = sum(N_mtDNA_mutations),
                  age_range = ifelse(all(is.na(Age_years)), "fetal",
                                     paste0(min(Age_years, na.rm = TRUE), "-",
                                            max(Age_years, na.rm = TRUE))),
                  n_female = sum(Sex == "Female", na.rm = TRUE),
                  n_male = sum(Sex == "Male", na.rm = TRUE),
                  .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

missing <- donor_metadata %>% filter(is.na(Sex) | (is.na(Age_years) & is.na(Gestation_pcw)))
if (nrow(missing)) {
  cat("\nDonors with incomplete age or sex, to be filled in by hand:\n")
  print(as.data.frame(missing[, c("Dataset","Donor_ID","Sex","Age_years")]), row.names = FALSE)
}

out <- paste0(root_dir, "/tables/donor_level_metadata.csv")
write_csv(donor_metadata, out)
cat("\nWritten to", out, "\n")

#-----------------------------------------------------------------------------------#
# config.R
#
# The only file you should need to edit to run this analysis on your own machine.
#
# Every analysis and figure script begins with:
#
#     source(here::here("config.R"))
#
# which sets the project paths, the shared plotting theme, and a helper for
# installing/loading packages. `here` locates the repository root from any
# working directory, so the scripts run whether you launch them from the project
# root, from RStudio, or from a subdirectory.
#-----------------------------------------------------------------------------------#

if (!requireNamespace("here", quietly = TRUE)) install.packages("here", repos = "https://cloud.r-project.org")

options(stringsAsFactors = FALSE)

#-----------------------------------------------------------------------------------#
# 1. SETTINGS YOU MAY NEED TO CHANGE
#-----------------------------------------------------------------------------------#

# Reference genome (GRCh37/hs37d5) FASTA, with its .fai index alongside.
# Needed only by the scripts that read trinucleotide context - the mutational
# signature profiles. Everything else runs without it.
genomeFile <- "~/R_work/reference_files/genome.fa"

#-----------------------------------------------------------------------------------#
# 2. PROJECT PATHS - derived automatically, no need to edit
#-----------------------------------------------------------------------------------#

root_dir   <- here::here()
data_dir   <- file.path(root_dir, "data")
tables_dir <- file.path(root_dir, "tables")

# Trailing slash: the scripts build paths as paste0(plots_dir, "Figure_01/x.pdf"),
# so these must end in a separator.
plots_dir   <- paste0(file.path(root_dir, "plots"), "/")
figures_dir <- plots_dir   # figures/ was merged into plots/; kept as an alias

# Plots that do not appear in the manuscript figures - exploratory output, and
# analyses produced during review.
rebuttal_figs_dir <- file.path(plots_dir, "additional_plots")
dir.create(rebuttal_figs_dir, showWarnings = FALSE, recursive = TRUE)
rebuttal_figs_dir <- paste0(rebuttal_figs_dir, "/")

# Frequently used data files
mito_data_file      <- file.path(data_dir, "mito_data.Rds")             # blood / normal haematopoiesis
nonblood_dir        <- file.path(data_dir, "nonblood")                  # per-cohort cross-tissue data
ref_file            <- file.path(data_dir, "Samples_metadata_ref.csv")  # blood individual metadata
nonblood_ref_file   <- file.path(data_dir, "metadata", "non_blood_metadata.xlsx")
colony_info_file    <- file.path(data_dir, "metadata", "colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt")
mito_cn_file        <- file.path(data_dir, "whole_genome_coverage_pileup_and_bedtools_annotated.csv")
mtdna_trinuc_freq_path <- file.path(data_dir, "mtDNA_trinuc_freqs_coding_dloop_heavy_light.Rds")
functions_file      <- file.path(data_dir, "mito_mutations_blood_functions.R")

genomeFile <- path.expand(genomeFile)

#-----------------------------------------------------------------------------------#
# 3. SHARED PLOTTING THEME
#
# Sized for the small multi-panel figures used throughout (5-8 pt text).
#-----------------------------------------------------------------------------------#

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  my_theme <- theme(text = element_text(family = "Helvetica"),
                    axis.text = element_text(size = 5),
                    axis.title = element_text(size = 7),
                    legend.text = element_text(size = 5),
                    legend.title = element_text(size = 7),
                    strip.text = element_text(size = 7),
                    legend.spacing = unit(1, "mm"),
                    legend.key.size = unit(5, "mm")) +
    theme(legend.key.height = unit(3, "mm"),
          legend.title = element_text(size = 8))

  # Larger variant for the R Markdown notebooks. Their figures are read on
  # screen rather than printed at panel size, so the manuscript's 5-8 pt text is
  # too small to be legible. Roughly 1.5x throughout.
  my_markdown_theme <- theme(text = element_text(family = "Helvetica"),
                             axis.text = element_text(size = 8),
                             axis.title = element_text(size = 10),
                             legend.text = element_text(size = 8),
                             legend.title = element_text(size = 10),
                             strip.text = element_text(size = 11),
                             legend.spacing = unit(1.2, "mm"),
                             legend.key.size = unit(7, "mm")) +
    theme(legend.key.height = unit(4, "mm"),
          legend.title = element_text(size = 10))

  # Default figure size for notebook chunks that do not set their own. Call
  # set_markdown_figure_defaults() from a notebook's setup chunk.
  set_markdown_figure_defaults <- function(width = 7, height = 4.5) {
    if (requireNamespace("knitr", quietly = TRUE))
      knitr::opts_chunk$set(fig.width = width, fig.height = height)
    invisible(NULL)
  }
}

#-----------------------------------------------------------------------------------#
# 4. NOTEBOOK PLOT SAVING
#
# The R Markdown notebooks display every plot inline in their rendered .html, so
# they do not need to write files to be useful - and the figure panels for the
# manuscript are produced by generate_figures/ instead. Saving from a notebook is
# therefore off by default.
#
# Set save_plots <- TRUE to have the notebooks also write their plots to disk.
# They go to plots/notebook_output/<notebook>/, a separate namespace, so they can
# never overwrite the figure panels in plots/Figure_NN/.
#-----------------------------------------------------------------------------------#

save_plots <- FALSE

#' Output directory for one notebook's plots
notebook_plots_dir <- function(notebook) paste0(plots_dir, "notebook_output/", notebook, "/")

#' ggsave(), honouring save_plots and creating the directory if needed
save_plot <- function(filename, ...) {
  if (!isTRUE(save_plots)) return(invisible(NULL))
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(filename = filename, ...)
}

#' gganimate::anim_save(), honouring save_plots
save_anim <- function(filename, ...) {
  if (!isTRUE(save_plots)) return(invisible(NULL))
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  gganimate::anim_save(filename = filename, ...)
}

#' Open a pdf device, honouring save_plots
#'
#' When saving is off the device is opened on the null file, so that the plotting
#' calls and the matching dev.off() that follow still work unchanged.
save_pdf <- function(file, ...) {
  if (!isTRUE(save_plots)) return(invisible(grDevices::pdf(file = nullfile(), ...)))
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(file = file, ...)
}

#-----------------------------------------------------------------------------------#
# 5. HELPERS
#-----------------------------------------------------------------------------------#

#' Load packages, installing any that are missing
#'
#' @param cran,bioc Character vectors of package names.
#' @param github Named character vector: names are package names, values are
#'   "user/repo" (e.g. c(treemut = "nangalialab/treemut")).
load_packages <- function(cran = NULL, bioc = NULL, github = NULL) {
  for (p in cran) if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  if (length(bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    for (p in bioc) if (!requireNamespace(p, quietly = TRUE))
      BiocManager::install(p, ask = FALSE, update = FALSE)
  }
  if (length(github)) {
    if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools", repos = "https://cloud.r-project.org")
    for (p in names(github)) if (!requireNamespace(p, quietly = TRUE))
      devtools::install_github(github[[p]])
  }
  invisible(lapply(c(cran, bioc, names(github)),
                   function(p) suppressPackageStartupMessages(
                     library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))))
}

#' Create a figure output directory and return its path
#'
#' e.g. fig_dir("Figure_06") -> "<root>/plots/Figure_06/"
fig_dir <- function(name) {
  d <- file.path(plots_dir, name)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  paste0(d, "/")
}

#' Source the project's shared analysis functions
source_project_functions <- function() source(functions_file)

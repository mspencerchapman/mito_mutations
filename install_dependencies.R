#-------------------------------------------------------------------------------
# install_dependencies.R
#
# Installs the R packages required by the analysis scripts in this repository.
# Run once, from the repository root:
#
#   Rscript install_dependencies.R
#
# Package versions used for the published analysis are listed in
# SESSIONINFO.md. The analysis was developed under R 4.1 on the Sanger farm and
# re-run under R 4.6; no version-specific behaviour is relied upon.
#-------------------------------------------------------------------------------

cran_packages <- c(
  # data manipulation / plotting
  "dplyr","tidyr","tibble","readr","readxl","stringr","forcats","reshape2",
  "ggplot2","RColorBrewer","scales","ggsci","ggridges","ggpubr","ggpmisc",
  "gridExtra","dichromat","gganimate","knitr",
  # phylogenetics
  "ape","phangorn","phylobase","phylosignal",
  # modelling / inference
  "lme4","lmerTest","abc","VGAM","lsa",
  # utilities
  "optparse","ids","devtools","remotes","BiocManager"
)

bioc_packages <- c(
  "deepSNV",         # shearwater mutation calling
  "ComplexHeatmap"
)

# Not on CRAN/Bioconductor - installed from GitHub.
# Sources below were read back from the installed packages' DESCRIPTION files.
# NB hdp: the original nicolaroberts/hdp no longer compiles under modern R - it
# uses the legacy Free()/Calloc() macros that R replaced with R_Free()/R_Calloc(),
# which current clang rejects. Use the Sanger fork, which is patched for this.
github_packages <- c(
  treemut   = "nangalialab/treemut",          # phylogeny mutation assignment
  rsimpop   = "nangalialab/rsimpop",          # clonal expansion simulation (Fig 6c)
  dndscv    = "im3sanger/dndscv",             # selection analysis
  BuenColors= "caleblareau/BuenColors",       # palettes
  mitovizR  = "robertopreste/mitovizR",       # mtDNA visualisation
  hdp       = "NickWilliamsSanger/hdp"        # mutational signature extraction (fork that builds on modern R)
)

# Returns the packages that are still missing after attempting installation, so
# the caller can report them. Without this the failures are silent and the
# script claims success while leaving dependencies uninstalled.
install_if_missing <- function(pkgs, installer) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (!length(missing)) { cat("  all present\n"); return(invisible(character(0))) }
  cat("  installing:", paste(missing, collapse = ", "), "\n")
  for (p in missing) try(installer(p))
  still_missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  invisible(still_missing)
}

cat("CRAN packages:\n")
failed <- install_if_missing(cran_packages, function(p) install.packages(p, repos = "https://cloud.r-project.org"))

cat("Bioconductor packages:\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
failed <- c(failed, install_if_missing(bioc_packages, function(p) BiocManager::install(p, ask = FALSE, update = FALSE)))

cat("GitHub packages:\n")
# remotes, not devtools: install_github() lives in remotes, and current devtools
# only Suggests it, so devtools alone errors with 'The package "remotes" is required'.
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
failed <- c(failed, install_if_missing(names(github_packages),
                                       function(p) remotes::install_github(github_packages[[p]])))

if (length(failed)) {
  cat("\nFAILED to install:", paste(failed, collapse = ", "), "\n")
  cat("Re-run this script, or install these manually, before running the analysis.\n")
  quit(status = 1)
}

cat("\nDone. Run sessionInfo() to record your versions.\n")

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
  "optparse","ids","devtools","BiocManager"
)

bioc_packages <- c(
  "deepSNV",         # shearwater mutation calling
  "ComplexHeatmap"
)

# Not on CRAN/Bioconductor - installed from GitHub.
# Sources below were read back from the installed packages' DESCRIPTION files,
# except hdp, which is not installed locally - confirm before relying on it.
github_packages <- c(
  treemut   = "nangalialab/treemut",          # phylogeny mutation assignment
  dndscv    = "im3sanger/dndscv",             # selection analysis
  BuenColors= "caleblareau/BuenColors",       # palettes
  mitovizR  = "robertopreste/mitovizR",       # mtDNA visualisation
  hdp       = "nicolaroberts/hdp"             # mutational signature extraction (source UNVERIFIED)
)

install_if_missing <- function(pkgs, installer) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (!length(missing)) { cat("  all present\n"); return(invisible(NULL)) }
  cat("  installing:", paste(missing, collapse = ", "), "\n")
  for (p in missing) try(installer(p))
}

cat("CRAN packages:\n")
install_if_missing(cran_packages, function(p) install.packages(p, repos = "https://cloud.r-project.org"))

cat("Bioconductor packages:\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
install_if_missing(bioc_packages, function(p) BiocManager::install(p, ask = FALSE, update = FALSE))

cat("GitHub packages:\n")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = "https://cloud.r-project.org")
install_if_missing(names(github_packages), function(p) devtools::install_github(github_packages[[p]]))

cat("\nDone. Run sessionInfo() to record your versions.\n")

#-----------------------------------------------------------------------------------#
# Generate_SourceData.R
#
# Source data workbooks for the main and Extended Data figures, written to
# source_data/ as one .xlsx per figure with one sheet per panel.
#
# Each figure script is sourced with ggsave() intercepted, so every panel it
# saves is captured together with the filename it was saved under - and those
# filenames already carry the panel letter (Fig1a., Fig2c., ...). The plotted
# values are then taken from the ggplot object and written out.
#
# Panels drawn with base graphics rather than ggplot - the phylogenies and
# mutation heatmaps of Fig. 3, Fig. 5a and Fig. 6b-d - have no single
# underlying table and are not included; the objects behind them are in the
# repository's data/ directory.
#-----------------------------------------------------------------------------------#

suppressMessages({library(dplyr); library(ggplot2); library(writexl); library(rlang)})
options(stringsAsFactors = FALSE)
source(here::here("config.R"))   #root_dir, plots_dir

out_dir <- paste0(root_dir, "/source_data/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#-----------------------------------------------------------------------------------#
### Capture every panel a figure script saves
#-----------------------------------------------------------------------------------#

.captured <- new.env(parent = emptyenv())

#Stand-in for ggsave: records the plot against the file it would have been
#written to, and writes nothing. Matches ggsave's argument order so that both
#ggsave(file, plot) and ggsave(filename=, plot=) are captured.
capture_ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  .captured$panels[[length(.captured$panels) + 1]] <-
    list(file = basename(filename), plot = plot)
  invisible(filename)
}

#Variables the plot actually maps: the plot-level aesthetics, every layer's
#aesthetics, and the faceting variables. Without this the sheets carry every
#column of the upstream data frame - 69 of them for Fig. 1a - rather than the
#handful of values a reader needs to redraw the panel.
mapped_vars <- function(p) {
  v <- character(0)
  grab <- function(m) if (!is.null(m)) for (q in m)
    v <<- c(v, tryCatch(all.vars(rlang::quo_get_expr(q)), error = function(e) character(0)))
  grab(p$mapping)
  for (lyr in p$layers) grab(lyr$mapping)
  fp <- p$facet$params
  for (k in c("facets", "rows", "cols"))
    if (!is.null(fp[[k]])) grab(fp[[k]])
  unique(v)
}

#The plotted values. Most panels carry their data on the plot; a few build it
#in the layers instead, so fall back to the first layer that has any.
panel_data <- function(p) {
  if (!inherits(p, "ggplot")) return(NULL)
  d <- tryCatch(p$data, error = function(e) NULL)
  if (is.null(d) || !is.data.frame(d) || !nrow(d)) {
    for (lyr in p$layers) {
      ld <- tryCatch(lyr$data, error = function(e) NULL)
      if (is.data.frame(ld) && nrow(ld)) { d <- ld; break }
    }
  }
  if (!is.data.frame(d) || !nrow(d)) return(NULL)
  keep <- intersect(mapped_vars(p), names(d))
  if (length(keep)) d <- d[, keep, drop = FALSE]
  #Excel cannot hold list-columns or S4; flatten anything that is not atomic
  d <- as.data.frame(d)
  for (j in names(d)) {
    if (!is.atomic(d[[j]])) d[[j]] <- vapply(d[[j]], function(x)
      paste(format(x), collapse = "; "), character(1))
    if (is.factor(d[[j]])) d[[j]] <- as.character(d[[j]])
  }
  d
}

#"Fig1a.mtDNA.copy.number...pdf" -> "Fig1a"; sheet names are capped at 31 chars
sheet_name <- function(file) {
  nm <- sub("\\.pdf$", "", file)
  panel <- regmatches(nm, regexpr("^(Fig|ExtDataFig)[0-9]+[a-z]?", nm))
  if (!length(panel)) panel <- substr(nm, 1, 20)
  rest <- sub("^(Fig|ExtDataFig)[0-9]+[a-z]?\\.?", "", nm)
  substr(gsub("[^A-Za-z0-9_ ]", "_", paste0(panel, if (nchar(rest)) paste0("_", rest))), 1, 31)
}

#kind is "main" or "ED"; the scripts name their panels FigNx. and ExtDataFigNx.
#respectively, which is what identifies a captured plot as a panel of a figure.
build <- function(fig_number, kind = "main") {
  stem   <- if (kind == "ED") "Generate_ExtData_Fig" else "Generate_Fig"
  prefix <- if (kind == "ED") "ExtDataFig" else "Fig"
  label  <- if (kind == "ED") paste("ED Fig", fig_number) else paste("Fig", fig_number)
  outnm  <- if (kind == "ED") paste0("SourceData_ExtendedDataFigure", fig_number, ".xlsx")
            else paste0("SourceData_Figure", fig_number, ".xlsx")
  script <- paste0(root_dir, "/generate_figures/", stem, fig_number, ".R")
  if (!file.exists(script)) { cat("  no script for", label, "\n"); return(invisible(NULL)) }
  .captured$panels <- list()

  e <- new.env(parent = globalenv())
  assign("ggsave", capture_ggsave, envir = e)
  #Base-graphics panels write themselves; send those devices to a temp file
  assign("pdf", function(file = NULL, ...) grDevices::pdf(tempfile(fileext = ".pdf"), ...), envir = e)
  suppressWarnings(suppressMessages(
    try(source(script, local = e, echo = FALSE), silent = TRUE)))
  while (!is.null(grDevices::dev.list())) grDevices::dev.off()

  sheets <- list(); skipped <- character(0)
  for (p in .captured$panels) {
    #Only panels of this figure, named FigNx. by the scripts
    if (!grepl(paste0("^", prefix, fig_number, "[a-z]"), p$file)) next
    d <- panel_data(p$plot)
    nm <- sheet_name(p$file)
    if (is.null(d)) { skipped <- c(skipped, p$file); next }
    if (nm %in% names(sheets)) nm <- substr(paste0(nm, "_2"), 1, 31)
    sheets[[nm]] <- d
  }

  if (!length(sheets)) { cat(sprintf("%s: no ggplot panels captured\n", label)); return(invisible(NULL)) }
  f <- paste0(out_dir, outnm)
  write_xlsx(sheets, path = f)
  cat(sprintf("%s: %d sheet(s) -> %s\n", label, length(sheets), basename(f)))
  for (nm in names(sheets)) cat(sprintf("    %-32s %d rows x %d cols\n", nm, nrow(sheets[[nm]]), ncol(sheets[[nm]])))
  if (length(skipped)) cat("    no tabular data for:", paste(skipped, collapse = ", "), "\n")
  invisible(f)
}

#-----------------------------------------------------------------------------------#
### Build them
#-----------------------------------------------------------------------------------#

cat("Writing source data to", out_dir, "\n\n--- Main figures ---\n")
for (n in 1:6) build(n, "main")
cat("\n--- Extended Data figures ---\n")
for (n in 1:9) build(n, "ED")
cat("\nDone. Panels drawn with base graphics rather than ggplot have no single\n")
cat("underlying table and are not included: Fig. 3, Fig. 5a, Fig. 6b-d, and\n")
cat("Extended Data Figs. 7 and 9, which are mutation heatmaps and phylogenies.\n")

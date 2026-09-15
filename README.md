# Mitochondrial mutation, drift and selection during human development and ageing

Code accompanying the manuscript *"Mitochondrial mutation, drift and selection
during human development and ageing"*.

This repository contains the mtDNA mutation-calling pipeline, the downstream
analysis scripts, and the scripts that generate each manuscript figure.

---

## Quick start

```bash
git clone <repository-url>
cd mito_mutations
Rscript install_dependencies.R     # one-off; see SESSIONINFO.md for versions
Rscript download_data.R            # fetch the processed data objects from Zenodo (~270 MB)
```

**Configuration.** Every script begins with `source(here::here("config.R"))`,
which sets the project paths, the shared plotting theme and a package helper.
`config.R` locates the repository root automatically, so there are no paths to
edit - with one exception. If you want to regenerate the mutational signature
profiles, set your reference genome in `config.R`:

```r
genomeFile <- "~/path/to/your/GRCh37/genome.fa"   # the only setting you may need to change
```

Everything else runs without it.

Figures can then be regenerated individually, e.g.:

```bash
Rscript generate_figures/Generate_Fig4.R
```

---

## System requirements

**Operating systems.** The analysis is plain R and has no OS-specific code.

| Tested on | R | Notes |
|---|---|---|
| macOS 26.6 (Apple silicon, `aarch64-apple-darwin23`) | 4.6.1 | environment for the re-run in this repository |
| Linux (x86_64) | 4.1 | environment of the original analysis |

Not tested on Windows. R >= 4.1 is expected to work; no version-specific
behaviour is relied upon.

**R packages.** 41 direct dependencies (33 CRAN, 2 Bioconductor, 6 GitHub),
which pull in 289 packages in total. They are listed with their sources in
[`install_dependencies.R`](install_dependencies.R); the exact versions used for
the published analysis are in [`SESSIONINFO.md`](SESSIONINFO.md).

**Non-R software.** None for the analysis or figure scripts. The upstream
variant-calling pipeline in `mtDNA_mutation_calling_pipeline/` additionally
needs `samtools`, `bedtools`, `java`/GATK and `perl`, and is not required to
reproduce any figure - the mutation calls it produces are included as data.

**Hardware.** No non-standard hardware. Runs on a normal desktop or laptop:

| | |
|---|---|
| RAM | 8 GB sufficient; peak observed 4.1 GB (largest figure script) |
| Disk | ~1 GB for the repository and data, plus ~1.2 GB for the R library |
| Cores | single-core; nothing here is parallelised |

The full ABC drift inference (`simulation_scripts_for_ABCs/`) is the one
compute-intensive step and was originally run on a compute cluster. Its outputs
are included in the repository, so the figures can be reproduced without
re-running it.

---

## Installation guide

```bash
git clone <repository-url>
cd mito_mutations
Rscript install_dependencies.R     # installs all R dependencies
```

The script installs only what is missing, reports anything that failed, and
exits non-zero if so.

**Typical install time.** ~2 minutes on the tested macOS machine (measured: 119 s
for a complete install into an empty R library, on a fast connection). Most of
that is download, since CRAN and Bioconductor ship pre-built binaries for macOS;
the six GitHub packages compile from source. On platforms without binary
packages - most Linux distributions - expect substantially longer, on the order
of 20-45 minutes, as every package is compiled.

To fetch the processed data objects (~270 MB, 16 files, hosted on Zenodo):

```bash
Rscript download_data.R
```

This is **not** needed for the demo below.

---

## Demo

Extended Data Fig. 12 is built entirely from data tracked in this repository, so
it runs immediately after installation with no data download.

```bash
Rscript generate_figures/Generate_ExtData_Fig12.R
```

**Expected output.** Eight PDF panels in `plots/Extended_Data_Figure_12/`,
named `*_rejection.pdf`: the posterior distributions of the mtDNA drift
parameter for each cohort (normal, MPN, MPN-non-coding, CML), under both the
sequential and individual ABC. These are the panels of Extended Data Fig. 12 in
the manuscript, and the copies already committed in that directory are the
reference - a fresh run reproduces them exactly (the only byte-level difference
is the creation date that PDFs embed). The same directory also holds the
equivalent panels from the neural-network ABC, for comparison.

**Expected run time.** ~1.5 seconds (measured: 1.4 s; peak memory 0.3 GB).

---

## Instructions for use

This repository is the analysis for a specific study rather than a general-purpose
tool, and most scripts assume this cohort's data structures and metadata. Two
components are directly reusable on other datasets:

**1. Drift inference by ABC.** Estimates the mtDNA drift parameter from a
phylogeny plus per-sample mutation VAFs, and is not specific to this cohort. See
**Analysis pipeline > Drift inference by ABC** below for the arguments. To apply
it to your own data you need a phylogeny (Newick), a matrix of mutation VAFs per
sample, and the cohort parameters - `units_per_year` (branch-length units: 1 for
a time-scaled tree) and the VAF threshold below which a mutation counts as
absent. These are set per cohort in `cohort_config` at the top of
`simulation_scripts_for_ABCs/Mitochondrial_drift_through_tree_ABC_SEQ_local.R`;
add an entry for your dataset rather than editing an existing one, as the values
are not interchangeable between cohorts.

**2. Mutation calling.** `mtDNA_mutation_calling_pipeline/` calls mtDNA variants
from BAMs with shearwater. It contains absolute paths to Sanger infrastructure
and LSF `bsub` submissions, so it needs adapting before running elsewhere (see
**Portability**).

The remaining analysis and figure scripts are study-specific: they can be read as
a record of how each result was produced, and re-run as-is to reproduce the
published figures, but are not designed to be pointed at other datasets.

**Reproducing the published figures.** Install, fetch the data, then run any
script in `generate_figures/`; each writes the panels for one manuscript figure
into `plots/Figure_NN/` or `plots/Extended_Data_Figure_NN/`. A typical script
takes seconds to about a minute (measured: 36 s for `Generate_Fig4.R`, the
heaviest, which runs a Wright-Fisher simulation). The narrative notebooks take
longer - a few minutes each, and up to ~50 minutes for `mtDNA_mut_phasing.Rmd`.
Scripts whose output is stochastic set a fixed seed, so repeated runs reproduce
the committed panels.

One known exception: `Generate_ExtData_Fig10.R` does not reproduce under the
Seurat version recorded in [`SESSIONINFO.md`](SESSIONINFO.md) (5.5.1).
`FindNeighbors()` returns an SNN graph in which every colony is a singleton, so
`FindClusters()` fails inside `GroupSingletons()` (`sample.int`: "invalid first
argument"). Setting `group.singletons=FALSE` avoids the error but returns a
single cluster containing every colony, so the panel cannot be reproduced by
forcing it through. The panel was originally produced in 2021-2022, under Seurat
4.x. Note that installing Seurat 4.x is not straightforward on a current system:
under R 4.6 its dependencies install but SeuratObject 4.1.4 and Seurat 4.4.0
themselves fail to build, so reproducing this panel needs an R and Seurat of that
period together - for example in a container - rather than just an older Seurat.
Extended Data Fig. 11 is a flow-sorting schematic and has no code.

---

## Repository structure

| Path | Contents |
|---|---|
| `data/` | Processed mutation calls, phylogenies, metadata and ABC posteriors (see **Data**) |
| `generate_figures/` | One script per manuscript figure - the entry point for reproducing results |
| `full_analysis_scripts/` | Full analyses behind each figure: compilation, drift, selection, signatures. `archive/` holds superseded versions, kept for provenance only |
| `simulation_scripts_for_ABCs/` | Forward simulations and approximate Bayesian computation for drift inference |
| `mtDNA_mutation_calling_pipeline/` | Upstream variant calling (shearwater, coverage, haplotype assignment) |
| `*.Rmd` | Narrative analysis notebooks with rendered `.html` output (see **Start here**) |
| `plots/`, `tables/` | Generated outputs. `plots/Figure_NN/` and `plots/Extended_Data_Figure_NN/` hold figure panels; `plots/additional_plots/` holds exploratory output that is not in the manuscript |

---

## Start here: the analysis notebooks

If you want to understand what the analysis does rather than regenerate a
specific figure, read the R Markdown notebooks first. They walk through the main
mutation analyses in order, with the reasoning written out alongside the code,
and each has a rendered `.html` you can read without running anything:

| Notebook | Covers |
|---|---|
| `mtDNA_mutations_blood.Rmd` | Coverage and copy number, mutation calling and filtering, burden with age, mutational signatures - the normal haematopoiesis dataset |
| `mtDNA_mutations_comparator_tissues.Rmd` | The same for the cross-tissue cohorts, plus tissue comparisons and heteroplasmic oocyte mutations |
| `Mitochondrial_drift_analysis.Rmd` | Drift inference from VAF distributions with age |
| `Nonblood_mtDNA_drift_analysis.Rmd` | Drift and homoplasmy across tissues |
| `mtDNA_mut_phasing.Rmd` | Phasing of mtDNA mutations |

Notebook figures use `my_markdown_theme` and larger default dimensions, since
they are read on screen rather than printed at panel size. The manuscript figure
scripts keep the smaller `my_theme`. Both are defined in `config.R`.

The notebooks display every plot inline, and do not write files by default -
the figure panels are produced by `generate_figures/` instead. To have a notebook
also save its plots, set `save_plots <- TRUE` in `config.R`; they are written to
`plots/notebook_output/<notebook>/`, kept separate so they cannot overwrite the
manuscript panels.

**These notebooks do not cover everything.** The ABC inference of drift rates
through phylogenies, the dN/dS selection analysis, the mutational signature
extraction itself, and the simulation work all live in the script folders below.
The notebooks are the readable entry point, not the complete record.

---

## Run order

Nothing needs to run in a fixed order to regenerate a *figure* - each
`generate_figures/` script is self-contained and reads the deposited data. The
ordering below matters only if you are rebuilding from scratch.

```
1. Rscript install_dependencies.R          # once
2. Rscript download_data.R                 # once, fetches the processed data
3. (optional) full_analysis_scripts/Compile_*.R
                                           # rebuilds the processed data objects
                                           # from raw calls - NOT needed if you
                                           # downloaded them in step 2
4. simulation_scripts_for_ABCs/Mitochondrial_drift_through_tree_ABC_SEQ_local.R
                                           # ABC drift inference, one run per
                                           # cohort; writes posterior tables
5. Rscript generate_figures/Generate_Fig1.R    ... Generate_Fig6.R
   Rscript generate_figures/Generate_ExtData_Fig1.R ... Generate_ExtData_Fig12.R
```

There is one script per manuscript figure: `Generate_Fig1.R` to
`Generate_Fig6.R`, and `Generate_ExtData_Fig1.R` to `Generate_ExtData_Fig12.R`
except `Fig11` - Extended Data Fig. 11 is a flow-sorting schematic with no code.

Only step 5 is needed to reproduce the figures from the deposited data. Steps 3
and 4 regenerate their own inputs and take considerably longer - the ABC is
roughly 35 minutes per cohort.

Two figure scripts depend on ABC output that is already included in the
repository (`data/Drift_ABC_clonal_expansions/`), so step 4 is only necessary if
you want to re-run the inference itself: `Generate_Fig6.R` (panel e) and
`Generate_ExtData_Fig12.R`.

---

## Analysis pipeline

The work runs in four stages. Stages 1-2 require the raw sequencing data and a
compute cluster; stages 3-4 run from the processed data in `data/`.

1. **Mutation calling** - `mtDNA_mutation_calling_pipeline/`
   Shearwater calling against matched normals, coverage summaries, contamination
   checks and germline haplotype assignment. Written for an LSF cluster and
   carries site-specific absolute paths (see **Portability**).

2. **Data compilation** - `full_analysis_scripts/Compile_*.R`
   Assembles per-cohort VAF matrices, phylogenies and metadata into the
   `mito_mutation_data_*.RDS` objects in `data/`.

3. **Analysis** - `full_analysis_scripts/`
   Mutation burden and signatures, drift, selection (dN/dS), lineage tracing.

4. **Figures** - `generate_figures/`
   One script per manuscript figure: `Generate_Fig1.R` - `Generate_Fig6.R` for the
   main figures, `Generate_ExtData_Fig*.R` for Extended Data figures. Outputs are
   written to `figures/Figure_NN/` and `figures/Extended_Data_Figure_NN/`, named by
   panel (e.g. `Fig4a.`, `ExtDataFig9c.`).

   Note that these scripts produce the individual **panels**, not the assembled
   figures: final composition, lettering and layout were done in Adobe Illustrator.
   Some figures are therefore one file per donor rather than per panel - Fig. 3, for
   example, writes one phylogeny PDF per donor, each of which is one of panels a-d.
   A few scripts also generate panels belonging to an Extended Data figure alongside
   their main figure (e.g. `Generate_Fig5.R` writes Extended Data Fig. 9).

### Drift inference by ABC

Effective mtDNA generation time is inferred by approximate Bayesian computation:
heteroplasmy is simulated forward through each donor's phylogeny under a
Wright-Fisher model, and simulations are accepted by distance to the observed
summary statistics.

```bash
cd simulation_scripts_for_ABCs

# sequential ABC - each mutation's posterior becomes the next mutation's prior
Rscript Mitochondrial_drift_through_tree_ABC_SEQ_local.R -t normal
Rscript Mitochondrial_drift_through_tree_ABC_SEQ_local.R -t MPN
Rscript Mitochondrial_drift_through_tree_ABC_SEQ_local.R -t MPN_nocoding

# individual ABC - every mutation fitted independently from the initial prior
Rscript Mitochondrial_drift_through_tree_ABC_SEQ_local.R -t normal -a individual
```

Key options: `-t` cohort, `-a` `sequential`|`individual`, `-i` simulations per
mutation (default 2e4, as published), `-j` single mutation index, `-c` cores,
`-s` seed. Each stage is cached, so an interrupted run resumes rather than
restarts.

Posteriors are written to `data/Drift_ABC_<cohort>/output/posterior_table_<mut>.Rds`.
`full_analysis_scripts/Plot_drift_seq_ABC_results.R` renders them, skipping any cohort not yet run.

The cohorts differ in three parameters that are **not interchangeable** - tree
time units, heteroplasmy detection threshold, and where the phylogeny is stored.
These are recorded in the `cohort_config` block at the top of the ABC script.

---

## Data

Data for this project lives in three places:

| What | Where | How to get it |
|---|---|---|
| Raw sequencing data | European Genome-phenome Archive (EGA) | Managed access - see the manuscript's Data Availability statement |
| Processed data objects (~270 MB) | Zenodo [10.5281/zenodo.22754723](https://doi.org/10.5281/zenodo.22754723) | `Rscript download_data.R` |
| Metadata, phylogenies, references, analysis products | This repository | included in the clone |

### Fetching the processed data

The bulk processed objects are too large to sit comfortably in git, so they are
deposited on Zenodo and downloaded on demand:

```bash
Rscript download_data.R            # fetch everything missing, verifying MD5s
Rscript download_data.R --list     # what is needed, and what is already present
Rscript download_data.R --only NW  # just one cohort, if that is all you need
Rscript download_data.R --check    # verify checksums of what you have
```

`data/zenodo_manifest.csv` is the authoritative list - filename, destination
path, size, MD5 and a description of each object. It is tracked in git, so the
expected contents are visible without downloading anything.

Each `mito_mutation_data_<cohort>.RDS` is a per-donor list containing VAF, depth
and shearwater matrices, the sample phylogeny (`tree` and the ultrametric
`tree.ultra`), germline assignments and copy-number-correlating mutations.
Cohort codes: `blood`, `NW` (blood MPN), `CML`, `lymph`, `KY` (lung organoid),
`HL` (colon), `SO` (colon IBD), `PR` (MUTYH mutant), `LM` (endometrium).

### Updating the deposited data

```bash
Rscript make_zenodo_manifest.R     # recompute sizes and checksums
bash prepare_zenodo_upload.sh      # stage the files flat, ready to upload
```

Then upload, publish a new Zenodo version, and set `ZENODO_RECORD_ID` in
`download_data.R`.

---

## Portability

Scripts in `mtDNA_mutation_calling_pipeline/` and some older analysis scripts
contain absolute paths to Sanger infrastructure (`/lustre/...`, `/nfs/...`) and
LSF `bsub` submissions. They are included for transparency and will need paths
adapted to run elsewhere. The analysis and figure scripts intended for re-use
resolve paths from `root_dir` instead.

---

## Citation

Please cite the manuscript. <!-- TODO: add full citation and DOI on acceptance -->

## License

The **code** in this repository is released under the
[MIT License](https://opensource.org/licenses/MIT), an OSI-approved licence -
see [LICENSE](LICENSE). You are free to use, modify and redistribute it,
including commercially, provided the copyright notice is retained.

The **processed data** deposited at
[Zenodo](https://doi.org/10.5281/zenodo.22754723), and the figures and
documentation here, are licensed under the
[Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) -
see [LICENSE-CC-BY-4.0](LICENSE-CC-BY-4.0). Creative Commons licences are not
intended for software, which is why the code carries a separate licence.

The R packages this code depends on (see `install_dependencies.R`) retain their
own licences; several are GPL-2 or GPL-3.

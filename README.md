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

This work is licensed under the
[Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
You are free to share and adapt the material for any purpose, including
commercially, provided you give appropriate credit. See [LICENSE](LICENSE).

The R packages this code depends on (see `install_dependencies.R`) retain their
own licences.

# Archived scripts

Superseded code, kept for provenance. **Nothing here needs to run** — each file
has been replaced by a current notebook or script, listed below. They are
retained so that the history of the analysis remains inspectable, not because
they are still part of the pipeline.

These scripts have **not** been updated to use `config.R`, and some still carry
old paths. Treat them as a record rather than as working code.

## Replaced by an R Markdown notebook

The notebooks in the repository root are the current version of these analyses.
The percentages are the proportion of each script's code that also appears in
its notebook, measured when the scripts were archived.

| Archived script | Superseded by | Overlap |
|---|---|---|
| `Nonblood_mutation_analysis_local_Apr2026.R` | `mtDNA_mutations_comparator_tissues.Rmd` | 98% |
| `Nonblood_mutation_analysis_local.R` | `mtDNA_mutations_comparator_tissues.Rmd` | 93% |
| `Nonblood_mutation_analysis.R` | `mtDNA_mutations_comparator_tissues.Rmd` | 90% |
| `Mitochondrial_drift_analysis.R` | `Mitochondrial_drift_analysis.Rmd` | 91% |
| `Mitochondrial_mut_analysis_Apr2026.R` | `mtDNA_mutations_blood.Rmd` | 87% |
| `mtDNA_mut_phasing.R` | `mtDNA_mut_phasing.Rmd` | 84% |
| `Mitochondrial_mut_analysis.R` | `mtDNA_mutations_blood.Rmd` | 81% |

The three `Nonblood_mutation_analysis` files are successive revisions of the
same analysis; `_local_Apr2026` was the last before it moved into the notebook.
The same applies to the two `Mitochondrial_mut_analysis` files.

## Replaced by another script

| Archived script | Superseded by |
|---|---|
| `Compile_nonblood_mito_data.R` | `../Compile_nonblood_mito_data_v2.R` |

## Never used

`mito_dNdS_simulation_NOT_USED/` — a simulation approach to mtDNA dN/dS that
did not end up in the manuscript.

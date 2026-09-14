# Project: Mitochondrial Drift Analysis

## Overview
This project analyzes variant allele frequency (VAF) data across tissues 
(blood, lymphoid, and other datasets) to study mitochondrial drift with age, 
using Poisson mixed-effects models and simulation-based inference.

## Code style
- Use tidyverse conventions (dplyr, ggplot2, magrittr pipes `%>%`)
- Prefer `glmer()` from lme4 for mixed models; Poisson family with log link 
  is standard here
- Continuous predictors (e.g. Age) should be scaled/centered before fitting 
  mixed models, to avoid eigenvalue/convergence warnings
- Keep random effects structures as simple as supportable by the data; 
  prefer `(1 | group)` or `(1 + x || group)` over full correlated random 
  slopes unless there's a clear justification
- Use `RColorBrewer` palettes for plots where relevant (e.g. "Spectral")
- Animations use `gganimate`; default renderer is `gifski_renderer()`

## File conventions
- Main analysis lives in `Mitochondrial_drift_analysis.Rmd`
- Plots are saved to a `plots_dir` variable, not hardcoded paths
- Use `n_near_homoplasmic_per_sample_FILT`-style descriptive object names 
  consistent with existing scripts

## What I want help with
- Adding roxygen-style or inline comments explaining statistical choices 
  (e.g. why a variable was scaled, why a random effect was simplified)
- Cleaning up repetitive plotting code (e.g. many similar ggplot chunks) 
  into reusable functions
- Flagging any model diagnostics that look off (convergence warnings, 
  overdispersion, VIF) and suggesting fixes consistent with prior fixes 
  in this project

## What to avoid
- Don't change variable/object names without asking, some are referenced 
  across multiple chunks
- Don't remove existing model specifications without flagging the reasoning 
  first

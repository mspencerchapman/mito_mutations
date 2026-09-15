# Session information

Package versions used for the re-run of the analysis in this repository,
recorded with `sessionInfo()` / `packageVersion()`. The original analysis was
developed under R 4.1 on Linux; it was re-run under the versions below.

Extended Data Fig. 10 does NOT reproduce under the Seurat version recorded below
(5.5.1) - see the note in generate_figures/Generate_ExtData_Fig10.R. That panel
was originally produced in 2021-2022, i.e. under Seurat 4.x (4.0 was released in
January 2021, 4.1 in early 2022; Seurat 5 reached CRAN in late 2023). The exact
patch version was not recorded. Seurat 4.x does not build under R 4.6, so that
panel needs an R and Seurat contemporary with the original analysis.

```
R version: 4.6.1
Platform:  aarch64-apple-darwin23
Running:   macOS Tahoe 26.6.2

dplyr          1.2.1
tidyr          1.3.2
tibble         3.3.1
readr          2.2.0
readxl         1.5.0
stringr        1.6.0
forcats        1.0.1
reshape2       1.4.5
ggplot2        4.0.3
RColorBrewer   1.1.3
scales         1.4.0
ggsci          5.2.0
ggridges       0.5.7
ggpubr         1.0.0
ggpmisc        1.0.0
gridExtra      2.3.1
dichromat      2.0.1
gganimate      1.0.11
knitr          1.51
ape            5.8.1
phangorn       2.12.1
phylobase      0.8.12
phylosignal    1.3.1
lme4           2.0.6
lmerTest       3.2.1
abc            2.2.2
VGAM           1.1.14
lsa            0.73.4
optparse       1.8.2
ids            1.0.1
devtools       2.5.2
remotes        2.5.0
BiocManager    1.30.27
deepSNV        1.58.0
ComplexHeatmap 2.28.0
treemut        1.2
rsimpop        2.4.0
dndscv         0.0.1.0
BuenColors     0.5.6
mitovizR       0.2.2
hdp            0.1.6
Seurat         5.5.1
```


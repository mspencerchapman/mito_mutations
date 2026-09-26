#-----------------------------------------------------------------------------------#
# Generate_Supplementary_Figs.R
#
# Figures that appear inside the Supplementary Notes rather than as main or
# Extended Data panels. Output goes to plots/Supplementary_Figures/.
#
# Registry of the supplementary figures and where each is produced:
#   Supp Fig 1  haplotype phylogenies (Note 1)          - mtDNA_mut_phasing.Rmd
#   Supp Fig 2  in vitro colony growth drift (Note 2)   - full_analysis_scripts/Drift_during_colony_growth.R
#   Supp Fig 3  HDP vs SigProfiler (Note 3)             - full_analysis_scripts/Mutational_signature_extraction_by_VAF.R
#   Supp Fig 4  heteroplasmic oocyte mutations (Note 6) - full_analysis_scripts/Heteroplasmic_oocyte_mutation_analysis.R
#   Supp Fig 5  drift from synonymous/non-coding muts   - THIS SCRIPT
#   Supp Fig 6  mature cell phenotyping (Note 11)        - gating strategy, not code-generated
#   Supp Fig 7  individual & sequential ABCs, 4 cohorts - THIS SCRIPT (was Extended Data Fig. 12)
#   Notes 10-11 flow cytometry gating                   - not code-generated
#-----------------------------------------------------------------------------------#

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(abc); library(lme4); library(dndscv)
})
options(stringsAsFactors = FALSE)
source(here::here("config.R"))   #root_dir, plots_dir, my_theme
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

supp_dir <- paste0(plots_dir,"Supplementary_Figures/")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

#-----------------------------------------------------------------------------------#
### Shared inputs
#-----------------------------------------------------------------------------------#
ref_df <- read.csv(paste0(root_dir,"/data/Samples_metadata_ref.csv")) %>% filter(Dataset!="Lymphocyte")
#Per-individual colours, matching the rest of the paper (see Generate_Fig4.R)
Individual_cols <- RColorBrewer::brewer.pal(12,"Paired")
names(Individual_cols) <- ref_df$Sample[order(ref_df$Age)]
mito_data <- readRDS(paste0(root_dir,"/data/mito_data.Rds"))
CN_correlating_muts <- readRDS(paste0(root_dir,"/data/CN_correlation.RDS"))
exclude_muts <- c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C",
                  "MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C")
mtDNA_CN <- 600
mutCN_cutoff <- 25

new_VAF_groups <- c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%",
                    "3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
VAF_categories_to_include_in_ABC <- new_VAF_groups[3:11]  #lowest two bins unreliable
vaf_bins <- 2^(-10:0)

nsamp_df <- bind_rows(Map(list=mito_data, Exp_ID=names(mito_data), function(list,Exp_ID)
  data.frame(exp_ID=Exp_ID, n_samp=length(colnames(list$matrices$SW))))) %>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))

#The ABC reference table. Empty VAF bins are stored as NA but mean zero mutations
#in that bin; left as NA, abc() returns an NA distance and the simulation can
#never be selected, which silently removes almost all low-generation simulations.
#See the equivalent repair in Generate_Fig4.R.
all_sumstats <- readRDS(paste0(root_dir,
  "/data/Drift_ABC_VAF_distribution/VAF_distribution_ABC_simulation_sumstats_combined.Rds"))
all_sumstats[,new_VAF_groups] <- lapply(all_sumstats[,new_VAF_groups],
                                        function(x) ifelse(is.na(x), 0, x))

#-----------------------------------------------------------------------------------#
### Generate SUPPLEMENTARY FIG. 5 ---------
### Drift inferred from synonymous and non-coding mutations only
#-----------------------------------------------------------------------------------#
# Sensitivity analysis for Supplementary Note 7: selection on protein-altering
# variants could in principle bias a drift estimate built on a neutral model, so
# the inference is repeated using only mutations that do not change the protein.
# Rejection sampling is used throughout, as in Fig. 4.

df_tidy <- bind_rows(Map(list=mito_data, exp_ID=names(mito_data), function(list,exp_ID){
  keep <- list$matrices$implied_mutCN > mutCN_cutoff |
    (matrix(!rownames(list$matrices$vaf) %in% CN_correlating_muts, ncol=1) %*%
     matrix(rep(1,ncol(list$matrices$vaf)), nrow=1))
  (list$matrices$vaf * list$matrices$SW * (list$matrices$ML_Sig=="N1") * keep) %>%
    as.data.frame() %>% tibble::rownames_to_column("mut_ref") %>%
    dplyr::select(-global) %>% gather("Sample","vaf",-mut_ref) %>%
    filter(!grepl("DEL|INS",mut_ref), !mut_ref %in% exclude_muts, vaf > 0) %>%
    mutate(exp_ID = exp_ID)
}))

#Annotate coding consequence with dndscv (ND6 excluded - light strand)
target_genes <- c("MT-CYB","MT-ND5","MT-ND2","MT-ND4","MT-ND1","MT-CO3",
                  "MT-ATP6","MT-ND3","MT-ATP8","MT-ND4L","MT-CO2","MT-CO1")
mtref_rda_path <- paste0(root_dir,"/data/mtref.rda")
annot <- dndscv(df_tidy %>% separate(mut_ref,c("chr","pos","ref","alt"),"_",remove=FALSE) %>%
                  dplyr::select("sampleID"=Sample,chr,pos,ref,alt),
                gene_list=target_genes, refdb=mtref_rda_path, numcode=2,
                max_coding_muts_per_sample=Inf, max_muts_per_gene_per_sample=Inf)$annotmuts %>%
  unite("mut_ref",chr,pos,ref,mut,sep="_") %>%
  dplyr::rename("Sample"=sampleID) %>% dplyr::select(Sample,mut_ref,impact)

df_annotated <- left_join(df_tidy, annot, by=c("Sample","mut_ref")) %>%
  replace_na(list(impact="Non-coding"))

#Summary statistics: mutations per sample in each VAF bin, counted once per donor
make_sumstats <- function(d) {
  bind_rows(lapply(names(mito_data), function(id)
    bind_rows(lapply(seq_along(vaf_bins), function(i) {
      lower <- if (i==1) 0 else vaf_bins[i-1]
      data.frame(exp_ID = id, VAF_range = new_VAF_groups[i],
                 abs_muts = d %>% filter(exp_ID==id, vaf>lower, vaf<=vaf_bins[i]) %>%
                   distinct(exp_ID,mut_ref) %>% nrow())
    })))) %>%
    left_join(nsamp_df, by="exp_ID") %>%
    mutate(muts_per_samp = abs_muts/n_samp) %>%
    dplyr::select(exp_ID,VAF_range,muts_per_samp) %>%
    pivot_wider(names_from="VAF_range", values_from="muts_per_samp") %>%
    mutate(across(where(is.numeric), ~replace_na(.x, 0)))
}

run_abc_posterior <- function(ss) {
  bind_rows(lapply(1:nrow(ss), function(i) {
    res <- abc(target = as.numeric(ss[i,VAF_categories_to_include_in_ABC]),
               param = all_sumstats[,1:2],
               sumstat = all_sumstats[,VAF_categories_to_include_in_ABC],
               tol = 0.05, transf = c("log","log"), method = "rejection")
    as.data.frame(res$unadj.values) %>% mutate(exp_ID = ss[[1]][i])
  })) %>%
    left_join(ref_df %>% dplyr::select(Sample,Age), by=c("exp_ID"="Sample")) %>%
    filter(!exp_ID %in% c("8pcw","18pcw"))   #foetal samples: WF assumptions do not hold
}

syn_posterior <- run_abc_posterior(make_sumstats(
  df_annotated %>% filter(impact %in% c("Synonymous","Non-coding"))))
#Matched control: the same pipeline on the complete mutation set, so that the two
#estimates differ only by the coding-consequence filter
all_posterior <- run_abc_posterior(make_sumstats(df_annotated))

report_drift <- function(post, tag) {
  m  <- lmer(total_generations ~ Age + (1|exp_ID), data=post)
  ml <- mtDNA_CN*365/m@beta[2]; ci <- mtDNA_CN*365/confint(m)["Age",]
  cat(sprintf("%-32s slope %.2f gen/yr | drift %.0f mitochondria days (95%% CI %.0f-%.0f)\n",
              tag, m@beta[2], ml, ci[2], ci[1]))
  list(model=m, ml=ml, ci=ci)
}
cat("\nSupplementary Fig. 5 - drift estimates (rejection ABC)\n")
syn_fit <- report_drift(syn_posterior, "synonymous + non-coding")
all_fit <- report_drift(all_posterior, "all mutations (control)")

#--- Supp Fig 5a: inferred mutation rate per donor -------------------------------
syn_rate <- syn_posterior %>% group_by(exp_ID) %>%
  summarise(Age = Age[1],
            med = median(muts_per_mitochondria_per_generation),
            lo  = quantile(muts_per_mitochondria_per_generation, 0.025),
            hi  = quantile(muts_per_mitochondria_per_generation, 0.975), .groups="drop")

SuppFig5a <- syn_rate %>%
  ggplot(aes(x=forcats::fct_reorder(exp_ID,Age), y=med, ymin=lo, ymax=hi))+
  geom_point(alpha=0.75,size=0.5)+
  geom_errorbar(width=0.3,alpha=0.5)+
  theme_classic()+
  #Three breaks only: five 6pt scientific labels collide at this panel width
  scale_y_continuous(breaks=c(1e-4,2e-4,3e-4),
                     labels=function(x) format(x, scientific=TRUE, digits=2))+
  labs(x="Individual", y="Mutations per\nmtDNA genome\nper generation")+
  my_theme+
  #As in Fig. 4e: after coord_flip() axis.title.y is the donor axis, which is
  #self-explanatory from the tick labels.
  theme(axis.title.y = element_blank())+
  coord_flip()
ggsave(paste0(supp_dir,"SuppFig5a.syn_noncoding_mutation_rate.pdf"), SuppFig5a, width=1.6, height=2)

#--- Supp Fig 5b: WF generations against age -------------------------------------
syn_band <- lmer_confidence_band(syn_fit$model, syn_posterior, predictor="Age")
SuppFig5b <- syn_posterior %>%
  mutate(exp_ID = factor(exp_ID, levels = ref_df$Sample[order(ref_df$Age)])) %>%
  ggplot(aes(x=Age,y=total_generations))+
  geom_ribbon(data=syn_band, aes(x=Age,ymin=lower,ymax=upper),
              inherit.aes=FALSE, fill="grey50", alpha=0.25)+
  geom_point(aes(col=exp_ID), alpha=0.05, size=0.25)+
  geom_abline(slope=syn_fit$model@beta[2], intercept=syn_fit$model@beta[1], linetype=1)+
  #Drop the first two colours: those belong to the foetal samples, excluded above
  scale_color_manual(values=Individual_cols[-c(1:2)])+
  scale_y_continuous(limits=c(0,1700))+
  theme_classic()+
  guides(colour=guide_legend(override.aes=list(alpha=1)))+
  labs(x="Age", y="Total WF generations\n(posterior distribution from ABC)", col="")+
  my_theme+
  theme(legend.key.size=unit(0.5,"mm"), legend.box.spacing=unit(0,"mm"))
ggsave(paste0(supp_dir,"SuppFig5b.syn_noncoding_generations_by_age.pdf"), SuppFig5b, width=2.6, height=2)

cat("Supplementary figures written to ", supp_dir, "\n", sep="")

#-----------------------------------------------------------------------------------#
### Generate SUPPLEMENTARY FIG. 7 ---------
### Individual and sequential ABCs for inference of mitochondrial drift
#-----------------------------------------------------------------------------------#
# a  normal blood
# b  MPN, including non-synonymous mutations
# c  MPN, non-synonymous mutations excluded
# d  CML
#
# For each cohort: the posterior for every mutation fitted independently from the
# initial prior ("individual"), and the sequential ABC in which each mutation's
# posterior becomes the next one's prior, refining yellow -> dark red. The final
# dark red distribution is the one shown in Fig. 6e.
#
# This was Extended Data Fig. 12; it moved to the Supplementary Information
# when the Extended Data figures were renumbered. The former two-panel Supp
# Fig 6 (MPN with and without non-synonymous mutations) is panels b and c here.
#
# Posteriors are produced by
#   simulation_scripts_for_ABCs/Mitochondrial_drift_through_tree_ABC_SEQ_local.R
# Cohorts with no posteriors on disk are skipped with a message.

suppressMessages({library(ggridges); library(RColorBrewer); library(readr)})

#Rejection posteriors, matching the Methods and the drift parameter quoted in the
#Results. Set to "neuralnet" to build from the regression-adjusted runs instead.
abc_method <- "rejection"
sfx <- if (abc_method=="rejection") "_rejection" else ""

#Shared prior, drawn once so the panels are directly comparable
set.seed(42)
prior_gt <- 10^runif(1e4, min=-1, max=2.7)

supp7_panels <- list(
  list(letter="a", cohort="normal",       muts="abc_muts_normal.csv"),
  list(letter="b", cohort="MPN",          muts="abc_muts_MPN.csv"),
  list(letter="c", cohort="MPN_nocoding", muts="abc_muts_MPN_nocoding.csv"),
  list(letter="d", cohort="CML",          muts="abc_muts_CML.csv"))
supp7_dir <- function(cohort, kind) paste0("Drift_ABC_clonal_expansions/Drift_ABC_",
                                           cohort,"_",kind,sfx)

load_posteriors <- function(dir, abc_muts) {
  d <- paste0(root_dir,"/data/",dir,"/output/")
  if (!dir.exists(d)) return(NULL)
  posts <- lapply(abc_muts$mut, function(mu) {
    f <- paste0(d,"posterior_table_",mu,".Rds")
    if (!file.exists(f)) return(NULL)
    data.frame(generation_time=readRDS(f)$generation_time, mut=mu)
  })
  if (all(sapply(posts,is.null))) return(NULL)
  bind_rows(c(list(data.frame(generation_time=prior_gt, mut="Prior")), posts)) %>%
    mutate(mut=factor(mut, levels=c("Prior", abc_muts$mut)))
}

#Sequential runs are ordered, so an ordered yellow->red ramp; individual runs are
#exchangeable, so a qualitative palette. scale>1 lets ridges overlap.
ridge_plot <- function(df, abc_type) {
  n <- length(levels(df$mut))
  cols <- if (abc_type=="sequential") colorRampPalette(brewer.pal(8,"YlOrRd"))(n)
          else colorRampPalette(brewer.pal(12,"Paired"))(n)
  names(cols) <- levels(df$mut)
  ggplot(df, aes(x=generation_time, y=mut, fill=mut))+
    geom_density_ridges(alpha=0.9, linewidth=0.3, col="black", scale=2)+
    scale_fill_manual(values=cols)+
    scale_x_log10(breaks=c(0.1,1,10,100), labels=c(0.1,1,10,100), limits=c(0.05,600))+
    theme_classic()+my_theme+
    theme(legend.position="none")+
    #The ridgeline axis is the mutation the posterior was fitted to (plus the
    #shared prior), not a count.
    labs(x="Generation time (days)", y="Mutation")
}

for (p in supp7_panels) {
  muts_file <- paste0(root_dir,"/simulation_scripts_for_ABCs/",p$muts)
  if (!file.exists(muts_file)) { cat("Supp Fig 7",p$letter,"- no manifest for",p$cohort,"- skipped\n"); next }
  abc_muts <- read_csv(muts_file, show_col_types=FALSE)
  ind <- load_posteriors(supp7_dir(p$cohort,"individual"), abc_muts)
  seq <- load_posteriors(supp7_dir(p$cohort,"sequential"), abc_muts)
  if (is.null(ind) && is.null(seq)) { cat("Supp Fig 7",p$letter,"- no posteriors for",p$cohort,"- skipped\n"); next }
  #Height scales with mutation count. The ED Fig. 12 version used
  #1.1 + 0.28*n for a full-page stack; these are tighter for the Supplementary.
  h <- 0.75 + 0.20*nrow(abc_muts)
  if (!is.null(ind)) ggsave(paste0(supp_dir,"SuppFig7",p$letter,".",p$cohort,"_individual_",abc_method,".pdf"),
                            ridge_plot(ind,"individual"), width=3.3, height=h)
  if (!is.null(seq)) ggsave(paste0(supp_dir,"SuppFig7",p$letter,".",p$cohort,"_sequential_",abc_method,".pdf"),
                            ridge_plot(seq,"sequential"), width=3.3, height=h)
  cat("Supp Fig 7",p$letter,"-",p$cohort,": individual",!is.null(ind)," sequential",!is.null(seq),"\n")
}

cat("Supplementary figures written to ", supp_dir, "\n", sep="")

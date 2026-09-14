#-----------------------------------------------------------------------------------#
# Generate_ExtData_Fig12.R
#
# Extended Data Fig. 12 | Individual and sequential approximate Bayesian
# computations (ABCs) for inference of mitochondrial drift.
#
#   a  normal blood
#   b  MPN - including non-synonymous mutations
#   c  MPN - non-synonymous mutations excluded
#   d  CML
#
# For each cohort: LEFT  = the posterior for each mutation fitted independently
#                          from the initial prior ("individual" ABC)
#                  RIGHT = the sequential ABC, where each mutation's posterior
#                          becomes the next mutation's prior. The gradual
#                          refinement runs yellow -> dark red, the dark red being
#                          the final result also shown in Fig. 6e.
#
# Posteriors are produced by
#   simulation_scripts_for_ABCs/Mitochondrial_drift_through_tree_ABC_SEQ_local.R
# Cohorts with no posteriors on disk are skipped with a message.
#-----------------------------------------------------------------------------------#

cran_packages=c("dplyr","ggplot2","ggridges","RColorBrewer","readr","gridExtra")
for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

options(stringsAsFactors = F)

#Set these file paths before running the script
source(here::here("config.R")) #sets root_dir, genomeFile, project paths and my_theme
plots_dir=paste0(root_dir,"/plots/")
dir.create(paste0(plots_dir,"Extended_Data_Figure_12"),showWarnings=FALSE,recursive=TRUE)
ed12_dir=paste0(plots_dir,"Extended_Data_Figure_12/")

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))+
  theme(legend.key.height=unit(3,"mm"),legend.title = element_text(size=8))

#The prior: log-uniform over 0.1-500 days, i.e. flat on the log10 x-axis used
#throughout. Drawn once and shared by every panel so they are directly comparable.
set.seed(42)
n_prior_draws=1e4
prior_gt<-10^runif(n_prior_draws,min=-1,max=2.7)

#Panel definitions: label, manifest, and the two ABC output directories.
#
#These point at the REJECTION posteriors, which is what the Methods describe
#("using the 'abc' R package, using the 'rejection' method") and what the drift
#parameter quoted in the Results is taken from. The neuralnet (regression
#adjusted) runs live in the directories without the _rejection suffix; set
#abc_method below to "neuralnet" to build the figure from those instead.
abc_method="rejection"
sfx<-if(abc_method=="rejection") "_rejection" else ""

panels<-list(
  list(letter="a",cohort="normal",       muts="abc_muts_normal.csv",       seq_dir=paste0("Drift_ABC_clonal_expansions/Drift_ABC_normal_sequential",sfx),     ind_dir=paste0("Drift_ABC_clonal_expansions/Drift_ABC_normal_individual",sfx)),
  list(letter="b",cohort="MPN",          muts="abc_muts_MPN.csv",          seq_dir=paste0("Drift_ABC_clonal_expansions/Drift_ABC_MPN_sequential",sfx),            ind_dir=paste0("Drift_ABC_clonal_expansions/Drift_ABC_MPN_individual",sfx)),
  list(letter="c",cohort="MPN_nocoding", muts="abc_muts_MPN_nocoding.csv", seq_dir=paste0("Drift_ABC_clonal_expansions/Drift_ABC_MPN_nocoding_sequential",sfx),   ind_dir=paste0("Drift_ABC_clonal_expansions/Drift_ABC_MPN_nocoding_individual",sfx)),
  list(letter="d",cohort="CML",          muts="abc_muts_CML.csv",          seq_dir=paste0("Drift_ABC_clonal_expansions/Drift_ABC_CML_sequential",sfx),            ind_dir=paste0("Drift_ABC_clonal_expansions/Drift_ABC_CML_individual",sfx))
)

#' Assemble prior + per-mutation posteriors for one ABC run
#'
#' Returns NULL if the run is not on disk, so the figure can be built from
#' whichever cohorts have been generated.
load_posteriors<-function(dir,abc_muts) {
  d<-paste0(root_dir,"/data/",dir,"/output/")
  if(!dir.exists(d)) return(NULL)
  posts<-lapply(abc_muts$mut,function(mu){
    f<-paste0(d,"posterior_table_",mu,".Rds")
    if(!file.exists(f)) return(NULL)
    data.frame(generation_time=readRDS(f)$generation_time,mut=mu)
  })
  if(all(sapply(posts,is.null))) return(NULL)
  dplyr::bind_rows(c(list(data.frame(generation_time=prior_gt,mut="Prior")),posts))%>%
    mutate(mut=factor(mut,levels=c("Prior",abc_muts$mut)))
}

#' Ridge plot of a set of posteriors
#'
#' Sequential runs are ordered (each posterior builds on the last), so they get an
#' ordered yellow->red ramp; individual runs are exchangeable, so they get a
#' qualitative palette. scale>1 lets the ridges overlap, to fit many mutations in
#' a short panel.
ridge_plot<-function(df,abc_type) {
  n<-length(levels(df$mut))
  cols<-if(abc_type=="sequential") colorRampPalette(RColorBrewer::brewer.pal(8,"YlOrRd"))(n)
        else colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(n)
  names(cols)<-levels(df$mut)
  ggplot(df,aes(x=generation_time,y=mut,fill=mut))+
    ggridges::geom_density_ridges(alpha=0.9,linewidth=0.3,col="black",scale=3)+
    scale_fill_manual(values=cols)+
    scale_x_log10(breaks=c(0.1,1,10,100),labels=c(0.1,1,10,100),limits=c(0.05,600))+
    theme_classic()+my_theme+
    theme(legend.position="none")+
    labs(x="Generation time (days)",y="Count")
}

for(p in panels) {
  muts_file<-paste0(root_dir,"/simulation_scripts_for_ABCs/",p$muts)
  if(!file.exists(muts_file)) {cat("ED Fig 12",p$letter,"- no manifest for",p$cohort,"- skipped\n"); next}
  abc_muts<-readr::read_csv(muts_file,show_col_types=FALSE)

  ind<-load_posteriors(p$ind_dir,abc_muts)
  seq<-load_posteriors(p$seq_dir,abc_muts)
  if(is.null(ind)&&is.null(seq)) {cat("ED Fig 12",p$letter,"- no posteriors for",p$cohort,"- skipped\n"); next}

  #Panel height scales with the number of mutations so the ridges stay legible
  h<-1.1+0.28*nrow(abc_muts)
  if(!is.null(ind)) ggsave(paste0(ed12_dir,"ExtDataFig12",p$letter,".",p$cohort,"_individual_",abc_method,".pdf"),ridge_plot(ind,"individual"),width=3.3,height=h)
  if(!is.null(seq)) ggsave(paste0(ed12_dir,"ExtDataFig12",p$letter,".",p$cohort,"_sequential_",abc_method,".pdf"),ridge_plot(seq,"sequential"),width=3.3,height=h)
  cat("ED Fig 12",p$letter,"-",p$cohort,": individual",!is.null(ind)," sequential",!is.null(seq),"\n")
}

cat("\nExtended Data Fig. 12 panels written to",ed12_dir,"\n")

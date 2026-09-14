#-------------------------------------------------------------------------------
# Plot_drift_seq_ABC_results.R
#
# Purpose: visualise the posterior distributions from the ABC (approximate
# Bayesian computation) inference of mitochondrial drift rates. This script does
# NOT run any inference - it reads pre-computed posterior tables written by the
# ABC pipeline and renders the figures.
#
# Two flavours of ABC are plotted for each setting:
#   "individual" - each mtDNA mutation is fitted independently from the same
#                  prior, so the panels are exchangeable replicates.
#   "sequential" - the posterior from one mutation becomes the prior for the
#                  next, so the panels are cumulative and should progressively
#                  narrow. This is why the two flavours get different colour
#                  ramps below (ordered YlOrRd vs. categorical Paired).
#
# Outputs (one pair per ABC run, plus one cross-setting comparison):
#   <abc_dir>/plots/seq_posts_<type>_<abc_type>.pdf
#   <abc_dir>/plots/seq_posts_ridges_<type>_<abc_type>.pdf
#   <comparison_plots_dir>/normal_disease_comparison_ridges_plot.pdf
#-------------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)

# The prior is drawn by Monte Carlo rather than evaluated analytically, so fix
# the seed to keep the plotted "Prior" band identical across panels and reruns.
set.seed(42)

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))+
  theme(legend.key.height=unit(3,"mm"),
        legend.title = element_text(size=8))

#-------------------------------------------------------------------------------
# Constants and paths
#-------------------------------------------------------------------------------

# The original ABC outputs lived on farm scratch and have been lost. Posteriors
# are now regenerated locally by
# simulation_scripts_for ABCs/Mitochondrial_drift_through_tree_ABC_MPN_SEQ_local.R,
# which writes the same <abc_dir>/output/posterior_table_<mut>.Rds layout under
# data/. Cohorts that have not yet been regenerated are skipped with a message
# rather than erroring, so this script runs against whatever is available.
study_dir<-paste0(path.expand("~/R_work/mito_mutations"),"/data/")

# Effective mtDNA copy number assumed by the forward simulations. Generation
# time is only identifiable jointly with copy number under the Wright-Fisher
# drift model, so posterior generation times are multiplied by this to express
# results on the (copy-number-adjusted) drift parameter scale.
mtDNA_CN_used_in_simulations=600

# Number of Monte Carlo draws used to represent the prior in the plots.
n_prior_draws=1e4

exp_df<-data.frame(dirs=paste0(study_dir,
                               c("Drift_ABC_normal_individual/",
                                 "Drift_ABC_normal_SEQ/",
                                 "Drift_ABC_MPN_individual/",
                                 "Drift_ABC_MPN/",
                                 "Drift_ABC_MPN_nocoding_individual/",
                                 "Drift_ABC_MPN_nocoding/",
                                 "Drift_ABC_CML_individual/",
                                 "Drift_ABC_CML/")),
           type=rep(c("normal","MPN","MPN_nocoding","CML"),each=2),
           abc_type=rep(c("individual","sequential"),times=4))

#-------------------------------------------------------------------------------
# The prior
#
# Generation time is given a log-uniform prior over 10^-1 to 10^2.7 days
# (~0.1 to 500 days): the plausible range spans several orders of magnitude and
# we have no reason to favour any scale within it, so uniform-on-the-log-scale
# is the appropriate weakly-informative choice (and is why every plot below uses
# scale_x_log10 - on that axis the prior is flat by construction).
# Starting VAF is uniform on (0, 1), i.e. uninformative.
#
# Drawn ONCE here and reused everywhere, so all panels show the same prior.
#-------------------------------------------------------------------------------

prior=data.frame(generation_time=10^runif(n_prior_draws,min=-1,max=2.7),
                 starting_vaf=runif(n_prior_draws))%>%
  mutate(mut="Prior")

#-------------------------------------------------------------------------------
# Helpers
#-------------------------------------------------------------------------------

#' Load the mutation manifest and posterior tables for a single ABC run
#'
#' @param abc_dir Directory of the ABC run; expected to contain an "abc_muts"
#'   csv (the mutations that were fitted, in the order they were fitted) and an
#'   "output/" subdirectory of posterior_table_<mut>.Rds files.
#' @return List with `abc_muts` (tibble) and `all_posts` (list of posteriors,
#'   in fitting order).
#' @details Uses absolute paths rather than setwd() so that the working
#'   directory is left untouched and the loop is not order-dependent.
load_abc_posteriors<-function(abc_dir) {
  muts_file<-list.files(path=abc_dir,pattern="abc_muts",full.names=TRUE)
  abc_muts<-readr::read_csv(muts_file)
  all_posts<-lapply(abc_muts$mut,function(mut) {
    readRDS(file.path(abc_dir,"output",paste0("posterior_table_",mut,".Rds")))
  })
  return(list(abc_muts=abc_muts,all_posts=all_posts))
}

#' Stack posteriors with the prior into a single long data frame for plotting
#'
#' @param all_posts List of posterior tables, in fitting order.
#' @param abc_muts Mutation manifest (must align with `all_posts`).
#' @param order_as_factor If TRUE, `order` is returned as a factor with levels
#'   0:n (0 = prior); if FALSE it is left numeric. Factor levels are needed when
#'   posteriors from different runs are later combined, so that the fitting
#'   order is not silently coerced or re-sorted.
#' @return Data frame with columns mut, order, starting_vaf, generation_time.
#'   `mut` is a factor with "Prior" as the first level, so it always plots at
#'   the top of the facet/ridge stack as the reference distribution.
build_comb_posts<-function(all_posts,abc_muts,order_as_factor=FALSE) {
  n_posts<-length(all_posts)
  make_order<-function(x) {if(order_as_factor) factor(x,levels=0:n_posts) else x}

  Map(post=all_posts,mut=abc_muts$mut,order=1:n_posts,function(post,mut,order) {
    as.data.frame(post)%>%
      dplyr::select(starting_vaf,generation_time)%>%
      mutate(mut=mut,order=make_order(order),.before=1)
  })%>%dplyr::bind_rows()%>%
    dplyr::bind_rows(prior%>%mutate(order=make_order(0)))%>%
    mutate(mut=factor(mut,levels=c("Prior",abc_muts$mut)))
}

#' Posterior summary of generation time: median and 95% credible interval
#'
#' @details Quantile-based (2.5%/97.5%) rather than an HPD interval, as the
#'   posteriors are skewed on the natural scale but roughly symmetric on the log
#'   scale on which they were inferred. Each summary is also rescaled by
#'   mtDNA_CN_used_in_simulations to give the drift parameter (see above).
summarise_posteriors<-function(comb_posts) {
  comb_posts%>%
    group_by(factor(order))%>%
    dplyr::summarise(lowerCI=quantile(generation_time,0.025),
                     median=median(generation_time),
                     upperCI=quantile(generation_time,0.975))%>%
    mutate(across(where(is.numeric), list(drift_param=function(x) x*mtDNA_CN_used_in_simulations)))
}

#-------------------------------------------------------------------------------
# Per-run plots: posterior of generation time for each mutation in turn
#-------------------------------------------------------------------------------

# Collect the summary tables so the numbers are available after the loop rather
# than being computed and discarded.
posterior_summaries<-list()

for(j in 1:nrow(exp_df)) {
  type=exp_df$type[j]
  abc_type=exp_df$abc_type[j]
  abc_dir=exp_df$dirs[j]

  plots_dir=paste0(abc_dir,"plots/")

  #Skip cohorts that have not been regenerated yet (see note on study_dir above)
  if(!dir.exists(paste0(abc_dir,"output"))) {
    cat("Skipping",type,abc_type,"- no posteriors at",abc_dir,"\n")
    next
  }
  dir.create(plots_dir,showWarnings=FALSE,recursive=TRUE)

  abc_data<-load_abc_posteriors(abc_dir)
  abc_muts<-abc_data$abc_muts
  all_posts<-abc_data$all_posts

  comb_posts<-build_comb_posts(all_posts,abc_muts)

  #Print the posterior parameter summary stats (i.e. median and 95% CI)
  posterior_summaries[[paste(type,abc_type,sep="_")]]<-summarise_posteriors(comb_posts)
  cat("\n---",type,abc_type,"---\n")
  print(posterior_summaries[[paste(type,abc_type,sep="_")]])

  seq_posts<-comb_posts%>%
    ggplot(aes(x=generation_time))+
    geom_histogram(bins=100,fill="pink",alpha=0.5,linewidth=0.1,col="gray")+
    scale_x_log10()+ #prior is flat on the log scale, so compare posteriors there
    facet_grid(rows=vars(mut))+
    theme_classic()+
    my_theme+
    theme(strip.text.y=element_text(angle=0))+
    labs(x="Generation time (days)",y="Count")

  ggsave(filename = paste0(plots_dir,"seq_posts_",type,"_",abc_type,".pdf"),seq_posts,width=3.5,height=5)

  # Sequential runs are ordered (each posterior builds on the last), so use a
  # sequential ramp; individual runs are exchangeable, so use a qualitative one.
  if(abc_type=="sequential") {
    mut_cols=colorRampPalette(RColorBrewer::brewer.pal(n=8,name="YlOrRd"))(length(all_posts)+1)
  } else if (abc_type=="individual") {
    mut_cols=colorRampPalette(RColorBrewer::brewer.pal(n=12,name="Paired"))(length(all_posts)+1)
  }

  names(mut_cols)<-levels(comb_posts$mut) #includes the "Prior" level
  seq_posts_ridges<-comb_posts%>%
    ggplot(aes(x=generation_time,y=mut,fill=mut))+
    ggridges::geom_density_ridges(alpha=0.9,linewidth=0.3,col="black",scale=3)+ #scale>1 lets ridges overlap, to fit many mutations in a small panel
    scale_fill_manual(values=mut_cols)+
    scale_x_log10(breaks=c(0.1,1,10,100),labels=c(0.1,1,10,100))+
    theme_classic()+
    my_theme+
    theme(strip.text.y=element_text(angle=0),legend.position = "none")+
    labs(x="Generation time (days)",y="Density")

  ggsave(filename = paste0(plots_dir,"seq_posts_ridges_",type,"_",abc_type,".pdf"),seq_posts_ridges,width=3.3,height=2)
}


#-------------------------------------------------------------------------------
# Plot the final posteriors for normal, MPN and CML drift rates against each other
#
# Only the LAST posterior of each sequential run is used: in a sequential ABC
# that posterior is conditioned on all mutations in the run, so it is the single
# summary of that setting. MPN_nocoding is deliberately excluded here - it is a
# sensitivity analysis of the MPN fit, not an independent setting.
#-------------------------------------------------------------------------------

plots_dir=paste0(study_dir,"comparison_plots/")
dir.create(plots_dir,showWarnings=FALSE,recursive=TRUE)

#Restrict to the settings whose sequential ABC has actually been regenerated.
comparison_settings<-c("normal","MPN","CML")
comparison_settings<-comparison_settings[sapply(comparison_settings,function(setting) {
  dir.exists(paste0(exp_df%>%filter(type==setting & abc_type=="sequential")%>%pull(dirs),"output"))
})]
cat("Cross-setting comparison using:",paste(comparison_settings,collapse=", "),"\n")

drift_rate_comparison<-lapply(comparison_settings, function(setting) {

  abc_dir<-exp_df%>%filter(type==setting & abc_type=="sequential")%>%pull(dirs)

  abc_data<-load_abc_posteriors(abc_dir)
  abc_muts<-abc_data$abc_muts
  all_posts<-abc_data$all_posts

  final_post_order<-nrow(abc_muts)

  comb_posts<-build_comb_posts(all_posts,abc_muts,order_as_factor=TRUE)

  final_post<-comb_posts%>%filter(order==final_post_order)%>%mutate(setting=setting)

  return(final_post)

})%>%dplyr::bind_rows()

normal_disease_comparison_ridges_plot<-drift_rate_comparison%>%
  dplyr::bind_rows(prior%>%mutate(order=factor(0),setting="prior"))%>%
  mutate(setting=factor(setting,levels=c("prior",comparison_settings)))%>% #prior first, so it reads as the reference row
  ggplot(aes(x=generation_time,y=setting,fill=setting))+
  ggridges::geom_density_ridges(linewidth=0.3)+
  scale_fill_brewer(palette = "RdPu")+
  scale_x_log10(breaks=c(0.1,1,10,100),labels=c(0.1,1,10,100))+
  theme_classic()+
  #ggridges::theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  my_theme+
  theme(strip.text.y=element_text(angle=0),legend.position = "none")+
  labs(x="Generation time (days)",y="Density")

ggsave(filename = paste0(plots_dir,"normal_disease_comparison_ridges_plot.pdf"),normal_disease_comparison_ridges_plot,width=2.5,height=2.2)

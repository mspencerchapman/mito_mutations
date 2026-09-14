#-----------------------------------------------------------------------------------#
# Generate_Fig6.R
#
# Fig. 6 | Understanding the implications of drift through simulation.
#
#   a  VAF distribution through life for a heteroplasmic mutation present in the
#      oocyte at 3.6% VAF (most lineages lose or fix it by old age)
#   b  The same mutation drifted through the KX004 phylogeny (illustrative)
#   c  Heteroplasmic drift through clonal expansions with contrasting dynamics
#      (2 yr = rapid growth, 35 yr = slow growth)
#   d  Heatmaps of the best shared marker mutations for the malignant clone in
#      one MPN and one CML individual        <-- SEE "PANEL D" BELOW, NOT YET WIRED UP
#   e  Posterior distributions of the inferred drift rate (Wright-Fisher
#      generation time) in CML, MPN and normal blood expanded clones
#
# Panels a-c are adapted from full_analysis_scripts/Implications_of_drift_analysis.R
# (which still carries the stale root_dir "~/R_work/mito_mutations_blood/").
# Panel e reads the ABC posteriors produced by
# simulation_scripts_for_ABCs/Mitochondrial_drift_through_tree_ABC_SEQ_local.R.
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# --------Load packages (and install if they are not installed yet)-------------------
#-----------------------------------------------------------------------------------#
cran_packages=c("devtools","ape","dplyr","tidyr","ggplot2","ggridges","RColorBrewer")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

# Panel c simulates clonal expansions under selection using rsimpop. This is the
# only figure script that needs it, so it is loaded lazily and panel c is skipped
# with a message if it is unavailable, rather than failing the whole script.
have_rsimpop <- requireNamespace("rsimpop", quietly = TRUE)
if(have_rsimpop) {library(rsimpop)} else {
  cat("NOTE: rsimpop is not installed - panel c will be skipped.\n",
      "      Install with: devtools::install_github('nangalialab/rsimpop')\n")
}

#-----------------------------------------------------------------------------------#
# ----------------------------------Set paths and import files------------------------
#-----------------------------------------------------------------------------------#

options(stringsAsFactors = F)

#Set these file paths before running the script
source(here::here("config.R")) #sets root_dir, genomeFile, project paths and my_theme
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R")) #supplies fisher_wright_drift, get_mito_mut_vaf_df, plot_tree, getTips, get_expanded_clade_nodes

#Set the key file paths using the root dir
plots_dir=paste0(root_dir,"/plots/")

#Create the figure output directories if they do not already exist
for(d in c("Figure_06")) dir.create(paste0(plots_dir,d),showWarnings=FALSE,recursive=TRUE)
fig6_dir=paste0(plots_dir,"Figure_06/")
dir.create(fig6_dir,showWarnings=FALSE,recursive=TRUE)

#Set the basic plotting theme for ggplot2
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

set.seed(42) #the panels below are stochastic simulations; fix the seed so the figure is reproducible

#-----------------------------------------------------------------------------------#
# Fig 6a | VAF distribution of a heteroplasmic oocyte mutation through life
#
# Drift is modelled in two phases because the rate is not constant through life:
#   development  - population size 450 x generation time 2.67 = 1,200 mito days
#                  (the drift parameter inferred for the 8pcw fetus)
#   post-natal   - population size 675 x generation time 23    = 15,500 mito days
#                  (the adult drift parameter)
# 274 days covers conception to birth; the remaining time is post-natal.
#-----------------------------------------------------------------------------------#

starting_oocyte_vaf=0.036
years_to_test=c(1,5,10,20,40,80)

sim_df<-lapply(years_to_test,function(years) {
  #Initial drift during development (here defined as pre-conception)
  vafs=sapply(1:5000,function(j) fisher_wright_drift(starting_oocyte_vaf,population_size=450,generation_time=2.67,total_time=274))
  #Subsequent drift during post-conception life
  final_vafs=Map(starting_vaf=vafs,years=years,f=function(starting_vaf,years) {
    vaf=fisher_wright_drift(starting_vaf,population_size=675,generation_time=23,total_time=((365*years) - 140))
    return(vaf)
  })
  return(data.frame(years=years,VAF=unlist(final_vafs)))
})%>%
  dplyr::bind_rows()

VAF.distribution.histogram<-sim_df%>%
  mutate(years=paste(years,"years"))%>%
  mutate(years=factor(years,levels=paste(years_to_test,"years")))%>%
  ggplot(aes(x=VAF))+
  geom_histogram(bins = 50,col="black",linewidth=0.2)+
  theme_classic()+
  facet_grid(rows=vars(years))+
  scale_fill_brewer(palette="Set3")+
  scale_y_log10()+ #the distribution becomes strongly bimodal (loss/fixation), so a log count axis keeps the middle visible
  my_theme+
  labs(x="VAF distribution",y="Count",fill="Time from\nconception\n(Years)")

ggsave(paste0(fig6_dir,"Fig6a.VAF_distribution_histogram.pdf"),VAF.distribution.histogram,height=3,width=2.5)

#-----------------------------------------------------------------------------------#
# Fig 6b | The same mutation drifted through a real adult phylogeny (KX004)
#
# Illustrative only: a single (adult) drift rate is used for simplicity, to show
# that independent drift down lineages fixes the mutation in several unrelated
# clades - i.e. shared heteroplasmy need not imply shared ancestry.
#-----------------------------------------------------------------------------------#

mito_data<-readRDS(paste0(root_dir,"/data/mito_data.Rds"))
exp_ID="KX004"
tree.ultra<-mito_data[[exp_ID]]$tree.ultra

vaf_df_het_oocyte=get_mito_mut_vaf_df(tree.ultra,
                                      node=tree.ultra$edge[1,1],
                                      starting_vaf=starting_oocyte_vaf,
                                      mito_cn=675,
                                      generation_time=23)
sample_vafs=vaf_df_het_oocyte$vaf[order(vaf_df_het_oocyte$node)][1:length(tree.ultra$tip.label)]
names(sample_vafs)<-tree.ultra$tip.label

pdf(file = paste0(fig6_dir,"Fig6b.example_phylogeny_hetmut_in_oocyte_simulation.pdf"),width = 7,height=2.5)
tree.ultra=plot_tree(tree.ultra,cex.label=0,bars=sample_vafs)
text(x = 0, y=-0.175*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",paste0("Max VAF: ",round(max(sample_vafs),digits=3)),pos = 4)
dev.off()

#-----------------------------------------------------------------------------------#
# Fig 6c | Drift through clonal expansions with contrasting dynamics
#
# Two expansions are simulated to the same final clone size but over different
# durations, by pairing the driver acquisition time with a fitness value: a long,
# gradual expansion (35 yr, lower fitness) and a rapid one (2 yr, higher fitness).
# The mtDNA copy number (1000) and generation time (20 days) are the highest
# density posterior estimates from the phylogeny-aware ABC, and the starting VAF
# in the MRCA is fixed at 0.2 for comparability across panels.
#-----------------------------------------------------------------------------------#

#Wrapper: detect the expanded clade, drop other tips, and drift from its MRCA
simulated_expansion_mito_vafs=function(tree,starting_vaf,mito_cn,generation_time,plot=T){
  expansion_node=get_expanded_clade_nodes(tree,min_clonal_fraction = 0.1)
  sub_tree<-drop.tip(tree,tree$tip.label[!tree$tip.label%in%c("s1",getTips(tree,expansion_node$nodes))])
  sub_tree$coords<-NULL
  starting_node<-sub_tree$edge[which(sub_tree$edge[,1]==sub_tree$edge[1,1]&!sub_tree$edge[,2]%in%1:length(sub_tree$tip.label)),2]
  vaf_df<-get_mito_mut_vaf_df(sub_tree = sub_tree,node=starting_node,starting_vaf=starting_vaf,
                              mito_cn=mito_cn,generation_time=generation_time)
  vaf_vec=vaf_df$vaf[order(vaf_df$node)][1:length(sub_tree$tip.label)]
  names(vaf_vec)=sub_tree$tip.label
  if(plot){
    sub_tree=plot_tree(sub_tree,cex.label=F,bars = vaf_vec)
    text(x = 0, y=-0.05*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",paste0("Max VAF: ",round(max(vaf_vec),digits=3)),pos = 4)
  }
  return(vaf_vec)
}

if(have_rsimpop) {
  #Each entry: years of expansion before sampling, and the driver fitness needed
  #to reach a comparable clone size in that time.
  #Parameters taken from full_analysis_scripts/Implications_of_drift_analysis.R.
  #Each scenario runs to 50 years and differs in when the driver is acquired, so
  #that the expansion window ranges from 35 years (slow, gradual growth) down to
  #2 years (rapid growth). Fitness rises as the window shortens, so that each
  #clone still reaches a comparable size by sampling.
  #The published Fig. 6c shows 35, 25, 15, 5 and 2 years; the original script also
  #simulates a 10 year scenario (driver at 40 yr, fitness 1.2), which is not shown.
  expansion_params<-list(list(years=35,nyears_driver_acquisition=15,fitness=0.3),
                         list(years=25,nyears_driver_acquisition=25,fitness=0.45),
                         list(years=15,nyears_driver_acquisition=35,fitness=0.6),
                         list(years=5, nyears_driver_acquisition=45,fitness=2.5),
                         list(years=2, nyears_driver_acquisition=48,fitness=5))

  #The driver clone can be lost to drift shortly after it is introduced, in which
  #case no expansion is present to simulate through - more likely for the short,
  #high-fitness scenarios. Retry until a clade reaches the 10% clonal fraction
  #that get_expanded_clade_nodes() looks for.
  max_attempts=20
  sim_tree_list.ultra<-lapply(expansion_params,function(p) {
    for(attempt in 1:max_attempts) {
      selsim=run_selection_sim(0.05,1/(2*190),
                               nyears_driver_acquisition = p$nyears_driver_acquisition,
                               target_pop_size = 5e4,nyears = 50,fitness=p$fitness)
      tree_m<-tryCatch({
        seltree100=get_subsampled_tree(get_tree_from_simpop(selsim),100)
        #Convert the tree to elapsed time, so branch lengths are comparable to the
        #molecular-time trees the drift model expects
        get_elapsed_time_tree(seltree100,mutrateperdivision=0.65,backgroundrate=16/365)
      },error=function(e) NULL)
      if(!is.null(tree_m) && length(get_expanded_clade_nodes(tree_m,min_clonal_fraction=0.1)$nodes)) {
        cat("Fig 6c:",p$years,"yr scenario - expansion found on attempt",attempt,"\n")
        return(tree_m)
      }
    }
    cat("Fig 6c:",p$years,"yr scenario - no expansion after",max_attempts,"attempts\n")
    return(NULL)
  })
  sim_tree_list.ultra<-Filter(Negate(is.null),sim_tree_list.ultra)
  names(sim_tree_list.ultra)<-sapply(expansion_params,function(p) p$years)[seq_along(sim_tree_list.ultra)]

  starting_vaf=0.2
  temp=Map(tree=sim_tree_list.ultra,years=names(sim_tree_list.ultra),f=function(tree,years) {
    #Skip scenarios where no clonal expansion was detected, rather than halting the
    #whole script: simulated_expansion_mito_vafs() calls get_expanded_clade_nodes()
    #and then getTips(), which errors on a zero-length node.
    if(!length(get_expanded_clade_nodes(tree,min_clonal_fraction=0.1)$nodes)) {
      cat("Fig 6c: no clonal expansion >=10% in the",years,"yr simulation - panel skipped\n")
      return(invisible(NULL))
    }
    fn<-paste0(fig6_dir,"Fig6c.simulation_plot_",starting_vaf,"_",years,"years.pdf")
    pdf(fn,width=2,height = 2.5)
    ok<-tryCatch({simulated_expansion_mito_vafs(tree,starting_vaf=starting_vaf,mito_cn=1000,generation_time=20);TRUE},
                 error=function(e){cat("Fig 6c:",years,"yr failed -",conditionMessage(e),"\n");FALSE})
    dev.off()
    if(!ok) unlink(fn) else cat("Fig 6c: wrote",years,"yr panel\n")
  })
} else {
  cat("Skipping Fig 6c (rsimpop not installed)\n")
}

#-----------------------------------------------------------------------------------#
# Fig 6d | Best shared marker mutations for the malignant clone
#
# Two example individuals: one MPN (PD9478) and one CML (PD51632). For each, the
# two best shared mtDNA marker mutations for the malignant clone are shown as a
# heatmap beneath the phylogeny, so that the marked samples can be read against
# the true clonal structure.
#
# VAFs below 3% are set to zero (rendered white): per the figure caption, samples
# with the mutation below that level are considered negative. The threshold is
# higher than the 1% used for the normal blood data (Fig. 5) because of the
# higher detection threshold in these datasets.
#-----------------------------------------------------------------------------------#

marker_mut_examples<-list(
  list(dataset="NW", exp_ID="PD9478",  muts=c("MT_9804_G_A","MT_11653_A_G")), #MPN
  list(dataset="CML",exp_ID="PD51632", muts=c("MT_919_A_G","MT_2256_T_C"))    #CML
)

vaf_negative_threshold=0.03 #samples below this are treated as negative (see caption)

#White for negative, then a YlOrRd ramp for increasing heteroplasmy - the same
#scheme as the shared-mutation heatmaps in Fig. 5, so the panels are comparable.
col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
names(col_scheme)<-seq(0,1,0.01)

for(ex in marker_mut_examples) {
  mito_data_ex<-readRDS(paste0(root_dir,"/data/nonblood/mito_mutation_data_",ex$dataset,".RDS"))
  list_ex<-mito_data_ex[[ex$exp_ID]]
  if(is.null(list_ex)) {cat("Fig 6d: no data for",ex$exp_ID,"- skipping\n"); next}

  tree.ultra_ex<-list_ex$tree.ultra
  #Apply the shearwater call filter, as elsewhere, so that unsupported calls do not
  #appear as low-level heteroplasmy
  vaf.mtx<-list_ex$matrices$vaf*list_ex$matrices$SW

  missing_muts<-ex$muts[!ex$muts%in%rownames(vaf.mtx)]
  if(length(missing_muts)) {cat("Fig 6d:",ex$exp_ID,"missing",paste(missing_muts,collapse=", "),"- skipping\n"); next}

  hm<-matrix(0,nrow=length(ex$muts),ncol=length(tree.ultra_ex$tip.label),
             dimnames=list(ex$muts,tree.ultra_ex$tip.label))
  for(i in 1:length(ex$muts)) {
    mut_vafs<-vaf.mtx[ex$muts[i],tree.ultra_ex$tip.label]
    mut_vafs[mut_vafs<vaf_negative_threshold]<-0
    mut_vafs<-round(mut_vafs,digits=2)
    hm[i,]<-col_scheme[as.character(mut_vafs)]
  }

  pdf(file=paste0(fig6_dir,"Fig6d.",ex$dataset,"_",ex$exp_ID,"_marker_mutations.pdf"),width=7,height=3)
  #plot_tree() returns the tree annotated with the plot geometry (ymax, coords);
  #add_mito_mut_heatmap() needs that ymax to place the heatmap rows, so the RETURN
  #VALUE must be passed on. Passing the original tree silently draws nothing (or
  #errors), because the stored objects carry no ymax.
  tree.ultra_plotted<-plot_tree(tree=tree.ultra_ex,cex.label=0,plot_axis=F,vspace.reserve=1.5)
  add_mito_mut_heatmap(tree=tree.ultra_plotted,heatmap=hm,border="gray",heatmap_bar_height=0.05,cex.label=0.5)
  dev.off()
  cat("Fig 6d: wrote",ex$dataset,ex$exp_ID,"\n")
}

#-----------------------------------------------------------------------------------#
# Fig 6e | Posterior drift rates in CML, MPN and normal blood expanded clones
#
# Uses the FINAL posterior of each cohort's sequential ABC - in a sequential ABC
# that posterior is conditioned on every mutation in the run, so it is the single
# summary for that setting. The prior is drawn once and shown as the reference row.
#-----------------------------------------------------------------------------------#

abc_dirs<-c(normal=paste0(root_dir,"/data/Drift_ABC_clonal_expansions/Drift_ABC_normal_sequential/"),
            MPN   =paste0(root_dir,"/data/Drift_ABC_clonal_expansions/Drift_ABC_MPN_sequential/"),
            CML   =paste0(root_dir,"/data/Drift_ABC_clonal_expansions/Drift_ABC_CML_sequential/"))

#Only include settings whose sequential ABC has been run
abc_dirs<-abc_dirs[sapply(abc_dirs,function(d) dir.exists(paste0(d,"output")))]
cat("Fig 6e using:",paste(names(abc_dirs),collapse=", "),"\n")

drift_rate_comparison<-Map(setting=names(abc_dirs),abc_dir=abc_dirs,f=function(setting,abc_dir) {
  muts_file<-list.files(path=abc_dir,pattern="abc_muts",full.names=TRUE)[1]
  abc_muts<-readr::read_csv(muts_file,show_col_types=FALSE)
  final_mut<-abc_muts$mut[nrow(abc_muts)] #the last mutation of the chain carries the accumulated posterior
  post<-readRDS(paste0(abc_dir,"output/posterior_table_",final_mut,".Rds"))
  data.frame(generation_time=post$generation_time,setting=setting)
})%>%dplyr::bind_rows()

#The prior: log-uniform over 0.1-500 days, i.e. flat on the log10 axis used below
prior<-data.frame(generation_time=10^runif(1e4,min=-1,max=2.7),setting="prior")

normal_disease_comparison_ridges_plot<-drift_rate_comparison%>%
  dplyr::bind_rows(prior)%>%
  mutate(setting=factor(setting,levels=c("prior",names(abc_dirs))))%>%
  ggplot(aes(x=generation_time,y=setting,fill=setting))+
  ggridges::geom_density_ridges(linewidth=0.3)+
  scale_fill_brewer(palette = "RdPu")+
  scale_x_log10(breaks=c(0.1,1,10,100),labels=c(0.1,1,10,100))+
  theme_classic()+
  my_theme+
  theme(strip.text.y=element_text(angle=0),legend.position = "none")+
  labs(x="Generation time (days)",y="Density")

ggsave(filename = paste0(fig6_dir,"Fig6e.normal_disease_comparison_ridges_plot.pdf"),
       normal_disease_comparison_ridges_plot,width=2.5,height=2.2)

cat("\nFigure 6 panels written to",fig6_dir,"\n")

#-------------------------------------------------------------------------------
# Mitochondrial_drift_through_tree_ABC_SEQ_local.R
#
# Local reconstruction of the sequential ABC that infers the mtDNA drift rate
# (effective generation time) from heteroplasmy levels across a clonal
# expansion with a known phylogeny. Handles both the "normal" (adult
# haematopoiesis), "MPN" and "MPN_nocoding" cohorts - see cohort_config below.
#
# Background: the original pipeline ran on the farm as three chained LSF stages
# (simulate -> extract sumstats in 200 array batches -> combine + abc), driven by
# controlling_sequential_ABC.sh. Those stage scripts and all of their outputs
# have been lost. The 200-way split was purely a parallelism device, not part of
# the method, so locally the three stages collapse into this single script.
#
# The statistical model is unchanged from the originals
# (Mitochondrial_drift_through_tree_simulations_MPNs_SEQUENTIAL_ABC.R for MPN,
# Mitochondrial_drift_through_tree_simulations.R for normal):
#   - same log-uniform prior on generation time (10^-1 to 10^2.7 days)
#   - same fixed mtDNA copy number of 600
#   - same 5 summary statistics
#   - same abc(tol=0.05, method="neuralnet")
# Deviations from those scripts are marked "DEVIATION:" below, and are either
# bug fixes or performance changes that do not alter the inference.
#
# Two ABC types are supported, matching the original analysis:
#   sequential - each mutation's posterior becomes the next mutation's prior, so
#                the estimate accumulates evidence across mutations (-a sequential)
#   individual - every mutation is fitted independently from the initial prior,
#                giving one free-standing posterior per mutation (-a individual)
#
# Usage:
#   Rscript Mitochondrial_drift_through_tree_ABC_SEQ_local.R -t normal -i 200 -j 1  # smoke test
#   Rscript Mitochondrial_drift_through_tree_ABC_SEQ_local.R -t normal              # full sequential run
#   Rscript Mitochondrial_drift_through_tree_ABC_SEQ_local.R -t normal -a individual
#-------------------------------------------------------------------------------

library(abc)
library(dplyr)
library(ggplot2)
library(ape)
library(optparse)
library(phylosignal)

option_list = list(
  make_option(c("-t", "--type"), action="store", default='MPN', type='character', help="cohort: 'normal' or 'MPN'"),
  make_option(c("-a", "--abc_type"), action="store", default='sequential', type='character', help="'sequential' (each posterior becomes the next prior) or 'individual' (every mutation fitted independently from the initial prior)"),
  make_option(c("-m", "--muts_df"), action="store", default='', type='character', help="csv with columns 'exp_ID' and 'mut', in the order they should be fitted; defaults to the cohort's manifest"),
  make_option(c("-j", "--j_index"), action="store", default=0, type='numeric', help="mutation index to run; 0 (default) runs all in sequence"),
  make_option(c("-i", "--iter"), action="store", default=2e4, type='numeric', help="simulations per mutation (2e4 = as published)"),
  make_option(c("-c", "--cores"), action="store", default=0, type='numeric', help="cores for the simulation/sumstat steps; 0 = detectCores()-1"),
  make_option(c("-s", "--seed"), action="store", default=42, type='numeric', help="RNG seed"),
  make_option(c("-o", "--out_dir"), action="store", default='', type='character', help="ABC directory; defaults to the cohort's directory under <root_dir>/data/")
  )

opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))
print(opt)

#-------------------------------------------------------------------------------
# Cohort configuration
#
# The normal and MPN analyses differ in three ways that are easy to miss and
# that silently corrupt the inference if crossed over:
#
# (1) units_per_year - the normal trees are scaled to MOLECULAR time (SNVs; the
#     KX004 tree is 1243.9 SNVs root-to-tip, i.e. 71.1 yr at 17.5 SNVs/yr),
#     whereas the MPN ultratrees are already scaled to YEARS (68.7 root-to-tip).
#     Using the wrong one misstates every generation count by a factor of 17.5.
# (2) absence_threshold - 0.02 for normal, raised to 0.03 for MPN because of the
#     higher detection threshold in the MPN data. Feeds n_absent and n_het.
# (3) where the tree lives - inside mito_data as $tree.ultra for normal, in a
#     separate ultratrees file for MPN.
#-------------------------------------------------------------------------------

cohort_config<-list(
  normal=list(mito_data_file="data/mito_data.Rds",
              trees_file=NA,            #trees are stored inside mito_data
              tree_slot="tree.ultra",
              units_per_year=17.5,      #trees in SNVs; ~17.5 SNVs/yr (Mitchell et al, 2022)
              absence_threshold=0.02,
              muts_file="abc_muts_normal.csv",
              out_subdir="Drift_ABC_normal_SEQ",
              out_subdir_individual="Drift_ABC_normal_individual"),
  MPN=list(mito_data_file="data/nonblood/mito_mutation_data_NW.RDS", #"NW" is the blood_MPN cohort
           trees_file="data/nonblood/MPN_ultratrees.RDS",
           tree_slot=NA,
           units_per_year=1,            #ultratrees already scaled to years
           absence_threshold=0.03,
           muts_file="abc_muts_MPN.csv",
           out_subdir="Drift_ABC_MPN",
           out_subdir_individual="Drift_ABC_MPN_individual"),
  #Sensitivity analysis: the same MPN data restricted to mutations with no
  #coding change, to check the drift estimate is not driven by selection on
  #protein-altering variants. Manifest recovered from
  #full_analysis_scripts/Nonblood_selection_analysis.R:251-252.
  MPN_nocoding=list(mito_data_file="data/nonblood/mito_mutation_data_NW.RDS",
                    trees_file="data/nonblood/MPN_ultratrees.RDS",
                    tree_slot=NA,
                    units_per_year=1,
                    absence_threshold=0.03,
                    muts_file="abc_muts_MPN_nocoding.csv",
                    out_subdir="Drift_ABC_MPN_nocoding",
                    out_subdir_individual="Drift_ABC_MPN_nocoding_individual")
)

if(!opt$t%in%names(cohort_config)) {stop("--type must be one of: ",paste(names(cohort_config),collapse=", "))}
config<-cohort_config[[opt$t]]
if(!opt$a%in%c("sequential","individual")) {stop("--abc_type must be 'sequential' or 'individual'")}
cat("Cohort:",opt$t,"| ABC type:",opt$a,"| units_per_year:",config$units_per_year,"| absence threshold:",config$absence_threshold,"\n")

#-------------------------------------------------------------------------------
# Paths - rooted at root_dir, following the convention in
# full_analysis_scripts/Nonblood_mtDNA_drift_analysis_local.R
#-------------------------------------------------------------------------------

root_dir="~/R_work/mito_mutations"
root_dir=path.expand(root_dir)

source(paste0(root_dir,"/data/mito_mutations_blood_functions.R")) #supplies find_latest_acquisition_node, getTips, add_ancestral_outgroup, plot_tree

out_subdir=if(opt$a=="individual") config$out_subdir_individual else config$out_subdir
abc_output_dir=if(nchar(opt$out_dir)) path.expand(opt$out_dir) else paste0(root_dir,"/data/",out_subdir,"/")
if(!grepl("/$",abc_output_dir)) {abc_output_dir<-paste0(abc_output_dir,"/")}

# Layout expected by Plot_drift_seq_ABC_results.R: an "abc_muts" csv in the ABC
# directory, posteriors under output/, figures under plots/.
dir.create(abc_output_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(paste0(abc_output_dir,"output"),showWarnings=FALSE)
dir.create(paste0(abc_output_dir,"plots"),showWarnings=FALSE)

# The simulations and the neuralnet regression are both stochastic and the
# published run was unseeded, so results will be close to but not identical to
# the original. L'Ecuyer-CMRG is required for reproducible parallel streams.
RNGkind("L'Ecuyer-CMRG")
set.seed(opt$seed)

n_cores=if(opt$cores>0) opt$cores else max(1,parallel::detectCores()-1)
cat("Using",n_cores,"cores\n")

#-------------------------------------------------------------------------------
# Data
#
# DEVIATION: the original re-derived germline calls, implied_mutCN column names
# and vaf.filt on every run (its lines 94-153). Those products are already
# stored in mito_mutation_data_NW.RDS, and - importantly - the ABC itself reads
# the RAW $matrices$vaf, never vaf.filt, so that whole block was vestigial here
# and is dropped. Nothing downstream changes.
#-------------------------------------------------------------------------------

mito_data<-readRDS(paste0(root_dir,"/",config$mito_data_file))
ultra_trees<-if(is.na(config$trees_file)) NULL else readRDS(paste0(root_dir,"/",config$trees_file))

muts_arg<-if(nchar(opt$m)) opt$m else config$muts_file
muts_file<-if(file.exists(muts_arg)) muts_arg else paste0(root_dir,"/simulation_scripts_for ABCs/",muts_arg)
abc_muts<-readr::read_csv(muts_file,show_col_types=FALSE)

# Keep a copy of the manifest alongside the posteriors, as the plotting script
# discovers the mutation order from it.
readr::write_csv(abc_muts,paste0(abc_output_dir,basename(muts_file)))

iter<-opt$iter
mtDNA_CN_used_in_simulations=600 #Copy number is fixed rather than inferred: it is only jointly identifiable with generation time, so fixing it lets the ABC hone in on the drift rate

#-------------------------------------------------------------------------------
# The drift simulator (unchanged from the original)
#
# Models independent Wright-Fisher drift down each lineage of a phylogeny:
# each "generation" is one round of binomial resampling of mito_cn genomes.
# Note the tree here is ultrametric and scaled to YEARS, so it is called with
# units_per_year=1 rather than the 17.5 SNVs/year molecular-clock default.
#-------------------------------------------------------------------------------

get_mito_mut_vaf_df=function(sub_tree,node,starting_vaf=0.5,units_per_year=17.5,mito_cn=1000,generation_time=1,vaf_df=NULL) {
  if(!is.numeric(node)|length(node)!=1) {stop("node must be a numeric of length=1")}
  if(is.null(vaf_df)){vaf_df<-data.frame(node=node,vaf=starting_vaf)}
  daughter_nodes<-sub_tree$edge[,2][which(sub_tree$edge[,1]==node)]
  for(daughter_node in daughter_nodes){

    #Get the branch length, and use to calculate the number of cell divisions the cell goes through
    branch_length<-sub_tree$edge.length[sub_tree$edge[,2]==daughter_node]
    curr_vaf<-vaf_df$vaf[vaf_df$node==node]
    ngen<-round((365*branch_length/units_per_year)/generation_time)
    if(ngen==0){ngen<-1}
    #This step performs the drift: several 'generations' of repeated binomial sampling with replacement
    for(i in 1:ngen) {curr_vaf<-rbinom(n=1,size=mito_cn,prob=curr_vaf)/mito_cn}

    #Add this data to the vaf_df
    vaf_df<-rbind(vaf_df,data.frame(node=daughter_node,vaf=curr_vaf))

    #Now perform function using the daughter as parent (iterative component of function)
    vaf_df<-get_mito_mut_vaf_df(sub_tree,node=daughter_node,starting_vaf=starting_vaf,mito_cn=mito_cn,generation_time=generation_time,vaf_df=vaf_df)
  }
  return(vaf_df)
}

#' Summary statistics of a set of tip VAFs
#'
#' The same five statistics are computed for the observed data and for every
#' simulation, so they must stay in one function. The absence threshold is
#' cohort-specific (see cohort_config).
#' Cmean is Abouheif's phylogenetic autocorrelation - the statistic that carries
#' the information about how drift is structured along the tree, as opposed to
#' the VAF distribution alone.
get_sumstats=function(sample_vafs,sub_tree.no.ancestral,absence_threshold=config$absence_threshold) {
  if(length(unique(sample_vafs))==1){
    Cmean=1 #A degenerate (all-identical) simulation is maximally autocorrelated by definition
  } else {
    sub_tree.4d<-phylobase::phylo4d(sub_tree.no.ancestral,tip.data=sample_vafs)
    #DEVIATION: methods="Cmean" and reps=1. The original called phyloSignal() with
    #its defaults, computing I, Lambda, K and K.star with 999 permutations each -
    #none of which are used, only $stat$Cmean. Cmean itself is deterministic given
    #tree + trait, so this is identical output for a fraction of the cost. This was
    #the step that forced the 200-job array on the farm.
    res<-phylosignal::phyloSignal(sub_tree.4d,methods="Cmean",reps=1)
    Cmean=res$stat$Cmean
  }
  return(data.frame(median_vaf=median(sample_vafs),
                    n_absent=sum(sample_vafs<absence_threshold),
                    n_homo=sum(sample_vafs>0.98),
                    n_het=sum(sample_vafs>absence_threshold&sample_vafs<0.98),
                    Cmean=Cmean))
}

#-------------------------------------------------------------------------------
# Run the sequential ABC
#
# DEVIATION: the original set j<-opt$j and then immediately overwrote it with
# `for(j in 1:6)`, so the -j argument the driver script passed was ignored. Here
# -j selects a single mutation and the default (0) runs the whole chain in order,
# which is what "sequential" requires anyway.
#-------------------------------------------------------------------------------

j_to_run<-if(opt$j==0) 1:nrow(abc_muts) else opt$j

for(j in j_to_run) {
  exp_ID=abc_muts$exp_ID[j]
  mut=abc_muts$mut[j]
  mut_RDS_file<-paste0(abc_output_dir,"abc.file.",exp_ID,".",mut,".Rds")
  mut_sim_file<-paste0(abc_output_dir,"sim.file.",exp_ID,".",mut,".Rds")
  mut_params_and_sumstats_file<-paste0(abc_output_dir,"params_and_sumstats.",exp_ID,".",mut,".Rds")
  posterior_table_file<-paste0(abc_output_dir,"output/posterior_table_",mut,".Rds")

  cat("\n================================================\n")
  cat(paste0("[",j,"/",nrow(abc_muts),"] ",exp_ID," ",mut,"\n"))
  cat("================================================\n")

  if(file.exists(posterior_table_file)) {
    cat("Posterior already exists - skipping (delete it to force a rerun)\n")
    next
  }

  #Set up the correct tree/ vaf matrix/ mut from the right individual.
  #For "normal" the ultrametric tree is stored inside mito_data; for MPN it
  #comes from the separate ultratrees file.
  tree.ultra<-if(is.null(ultra_trees)) mito_data[[exp_ID]][[config$tree_slot]] else ultra_trees[[exp_ID]]
  if(is.null(tree.ultra)) {stop("No tree found for ",exp_ID)}
  vaf<-mito_data[[exp_ID]]$matrices$vaf

  tree.ultra<-keep.tip(tree.ultra,tip = tree.ultra$tip.label[tree.ultra$tip.label%in%colnames(vaf)])

  if(!"Ancestral"%in%tree.ultra$tip.label) {tree.ultra<-add_ancestral_outgroup(tree=tree.ultra,outgroup_name="Ancestral")}
  tree.ultra$node.label<-tree.ultra$edge[,2][!tree.ultra$edge[,2]%in%1:length(tree.ultra$tip.label)]

  #Create the subtree: restrict to the clade descended from the latest node at
  #which the mutation could have been acquired, as drift is only modelled from
  #the point the heteroplasmy was present.
  pos_samples<-names(vaf)[which(vaf[mut,]>0.05)]
  #DEVIATION: restrict to samples that are actually tips of the (already pruned)
  #tree. A few mutations - e.g. PD9478 MT_11653_A_G and PD6629 MT_14168_T_C in
  #the MPN_nocoding set - are called in a sample that has no tip in the
  #phylogeny, and find_latest_acquisition_node then returns a length-zero node
  #and the run dies. Such a sample cannot be placed on the tree, is not
  #simulated, and does not enter the observed summary statistics (which are
  #taken over sub_tree.no.ancestral$tip.label), so dropping it here is what the
  #rest of the script already assumes. No effect where every positive sample is
  #a tip, which is the case for all of the normal and main MPN mutations.
  n_pos_all<-length(pos_samples)
  pos_samples<-pos_samples[pos_samples%in%tree.ultra$tip.label]
  if(length(pos_samples)<n_pos_all) {cat("NOTE:",n_pos_all-length(pos_samples),"of",n_pos_all,"positive samples are not tips of the tree and were dropped\n")}
  if(!length(pos_samples)) {stop("No positive samples remain on the tree for ",mut)}
  latest_acquisition_node=find_latest_acquisition_node(tree.ultra,pos_samples)
  sub_tree=drop.tip(tree.ultra,tree.ultra$tip.label[!tree.ultra$tip.label%in%c("Ancestral",getTips(tree.ultra,latest_acquisition_node))],trim.internal = T)
  sub_tree$coords<-NULL

  #Find the 'root' of the mutation tree
  starting_node<-sub_tree$edge[which(sub_tree$edge[,1]==sub_tree$edge[1,1]&!sub_tree$edge[,2]%in%1:length(sub_tree$tip.label)),2]
  sub_tree.no.ancestral<-drop.tip(sub_tree,"Ancestral")

  cat("Subtree:",length(sub_tree.no.ancestral$tip.label),"tips |",length(pos_samples),"samples with VAF>0.05\n")

  #Diagnostic plot of the tree, the observed VAFs and the phylogenetic correlogram.
  #Wrapped in try() so that a plotting failure cannot lose a long simulation run.
  try({
    pdf(paste0(abc_output_dir,"plots/",exp_ID,"_",mut,".pdf"),width=8,height=4)
    sub_tree_plotted=plot_tree(sub_tree,cex.label=0,bars = vaf[mut,])
    text(x = 0, y=-0.05*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",paste0("Max VAF: ",round(max(vaf[mut,]),digits = 3)),pos = 4)
    text(x = 0, y=0.95*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",mut,pos = 4)
    sub_tree.4d<-phylobase::phylo4d(sub_tree.no.ancestral,tip.data=t(vaf[mut,sub_tree.no.ancestral$tip.label]))
    barplot.phylo4d(sub_tree.4d)
    mut.crlg <- phyloCorrelogram(sub_tree.4d, trait = mut)
    plot(mut.crlg)
    dev.off()
  },silent=FALSE)

  #Observed summary statistics - computed with the same function as the simulations
  data_sample_vafs<-as.numeric(vaf[mut,sub_tree.no.ancestral$tip.label])
  sumstats.data<-get_sumstats(data_sample_vafs,sub_tree.no.ancestral)
  cat("Observed sumstats: ");print(sumstats.data)

  #-----------------------------------------------------------------------------
  # Prior
  #-----------------------------------------------------------------------------
  if(opt$a=="individual"|j==1) {
    #Every mutation of an "individual" ABC is fitted independently, so it always
    #starts from the initial log-uniform prior; in a sequential ABC only the
    #first mutation does.
    prior_table<-data.frame(starting_vaf=runif(n=iter),
                            mito_CN=rep(mtDNA_CN_used_in_simulations,iter),
                            log_generation_time=runif(n=iter,min=-1,max=2.7))%>%
      dplyr::mutate(generation_time=10^log_generation_time)
  } else {
    #Otherwise the previous mutation's posterior becomes this mutation's prior.
    #DEVIATION: the original mutate()d the posterior table in place, which only
    #works if it happens to have exactly `iter` rows - it has tol*iter. Building a
    #fresh table of `iter` rows by resampling is the same operation without that
    #constraint. starting_vaf is redrawn because it is a per-mutation nuisance
    #parameter, not shared across mutations; only the drift rate carries over.
    prev_posterior_file<-paste0(abc_output_dir,"output/posterior_table_",abc_muts$mut[j-1],".Rds")
    if(!file.exists(prev_posterior_file)) {stop("Missing posterior for the previous mutation: ",prev_posterior_file)}
    posterior_table<-readRDS(prev_posterior_file)
    prior_table<-data.frame(starting_vaf=runif(n=iter),
                            mito_CN=rep(mtDNA_CN_used_in_simulations,iter),
                            log_generation_time=sample(posterior_table$log_generation_time,size=iter,replace=TRUE))%>%
      dplyr::mutate(generation_time=10^log_generation_time)
  }

  #-----------------------------------------------------------------------------
  # Simulate
  #-----------------------------------------------------------------------------
  if(file.exists(mut_sim_file)){
    cat("Loading existing simulations\n")
    sim_out_list<-readRDS(mut_sim_file)
  } else {
    cat("Running",iter,"simulations\n")
    t0<-Sys.time()
    sim_out_list<-parallel::mclapply(1:iter,function(i) {
      #Get the starting parameters from the prior table
      starting_vaf=prior_table$starting_vaf[i]
      mito_cn=prior_table$mito_CN[i] #DEVIATION (bug fix): the original read prior_table$mito_cn, but the column is mito_CN. R's $ is case-sensitive, so this was NULL and copy number silently dropped out of the params vector, which would then be rejected by the transf argument to abc().
      generation_time=prior_table$generation_time[i]

      vaf_df<-get_mito_mut_vaf_df(sub_tree,node=starting_node,starting_vaf=starting_vaf,units_per_year=config$units_per_year,mito_cn=mito_cn,generation_time=generation_time)
      sample_vafs=vaf_df$vaf[order(vaf_df$node)][1:length(sub_tree.no.ancestral$tip.label)]
      names(sample_vafs)<-sub_tree.no.ancestral$tip.label

      return(list(params=c(starting_vaf=starting_vaf,mito_cn=mito_cn,generation_time=generation_time),
                  sample_vafs=sample_vafs))
    },mc.cores=n_cores)
    cat("  done in",round(difftime(Sys.time(),t0,units="mins"),1),"mins\n")
    saveRDS(sim_out_list,file=mut_sim_file)
  }

  #-----------------------------------------------------------------------------
  # Summary statistics
  #-----------------------------------------------------------------------------
  if(file.exists(mut_params_and_sumstats_file)) {
    cat("Loading existing params/sumstats\n")
    params_and_sumstats<-readRDS(mut_params_and_sumstats_file)
  } else {
    cat("Extracting summary statistics\n")
    t0<-Sys.time()
    params=dplyr::bind_rows(lapply(sim_out_list,function(list) return(list$params)))
    sumstats=dplyr::bind_rows(parallel::mclapply(sim_out_list,function(list) {
      get_sumstats(list$sample_vafs,sub_tree.no.ancestral)
    },mc.cores=n_cores))
    cat("  done in",round(difftime(Sys.time(),t0,units="mins"),1),"mins\n")
    params_and_sumstats=list(params=params,sumstats=sumstats)
    saveRDS(params_and_sumstats,file=mut_params_and_sumstats_file)
  }

  #-----------------------------------------------------------------------------
  # ABC
  #
  # mito_CN is fixed at 600, so it has zero variance and cannot be regressed on;
  # it is dropped here and the matching entry removed from transf. The original
  # passed it, but only because the case bug above had already stripped it from
  # params - so the effective model is unchanged.
  #-----------------------------------------------------------------------------
  cat("Running the ABC\n")
  abc_params<-params_and_sumstats$params[,c("starting_vaf","generation_time")]
  abc.nn<-abc(target = as.numeric(sumstats.data),
              param = abc_params,
              sumstat = params_and_sumstats$sumstats,
              tol=0.05,
              transf = c("log","log"),
              method = "neuralnet")
  saveRDS(abc.nn,file=mut_RDS_file)

  #-----------------------------------------------------------------------------
  # Posterior table
  #
  # This is the file the original _combine_sumstats.R wrote and that
  # Plot_drift_seq_ABC_results.R reads. Its format is pinned by how the original
  # script read it back in as a prior: it needs log_generation_time, plus
  # starting_vaf and generation_time for plotting.
  #-----------------------------------------------------------------------------
  posterior_table<-data.frame(starting_vaf=abc.nn$adj.values[,"starting_vaf"],
                              generation_time=abc.nn$adj.values[,"generation_time"])%>%
    dplyr::mutate(log_generation_time=log10(generation_time))
  saveRDS(posterior_table,file=posterior_table_file)

  cat("Posterior generation time (days): median",round(median(posterior_table$generation_time),3),
      " 95% CI",paste(round(quantile(posterior_table$generation_time,c(0.025,0.975)),3),collapse=" - "),"\n")
  cat("Wrote",posterior_table_file,"\n")
}

cat("\nDone.\n")

#Mitochondrial simulation
#This script is quite slow, therefore it is run as a program for each mutation at a time
#To run script do e.g.
#./Mitochondrial_drift_simulation.R -j 1

#Iterative modelling function to calculate the VAFs of the samples based on the:
#(1) starting VAF
#(2) mitochondrial copy number
#(3) frequency of cell divisions relative to somatic SNVs

library(dplyr)
library(ape)
library(optparse)
library(phylosignal)

option_list = list(
  make_option(c("-j", "--j_index"), action="store", default='1', type='numeric', help="mutation index within the dataframe"),
  make_option(c("-b", "--batch"), action="store", default='1', type='numeric', help="batch index")
  )
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))
print(opt)

###################################

#Source functions needed for the script
my_working_directory<-"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood"
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#############
root_dir="/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood"
############


#========================================#
# Import data ####
#========================================#

mito_cn=read.csv("/lustre/scratch125/casm/team268im/al28/mtDNA/whole_genome_coverage_pileup_and_bedtools_annotated.csv",header=T)

#Import metadata relating to all individuals studied
ref_df<-readxl::read_excel("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/non_blood_metadata.xlsx")

exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C") #These are the ones to exclude for the endometrial analysis

MPN_ultra_trees<-readRDS("MPN_ultratrees.RDS")

all_cohorts=c("NW")
all_mito_datasets<-lapply(all_cohorts,function(dataset) {
  if(dataset=="blood") {
    dataset_mito_data<-readRDS("../mito_data.Rds")
    CN_correlating_muts<-readRDS("../CN_correlation.RDS")
  } else {
    CN_correlating_muts=readRDS(paste0("CN_correlation_",dataset,".RDS"))
    dataset_mito_data=readRDS(paste0("mito_mutation_data_",dataset,".RDS"))
  }
  
  rho_cut_off=0
  
  #Deal with the germline mutations (if not done already)
  dataset_mito_data<-Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) {
    if(is.null(list$germline_muts)) {
      list$germline_muts<-identify_germline(matrices=list$matrices,threshold=0.9)
      cat(list$germline_muts,sep="\n")
      list$matrices<-reverse_germline(matrices=list$matrices,threshold=0.9)
    }
    return(list)
  })
  
  #Sort out the names of the mutCN column
  dataset_mito_data<-lapply(dataset_mito_data,function(list) {
    colnames(list$matrices$implied_mutCN)<-stringr::str_split(colnames(list$matrices$implied_mutCN),pattern = '\\.\\.\\.',simplify=T)[,1]
    return(list)
  })
  
  #Add the CN correlating muts vector to each patient (slightly inefficient but programatically more straight-forward)
  dataset_mito_data<-lapply(dataset_mito_data,function(list) {
    list$CN_correlating_muts<-CN_correlating_muts
    return(list)
  })
  
  
  cat("Adding the filtered VAF matrix.",sep="\n")
  #Annotate specific mutations with their most likely signature using the 'sig_ref' dataframe
  dataset_mito_data<-Map(list=dataset_mito_data,this_exp_ID=names(dataset_mito_data),function(list,this_exp_ID) {
    cat(this_exp_ID,sep="\n")
    
    if(dataset=="lymph") {list$tree$tip.label<-unique(list$sample_shearwater_calls$sampleID)}
    
    mutCN_cutoff=25 #If the mitochondrial copy number is over 25, retain mutation even if is in the "CN correlating muts" list
    CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
    vaf.filt<-(list$matrices$vaf[,list$tree$tip.label]*list$matrices$SW[,list$tree$tip.label][,list$tree$tip.label]*CN_correlating_mut_removal_mat[,list$tree$tip.label])[!grepl("DEL|INS",rownames(list$matrices$vaf))&
                                                                                                                                                                            list$rho_vals>rho_cut_off&
                                                                                                                                                                            !rownames(list$matrices$vaf)%in%exclude_muts,]
    
    list$matrices$CN_correlating_mut_removal_mat<-CN_correlating_mut_removal_mat
    list$matrices$vaf.filt<-vaf.filt
    return(list)
  })
  
  if(dataset=="lymph") {
    dataset_mito_data<-Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) {
      if(!exp_ID=="TX001") {list$het_oocyte_muts<-c()}
      return(list)
    })
  }
  return(dataset_mito_data)
})
names(all_mito_datasets)<-all_cohorts
mito_data<-all_mito_datasets$NW

#========================================#
# Select mutations for ABC ####
#========================================#
#One mutation is selected based on the "-j" option provided when executing the script
abc_muts<-data.frame(exp_ID=c("PD9478","PD5179","PD5847","PD5182","PD6646","PD4781"),
                     mut=c("MT_9804_G_A","MT_14111_T_C","MT_5136_T_C","MT_14452_A_G","MT_2270_A_G","MT_12565_T_C"))

abc_output_dir=paste0(root_dir,"/Drift_ABC_time_tree_MPN/")

j<-opt$j
batch_number<-opt$b
batch_size=100
exp_ID=abc_muts$exp_ID[j]
mut=abc_muts$mut[j]
mut_sim_file<-paste0(abc_output_dir,"sim.file.",exp_ID,".",mut,".Rds")
mut_params_and_sumstats_file<-paste0(abc_output_dir,"params_and_sumstats.","batch",batch_number,"_",exp_ID,".",mut,".Rds")

cat(paste0(exp_ID,'\n'))
cat(paste0(mut,'\n'))
cat(paste("Batch",batch_number),sep="\n")

#========================================#
# Get the sub_tree used for the simulations from the data ####
#========================================#
tree.ultra<-MPN_ultra_trees[[exp_ID]]
vaf<-mito_data[[exp_ID]]$matrices$vaf

tree.ultra<-keep.tip(tree.ultra,tip = tree.ultra$tip.label[tree.ultra$tip.label%in%colnames(vaf)])

if(!"Ancestral"%in%tree.ultra$tip.label) {tree.ultra<-add_ancestral_outgroup(tree=tree.ultra,outgroup_name="Ancestral")}
tree.ultra$node.label<-tree.ultra$edge[,2][!tree.ultra$edge[,2]%in%1:length(tree.ultra$tip.label)]

#Create the subtree
pos_samples<-names(vaf)[which(vaf[mut,]>0.05)]
latest_acquisition_node=find_latest_acquisition_node(tree.ultra,pos_samples)
sub_tree=drop.tip(tree.ultra,tree.ultra$tip.label[!tree.ultra$tip.label%in%c("Ancestral",getTips(tree.ultra,latest_acquisition_node))],trim.internal = T)
sub_tree$coords<-NULL

#Find the 'root' of the mutation tree
starting_node<-sub_tree$edge[which(sub_tree$edge[,1]==sub_tree$edge[1,1]&!sub_tree$edge[,2]%in%1:length(sub_tree$tip.label)),2]
sub_tree.no.ancestral<-drop.tip(sub_tree,"Ancestral")


#========================================#
# Import the simulation data & extract summary statistics ####
#========================================#

if(file.exists(mut_sim_file)){
  sim_out_list<-readRDS(mut_sim_file)
  iter=length(sim_out_list)
} else {
  cat("No simulation file found in expected location\n")
}
params=dplyr::bind_rows(lapply(sim_out_list,function(list) return(list$params)))

cat("Extracting the summary statistics from the simulations\n")
starting_idx=1+(batch_size*(batch_number-1))
ending_idx=batch_number*batch_size

cat(paste("Starting index:",starting_idx),sep="\n")
cat(paste("Ending index:",ending_idx),sep="\n")

if(starting_idx>iter) {stop(cat('Batch number is too high - no further simulations'))}
if(ending_idx>iter) {ending_idx<-iter}

sumstats=dplyr::bind_rows(Map(list=sim_out_list[starting_idx:ending_idx],i=starting_idx:ending_idx,function(list,i) {
  if(i%%5==0){cat(i,sep="\n")}
  
  sample_vafs=list$sample_vafs
  if(length(unique(sample_vafs))==1){
    Cmean=1
    sigma=NA
  } else {
    sub_tree.4d<-phylobase::phylo4d(sub_tree.no.ancestral,tip.data=sample_vafs)
    res<-phyloSignal(sub_tree.4d)
    mut.crlg <- phyloCorrelogram(sub_tree.4d,trait="dt")
    Cmean=res$stat$Cmean
    sigma=mut.crlg$sigma
  }
  
  return(data.frame(median_vaf=median(sample_vafs),
                    n_absent=sum(sample_vafs<0.03),
                    n_homo=sum(sample_vafs>0.98),
                    n_het=sum(sample_vafs>0.03&sample_vafs<0.98),
                    Cmean=Cmean,
                    sigma=sigma))
}))
params_and_sumstats=list(params=params,sumstats=sumstats)
saveRDS(params_and_sumstats,file=mut_params_and_sumstats_file)



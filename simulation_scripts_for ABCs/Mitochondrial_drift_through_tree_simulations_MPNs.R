#Mitochondrial simulation
#This script is quite slow, therefore it is run as a program for each mutation at a time
#To run script do e.g.
#./Mitochondrial_drift_simulation.R -j 1

#Iterative modelling function to calculate the VAFs of the samples based on the:
#(1) starting VAF
#(2) mitochondrial copy number
#(3) frequency of cell divisions relative to somatic SNVs

library(abc)
library(dplyr)
library(ggplot2)
library(ape)
library(optparse)
library(phylosignal)

option_list = list(
  make_option(c("-j", "--j_index"), action="store", default='1', type='numeric', help="mutation index within the dataframe")
  )
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))
print(opt)


# R_scripts_dir = "~/R_work/my_functions/"
# treemut_dir="~/R_work/treemut/"
# R_scripts=list.files(R_scripts_dir,pattern = ".R",full.names = T)
# sapply(R_scripts[-2],source)

###################################

#Source functions needed for the script
my_working_directory<-"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood"
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)
genomeFile="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"
root_dir="/lustre/scratch125/casm/team268im/al28/mtDNA"

#Import the mitochondrial copy number data
mito_cn=read.csv("/lustre/scratch125/casm/team268im/al28/mtDNA/whole_genome_coverage_pileup_and_bedtools_annotated.csv",header=T)

translate=data.frame(dataset=c("KY","HL","SO","PR","LM","NW","lymph"),
                     al_ref=c("lung_organoid","colon","colon_ibd","muty_mutant","endometrium","blood_MPN","immune"))

muty_samples=c("PD44887","PD44888","PD44889","PD44890","PD44891")

#Import metadata relating to all indviduals studied
ref_df<-readxl::read_excel("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/non_blood_metadata.xlsx")

#Import colony-level data for the lymphoid dataset as there are multiple different cell types
colony_info<-read.delim("colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt")
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C") #These are the ones to exclude for the endometrial analysis

MPN_ultra_trees<-readRDS("MPN_ultratrees.RDS")

##This is the key iterative simulation function for this simulation script. It takes:
#(1) a tree structure ("sub_tree")
#(2) the node number of the tree root, from which to start the drift
#(2) a starting VAF of a theoretical mitochondrial mutation in the root of the tree
#(3) population size and generation time parameters that define the Wright-Fisher drift rate
#It then models independent drift down each lineage and outputs a data frame the of the mutation VAF at each node /tip
#Note that it assumes that the tree is scale to "molecular time" (i.e. SNVs) and that the mutation rate is constant at ~17.5 SNVs per year (Mitchell et al, 2022)
#However, if the tree is already scaled to actual time, can change the 'units_per_year' metric to 1 (if scaled to years) or 365 (if scaled to days)
get_mito_mut_vaf_df=function(sub_tree,node,starting_vaf=0.5,units_per_year=17.5,mito_cn=1000,generation_time=1,vaf_df=NULL) {
  if(!is.numeric(node)|length(node)!=1) {stop("node must be a numeric of length=1")}
  if(is.null(vaf_df)){vaf_df<-data.frame(node=node,vaf=starting_vaf)}
  daughter_nodes<-sub_tree$edge[,2][which(sub_tree$edge[,1]==node)]
  for(daughter_node in daughter_nodes){
    
    #Get the branch length, and use to calculate the number of cell divisions the cell goes through (combined with the 'ngen_per_mut' variable)
    branch_length<-sub_tree$edge.length[sub_tree$edge[,2]==daughter_node]
    curr_vaf<-vaf_df$vaf[vaf_df$node==node]
    ngen<-round((365*branch_length/units_per_year)/generation_time) #Mutation rate is ~17.5 mutations per year (see Mitchell et al)
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

##CHOOSE DATASET AND START ANALYSIS
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




#############
root_dir="/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood"
mito_data<-all_mito_datasets$NW
############


#One mutation is selected based on the "-n" option provided when executing the script
#####-----------Select mutations for ABC-----------------
abc_muts<-data.frame(exp_ID=c("PD9478","PD5179","PD5847","PD5182","PD6646","PD4781"),
                     mut=c("MT_9804_G_A","MT_14111_T_C","MT_5136_T_C","MT_14452_A_G","MT_2270_A_G","MT_12565_T_C"))


########
abc_output_dir=paste0(root_dir,"/Drift_ABC_time_tree_MPN/")
########

j<-opt$j

for(j in 1:6) {
  exp_ID=abc_muts$exp_ID[j]
  mut=abc_muts$mut[j]
  mut_RDS_file<-paste0(abc_output_dir,"abc.file.",exp_ID,".",mut,".Rds")
  mut_sim_file<-paste0(abc_output_dir,"sim.file.",exp_ID,".",mut,".Rds")
  mut_params_and_sumstats_file<-paste0(abc_output_dir,"params_and_sumstats.",exp_ID,".",mut,".Rds")
  
  cat(paste0(exp_ID,'\n'))
  cat(paste0(mut,'\n'))
  
  force_rerun=T
  if(file.exists(mut_RDS_file)&!force_rerun){
    next
  } else {
    cat("The RDS file for this mutation does not currently exist\n")
    
    #Set up the correct tree/ vaf matrix/ mut from the right individual
    tree.ultra<-MPN_ultra_trees[[exp_ID]]
    vaf<-mito_data[[exp_ID]]$matrices$vaf
    
    tree.ultra<-keep.tip(tree.ultra,tip = tree.ultra$tip.label[tree.ultra$tip.label%in%colnames(vaf)])
    
    if(!"Ancestral"%in%tree.ultra$tip.label) {tree.ultra<-add_ancestral_outgroup(tree=tree.ultra,outgroup_name="Ancestral")}
    tree.ultra$node.label<-tree.ultra$edge[,2][!tree.ultra$edge[,2]%in%1:length(tree.ultra$tip.label)]
    iter<-2e4
    
    #Create the subtree
    pos_samples<-names(vaf)[which(vaf[mut,]>0.05)]
    latest_acquisition_node=find_latest_acquisition_node(tree.ultra,pos_samples)
    sub_tree=drop.tip(tree.ultra,tree.ultra$tip.label[!tree.ultra$tip.label%in%c("Ancestral",getTips(tree.ultra,latest_acquisition_node))],trim.internal = T)
    sub_tree$coords<-NULL
    pdf(paste0(abc_output_dir,exp_ID,"_",mut,".pdf"),width=8,height=4)
    sub_tree=plot_tree(sub_tree,cex.label=0,bars = vaf[mut,])
    text(x = 0, y=-0.05*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",paste0("Max VAF: ",round(max(vaf[mut,]),digits = 3)),pos = 4)
    text(x = 0, y=0.95*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",mut,pos = 4)
    
    #Find the 'root' of the mutation tree
    starting_node<-sub_tree$edge[which(sub_tree$edge[,1]==sub_tree$edge[1,1]&!sub_tree$edge[,2]%in%1:length(sub_tree$tip.label)),2]
    sub_tree.no.ancestral<-drop.tip(sub_tree,"Ancestral")
    sub_tree.4d<-phylobase::phylo4d(sub_tree.no.ancestral,tip.data=t(vaf[mut,sub_tree.no.ancestral$tip.label]))
    barplot.phylo4d(sub_tree.4d)
    res<-phyloSignal(sub_tree.4d)
    mut.crlg <- phyloCorrelogram(sub_tree.4d, trait = mut)
    mut.crlg$sigma
    plot(mut.crlg)
    dev.off()
    
    data_sample_vafs<-as.numeric(vaf[mut,sub_tree.no.ancestral$tip.label])
    sumstats.data<-c(median_vaf=median(data_sample_vafs),
                     n_absent=sum(data_sample_vafs<0.03), #Set as 0.03 due to increased threshold in the MPN data
                     n_homo=sum(data_sample_vafs>0.98),
                     n_het=sum(data_sample_vafs>0.03&data_sample_vafs<0.98),
                     Cmean=res$stat$Cmean,
                     sigma=mut.crlg$sigma)
    
    if(file.exists(mut_sim_file)){
      sim_out_list<-readRDS(mut_sim_file)
    } else {
      
      cat("No existing simulation file, starting the simulations\n")
      sim_out_list<-lapply(1:iter,function(i) {
        if(i%%10==0){print(i)}
        
        #Set the starting parameters
        starting_vaf=runif(1)
        mito_cn=600 #run with fixed population size of 600 -> makes it easier to hone in on the correct drift rate
        log_generation_time=runif(1,-1,2.7) #Log of the generation time (in days), therefore prior goes between 0.1 days and 500 days
        generation_time=10^log_generation_time
        
        vaf_df<-get_mito_mut_vaf_df(sub_tree,node=starting_node,starting_vaf=starting_vaf,units_per_year=1,mito_cn=mito_cn,generation_time=generation_time)
        sample_vafs=vaf_df$vaf[order(vaf_df$node)][1:length(sub_tree.no.ancestral$tip.label)]
        vaf_mat<-matrix(sample_vafs,dimnames=list(sub_tree.no.ancestral$tip.label,"vafs"))
        names(sample_vafs)<-sub_tree.no.ancestral$tip.label
        
        return(list(params=c(starting_vaf=starting_vaf,mito_cn=mito_cn,generation_time=generation_time),
                    sample_vafs=sample_vafs))
      })
      cat("Saving the simulations list\n")
      saveRDS(sim_out_list,file=mut_sim_file)
    }
    
    if(file.exists(mut_params_and_sumstats_file)) {
      params_and_sumstats<-readRDS(mut_params_and_sumstats_file)
    } else {
      params=dplyr::bind_rows(lapply(sim_out_list,function(list) return(list$params)))
      
      cat("Extracting the summary statistics from the simulations\n")
      sumstats=dplyr::bind_rows(Map(list=sim_out_list,i=1:iter,function(list,i) {
        if(i%%100==0){cat(i,sep="\n")}
        
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
      
    }
    
    cat("Saving the abc object\n")
    abc.nn<-abc(target = sumstats.data[1:5],param = params_and_sumstats$params,sumstat = params_and_sumstats$sumstats[,1:5],tol=0.05,transf = c("log","none","log"),method = "neuralnet")
    saveRDS(abc.nn,file=mut_RDS_file)
  }
}

#Pulling together the abc files and multiply together the posteriors
library(abc)
library(dplyr)
library(ggplot2)
my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"),
                strip.text.x = element_text(size=6))

abc_output_files<-list.files(path = abc_output_dir,pattern="abc.file.",full.names = F)
abc_values<-lapply(abc_output_files[-3],function(file){
  abc.nn<-readRDS(paste0(abc_output_dir,file))
  exp_ID<-stringr::str_split(file,pattern="\\.",simplify = T)[,3]
  mut<-stringr::str_split(file,pattern="\\.",simplify = T)[,4]
  return(data.frame(exp_ID=exp_ID,mut=mut,adj.values=log10(abc.nn$adj.values[,3]),unadj.values=log10(abc.nn$unadj.values[,3])))
})

#Group posterior results into "bins" of width 0.1 (on a log scale)
abc_method="nn"
if(abc_method=="rej"){
  col="unadj.values"
} else if(abc_method=="nn") {
  col="adj.values"
}
bin_width=0.1
bins=seq(-1,2.7-bin_width,bin_width)
all_post<-lapply(abc_values,function(df){
  n_per_bin<-sapply(bins,function(ll) sum(df[[col]]>ll&df[[col]]<ll+bin_width))
  names(n_per_bin)<-bins
  return(n_per_bin)
})%>%
  dplyr::bind_rows()

#Now multiply these posterior distributions across all mutations & normalize
#This gives a final posterior distribution for the drift that it most consistent with ALL individual mutations
comb_post<-apply(all_post,2,prod) #Multiply
comb_post<-comb_post/sum(comb_post) #Normalize

figures_dir=paste0(abc_output_dir,"plots/")
individual_mutation_posteriors<-abc_values%>%
  dplyr::bind_rows()%>%
  mutate(generation_time=10^get(col))%>%
  ggplot(aes(x=generation_time))+
  geom_histogram()+
  facet_wrap(~mut,nrow = 2)+
  theme_bw()+
  scale_x_log10(limits=c(0.1,500))+
  labs(x="Generation time (days)",y="Posterior distribution")+
  my_theme

ggsave(filename = paste0(figures_dir,"MPN_Drift_ABC_individual_posteriors_",abc_method,".pdf"),individual_mutation_posteriors,width = 7,height=2.5)

combined_posterior<-data.frame(log_generation_time=as.numeric(names(comb_post)),
                               posterior_probability=comb_post)%>%
  mutate(generation_time=10^(log_generation_time+0.5*bin_width))%>%
  ggplot(aes(x=generation_time,
             y=posterior_probability))+
  geom_bar(col="black",fill="lightblue",stat="identity")+
  geom_point(size=0.3)+
  scale_x_log10()+
  geom_line(linewidth=0.25)+
  theme_bw()+
  labs(x="Generation time (days)",y="Posterior distribution")+
  my_theme

ggsave(filename = paste0(figures_dir,"MPN_Drift_ABC_combined_posteriors_",abc_method,".pdf"),combined_posterior,width = 2,height=1.5)

#This cumulative probability plot helps to estimate the 95% posterior interval
posterior_probability_cumsum<-data.frame(log_generation_time=as.numeric(names(comb_post)),
                                         posterior_probability=comb_post)%>%
  mutate(log_generation_time=log_generation_time+0.5*bin_width)%>%
  mutate(cumsum=cumsum(posterior_probability))%>%
  ggplot(aes(x=log_generation_time,y=cumsum))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = c(0.025,0.975),linetype=2)+
  scale_x_continuous(breaks=seq(-1,2.7,by=0.05))+
  labs(x="Log10 effective generation time (days)",y="Cumulative sum of posterior probability")+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90))

#Using this plot estimate posterior intervals as:
#0.15 - 0.55, max posterior density = 0.3

#What are these as drift parameters?
600*10^(c(0.05,0.5,1.2))
600*10^(c(0.05,0.3,0.55))

ggsave(filename = paste0(figures_dir,"MPN_Drift_ABC_posterior_probability_cumsum_",abc_method,".pdf"),posterior_probability_cumsum,width = 5,height=5)



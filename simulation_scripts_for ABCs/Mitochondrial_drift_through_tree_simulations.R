#!/software/R-4.1.0/bin/Rscript

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

option_list = list(
  make_option(c("-j", "--j_index"), action="store", default='1', type='numeric', help="mutation index within the dataframe"),
  )
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))
print(opt)

##This is the key iterative simulation function for this simulation script. It takes:
#(1) a tree structure ("sub_tree")
#(2) the node number of the tree root, from which to start the drift
#(2) a starting VAF of a theoretical mitochondrial mutation in the root of the tree
#(3) population size and generation time parameters that define the Wright-Fisher drift rate
#It then models independent drift down each lineage and outputs a data frame the of the mutation VAF at each node /tip
#Note that it assumes that the tree is scale to "molecular time" (i.e. SNVs) and that the mutation rate is constant at ~17.5 SNVs per year (Mitchell et al, 2022)
get_mito_mut_vaf_df=function(sub_tree,node,starting_vaf=0.5,mito_cn=1000,generation_time=1,vaf_df=NULL) {
  if(is.null(vaf_df)){vaf_df<-data.frame(node=node,vaf=starting_vaf)}
  daughter_nodes<-sub_tree$edge[,2][which(sub_tree$edge[,1]==node)]
  for(daughter_node in daughter_nodes){
    
    #Get the branch length, and use to calculate the number of cell divisions the cell goes through (combined with the 'ngen_per_mut' variable)
    branch_length<-sub_tree$edge.length[sub_tree$edge[,2]==daughter_node]
    curr_vaf<-vaf_df$vaf[vaf_df$node==node]
    ngen<-round((365*branch_length/17.5)/generation_time) #Mutation rate is ~17.5 mutations per year (see Mitchell et al)
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

R_scripts_dir = "~/R_work/my_functions/"
treemut_dir="~/R_work/treemut/"
R_scripts=list.files(R_scripts_dir,pattern = ".R",full.names = T)
sapply(R_scripts[-2],source)

root_dir="~/R_work/mito_mutations_blood/"
mito_data_file=paste0(root_dir,"data/mito_data.Rds")
mito_data<-readRDS(mito_data_file)

#One mutation is selected based on the "-n" option provided when executing the script
#####-----------Select mutations for ABC-----------------
abc_muts<-data.frame(exp_ID=c("KX008","KX003","KX004","KX004","KX004","KX004","KX004","KX007","KX003","KX004"),
                     mut=c("MT_13552_G_A","MT_9151_A_G","MT_4232_T_C","MT_4487_A_G","MT_4965_A_G","MT_6379_T_C","MT_11790_T_C","MT_6180_G_A","MT_12396_T_C","MT_14035_T_C"))
abc_output_dir=paste0(root_dir,"data/Drift_ABC2/")

j<-opt$j
exp_ID=abc_muts$exp_ID[j]
mut=abc_muts$mut[j]
mut_RDS_file<-paste0(abc_output_dir,"abc.file.",exp_ID,".",mut,".Rds")
mut_sim_file<-paste0(abc_output_dir,"sim.file.",exp_ID,".",mut,".Rds")
mut_params_and_sumstats_file<-paste0(abc_output_dir,"params_and_sumstats.",exp_ID,".",mut,".Rds")

cat(paste0(exp_ID,'\n'))
cat(paste0(mut,'\n'))

if(file.exists(mut_RDS_file)){
  next
} else {
  cat("The RDS file for this mutation does not currently exist\n")
  
  #Set up the correct tree/ vaf matrix/ mut from the right individual
  tree.ultra<-mito_data[[exp_ID]]$tree.ultra
  tree.ultra$node.label<-tree.ultra$edge[,2][!tree.ultra$edge[,2]%in%1:length(tree.ultra$tip.label)]
  vaf<-mito_data[[exp_ID]]$matrices$vaf
  iter<-2e4
  
  #Create the subtree
  pos_samples<-names(vaf)[which(vaf[mut,]>0.05)]
  latest_acquisition_node=find_latest_acquisition_node(tree.ultra,pos_samples)
  sub_tree=drop.tip(tree.ultra,tree.ultra$tip.label[!tree.ultra$tip.label%in%c("Ancestral",getTips(tree.ultra,latest_acquisition_node))],trim.internal = T)
  sub_tree$coords<-NULL
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
  
  data_sample_vafs<-as.numeric(vaf[mut,sub_tree.no.ancestral$tip.label])
  sumstats.data<-c(median_vaf=median(data_sample_vafs),
                   n_absent=sum(data_sample_vafs<0.02),
                   n_homo=sum(data_sample_vafs>0.98),
                   n_het=sum(data_sample_vafs>0.02&data_sample_vafs<0.98),
                   Cmean=res$stat$Cmean,
                   sigma=mut.crlg$sigma)
  
  if(file.exists(mut_sim_file)){
    sim_out_list<-readRDS(mut_sim_file)
  } else {
    
    cat("No existing simulation file, starting the simulations\n")
    sim_out_list<-lapply(1:iter,function(i) {
      if(i%%500==0){print(i)}
      
      #Set the starting parameters
      starting_vaf=runif(1)
      mito_cn=675 #run with fixed population size of 675 -> makes it easier to hone in on the correct drift rate
      log_generation_time=runif(1,-1,2.7)
      generation_time=10^log_generation_time
      
      vaf_df<-get_mito_mut_vaf_df(sub_tree,node=starting_node,starting_vaf=starting_vaf,mito_cn=mito_cn,generation_time=generation_time)
      sample_vafs=vaf_df$vaf[order(vaf_df$node)][1:length(sub_tree.no.ancestral$tip.label)]
      vaf_mat<-matrix(sample_vafs,dimnames=list(sub_tree.no.ancestral$tip.label,"vafs"))
      names(sample_vafs)<-sub_tree.no.ancestral$tip.label
      
      return(list(params=c(starting_vaf=starting_vaf,mito_cn=mito_cn,generation_time=generation_time),
                  sample_vafs=sample_vafs))
    })
    cat("Saving the simulations list\n")
    saveRDS(sim_out_list,file=mut_sim_file)
  }
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
                      n_absent=sum(sample_vafs<0.02),
                      n_homo=sum(sample_vafs>0.98),
                      n_het=sum(sample_vafs>0.02&sample_vafs<0.98),
                      Cmean=Cmean,
                      sigma=sigma))
  }))
  params_and_sumstats=list(params=params,sumstats=sumstats)
  saveRDS(params_and_sumstats,file=mut_params_and_sumstats_file)
  
  cat("Saving the abc object\n")
  abc.nn<-abc(target = sumstats.data,param = params,sumstat = sumstats,tol=0.05,transf = c("none","none","log"),method = "neuralnet")
  saveRDS(abc.nn,file=mut_RDS_file)
}


library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)

root_dir="~/R_work/mito_mutations_blood/"
figures_dir=paste0(root_dir,"figures/")
source(paste0(root_dir,"data/mito_mutations_blood_functions.R"))


###-----------------------------DRIFT SIMULATIONS OF HETEROPLASMIC OOCYTE MUTATIONS-------------------------------
#Simulate drift of heteroplasmic oocyte mutation, approximate here as occurring at two different rates:
#(1) the higher drift rate apparent during early development (here defined as pre-conception)
#(2) the slower drift rate occurring post development (here defined as post-conception)
#Assess the VAF distribution across a population of cells at different ages of the individual
years_to_test=c(1,5,10,20,40,80)
sim_df<-lapply(years_to_test,function(years) {
  #Initial drift during development (here defined as pre-conception)
  #Pop size (450) x generation time (2.67) = 1200 (the development drift parameter)
  vafs=sapply(1:5000,function(j) fisher_wright_drift(0.036,population_size=450,generation_time=2.67,total_time=274))
  #Subsequent drift during post-conception life
  #Pop size (675) x generation time (23) = 15,500 (the adult drift parameter)
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
  geom_histogram(bins = 50,col="black",size=0.2)+
  theme_classic()+
  facet_grid(rows=vars(years))+
  scale_fill_brewer(palette="Set3")+
  scale_y_log10()+
  labs(x="VAF distribution",y="Count",fill="Time from\nconception\n(Years)")+
  my_theme+
  theme(strip.text.y = element_text(size=5))
ggsave(paste0(figures_dir,"Figure_06/VAF_distribution_histogram.pdf"),VAF.distribution.histogram,height=3,width=2.5)

#Now simulate the final distribution of mutation VAFs through the phylogeny using one of the adult phylogenies (KX004) as an example
#Here we just use a single drift rate (the adult rate) for simplicity
mito_data<-readRDS(paste0(root_dir,"data/mito_data.Rds"))
exp_ID="KX004"
vaf_df_het_oocyte=get_mito_mut_vaf_df(mito_data[[exp_ID]]$tree.ultra,node=mito_data[[exp_ID]]$tree.ultra$edge[1,1],starting_vaf = 0.036,mito_cn = 675,generation_time = 18)
sample_vafs=vaf_df_het_oocyte$vaf[order(vaf_df_het_oocyte$node)][1:length(mito_data[[exp_ID]]$tree.ultra$tip.label)]
names(sample_vafs)<-mito_data[[exp_ID]]$tree.ultra$tip.label

pdf(file = paste0(figures_dir,"Figure_06/example_phylogeny_hetmut_in_oocyte_simulation.pdf"),width = 7,height=2.5)
mito_data[[exp_ID]]$tree.ultra=plot_tree(mito_data[[exp_ID]]$tree.ultra,cex.label=0,bars=sample_vafs)
text(x = 0, y=-0.175*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",paste0("Max VAF: ",max(sample_vafs)),pos = 4)
dev.off()


###-----------------------------DRIFT SIMULATIONS THROUGH CLONAL EXPANSIONS WITH DIFFERENT DYNAMICS-------------------------------
#The clonal expansions in normal ageing blood tend to be long slow gradual expansions over several decades.
#Our data tells us how mitochondrial mutations behave as potential lineage markers in this context.
#However, other expansions (e.g. malignancies) may occur over much shorter time periods.
#This script is therefore to test if mitochondrial mutations may act as better
#lineage markers over shorter time periods.
#Note that implicit in these simulations is that the drift rate is fairly similar in these different
#contexts.  This may not be true as drift is likely to increase during more rapid cell division.
#Therefore these simulations represent a "best case scenario" for the behaviour of mitochondrial
#mutations as linage tracing markers during faster clonal expansions.

rsimpop2_dir="~/R_work/rsimpop2"
setwd(rsimpop2_dir)
dyn.load("src/rsimpop_interface.so")
source("R/sim_pop_v2.R")
source("R/plot_tree.R")
source("R/plot_tree_annots.R")
source("R/sim_pop_v2.R")
source("R/wrapped_sims.R")
setwd(root_dir)

n_sampled=30 #The number of simulated single cells picked from the clonal expansion

##This is the key iterative simulation function to simulate drift through the tree. It takes:
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
    ngen<-round((365*branch_length/17.5)/generation_time)
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

###-----------------------------Use a simulation of selection-------------------------------
#This is a wrapper function for rsimpop that simulations a phylogeny, with introduction of a driver mutation at a specific time point
run_selection_sim=function (initial_division_rate,
                            final_division_rate,
                            target_pop_size = 1e+05,
                            nyears_driver_acquisition = 15,
                            nyears = 40,
                            fitness = 0.2) {
    cfg = list(compartment = data.frame(val = c(0, 1), rate = c(-1,initial_division_rate), popsize = c(1, target_pop_size)),
               info = data.frame(val = c(0, 1, 1), fitness = c(0, 0,fitness), population = rep(0, 3)))
    params = list(n_sim_days = 50 * 365, b_stop_at_pop_size = 1, b_stop_if_empty = 0)
    growthphase = sim_pop(NULL, params = params, cfg)
    hscDivTime = 1/(2 * final_division_rate)
    tree0 = get_tree_from_simpop(growthphase)
    cfg$compartment$rate[2] = final_division_rate
    cfg$compartment$popsize[2] = target_pop_size
    years = nyears
    params[["n_sim_days"]] = nyears_driver_acquisition * 365
    mutPerYear = 15
    params[["b_stop_at_pop_size"]] = 0
    adult1 = sim_pop(tree0, params = params, cfg)
    adult1 = combine_simpops(growthphase, adult1)
    tree1 = get_tree_from_simpop(adult1)
    params[["n_sim_days"]] = nyears * 365
    params[["b_stop_if_empty"]] = 1
    dc = 0
    tries = 0
    tree1_tmp = tree1
    while (dc/target_pop_size < 0.05 && tries < 100) {
        cat("No driver found: tries=", tries, "\n")
        idx = sample(length(tree1$tip.label) - 1, 1) + 1
        celltype = rep(NA, length(tree1$tip.label))
        celltype[idx] = -1
        tree1_tmp = assign_celltype(tree1, celltype, tree1$cfg)
        adult2 = sim_pop(tree1_tmp, params = params, tree1_tmp$cfg)
        dc = adult2$cfg$info$population[3]
        tries = tries + 1
    }
    adult2 = combine_simpops(adult1, adult2)
    fulltree = get_tree_from_simpop(adult2)
    return(fulltree)
}

####-------------SIMULATE THE CLONAL EXPANSIONS---------------

##Early driver acquisition - 35 years from acquisition to sampling
#Introduce driver after 15 years, and simulation runs to 50 years i.e. 35 years of gradual expansion
#NB. the relative increased fitness of the driver clone is fairly slight (0.3) to reflect the gradual expansion dynamics
selsim=run_selection_sim(0.05,1/(2*190),nyears_driver_acquisition = 15,target_pop_size = 5e4,nyears = 50,fitness=0.3)
seltree=get_tree_from_simpop(selsim)
seltree100=get_subsampled_tree(seltree,100)
tree=plot_tree(seltree100,cex.label=0)
tree_m=get_elapsed_time_tree(seltree100,mutrateperdivision=0.65,backgroundrate=16/365)
tree_35=plot_tree(tree_m,cex.label=0)

##Early mid driver acquisition - 25 years from acquisition to sampling
#Introduce driver after 25 years, and simulation runs to 50 years i.e. 25 years of  expansion
#NB. the relative increased fitness of the driver clone is increased (0.45) to reflect the gradual expansion dynamics
selsim=run_selection_sim(0.05,1/(2*190),nyears_driver_acquisition = 25,target_pop_size = 5e4,nyears = 50,fitness=0.45)
seltree=get_tree_from_simpop(selsim)
seltree100=get_subsampled_tree(seltree,100)
tree=plot_tree(seltree100,cex.label=0)
tree_m=get_elapsed_time_tree(seltree100,mutrateperdivision=0.65,backgroundrate=16/365)
tree_25=plot_tree(tree_m,cex.label=0)

##Late mid driver acquisition - 15 years from acquisition to sampling
#Introduce driver after 35 years, and simulation runs to 50 years i.e. 15 years of  expansion
#NB. the relative increased fitness of the driver clone is further increased (0.6) 
selsim=run_selection_sim(0.05,1/(2*190),nyears_driver_acquisition = 35,target_pop_size = 5e4,nyears = 50,fitness=0.6)
seltree=get_tree_from_simpop(selsim)
seltree100=get_subsampled_tree(seltree,100)
tree=plot_tree(seltree100,cex.label=0)
tree_m=get_elapsed_time_tree(seltree100,mutrateperdivision=0.65,backgroundrate=16/365)
tree_15=plot_tree(tree_m,cex.label=0)

##Late mid driver acquisition - 10 years from acquisition to sampling
#Introduce driver after 35 years, and simulation runs to 50 years i.e. 15 years of  expansion
#NB. the relative increased fitness of the driver clone is further increased (0.6) 
selsim=run_selection_sim(0.05,1/(2*190),nyears_driver_acquisition = 40,target_pop_size = 5e4,nyears = 50,fitness=1.2)
seltree=get_tree_from_simpop(selsim)
seltree100=get_subsampled_tree(seltree,100)
tree=plot_tree(seltree100,cex.label=0)
tree_m=get_elapsed_time_tree(seltree100,mutrateperdivision=0.65,backgroundrate=16/365)
tree_10=plot_tree(tree_m,cex.label=0)

##Late driver acquisition - 5 years from acquisition to sampling
#Introduce driver after 45 years, and simulation runs to 50 years i.e. 5 years of  rapid expansion
#NB. the relative increased fitness of the driver clone is very high (2.5) to reflect the rapid expansion dynamics
selsim=run_selection_sim(0.05,1/(2*190),nyears_driver_acquisition = 45,target_pop_size = 5e4,nyears = 50,fitness=2.5)
seltree=get_tree_from_simpop(selsim)
seltree100=get_subsampled_tree(seltree,100)
tree=plot_tree(seltree100,cex.label=0)
tree_m=get_elapsed_time_tree(seltree100,mutrateperdivision=0.65,backgroundrate=16/365)
tree_5=plot_tree(tree_m,cex.label=0)

##Very late driver acquisition - 2 years from acquisition to sampling
#Introduce driver after 48 years, and simulation runs to 50 years i.e. 2 years of  very rapid expansion
#NB. the relative increased fitness of the driver clone is even higher (5) to reflect the rapid expansion dynamics
selsim=run_selection_sim(0.05,1/(2*190),nyears_driver_acquisition = 48,target_pop_size = 5e4,nyears = 50,fitness=5)
seltree=get_tree_from_simpop(selsim)
seltree100=get_subsampled_tree(seltree,100)
tree=plot_tree(seltree100,cex.label=0)
tree_m=get_elapsed_time_tree(seltree100,mutrateperdivision=0.65,backgroundrate=16/365)
tree_2=plot_tree(tree_m,cex.label=0)


####-------------SIMULATE THE DRIFT OF MITOCHONDRIAL MUTATIONS---------------

#Wrapper function to simulate mitochondrial drift through the clonal expansion
#This detects the clonal expansion and drops the other tips
#Then simulates the drift from the MRCA of the expansion
simulated_expansion_mito_vafs=function(tree,starting_vaf,mito_cn,generation_time,plot=T){
  expansion_node=get_expanded_clade_nodes(tree,min_clonal_fraction = 0.1)
  sub_tree<-drop.tip(tree,tree$tip.label[!tree$tip.label%in%c("s1",getTips(tree,expansion_node$nodes))])
  sub_tree$coords<-NULL
  starting_node<-sub_tree$edge[which(sub_tree$edge[,1]==sub_tree$edge[1,1]&!sub_tree$edge[,2]%in%1:length(sub_tree$tip.label)),2]
  vaf_df<-get_mito_mut_vaf_df(sub_tree = sub_tree,
                              node=starting_node,
                              starting_vaf=starting_vaf,
                              mito_cn=mito_cn,
                              generation_time=generation_time)
  vaf_vec=vaf_df$vaf[order(vaf_df$node)][1:length(sub_tree$tip.label)]
  names(vaf_vec)=sub_tree$tip.label
  if(plot){
    sub_tree=plot_tree(sub_tree,cex.label=F,bars = vaf_vec)
    text(x = 0, y=-0.05*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",paste0("Max VAF: ",round(max(vaf_vec),digits = 3)),pos = 4)
  }
  return(vaf_vec)
}

#Calculate what proportion of clade samples are marked according to the starting vaf and years since origin of expansion
starting_vafs=c(0.01,0.05,0.1,0.2,0.4,0.8,0.99) #starting_vafs=c(0.01,0.05)

#Make a list of the different trees, and then make ultrametric versions
sim_tree_list=list(years_35=tree_35,years_25=tree_25,years_15=tree_15,years_10=tree_10,years_5=tree_5,years_2=tree_2)
sim_tree_list.ultra=lapply(sim_tree_list,function(tree) {
  tree.ultra<-make.ultrametric.tree(tree)
  tree.ultra$coords<-NULL
  tree.ultra$edge.length[is.infinite(tree.ultra$edge.length)]<-0
  tree.ultra$edge.length<-tree.ultra$edge.length*mean(get_mut_burden(tree))
  return(tree.ultra)
})

#Now simulate drift from different starting VAFs in the MRCA for each tree
#Repeat the simulation 100 times for each starting VAF and tree
all_years=c(35,25,15,10,5,2)
vaf_detection_threshold=0.01 #vaf_detection_threshold=0

all_trees_df<-Map(tree=sim_tree_list,years=all_years,function(tree,years) {
  cat(paste0(years," years of expansion tree\n"))
  tree_df<-lapply(starting_vafs,function(starting_vaf) {
    cat(paste0(starting_vaf,"\n"))
    vaf_list<-lapply(1:100,function(i) simulated_expansion_mito_vafs(tree,starting_vaf=starting_vaf,mito_cn=1000,generation_time=20,plot=F))
    prop_marked=sapply(vaf_list,function(vaf_vec) sum(vaf_vec>vaf_detection_threshold)/length(vaf_vec))
    return(data.frame(years=years,starting_vaf=starting_vaf,prop_marked=prop_marked))
  })%>%dplyr::bind_rows()
  return(tree_df)
})%>%dplyr::bind_rows()


##-------------------VISUALIZE THE RESULTS-------------------

#Plot the expanded clades with simulated vafs, with starting VAF of 0.2 in MRCA (for illustrative purposes)
temp=Map(tree=sim_tree_list.ultra,years=all_years,f=function(tree,years) {
  starting_vaf=0.2
  pdf(paste0(figures_dir,"Figure_06/simulation_plot_",starting_vaf,"_",years,"years.pdf"),width=2,height = 2.5)
  simulated_expansion_mito_vafs(tree,starting_vaf=starting_vaf,mito_cn=1000,generation_time=20)
  dev.off()
})

#Plot the proportion of clonal expansion samples that would have the original marker mutation at >1% VAF in different contexts
sim_prop_marked<-all_trees_df%>%
  mutate(years=factor(paste(years,"years"),levels=paste(rev(all_years),"years")))%>%
  ggplot(aes(x=factor(starting_vaf),y=prop_marked))+
  geom_boxplot()+
  geom_jitter(size=0.5,alpha=0.2,width=0.3,height=0.01)+
  facet_grid(cols=vars(years))+
  theme_bw()+
  my_theme+
  labs(x="Mitochondrial mutation VAF in MRCA",y=paste0("Proportion of expansion with mutation >",100*vaf_detection_threshold,"% VAF"))
ggsave(filename = paste0(figures_dir,"Figure_06/sim_marked_prop.pdf"),plot = sim_prop_marked,height=2.5,width=7)




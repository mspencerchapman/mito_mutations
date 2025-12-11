library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(phylosignal)
options(stringsAsFactors = F)

#Set these file paths before running the script
root_dir="~/R_work/mito_mutations_blood"
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"/data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
figures_dir=paste0(root_dir,"/figures/")

#Set the plotting theme for ggplot2
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

#Now import the mitochondrial mutation data
mito_data_file=paste0(root_dir,"/data/mito_data.Rds")
mito_data<-readRDS(mito_data_file)
CN_correlating_muts<-readRDS(paste0(root_dir,"/data/CN_correlation.RDS"))
exclude_muts<-c("MT_16183_A_C",CN_correlating_muts) #The "MT_16183_A_C mut is a germline mutation in the 8pcw foetus, but not recognised as such due to a very common deletion result from the homopolymer that results

##----------------------LOOK FOR LOW LEVEL HETEROPLASMIC MUTATIONS-------------------

#Analysis only in the foetal/ cord blood phylogenies
young_IDs=c("8pcw","18pcw","CB001","CB002")

#Add a "stats_df" object that gives the phylogenetic signal/ beta-binomial values for each shared mutation
mito_data_young<-Map(list=mito_data[young_IDs],exp_ID=young_IDs,function(list,exp_ID) {
  print(exp_ID)
  
  #If phylogenetic signal data is already embedded, use this. Otherwise just add for shared mutations.
  if(is.null(list$phylosignal)) {
    list$stats_df<-data.frame(exp_ID=exp_ID,
                                   rho_val=list$rho_vals,
                                   n_pos=rowSums(list$matrices$SW),
                                   mean_vaf=rowMeans(list$matrices$vaf),
                                   max_vaf=apply(list$matrices$vaf,1,max))%>%
      tibble::rownames_to_column(var="mut_ref")
    
    print(nrow(list$stats_df))
    
    tree.noancestral<-drop.tip(list$tree.ultra,"Ancestral")
    tree.4d<-phylobase::phylo4d(tree.noancestral,tip.data=t(list$matrices$vaf[list$stats_df$mut_ref,tree.noancestral$tip.label]))
    phylosignal<-phyloSignal(tree.4d)
    
    list$stats_df$Cmean_pvalue<-phylosignal$pvalue$Cmean
  } else {
    list$stats_df<-data.frame(exp_ID=exp_ID,
                                   Cmean_pvalue=list$phylosignal$pvalue$Cmean,
                                   rho_val=list$rho_vals,
                                   n_pos=rowSums(list$matrices$SW),
                                   mean_vaf=rowMeans(list$matrices$vaf),
                                   max_vaf=apply(list$matrices$vaf,1,max))%>%
      tibble::rownames_to_column(var="mut_ref")%>%
      dplyr::filter(n_pos>1)
  }
  return(list)
})

#Now pull out the shared mutations for which the MRCA is the tree root and visualize
resave_plots=F
het_oocyte_muts<-lapply(young_IDs,function(Exp_ID) {
  low_level_het_muts<-mito_data_young[[Exp_ID]]$stats_df%>%
    dplyr::bind_rows()%>%
    dplyr::filter(n_pos>=2 & !grepl("INS|DEL",mut_ref))%>%
    mutate(Cmean_qvalue=p.adjust(Cmean_pvalue,method = "BH"))%>%
    dplyr::filter(!mut_ref%in%exclude_muts)%>%
    dplyr::filter(max_vaf>0.01)%>%
    dplyr::filter(!(mean_vaf>0.003&Cmean_pvalue>0.01&rho_val<=5e-3))%>% #Final filter to take out a few remaining artefacts that are present at fairly high global vaf but are not phylocorrelated and have low dispersion
    filter(exp_ID==Exp_ID)%>%
    pull(mut_ref)
  print(low_level_het_muts)
  
  #find MRCA of pos samples for each mutation
  define_pos=0.005
  tree.no.ancestral<-drop.tip(mito_data[[Exp_ID]]$tree.ultra,"Ancestral")
  mut_MRCA<-sapply(low_level_het_muts,function(mut) {
    pos_samples=colnames((mito_data[[Exp_ID]]$matrices$vaf*mito_data[[Exp_ID]]$matrices$SW)[mut,(mito_data[[Exp_ID]]$matrices$vaf*mito_data[[Exp_ID]]$matrices$SW)[mut,]>define_pos])
    pos_samples<-pos_samples[pos_samples%in%tree.no.ancestral$tip.label]
    if(length(pos_samples)>1){
      MRCA_node<-find_latest_acquisition_node(tree.no.ancestral,pos_samples)
      return(MRCA_node)
    } else {
      return("fail")
    }
  })
  
  tree_root=ifelse(Exp_ID=="18pcw",236,getRoot(tree.no.ancestral))
  MRCA_is_root_muts<-low_level_het_muts[mut_MRCA==tree_root]
  return(MRCA_is_root_muts)
  
  if(resave_plots) {
    pdf(file = paste0(figures_dir,"Figure_03/Low_level_hetmut_",Exp_ID,"_phylos.pdf"),width = 18,height=5)
    par(mfrow=c(2,6))
    temp=lapply(MRCA_is_root_muts,function(mut){
      plot_tree(mito_data[[Exp_ID]]$tree.ultra,cex.label = 0,bars = mito_data[[Exp_ID]]$matrices$vaf[mut,],plot_axis = T,title = mut)
      text(x = 0, y=-0.05*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",paste0("Max VAF: ",round(max(mito_data[[Exp_ID]]$matrices$vaf[mut,]),digits = 3),"; Mean VAF: ",round(mito_data_young[[Exp_ID]]$stats_df$mean_vaf[mito_data_young[[Exp_ID]]$stats_df$mut==mut],digits=3),"; Cmean pval = ",mito_data_young[[Exp_ID]]$stats_df$Cmean_pval[mito_data_young[[Exp_ID]]$stats_df$mut==mut]),pos = 4)
    })
    dev.off()
  }
})
names(het_oocyte_muts)<-young_IDs
saveRDS(het_oocyte_muts,paste0(root_dir,"/data/het_oocyte_muts.Rds"))

##----------------------MODELLING ALLELIC DRIFT IN EMBRYOGENESIS----------------

##THIS IS DONE USING THE HETEROPLASMIC MUTATIONS FROM THE 8 PCW FOETUS

#Review the VAF distributions of the heteroplasmic oocyte mutations
as.data.frame(mito_data$`8pcw`$matrices$vaf[het_oocyte_muts[[1]],])%>%
  tibble::rownames_to_column(var="mut_ref")%>%
  gather(-mut_ref,key="Sample",value="vaf")%>%
  ggplot(aes(x=vaf))+
  geom_histogram()+
  theme_classic()+
  facet_grid(rows=vars(mut_ref))

#Define the WF drift function used in the simulations
wright_fisher_drift=function(starting_vaf,population_size,generation_time,total_time){
  curr_vaf=starting_vaf
  total_generations=round(total_time/generation_time)
  for(i in 1:total_generations) {
    curr_vaf<-rbinom(1,size = population_size,prob = curr_vaf)/population_size
  }
  return(curr_vaf)
}

#Set up empty lists to save results of the simulations and abc results
all_params_and_sumstats=list()
abc_res=list()
all_targets=list()
nsim=1e4

#Define the prior distribution for the generation times - this is on a log scale, so is between 0.1 and 500 days
#Have the same vector of generation times for all mutation simulations
log_gen_times<-runif(nsim,-1,2.7)

#Set up the matrices
for(j in 1:length(het_oocyte_muts[[1]])){
  mut=het_oocyte_muts[[1]][j]
  print(mut)
  params_and_sumstats_file=paste0(root_dir,"/data/Drift_ABC_foetal/params_and_sumstats_",mut,".RDS")
  param.names=c("starting_vaf","population_size","log_generation_time")
  sumstat.names=c("mean_vaf","min_vaf","max_vaf","quant_0.05","quant_0.5","quant_0.95","sd_vaf")
  
  #Get the summary statistics of the data - the "target" for the ABC
  vaf_vec=as.numeric(mito_data$`8pcw`$matrices$vaf[mut,mito_data$`8pcw`$tree$tip.label])
  target=c(mean(vaf_vec),min(vaf_vec), max(vaf_vec),quantile(vaf_vec,0.05),quantile(vaf_vec,0.5),quantile(vaf_vec,0.95),sd(vaf_vec))
  all_targets[[j]]<-target
  
  if(file.exists(params_and_sumstats_file)){
    print("Importing existing params and sumstats file")
    params_and_sumstats<-readRDS(params_and_sumstats_file)
    params=params_and_sumstats$params
    sumstats=params_and_sumstats$sumstats
    
  } else {
    print("No pre-existing params and sumstats file found. Starting simulations")
    
    #Set up empty matrices for saving the parameters ("params") and summary statistics ("sumstats) for each simulation
    params=matrix(NA,nrow=nsim,ncol=length(param.names),dimnames=list(1:nsim,param.names))
    sumstats=matrix(NA,nrow=nsim,ncol=length(sumstat.names),dimnames=list(1:nsim,sumstat.names))
    
    #Draw the parameters from the priors
    params[,"starting_vaf"]<-mean(vaf_vec) #Starting VAF set as the mean VAF across samples. This is the most likely starting VAF and makes the ABC most efficient.
    params[,"population_size"]<-450
    params[,"log_generation_time"]<-log_gen_times
    
    #Perform the actual simulations using the WF drift function, saving the summary statistics at the end of each simulation
    for(i in 1:nsim){
      if(i%%100==0){print(i)}
      generation_time=10^params[,"log_generation_time"][i]
      final_vaf_vec<-sapply(1:277,function(j) wright_fisher_drift(params[,"starting_vaf"][i],population_size=params[,"population_size"][i],generation_time=generation_time,total_time=56))
      sumstats[i,]<-c(mean(final_vaf_vec),min(final_vaf_vec),max(final_vaf_vec),quantile(final_vaf_vec,0.05),quantile(final_vaf_vec,0.5),quantile(final_vaf_vec,0.95),sd(final_vaf_vec))
    }
    
    #Combine and save the params and sumstats objects for this mutation
    params_and_sumstats=list(params=params,sumstats=sumstats)
    saveRDS(object=params_and_sumstats,file = params_and_sumstats_file)
  }
  
  #Combine these objects into a single list
  all_params_and_sumstats[[j]]<-params_and_sumstats
  
  #Perform ABC for the individual mutations - for comparison with the combined ABC (to check that there are no major outliers)
  print("Performing the ABC inference")
  abc.nn<-abc::abc(target=target,param=params,sumstat=sumstats,tol=0.1,transf = c("none","log","none"),method="neuralnet")
  abc_res[[j]]<-abc.nn
  print("Saving the ABC output")
  saveRDS(abc.nn,file=paste0(root_dir,"/data/Drift_ABC_foetal/abc.nn.file.",mut,".RDS"))
}

#Do combined abc using sumstats from all mutations
#The params (generation times) are matched for all simulations
all_sumstats<-Map(list=all_params_and_sumstats,mut=het_oocyte_muts[[1]],f=function(list,mut) {
  c_names<-colnames(list$sumstats)
  colnames(list$sumstats)<-paste(mut,c_names,sep="_")
  return(list$sumstats)
})

all_sumstats_mat<-Reduce(cbind,all_sumstats)
all_params<-all_params_and_sumstats[[1]]$params[,3,drop=F]
all_targets_vec=Reduce(c,all_targets)
abc.comb.nn<-abc::abc(target=all_targets_vec,param=all_params,sumstat=all_sumstats_mat,tol=0.1,transf = c("none"),method="neuralnet")

foetal_gentime_posterior<-10^(quantile(abc.comb.nn$adj.values,c(0.025,0.5,0.975)))
foetal_gentime_posterior_plot<-data.frame(type="posterior",log_gen_time=abc.comb.nn$adj.values)%>%
  rbind(data.frame(type="prior",log_gen_time=as.numeric(all_params)))%>%
  mutate(gen_time=10^log_gen_time)%>%
  ggplot(aes(x=gen_time,col=type))+
  stat_density(geom="line",position="identity")+
  labs(x="WF generation time (days)",y="Density",col="")+
  theme_bw()+
  scale_x_log10(limits=c(0.1,500))+
  scale_y_continuous(limits=c(0,20))+
  my_theme+
  theme(legend.position = "top",legend.box.spacing = unit(0,"mm"))+
  annotate(geom="text",
           x=10,
           y=17,
           size=2,
           label=paste0("WF generation time = ",round(foetal_gentime_posterior[2],digits = 2)," days\n (95% PI: ",round(foetal_gentime_posterior[1],digits = 2),"-",round(foetal_gentime_posterior[3],digits = 2),")"))

ggsave(filename = paste0(figures_dir,"Figure_04/Foetal_hetmuts_gen_posterior.pdf"),foetal_gentime_posterior_plot,width = 1.75,height=2)

##Review the individual posterior distributions
gen_time_post_df<-Map(abc=abc_res,mut=het_oocyte_muts[[1]],f=function(abc,mut){
  return(data.frame(mut=mut,log_gen_time=abc$unadj.values[,3]))
})%>%
  dplyr::bind_rows()

gen_time_post_df%>%
  mutate(gen_time=10^log_gen_time)%>%
  ggplot(aes(x=gen_time))+
  #geom_histogram()+
  geom_density()+
  scale_x_log10()+
  facet_grid(rows=vars(mut))+
  theme_bw()+
  labs(x="WF generation time (days)",y="Count")+
  my_theme


#-----------------------------------------------------------------------------------#
# Generate_ExtData_Fig9.R
#
# Extended Data Fig. 9 - clonal marking of expanded clades by mtDNA mutations.
#
#   a  clonal fraction of expanded clades per individual
#   b  correlation between clade MRCA and the proportion of samples marked
#   c  phylogenetic signal against the number of samples sharing a mutation
#   d  phylogenetic signal against global VAF
#   e  as d, binned
#
# Panels were previously produced by Generate_Fig5.R; this script generates them
# on their own, so each manuscript figure has one script. The data preparation is
# therefore shared with Generate_Fig5.R by design.
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# --------Load packages (and install if they are not installed yet)-------------------
#-----------------------------------------------------------------------------------#
cran_packages=c("devtools","ape","stringr","dplyr","tidyr","ggplot2","gridExtra","phylosignal")
bioconductor_packages=c("MutationalPatterns","BSgenome","BSgenome.Hsapiens.UCSC.hg19","TxDb.Hsapiens.UCSC.hg19.knownGene")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if (!require("BiocManager", quietly = T, warn.conflicts = F))
  install.packages("BiocManager")
for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

if(!require("dndscv", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("im3sanger/dndscv")
  library("dndscv",character.only=T,quietly = T, warn.conflicts = F)
}

#-----------------------------------------------------------------------------------#
# ----------------------------------Set paths and import files------------------------
#-----------------------------------------------------------------------------------#

options(stringsAsFactors = F)

#Set these file paths before running the script
#genomeFile is set in config.R
source(here::here("config.R")) #sets root_dir, genomeFile, project paths and my_theme
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"/data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
#plots_dir, rebuttal_figs_dir and my_theme all come from config.R

#Create the figure output directories if they do not already exist
for(d in c("Extended_Data_Figure_09")) dir.create(paste0(plots_dir,d),showWarnings=FALSE,recursive=TRUE)

#Read in the mitochondrial copy number data
mito_cn=read.csv(paste0(root_dir,"/data/whole_genome_coverage_pileup_and_bedtools_annotated.csv"),header=T)

#Combine the sample level metadata info for the adult and foetal blood samples
sample_level_metadata_EM<-read.csv(paste0(root_dir,"/data/EM_sample_level_metadata.csv"))%>%dplyr::rename("exp_ID"="donor_id","Sample"="PDID","Cell_type"="cell_type","coverage"="mean_depth")
sample_level_metadata_foetal<-read.csv(paste0(root_dir,"/data/foetal_sample_level_metadata.csv"))%>%mutate(Sample=paste(Donor_ID,"hum",sep="_"))%>%dplyr::select(-Donor_ID,-Percentage)
sample_level_metadata<-dplyr::bind_rows(sample_level_metadata_EM,sample_level_metadata_foetal)

phenotype_data_EM<-read.csv(paste0(root_dir,"/data/Summary_pheno_pdid.csv"),stringsAsFactors = F)%>%
  dplyr::rename("exp_ID"=Donor_ID,"Sample"=PDID)

#Import individual level metadata for the adult and foetal blood samples
ref_df=read.csv(ref_file)%>%filter(Dataset!="Lymphocyte")
Individual_cols=RColorBrewer::brewer.pal(12,"Paired")
names(Individual_cols)<-ref_df$Sample[order(ref_df$Age)]

mito_cn<-mito_cn%>%
  left_join(sample_level_metadata,by="Sample",relationship="many-to-many")%>%
  left_join(ref_df%>%dplyr::rename("exp_ID"=Sample))
mito_cn$exp_ID<-factor(mito_cn$exp_ID,levels=ref_df$Sample[order(ref_df$Age)]) #Make the exp_ID a factor, with levels increasing by individual age
mito_cn<-left_join(mito_cn,phenotype_data_EM,relationship="many-to-many")

#Now import the mitochondrial mutation data
mito_data_file=paste0(root_dir,"/data/mito_data.Rds")
mito_data<-readRDS(mito_data_file)
CN_correlating_muts<-readRDS(paste0(root_dir,"/data/CN_correlation.RDS"))

#Define the 'black listed' mutation set - those with recurent artefacts despite the Shearwater filtering
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C")

#Define the "old individuals" used to assess mitochondrial mutations as lineage tracing markers
old_individuals=c("KX003","KX004","KX007","KX008")

#-----------------------------------------------------------------------------------#
# --------MEASURE 'CLONAL MARKING' OF EXPANDED CLADES BY MITOCHONDRIAL MUTATIONS-----
#-----------------------------------------------------------------------------------#

#Define function used to recognise the clonal expansions
get_expanded_clade_nodes=function(tree,height_cut_off=100,min_clonal_fraction=0.02,min_samples=1){
  nodeheights=nodeHeights(tree)
  
  #This pulls out nodes that fulfill on the criteria: branches cross the cut-off & contain the minimum proportion of samples
  nodes=tree$edge[,2][nodeheights[,1] < height_cut_off &
                        !nodeheights[,2] < height_cut_off &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))/length(tree$tip.label)})>min_clonal_fraction &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))})>=min_samples]
  df=data.frame(nodes=nodes,n_samples=sapply(nodes,function(node) {length(getTips(tree,node))}),MRCA_time=sapply(nodes,function(node) {nodeheight(tree,node)}),clonal_fraction=sapply(nodes,function(node) {length(getTips(tree,node))/length(tree$tip.label)}))
  return(df)
}
#-----------------------------------------------------------------------------------#
### Generate EXTENDED DATA FIG. 9A ---------
#-----------------------------------------------------------------------------------#

expanded.clades.plot<-dplyr::bind_rows(Map(list=mito_data[old_individuals],exp_ID=old_individuals,function(list,exp_ID){
  exp_nodes<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = 100,min_clonal_fraction=0.01)
  exp_nodes$exp_ID<-exp_ID
  return(exp_nodes)
}))%>%
  arrange(desc(n_samples))%>%
  ggplot(aes(x=factor(exp_ID,levels = c("KX007","KX008","KX004","KX003")),y=clonal_fraction,fill=n_samples))+
  geom_bar(stat="identity",position="stack",col="black")+
  scale_fill_gradient(low="lightgrey",high="darkred")+
  labs(x="Individual",y="Clonal fraction",fill="Number of\n samples\n in clone")+
  theme_bw()+
  my_theme
ggsave(filename=paste0(plots_dir,"Extended_Data_Figure_09/ExtDataFig9a.expanded_clades_plot.pdf"),expanded.clades.plot,width=3,height=2.5)

expanded_clades_df<-Map(list=mito_data[old_individuals],exp_ID=old_individuals,function(list,exp_ID){
  cat(paste0(exp_ID,"\n"))
  
  mutCN_cutoff=25 #If the mutant mitochondrial copy number is over 25, retain mutation even if is in the "CN correlating muts" list, this number is set empirically.
  CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
  vaf.filt<-(list$matrices$vaf*list$matrices$SW*(list$matrices$ML_Sig=="N1")*CN_correlating_mut_removal_mat)
  
  #Review how many expanded clades have reliable mitochondrial marker mutations
  marker_mut_cutoff<-0.01
  pos_mut_cutoff<-0.01
  exp_nodes<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = 100,min_clonal_fraction=0.01)
  full_df<-dplyr::bind_cols(data.frame(exp_ID=rep(exp_ID,nrow(exp_nodes))),
                            exp_nodes,
                            data.frame(marker_mut_cutoff=rep(marker_mut_cutoff,nrow(exp_nodes))),
                            exp_nodes_muts<-dplyr::bind_rows(lapply(exp_nodes$nodes,function(node) {
                              node_samples=getTips(list$tree.ultra,node)
                              if(any(vaf.filt[,node_samples]>marker_mut_cutoff)){
                                node_homo_muts<-names(rowSums(vaf.filt[,node_samples,drop=F]>marker_mut_cutoff)[rowSums(vaf.filt[,node_samples,drop=F]>marker_mut_cutoff)>0])
                                node_homo_muts<-node_homo_muts[!is.na(node_homo_muts)]
                                
                                pos_samples_per_mut<-sapply(node_homo_muts,function(mut){
                                  n_samples<-sum(vaf.filt[mut,node_samples]>pos_mut_cutoff)
                                  return(n_samples)
                                })
                                
                                mean_het_of_pos<-sapply(node_homo_muts,function(mut){
                                  pos_samples<-node_samples[vaf.filt[mut,node_samples]>pos_mut_cutoff]
                                  return(mean(as.numeric(vaf.filt[mut,pos_samples])))
                                })
                                
                                return(data.frame(nmuts=length(node_homo_muts),
                                                  homo_muts=paste0(node_homo_muts,collapse=","),
                                                  BMM=node_homo_muts[which.max(pos_samples_per_mut)],
                                                  pos_samples_per_mut=paste0(pos_samples_per_mut,collapse=","),
                                                  max_pos_samples=max(pos_samples_per_mut),
                                                  max_pos_prop=max(pos_samples_per_mut)/length(node_samples),
                                                  mean_heteroplasmy=mean(as.numeric(vaf.filt[node_homo_muts[which.max(pos_samples_per_mut)],node_samples])),
                                                  mean_het_of_pos=mean_het_of_pos[which.max(pos_samples_per_mut)]))
                              } else {
                                return(data.frame(nmuts=0,homo_muts=NA,pos_samples_per_mut=NA,max_pos_samples=NA,max_pos_prop=NA,mean_heteroplasmy=NA))
                              }
                            })))
  return(dplyr::select(full_df,-homo_muts,-pos_samples_per_mut))
})%>%dplyr::bind_rows()
#-----------------------------------------------------------------------------------#
### Generate EXTENDED DATA FIG. 9B ---------
#-----------------------------------------------------------------------------------#

#Show correlation of lineage marker with the time of the most recent common ancestor of the clone (MRCA)
MRCA.prop.correlation<-expanded_clades_df%>%
  ggplot(aes(x=MRCA_time,y=max_pos_prop,col=mean_heteroplasmy))+
  geom_point(aes(size=clonal_fraction),alpha=0.75)+
  scale_x_continuous(limits=c(0,NA))+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))+
  labs(x="Molecular time of clade's MRCA",
       y="Proportion of clade\nwith best lineage marker ",
       col="Mean\nheteroplasmy",
       size="Clade size\n(clonal fraction)")+
  theme_bw()+
  geom_smooth(col="black",size=0.6,method="lm")+
  my_theme
summary(lm(max_pos_prop~MRCA_time,data=expanded_clades_df))
ggsave(filename=paste0(plots_dir,"Extended_Data_Figure_09/ExtDataFig9b.MRCA_prop_correlation_plot.pdf"),MRCA.prop.correlation,width=4,height=2.5)
#-----------------------------------------------------------------------------------#
#----PULL OUT THE SHARED MUTATIONS - embed as an additional object within the mito_data list-----
#-----------------------------------------------------------------------------------#

#Find those mitochondrial mutations that are present in more than one sample at the specified cut_off
#Calculate the "phylogenetic signal" of each of these mutations i.e. the degree to which they follow the phylogeny
mito_data<-Map(list=mito_data,exp_ID=names(mito_data), function(list,exp_ID) {
  cat(exp_ID,sep="\n")
  mut_vaf_cutoff=0.01
  
  vaf.filt<-list$matrices$vaf.filt
  tree.ultra<-list$tree.ultra
  
  shared_muts<-rownames(vaf.filt)[rowSums(vaf.filt>mut_vaf_cutoff,na.rm = T)>1]
  n_pos<-rowSums(vaf.filt[shared_muts,]>mut_vaf_cutoff,na.rm = T)
  print(paste("There are",length(shared_muts),"shared mutations"))
  
  #Test shared muts with phylosignal
  if("phylosignal"%in%names(list)) {
    res_cor<-list$phylosignal
  } else {
    tree4d<-phylobase::phylo4d(drop.tip(tree.ultra,"Ancestral"),tip.data=t(vaf.filt[shared_muts,list$tree$tip.label]))
    res_cor=phyloSignal(tree4d)
  }
  
  #Add phylosignal info onto the dataframe
  shared_muts_df<-data.frame(mut=shared_muts,n_pos=n_pos)%>%
    tidyr::separate(mut,into=c("Chrom","Pos","Ref","Alt"),sep="_",remove=F)%>%
    mutate(Pos=as.numeric(Pos))%>%
    mutate(global_VAF=list$matrices$vaf[mut,"global"])%>%
    mutate(lambda_pval=res_cor$pval[mut,"Lambda"],Cmean_pval=res_cor$pval[mut,"Cmean"])%>%
    mutate(lambda=res_cor$stat[mut,"Lambda"],Cmean=res_cor$stat[mut,"Cmean"])%>%
    arrange(Pos)
  list$shared_muts_df<-shared_muts_df
  return(list)
})

#Save the file with the new shared mutations data
#Only write back when the cached result is actually new. This keeps the
#deposited data file byte-stable: without the guard, simply running this
#script changes mito_data.Rds and so invalidates its published checksum
#(see data/zenodo_manifest.csv).
if(!all(sapply(mito_data,function(list) !is.null(list$shared_muts_df)))) {
  saveRDS(mito_data,file=mito_data_file)
} else {
  cat("shared_muts_df already present for every donor - mito_data.Rds left unchanged\n")
}

#Combine this info into a single data frame
shared_muts_old_combined<-dplyr::bind_rows(Map(exp_ID=names(mito_data),list=mito_data,function(exp_ID,list) {
  bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
  list$shared_muts_df%>%
    left_join(bb_df)%>%
    mutate(exp_ID=exp_ID)
}))%>%mutate(signif=Cmean_pval<0.05)%>%
  dplyr::filter(exp_ID%in%old_individuals & !mut%in%CN_correlating_muts)
#-----------------------------------------------------------------------------------#
### Generate EXTENDED DATA FIG. 9C ---------
#-----------------------------------------------------------------------------------#

phylosignal.by.nsamples<-shared_muts_old_combined%>%
  group_by(exp_ID,signif,n_pos)%>%
  summarise(n=n())%>%
  mutate(n_pos_limited=ifelse(n_pos>=10,"≥10",as.character(n_pos)))%>%
  mutate(n_pos_limited=factor(n_pos_limited,levels=c(as.character(2:10),"≥10")))%>%
  ggplot(aes(x=n_pos_limited,y=n,fill=signif))+
  geom_bar(stat="identity")+
  theme_bw()+
  labs(x=str_wrap("Number of samples sharing mutation (VAF > 1%)",width=30),
       y="Count",
       fill="Significant \nphylogenetic \nsignal")+
  my_theme+theme(legend.margin = margin(t=0.1,unit="mm"))

ggsave(filename=paste0(plots_dir,"Extended_Data_Figure_09/ExtDataFig9c.Phylosignal_by_nsamples.pdf"),phylosignal.by.nsamples,width=2,height=2)

#-----------------------------------------------------------------------------------#
### Generate EXTENDED DATA FIG. 9D ---------
#-----------------------------------------------------------------------------------#

phylosignal.by.globalvaf<-shared_muts_old_combined%>%
  ggplot(aes(x=global_VAF,y=Cmean_pval,col=signif,size=n_pos))+
  geom_point(alpha=0.5)+
  scale_x_log10(labels=scales::label_number(accuracy = 0.001))+
  scale_y_log10()+
  geom_vline(xintercept = 0.005,linetype=2)+
  theme_bw()+
  labs(x="Global VAF",
       y="Cmean p-value",
       size="No of positive\nsamples",
       col="Significant\nphylogenetic\nsignal")+
  my_theme+theme(legend.margin = margin(t=0.1,unit="mm"))
ggsave(filename=paste0(plots_dir,"Extended_Data_Figure_09/ExtDataFig9d.Phylosignal_by_globalvaf.pdf"),phylosignal.by.globalvaf,width=3.5,height=2)

phylosignal.by.globalvaf.binned<-shared_muts_old_combined%>%
  mutate(bin=ifelse(global_VAF<0.005,"<0.5%",ifelse(global_VAF>0.01,">1%","0.5-1%")))%>%
  group_by(bin,signif)%>%
  summarise(n=n())%>%
  ggplot(aes(x=factor(bin,levels=c("<0.5%","0.5-1%",">1%")),y=n,fill=signif))+
  geom_bar(stat="identity",position="stack")+
  theme_bw()+
  labs(x="Global VAF group",y="Count",fill="Significant\nphylogenetic\nsignal")+
  my_theme+theme(legend.margin = margin(t=0.1,unit="mm"),axis.text.x = element_text(angle = 90))
ggsave(filename=paste0(plots_dir,"Extended_Data_Figure_09/ExtDataFig9e.Phylosignal_by_globalvaf_binned.pdf"),phylosignal.by.globalvaf.binned,width=1.5,height=2)

#-----------------------------------------------------------------------------------#
# Generate_ExtData_Fig8.R
#
# Extended Data Fig. 8 - mtDNA mutation heatmaps on the phylogenies of the older
# individuals, for mutations with and without significant phylogenetic signal.
# One pair of panels per individual (KX003, KX007, KX008); the equivalent panels
# for KX004 are Fig. 5a and are produced by Generate_Fig5.R.
#
# Panels were previously produced by Generate_Fig5.R, which wrote KX004 to
# Figure_05/ and the rest here via an ifelse(). This script generates Extended
# Data Fig. 8 on its own, so each manuscript figure has one script. The data
# preparation is therefore shared with Generate_Fig5.R by design.
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
for(d in c("Extended_Data_Figure_08")) dir.create(paste0(plots_dir,d),showWarnings=FALSE,recursive=TRUE)

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
### Generate EXTENDED DATA FIG. 8 ---------
#-----------------------------------------------------------------------------------#

#Only the individuals whose panels form Extended Data Fig. 8; KX004 is Fig. 5a.
ed8_individuals<-old_individuals[old_individuals!="KX004"]

#Generate the mutations with/ without phylosignal separately
#These can then been combined in illustrator/ inkscape

#Visualize the mutations that show significant phylogenetic signal
col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
temp=Map(list=mito_data[ed8_individuals],exp_ID=ed8_individuals,f=function(list,exp_ID){
  bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
  plot_muts<-list$shared_muts_df%>%
    filter(Cmean_pval<0.05)%>%
    pull(mut)
  
  vaf.mtx<-list$matrices$vaf*list$matrices$SW
  
  names(col_scheme)<-seq(0,1,0.01)
  hm<-matrix(0,nrow=length(plot_muts),ncol=length(list$tree$tip.label),dimnames = list(plot_muts,list$tree$tip.label))
  for(i in 1:length(plot_muts)) {
    mut<-plot_muts[i]
    mut_vafs<-vaf.mtx[mut,list$tree$tip.label]
    mut_vafs[mut_vafs<0.01]<-0
    mut_vafs<-round(mut_vafs,digits=2)
    hm[i,]<-col_scheme[as.character(mut_vafs)]
  }
  plot_muts.clustered<-hclust(dist(vaf.mtx[plot_muts,]))
  par(mfrow=c(1,1))
  fig<-"Extended_Data_Figure_08/"
  pdf(file = paste0(plots_dir,fig,exp_ID,"_phylosignal.pdf"),width = 7,height=4)
  list$tree.ultra=plot_tree(tree = list$tree.ultra,cex.label = 0,plot_axis=F,vspace.reserve = 3.1)
  add_mito_mut_heatmap(tree=list$tree.ultra,heatmap=hm[plot_muts.clustered$order,],border="gray",heatmap_bar_height=0.05,cex.label = 0.25)
  dev.off()
})

#Visualize the mutations that do not show significant phylogenetic signal
temp=Map(list=mito_data[ed8_individuals],exp_ID=ed8_individuals,f=function(list,exp_ID){
  bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
  plot_muts<-list$shared_muts_df%>%
    filter(Cmean_pval>0.05)%>%
    pull(mut)
  
  vaf.mtx<-list$matrices$vaf*list$matrices$SW
  
  names(col_scheme)<-seq(0,1,0.01)
  hm<-matrix(0,nrow=length(plot_muts),ncol=length(list$tree$tip.label),dimnames = list(plot_muts,list$tree$tip.label))
  for(i in 1:length(plot_muts)) {
    mut<-plot_muts[i]
    mut_vafs<-vaf.mtx[mut,list$tree$tip.label]
    mut_vafs[mut_vafs<0.01]<-0
    mut_vafs<-round(mut_vafs,digits=2)
    hm[i,]<-col_scheme[as.character(mut_vafs)]
  }
  plot_muts.clustered<-hclust(dist(vaf.mtx[plot_muts,]))
  par(mfrow=c(1,1))
  fig<-"Extended_Data_Figure_08/"
  pdf(file = paste0(plots_dir,fig,exp_ID,"_no_phylosignal.pdf"),width = 7,height=4)
  list$tree.ultra=plot_tree(tree = list$tree.ultra,cex.label = 0,plot_axis=F,vspace.reserve = 3)
  add_mito_mut_heatmap(tree=list$tree.ultra,heatmap=hm[plot_muts.clustered$order,],border="gray",heatmap_bar_height=0.05,cex.label = 0.25)
  dev.off()
})

#Plot the scale legend for the VAF colour scheme
pdf(file=paste0(plots_dir,"Figure_05/Heatmap_scale_bar.pdf"),width=2,height=5)
par(mfrow=c(1,1))
autoimage::legend.scale(
  c(0,1),
  col = col_scheme,
  horizontal = F
)
dev.off()

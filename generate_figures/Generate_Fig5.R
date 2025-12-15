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
genomeFile="~/Documents/Reference_files/genome.fa" #This should be the hg37 genome file
root_dir="~/R_work/mito_mutations"
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"/data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
plots_dir=paste0(root_dir,"/plots/")
rebuttal_figs_dir=paste0(root_dir,"/rebuttal_plots/")

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
### Generate SUPPLEMENTARY FIG. 9A ---------
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
ggsave(filename=paste0(plots_dir,"Supp_Figure_09/SuppFig9a.expanded_clades_plot.pdf"),expanded.clades.plot,width=3,height=2.5)

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
### Generate FIG. 5B ---------
#-----------------------------------------------------------------------------------#

node_factor_levels=expanded_clades_df%>%arrange(n_samples)%>%mutate(levels=str_c(exp_ID,nodes,sep = "_"))%>%pull(levels)
expanded.clades.marking.plot<-expanded_clades_df%>%
  mutate(n_neg_samples=n_samples-max_pos_samples)%>%
  dplyr::select(exp_ID,nodes,n_samples,max_pos_samples,n_neg_samples)%>%
  mutate(max_pos_samples=ifelse(max_pos_samples==1,0,max_pos_samples))%>%
  mutate(n_neg_samples=n_samples-max_pos_samples)%>%
  gather(-exp_ID,-nodes,-n_samples,key="Pos_or_neg",value="n_samples")%>%
  mutate(Pos_or_neg=ifelse(Pos_or_neg=="n_neg_samples","Absent","Present"))%>%
  mutate(Pos_or_neg=factor(Pos_or_neg,levels=c("Present","Absent")))%>%
  mutate(levels=str_c(exp_ID,nodes,sep = "_"))%>%
  ggplot(aes(x=factor(levels,levels=node_factor_levels),y=n_samples,fill=Pos_or_neg))+
  geom_bar(position="stack",stat="identity",col="black",linewidth=0.15)+
  scale_fill_brewer(palette="Set2")+
  facet_grid(~exp_ID,scales="free",space = "free")+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  labs(x="Clonal expansion",y=str_wrap("Number of samples within expansion",width=15),fill=str_wrap("Best mitochondrial marker mutation",width=10))+
  my_theme
ggsave(filename=paste0(plots_dir,"Figure_05/Fig5b.expanded_clades_marking_plot.pdf"),expanded.clades.marking.plot,width=7,height=2.3)

#-----------------------------------------------------------------------------------#
### Generate SUPPLEMENTARY FIG. 9B ---------
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
ggsave(filename=paste0(plots_dir,"Supp_Figure_09/SuppFig9b.MRCA_prop_correlation_plot.pdf"),MRCA.prop.correlation,width=4,height=2.5)

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
saveRDS(mito_data,file=mito_data_file)

#Combine this info into a single data frame
shared_muts_old_combined<-dplyr::bind_rows(Map(exp_ID=names(mito_data),list=mito_data,function(exp_ID,list) {
  bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
  list$shared_muts_df%>%
    left_join(bb_df)%>%
    mutate(exp_ID=exp_ID)
}))%>%mutate(signif=Cmean_pval<0.05)%>%
  dplyr::filter(exp_ID%in%old_individuals & !mut%in%CN_correlating_muts)

#-----------------------------------------------------------------------------------#
### Generate SUPPLEMENTARY FIG. 9C ---------
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

ggsave(filename=paste0(plots_dir,"Supp_Figure_09/SuppFig9c.Phylosignal_by_nsamples.pdf"),phylosignal.by.nsamples,width=2,height=2)

#-----------------------------------------------------------------------------------#
### Generate SUPPLEMENTARY FIG. 9D ---------
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
ggsave(filename=paste0(plots_dir,"Supp_Figure_09/SuppFig9d.Phylosignal_by_globalvaf.pdf"),phylosignal.by.globalvaf,width=3.5,height=2)

phylosignal.by.globalvaf.binned<-shared_muts_old_combined%>%
  mutate(bin=ifelse(global_VAF<0.005,"<0.5%",ifelse(global_VAF>0.01,">1%","0.5-1%")))%>%
  group_by(bin,signif)%>%
  summarise(n=n())%>%
  ggplot(aes(x=factor(bin,levels=c("<0.5%","0.5-1%",">1%")),y=n,fill=signif))+
  geom_bar(stat="identity",position="stack")+
  theme_bw()+
  labs(x="Global VAF group",y="Count",fill="Significant\nphylogenetic\nsignal")+
  my_theme+theme(legend.margin = margin(t=0.1,unit="mm"),axis.text.x = element_text(angle = 90))
ggsave(filename=paste0(plots_dir,"Supp_Figure_09/SuppFig9e.Phylosignal_by_globalvaf_binned.pdf"),phylosignal.by.globalvaf.binned,width=1.5,height=2)

#-----------------------------------------------------------------------------------#
### Generate FIG. 5A/ Supplementary FIG 8 ---------
#-----------------------------------------------------------------------------------#

#Generate the mutations with/ without phylosignal separately
#These can then been combined in illustrator/ inkscape

#Visualize the mutations that show significant phylogenetic signal
col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
temp=Map(list=mito_data[old_individuals],exp_ID=old_individuals,f=function(list,exp_ID){
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
  fig<-ifelse(exp_ID=="KX004","Figure_05/","Supp_Figure_08/")
  pdf(file = paste0(plots_dir,fig,exp_ID,"_phylosignal.pdf"),width = 7,height=4)
  plot_tree(tree = list$tree.ultra,cex.label = 0,plot_axis=F,vspace.reserve = 3.1)
  add_mito_mut_heatmap(tree=list$tree.ultra,heatmap=hm[plot_muts.clustered$order,],border="gray",heatmap_bar_height=0.05,cex.label = 0.25)
  dev.off()
})

#Visualize the mutations that do not show significant phylogenetic signal
temp=Map(list=mito_data[old_individuals],exp_ID=old_individuals,f=function(list,exp_ID){
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
  fig<-ifelse(exp_ID=="KX004","Figure_05/","Supp_Figure_08/")
  pdf(file = paste0(plots_dir,fig,exp_ID,"_no_phylosignal.pdf"),width = 7,height=4)
  plot_tree(tree = list$tree.ultra,cex.label = 0,plot_axis=F,vspace.reserve = 3)
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

#-----------------------------------------------------------------------------------#
### Generate FIG. 5C
#-----------------------------------------------------------------------------------#

## Review the Mitochondrial Mutations with patchy marking of subclones
selected_muts=data.frame(exp_ID=c("KX003","KX003","KX004","KX003","KX004","KX004","KX004","KX004"),
                         mut=c("MT_9151_A_G","MT_8610_T_C","MT_11790_T_C","MT_8157_T_C","MT_6379_T_C","MT_3332_T_C","MT_4965_A_G","MT_4232_T_C"))

#The colour scheme relates to the coverage of the mutation position in each sample
bins=seq(0,10000,50)
length(bins[bins<=1000])
coverage_col_scheme=c(colorRampPalette(brewer.pal(n=8,name="YlGnBu"))(length(bins[bins<=1000])),rep("#0C2C84",length(bins)-length(bins[bins<=1000])))
names(coverage_col_scheme)<-bins

#Define the custom function to add the coverage heatmaps
add_coverage_heatmap=function(tree,heatmap,heatvals=NULL,border="white",cex.label=2){
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.01*ymax
  gap=tree$vspace.reserve/dim(heatmap)[1]
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-0.025*ymax
    #bot=top-(0.05/dim(heatmap)[1])*ymax
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border,lwd = 0.25)
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      text(labels[i],x=-0.5,y=0.5*(top+bot),pos = 2,cex = cex.label)
    }
    top=bot
  }
  tree
}

pdf(file=paste0(plots_dir,"Figure_05/Coverage_scale_bar.pdf"),width=2,height=5)
par(mfrow=c(1,1))
autoimage::legend.scale(
  c(0,1000),
  col = coverage_col_scheme[1:length(bins[bins<=1000])],
  horizontal = F,
  axis.args=list(at=seq(100,1000,100),labels=c(seq(100,900,100),">1000"))
)
dev.off()

marker_mut_cutoff=0.01
for(i in 1:nrow(selected_muts)){
  exp_ID<-selected_muts$exp_ID[i]
  mut<-selected_muts$mut[i]
  
  list<-mito_data[[exp_ID]]
  vaf.mtx<-list$matrices$vaf*list$matrices$SW
  vaf.mtx<-vaf.mtx[,-which(colnames(vaf.mtx)=="global")]
  tree.ultra<-list$tree.ultra
  
  pos_samples<-names(vaf.mtx)[which(vaf.mtx[mut,]>marker_mut_cutoff)]
  print(pos_samples)
  
  latest_acquisition_node=find_latest_acquisition_node(tree.ultra,pos_samples)
  
  sub_tree=drop.tip(tree.ultra,tree.ultra$tip.label[!tree.ultra$tip.label%in%c("Ancestral",getTips(tree.ultra,latest_acquisition_node))],trim.internal = T)
  sub_tree$coords<-NULL
  pdf(file=paste0(plots_dir,"Figure_05/Fig5c.",exp_ID,"_",mut,"_sub_tree_with_coverage.pdf"),width=2,height=2.5)
  
  mut_string<-stringr::str_split(mut,pattern="_",simplify=T)
  sub_tree=plot_tree(sub_tree,cex.label=0,bars = vaf.mtx[mut,],vspace.reserve = 5,cex.axis=0.5,title = paste0("MT ",mut_string[2],": ",mut_string[3],">",mut_string[4]))
  text(x = 0, y=-0.25*par()[['yaxp']][2],cex = 0.5,font=3,col="#00000095",paste0("Max VAF: ",round(max(vaf.mtx[mut,]),digits = 3)),pos = 4)
  text(x = 0, y=1.1*par()[['yaxp']][2],cex = 0.75,font=3,col="#00000095",mut,pos = 4)
  
  #Plot heatmap of coverage at the mutation site (show that zero VAF is not due to low coverage)
  hm<-matrix(0,nrow=1,ncol=length(sub_tree$tip.label),dimnames = list("Coverage",sub_tree$tip.label))
  mut_coverage<-list$matrices$NR[mut,sub_tree$tip.label[which(sub_tree$tip.label!="Ancestral")]]
  mut_coverage_rounded<-round(mut_coverage,digits = -2)
  mut_coverage_rounded[mut_coverage_rounded>10000]<-10000
  hm[1,]<-c(coverage_col_scheme[as.character(mut_coverage_rounded)],NA)
  add_coverage_heatmap(tree=sub_tree,heatmap=hm,border="gray",cex.label = 0.4)
  dev.off()
}

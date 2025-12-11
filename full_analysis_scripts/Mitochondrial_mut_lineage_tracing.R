library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(phylosignal)
options(stringsAsFactors = F)

#Set these file paths before running the script
genomeFile="~/R_work/reference_files/genome.fa"
root_dir="~/R_work/mito_mutations_blood/"
source(paste0(root_dir,"data/mito_mutations_blood_functions.R"))

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"data/Samples_metadata_ref.csv")
figures_dir=paste0(root_dir,"figures/")

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
mito_data_file=paste0(root_dir,"data/mito_data.Rds")
mito_data<-readRDS(mito_data_file)
CN_correlating_muts<-readRDS(paste0(root_dir,"data/CN_correlation.RDS"))
CN_correlating_muts<-CN_correlating_muts[-which(CN_correlating_muts=="MT_16519_T_T")]

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
ggsave(filename=paste0(figures_dir,"Supp_Figure_06/expanded_clades_plot.pdf"),expanded.clades.plot,width=3,height=2.5)

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
  geom_bar(position="stack",stat="identity",col="black",size=0.15)+
  scale_fill_brewer(palette="Set2")+
  facet_grid(~exp_ID,scales="free",space = "free")+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  labs(x="Clonal expansion",y=str_wrap("Number of samples within expansion",width=15),fill=str_wrap("Best mitochondrial marker mutation",width=10))+
  my_theme
ggsave(filename=paste0(figures_dir,"Figure_05/expanded_clades_marking_plot.pdf"),expanded.clades.marking.plot,width=7,height=2.3)

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
ggsave(filename=paste0(figures_dir,"Supp_Figure_06/MRCA_prop_correlation_plot.pdf"),MRCA.prop.correlation,width=4,height=2.5)

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

shared_muts_old_combined%>%
  ggplot(aes(x=n_pos,y=global_VAF,col=signif))+
  geom_jitter(width=0.02,alpha=0.5)+
  scale_x_log10()+
  scale_y_log10()+
  theme_classic()+
  my_theme

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

ggsave(filename=paste0(figures_dir,"Supp_Figure_06/Phylosignal_by_nsamples.pdf"),phylosignal.by.nsamples,width=2,height=2)

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
ggsave(filename=paste0(figures_dir,"Supp_Figure_06/Phylosignal_by_globalvaf.pdf"),phylosignal.by.globalvaf,width=3.5,height=2)

phylosignal.by.globalvaf.binned<-shared_muts_old_combined%>%
  mutate(bin=ifelse(global_VAF<0.005,"<0.5%",ifelse(global_VAF>0.01,">1%","0.5-1%")))%>%
  group_by(bin,signif)%>%
  summarise(n=n())%>%
  ggplot(aes(x=factor(bin,levels=c("<0.5%","0.5-1%",">1%")),y=n,fill=signif))+
  geom_bar(stat="identity",position="stack")+
  theme_bw()+
  labs(x="Global VAF group",y="Count",fill="Significant\nphylogenetic\nsignal")+
  my_theme+theme(legend.margin = margin(t=0.1,unit="mm"),axis.text.x = element_text(angle = 90))
ggsave(filename=paste0(figures_dir,"Supp_Figure_06/Phylosignal_by_globalvaf_binned.pdf"),phylosignal.by.globalvaf.binned,width=1.5,height=2)

#Visualize the mutations that show significant phylogenetic signal
par(mfrow=c(2,2))
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
  fig<-ifelse(exp_ID=="KX004","Figure_05/","Supp_Figure_05/")
  pdf(file = paste0(figures_dir,fig,exp_ID,"_phylosignal.pdf"),width = 7,height=4)
  plot_tree(tree = list$tree.ultra,cex.label = 0,plot_axis=F,vspace.reserve = 3.1)
  add_mito_mut_heatmap(tree=list$tree.ultra,heatmap=hm[plot_muts.clustered$order,],border="gray",heatmap_bar_height=0.05,cex.label = 0.25)
  dev.off()
})

#Visualize the mutations that do not show significant phylogenetic signal
par(mfrow=c(2,2))
col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
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
  fig<-ifelse(exp_ID=="KX004","Figure_05/","Supp_Figure_05/")
  pdf(file = paste0(figures_dir,fig,exp_ID,"_no_phylosignal.pdf"),width = 7,height=4)
  plot_tree(tree = list$tree.ultra,cex.label = 0,plot_axis=F,vspace.reserve = 3)
  add_mito_mut_heatmap(tree=list$tree.ultra,heatmap=hm[plot_muts.clustered$order,],border="gray",heatmap_bar_height=0.05,cex.label = 0.25)
  dev.off()
})


#Plot the scale legend for the VAF colour scheme
pdf(file=paste0(figures_dir,"Supp_Figure_05/Heatmap_scale_bar.pdf"),width=2,height=5)
par(mfrow=c(1,1))
autoimage::legend.scale(
  c(0,1),
  col = col_scheme,
  horizontal = F
)
dev.off()

##----------------Review the Mitochondrial Mutations with patchy marking of subclones--------------------------
selected_muts=data.frame(exp_ID=c("KX003","KX003","KX004","KX003","KX004","KX004","KX004","KX004"),
                         mut=c("MT_9151_A_G","MT_8610_T_C","MT_11790_T_C","MT_8157_T_C","MT_6379_T_C","MT_3332_T_C","MT_4965_A_G","MT_4232_T_C"))

bins=seq(0,10000,50)
length(bins[bins<=1000])
coverage_col_scheme=c(colorRampPalette(brewer.pal(n=8,name="YlGnBu"))(length(bins[bins<=1000])),rep("#0C2C84",length(bins)-length(bins[bins<=1000])))
names(coverage_col_scheme)<-bins

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

pdf(file=paste0(figures_dir,"Figure_05/Coverage_scale_bar.pdf"),width=2,height=5)
par(mfrow=c(1,1))
autoimage::legend.scale(
  c(0,1000),
  col = coverage_col_scheme[1:length(bins[bins<=1000])],
  horizontal = F,
  axis.args=list(at=seq(100,1000,100),labels=c(seq(100,900,100),">10000"))
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
  pdf(file=paste0(figures_dir,"Figure_05/",exp_ID,"_",mut,"_sub_tree_with_coverage.pdf"),width=2,height=2.5)
  sub_tree=plot_tree(sub_tree,cex.label=0,bars = vaf.mtx[mut,],vspace.reserve = 5,cex.axis=0.5)
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

#-----------------------------------------------------------------------------------#
##----------------INFERRED CLONES----------------------
#-----------------------------------------------------------------------------------#

library(Seurat)

#Define function to get clusters
seuratSNN <- function(matSVD, resolution = 1, k.param = 10){ 
  set.seed(1)
  rownames(matSVD) <- make.unique(rownames(matSVD))
  obj <- FindNeighbors(matSVD, k.param = k.param, annoy.metric = "cosine")
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

output.dir <- paste0(root_dir,"data/mito_mut_clones")
dir.create(output.dir,showWarnings = F)

# store the information for the heatmap here
df.list <- list()
matrix.list <- list()

# iterate over multiple patients here (if there are multiple)
patient.ids<-c("KX003","KX004","KX007","KX008")
for (i in 1:length(patient.ids)){
  
  # read in one of the mgatk SE objects
  patient.tmp <- patient.ids[i]
  cat(paste0("\nProcessing patient...", patient.tmp, "\n"))
  vaf.mtx <- mito_data[[patient.tmp]]$matrices$vaf.filt
  
  #Remove mutations present in all samples & keep only those present at global VAF >0.5%
  vaf.mtx<-vaf.mtx[mito_data[[patient.tmp]]$matrices$vaf$global>0.005 &
                     !grepl("INS|DEL",rownames(mito_data[[patient.tmp]]$matrices$vaf)),]
  if(any(colnames(vaf.mtx)=="global")) {
    vaf.mtx<-vaf.mtx[,-which(colnames(vaf.mtx)=="global")]
  }
  
  
  # filter them according to a certain mean coverage
  # coverage.data <- coverage.data %>%  dplyr::filter(mean_cov >= 10)
  
  # remove variants which are not present at a VAF > 1% at least once
  mutations <- rownames(vaf.mtx)
  filtered.mutations <- mutations[rowSums(vaf.mtx > 0.01) > 1]
  vaf.filtered.mtx <- vaf.mtx[filtered.mutations,]
  
  #-----------------------------------------------------------------------------------#
  ## DEFINE CLONES BASED ON THE MITOCHONDRIAL MUTATIONS -------------------------------
  #-----------------------------------------------------------------------------------#
  
  # get clusters with cluster resolution 10 and knn 50
  clusters <- seuratSNN((t(vaf.filtered.mtx)),resolution= 10,k.param=5)
  clusters <- str_pad(clusters, 2, pad = "0")
  
  # assign colours to clusters
  names_clusters <- unique(clusters)
  cluster_cols<- c("lightgray","#1f77b4","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf","#ad494a","#e7ba52","#8ca252","#756bb1","#636363","#aec7e8", brewer.pal(12, "Paired"))
  vec_go <- cluster_cols[1:length(names_clusters)]
  names(vec_go) <- sort(names_clusters)
  
  #Change the colour of the biggest cluster to light grey - this is generally the "NULL" cluster with no informative mutations
  biggest_cluster<-names(table(clusters))[which.max(table(clusters))]
  original_cluster_0<-which(clusters=="00")
  new_cluster_0<-which(as.character(clusters)==biggest_cluster)
  clusters[new_cluster_0]<-"00"
  clusters[original_cluster_0]<-biggest_cluster
  
  # Make data.frame for cluster_id and sample relationship
  df <- data.frame(
    sample_id = colnames(vaf.filtered.mtx), 
    cluster_id = as.character(clusters), 
    stringsAsFactors = F
  ) %>% arrange(clusters)
  
  #Work out which cluster is the 'no shared muts' cluster and give a bland colour
  mean_cluster_muts<-sapply(unique(clusters),function(no){
    cluster_samples<-df%>%dplyr::filter(cluster_id==no)%>%pull(sample_id)
    sum((vaf.filtered.mtx[,cluster_samples]>0.01))/length(cluster_samples)
  })
  min_mut_cluster<-unique(clusters)[which.min(mean_cluster_muts)]
  vec_go[min_mut_cluster]<-"lightgrey"
  df <- df[!duplicated(df),]
  
  # store the df
  df.list[[i]] <- df
  write.table(df, paste0(output.dir, "/", patient.tmp, "_mtdna_clone_assignment.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  
  #-----------------------------------------------------------------------------------#
  ## Create a heatmap -------------------------------
  #-----------------------------------------------------------------------------------#
  
  # make a heatmap annotation
  ha_col <- HeatmapAnnotation(Cluster = as.character(df$cluster_id),
                              col = list(Cluster = vec_go))
  
  # replace variants below 1% VAF
  afp <- vaf.filtered.mtx 
  aftree <- afp
  #afp[afp < 0.01] <- 0
  afp[afp > 0.1] <- 0.1
  
  # order the mutations according to abundance
  mean_vaf <- rowMeans(afp)
  mut_order <- names(mean_vaf[order(mean_vaf, decreasing = T)])
  ordered_afp <- afp[mut_order,]
  colnames(ordered_afp) <- paste0(df$sample_id, "_", colnames(ordered_afp))
  
  # store the other matrix as well
  matrix.list[[i]] <- data.matrix(ordered_afp)
  
  # cluster hierarchically within each group
  # here we define the groups we want to cluster according to
  groups <- unique(clusters)[order(unique(clusters))]
  ordered_names <- c()
  
  # perform hierarchical clustering within each cluster
  k <- 1
  for (k in 1:length(groups)){
    
    # which group?
    print(groups[k])
    
    # get the barcodes from the respective group
    barcodes <- df[grep(groups[k], df$cluster_id), "sample_id"]
    
    if (length(barcodes) > 1){
      # iterate over all chrs and subset the matrix
      matrix <- vaf.filtered.mtx[,colnames(vaf.filtered.mtx) %in% barcodes]
      avgd <- colMeans(matrix)
      
      # Order cells with hierarchical clustering
      dist.centered.matrix <- stats::dist(as.matrix(avgd), method = "euclidean")
      hc <- hclust(dist.centered.matrix, method = "ward.D2")
      
      # make a vector with the right order of barcodes per group
      ordered_names <- c(ordered_names, hc$labels[hc$order]) 
      
    } else {
      next
    }
    
  }
  
  hm <- Heatmap((data.matrix(afp)[,as.character(ordered_names)]),  # var_order
                col=as.character(BuenColors::jdb_palette("solar_rojos",type="continuous")),
                show_row_names = FALSE, 
                top_annotation=ha_col,
                cluster_columns = FALSE,
                name = "AF",use_raster = FALSE,
                row_names_gp = gpar(fontsize = 10),
                cluster_rows = TRUE, 
                show_column_names = FALSE)
  
  # save heatmaps
  pdf(paste0(figures_dir, "Supp_Figure_08/", patient.tmp, "_mito_mutation_heatmap.pdf"), width=15, height=8)
  print(hm) 
  dev.off()
  
  #-----------------------------------------------------------------------------------#
  ## CREATE "TREES" AKA DENDROGRAMS -------------------------------
  #-----------------------------------------------------------------------------------#

  # Get group means 
  matty <- sapply(names(vec_go), function(cluster){
    cells <- df %>% dplyr::filter(cluster_id == cluster) %>% pull(sample_id) %>% as.character()
    Matrix::rowMeans(sqrt(afp[,cells]))
  })
  
  if(length(groups) > 2){
    
    # Do cosine distance; note that we used sqrt transformation already when creating the pseudo bulk-cell matrix
    mito.hc <- hclust(dist(lsa::cosine((matty))))
    plot(mito.hc)
    
    pdf(paste0(figures_dir,"Supp_Figure_08/", patient.tmp, "_hierarchical_tree.pdf"), width = 5, height = 5)
    print(plot(mito.hc))
    dev.off()
    
  } else {
    next
  }
}

#Import clustered clones to overlay onto true tree
par(mfrow=c(2,2))
cluster_cols<- c("lightgray","#1f77b4","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b","#e377c2","#7f7f7f",
                 "#bcbd22","#17becf","#ad494a","#e7ba52","#8ca252","#756bb1","#636363","#aec7e8", brewer.pal(12, "Paired"))
length(cluster_cols)
mito_data=Map(list=mito_data[old_individuals],exp_ID=old_individuals,function(list,exp_ID){
  exp_clones<-read.delim(paste0(root_dir,"data/mito_mut_clones/",exp_ID,"_mtdna_clone_assignment.txt"))
  n_clones<-length(unique(exp_clones$cluster_id))
  exp_cluster_cols<-cluster_cols[1:n_clones]
  names(exp_cluster_cols)<-unique(exp_clones$cluster_id)
  
  #Change the colour of the biggest cluster to light grey
  biggest_cluster<-names(table(exp_clones$cluster_id))[which.max(table(exp_clones$cluster_id))]
  original_cluster_0<-which(exp_clones$cluster_id==0)
  new_cluster_0<-which(as.character(exp_clones$cluster_id)==biggest_cluster)
  exp_clones$cluster_id[new_cluster_0]<-0
  exp_clones$cluster_id[original_cluster_0]<-as.integer(biggest_cluster)
  
  #Plot the tree annotated with the clusters
  par(mfrow=c(1,1))
  pdf(file=paste0(figures_dir,"Supp_Figure_08/Inferred_clones_",exp_ID,".pdf"),width = 7,height=3)
  plot_tree(tree = list$tree.ultra,cex.label = 0)
  clone_hm<-matrix(NA,nrow=1,ncol=length(list$tree.ultra$tip.label),dimnames = list("Clones",list$tree.ultra$tip.label))
  for(i in 1:nrow(exp_clones)){
    if(exp_clones$sample_id[i]%in%list$tree.ultra$tip.label){
      clone_hm[1,exp_clones$sample_id[i]]<-exp_cluster_cols[as.character(exp_clones$cluster_id[i])]
    }
  }
  
  add_heatmap(tree=list$tree.ultra,heatmap=clone_hm,border="gray",cex.label = 1)
  legend("topleft", inset=c(.01,.01), title="Clone no.",
         names(exp_cluster_cols), fill=exp_cluster_cols, horiz=F, cex=0.7,ncol=1)
  dev.off()
  list$clusters<-exp_clones
  list$cluster_cols<-exp_cluster_cols
  list$cluster_hm<-clone_hm
  return(list)
})

#Plot the clone assignments of the expanded clades
expanded_clades_cluster_assignments<-lapply(1:nrow(expanded_clades_df),function(i) {
  exp_ID<-expanded_clades_df$exp_ID[i]
  node<-expanded_clades_df$nodes[i]
  Samples<-getTips(mito_data[[exp_ID]]$tree.ultra,node = node)
  df_node<-data.frame(exp_ID=exp_ID,node=node,sample_id=Samples)%>%
    left_join(mito_data[[exp_ID]]$clusters)
  return(df_node)
})%>%
  dplyr::bind_rows()%>%
  mutate(node=paste(exp_ID,node,sep="_"))

#Add a proportions variable to adjust for clone size
expanded_clades_cluster_assignments$prop<-sapply(1:nrow(expanded_clades_cluster_assignments),function(i){
  clone_total_n=sum(expanded_clades_cluster_assignments$node==expanded_clades_cluster_assignments$node[i])
  return(1/clone_total_n)
})
n_clones=1+max(expanded_clades_cluster_assignments$cluster_id)
expansion.assignment.proportions.plot<-expanded_clades_cluster_assignments%>%
  ggplot(aes(x=factor(node,levels=expanded_clades_df%>%arrange(n_samples)%>%mutate(nodes=paste(exp_ID,nodes,sep="_"))%>%pull(nodes)),y=prop,fill=factor(cluster_id,levels=0:(n_clones-1))))+
  geom_bar(stat="identity",position="stack",col="black",size=0.05,width = 0.7)+
  scale_fill_manual(values = cluster_cols[1:n_clones],drop=F)+
  facet_grid(cols=vars(factor(exp_ID,levels=c("KX004","KX003","KX007","KX008"))),scales="free",space = "free")+
  theme_bw()+
  theme(axis.text.x = element_blank(),legend.key.height = unit(1.5,"mm"),legend.key.width = unit(5,"mm"))+
  my_theme+
  labs(fill="Clone\nassignments",
       x="Clonal expansions",
       y="Clone proportions")
ggsave(filename=paste0(figures_dir,"Supp_Figure_08/expansion_clone_assignment_props_plot.pdf"),expansion.assignment.proportions.plot,width=7,height=2)

#Look at 'inappropriate assignment' to clones i.e. singleton cells assigned to a cluster
singleton_samples<-Map(list=mito_data[old_individuals],Exp_ID=old_individuals,function(list,Exp_ID){
  #Get the 'expanded clades' - defined loosely as any clade with a MRCA with another sample > 100 mutations
  exp_expanded_clades<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off=100,min_clonal_fraction=0.001+(1/length(list$tree.ultra$tip.label)))
  #Get the list of samples within any of those clades
  expanded_clade_samples<-unlist(lapply(exp_expanded_clades$nodes,function(node) {getTips(tree=list$tree.ultra,node=node)}))
  #Return the list of samples that aren't in any of the expanded clades i.e. the 'singletons'
  return(list$tree.ultra$tip.label[!list$tree.ultra$tip.label%in%c(expanded_clade_samples,"Ancestral")])
})

#Summarise how many of these singletons there are, and how many have been assigned to a cluster that isn't '0' (the 'default' cluster with no shared mutations)
dplyr::bind_rows(Map(list=mito_data,exp_ID=names(mito_data),f=function(list,exp_ID) cbind(list$clusters,exp_ID)))%>%
  dplyr::filter(sample_id%in%unlist(singleton_samples))%>%
  #group_by(exp_ID)%>%
  summarise(n=n(),n0=sum(cluster_id==0),nAssigned=sum(cluster_id!=0),prop_assigned=sum(cluster_id!=0)/n())

#Visualize this assignment
singleton.assignment.plot<-dplyr::bind_rows(Map(list=mito_data,exp_ID=names(mito_data),f=function(list,exp_ID) cbind(list$clusters,exp_ID)))%>%
  dplyr::filter(sample_id%in%unlist(singleton_samples))%>%
  ggplot(aes(x=factor(exp_ID,levels=c("KX004","KX003","KX007","KX008")),y=1,fill=factor(cluster_id,levels=0:(n_clones-1))))+
  geom_bar(stat="identity",col="black",size=0.1,position="stack")+
  scale_fill_manual(values = cluster_cols,drop=F)+
  theme_bw()+
  labs(fill="Clone\nassignments",
       x="Singleton colony assignments",
       y="Count")+
  #coord_flip()+
  theme(axis.text.x=element_text(angle=90),legend.position ="right",legend.direction = "vertical",legend.key.height = unit(1.5,"mm"),legend.key.width = unit(5,"mm"))+
  my_theme
ggsave(filename=paste0(figures_dir,"Supp_Figure_08/singleton_clone_assignment_plot.pdf"),singleton.assignment.plot,width=3,height=2.5)


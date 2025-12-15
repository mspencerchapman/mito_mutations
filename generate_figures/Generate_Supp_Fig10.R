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


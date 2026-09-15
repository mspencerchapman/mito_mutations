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
for(d in c("Extended_Data_Figure_10")) dir.create(paste0(plots_dir,d),showWarnings=FALSE,recursive=TRUE)

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

#Neither ComplexHeatmap (the heatmap) nor igraph (the clustering above) is
#loaded by the package block at the top of the script.
library(ComplexHeatmap)
library(igraph)

#Define function to get clusters
#Shared-nearest-neighbour clustering of colonies by their mtDNA mutation
#profiles, following the approach of Lareau et al. as described in the methods.
#
#The parameters are those stated there, and are deliberately kept:
#  - cosine distance between colonies
#  - k.param = 5   (neighbours per colony, including itself, as in Seurat)
#  - resolution = 10 (Lareau et al. used 1, which merged clusters inappropriately)
#  - the RAW VAF matrix, not its square root
#  - Seurat's prune.SNN = 1/15 applied to the Jaccard weights
#and the mutations are filtered upstream exactly as the methods state: present at
#>1% VAF in more than one sample, and at a global VAF > 0.5%.
#
#This previously called Seurat's FindNeighbors()/FindClusters(). It no longer
#does, for two reasons. Most colonies carry no mutation above threshold, so their
#profile is the zero vector and cosine distance between them is undefined (0/0);
#under current Seurat the neighbour search degenerates to a fully connected graph
#and no structure is recovered at all. The two operations actually needed - a
#cosine kNN/SNN graph and Louvain - are therefore implemented directly, and the
#uninformative colonies are assigned explicitly to the "NULL" cluster that the
#colouring below already expects, rather than being left to undefined behaviour.
#
#IMPORTANT: the clustering is sensitive to k.param and resolution, and igraph's
#resolution is not on the same scale as Seurat's, because the two use different
#Louvain implementations. Measured against the committed assignments, on the
#colonies carrying informative mutations, this reproduces the original partition
#with an adjusted Rand index of roughly 0.5-1.0 depending on the individual, and
#no single resolution is best for all four. It is therefore an approximation, and
#the committed assignments in data/mito_mut_clones/ - the record of the published
#clustering - take precedence wherever they exist (see below).
snn_clusters <- function(mat, resolution = 1, k.param = 10) {
  rownames(mat) <- make.unique(rownames(mat))
  #A few VAFs are NA in the source matrices. An NA means the site was not
  #measured in that colony, which carries no evidence of a shared mutation, so
  #treat it as absent. Left as NA it would make the colony's norm NA and so
  #propagate NA cluster labels downstream.
  mat[is.na(mat)] <- 0
  informative <- sqrt(rowSums(mat^2)) > 0
  clusters <- rep("NULL", nrow(mat))
  names(clusters) <- rownames(mat)
  if (sum(informative) < 3) return(clusters)

  x  <- mat[informative, , drop = FALSE]
  xn <- x / sqrt(rowSums(x^2))
  cs <- xn %*% t(xn)                       #cosine similarity between colonies
  k  <- min(k.param, nrow(x) - 1)

  #k nearest neighbours of each colony
  adj <- matrix(0, nrow(x), nrow(x))
  for (i in seq_len(nrow(x))) adj[i, order(-cs[i, ])[seq_len(k + 1)]] <- 1
  diag(adj) <- 0

  #SNN weight = Jaccard overlap of the two neighbour sets, pruned as Seurat does
  inter <- adj %*% t(adj)
  deg   <- rowSums(adj)
  snn   <- inter / pmax(outer(deg, deg, "+") - inter, 1)
  snn[snn < 1/15] <- 0

  g  <- igraph::graph_from_adjacency_matrix(snn, mode = "undirected",
                                            weighted = TRUE, diag = FALSE)
  cl <- igraph::cluster_louvain(g, resolution = resolution)
  memb <- as.character(igraph::membership(cl))

  #Merge single-colony clusters into whichever cluster they are most strongly
  #connected to, as Seurat's GroupSingletons() did. The code downstream orders
  #colonies within each cluster by hierarchical clustering and only handles
  #clusters of two or more, so singletons would otherwise be dropped from the
  #heatmap and its annotation would no longer match. A singleton with no
  #connectivity at all stays put.
  repeat {
    sizes <- table(memb)
    singles <- names(sizes)[sizes == 1]
    if (!length(singles)) break
    moved <- FALSE
    for (sg in singles) {
      i <- which(memb == sg)
      if (length(i) != 1) next   #already reassigned earlier in this pass
      others <- memb != sg
      if (!any(others)) next
      conn <- tapply(snn[i, others], memb[others], sum)
      if (!length(conn) || max(conn) <= 0) next
      memb[i] <- names(conn)[which.max(conn)]
      moved <- TRUE
    }
    if (!moved) break
  }

  clusters[informative] <- memb
  clusters
}

output.dir <- paste0(root_dir,"/data/mito_mut_clones")
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
  #na.rm: without it a mutation with any NA VAF makes rowSums() NA, so an
  #all-NA row was selected into the matrix below (KX003 and KX007).
  filtered.mutations <- mutations[rowSums(vaf.mtx > 0.01, na.rm = TRUE) > 1]
  vaf.filtered.mtx <- vaf.mtx[filtered.mutations,]
  
  #-----------------------------------------------------------------------------------#
  ## DEFINE CLONES BASED ON THE MITOCHONDRIAL MUTATIONS -------------------------------
  #-----------------------------------------------------------------------------------#
  
  #The clone assignments in data/mito_mut_clones/ are the record of the clustering
  #behind the published figure, so use them when they exist: the heatmap, the
  #dendrogram of cluster means and the clone overlays are then all derived from
  #one clustering, and reproduce the published panels. snn_clusters() below is
  #the fallback for regenerating them from scratch - delete the files to do so.
  clone_assignment_file<-paste0(output.dir, "/", patient.tmp, "_mtdna_clone_assignment.txt")

  if(file.exists(clone_assignment_file)) {
    cat("  using the committed clone assignments for", patient.tmp, "\n")
    committed<-read.delim(clone_assignment_file, stringsAsFactors = FALSE)
    clusters<-str_pad(as.character(committed$cluster_id), 2, pad = "0")
    names(clusters)<-committed$sample_id
    clusters<-clusters[colnames(vaf.filtered.mtx)]
    #already the final labels, so no relabelling below

  } else {
    # get clusters with cluster resolution 10 and knn 50
    clusters <- snn_clusters(t(vaf.filtered.mtx), resolution = 10, k.param = 5)
    clusters <- str_pad(clusters, 2, pad = "0")

    #Change the label of the biggest cluster to "00" - the colonies with no
    #informative mutations - so it takes the light grey heading the palette.
    biggest_cluster<-names(table(clusters))[which.max(table(clusters))]
    original_cluster_0<-which(clusters=="00")
    new_cluster_0<-which(as.character(clusters)==biggest_cluster)
    clusters[new_cluster_0]<-"00"
    clusters[original_cluster_0]<-biggest_cluster
  }

  #NB the colours are keyed by cluster label and are assigned after the swap
  #above, not before: the matrix of cluster means further down is built from
  #names(vec_go), so names taken before the swap would refer to a label that no
  #longer exists and yield an all-NA column.
  # assign colours to clusters
  names_clusters <- unique(clusters)
  cluster_cols<- c("lightgray","#1f77b4","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf","#ad494a","#e7ba52","#8ca252","#756bb1","#636363","#aec7e8", brewer.pal(12, "Paired"))
  vec_go <- cluster_cols[1:length(names_clusters)]
  names(vec_go) <- sort(names_clusters)
  
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
  #Only write when absent. These clone assignments are tracked input data that
  #other scripts read back, and this script does not cluster correctly against
  #the clustering has been reimplemented (see snn_clusters above), so a run here
  #would replace the committed assignments - produced by the original Seurat
  #clustering - with ones from the new implementation. Delete the file to recompute.
  clone_assignment_file<-paste0(output.dir, "/", patient.tmp, "_mtdna_clone_assignment.txt")
  if(!file.exists(clone_assignment_file)) {
    write.table(df, clone_assignment_file, col.names = T, row.names = F, quote = F, sep = "\t")
  } else {
    cat("  clone assignments already exist for", patient.tmp, "- not overwriting\n")
  }
  
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
  pdf(paste0(plots_dir, "Extended_Data_Figure_10/", patient.tmp, "_mito_mutation_heatmap.pdf"), width=15, height=8)
  print(hm) 
  dev.off()
  
  #-----------------------------------------------------------------------------------#
  ## CREATE "TREES" AKA DENDROGRAMS -------------------------------
  #-----------------------------------------------------------------------------------#

  # Get group means 
  matty <- sapply(names(vec_go), function(cluster){
    cells <- df %>% dplyr::filter(cluster_id == cluster) %>% pull(sample_id) %>% as.character()
    #na.rm: a few VAFs are NA in the source matrices (9 for KX003, 2 for KX007),
    #which would otherwise make the whole cluster mean NA and the cosine
    #similarity below undefined.
    Matrix::rowMeans(sqrt(afp[,cells]), na.rm = TRUE)
  })
  
  #The "NULL" cluster holds the colonies with no mutation above threshold, so its
  #mean profile is all zero and cosine similarity against it is undefined (0/0).
  #A dendrogram of mutation profiles has no meaningful position for a cluster
  #with no mutations, so drop any all-zero column before clustering.
  matty <- matty[, colSums(abs(matty)) > 0, drop = FALSE]

  if(ncol(matty) > 2){
    
    # Do cosine distance; note that we used sqrt transformation already when creating the pseudo bulk-cell matrix
    mito.hc <- hclust(dist(lsa::cosine((matty))))
    plot(mito.hc)
    
    pdf(paste0(plots_dir,"Extended_Data_Figure_10/", patient.tmp, "_hierarchical_tree.pdf"), width = 5, height = 5)
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
  exp_clones<-read.delim(paste0(root_dir,"/data/mito_mut_clones/",exp_ID,"_mtdna_clone_assignment.txt"))
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
  pdf(file=paste0(plots_dir,"Extended_Data_Figure_10/Inferred_clones_",exp_ID,".pdf"),width = 7,height=3)
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
#The expanded clades, as built for Fig. 5 and Extended Data Fig. 9. This script
#used expanded_clades_df without defining it, so it only ever ran in a session
#where one of those scripts had already been sourced.
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
  geom_bar(stat="identity",position="stack",col="black",linewidth=0.05,width = 0.7)+
  scale_fill_manual(values = cluster_cols[1:n_clones],drop=F)+
  facet_grid(cols=vars(factor(exp_ID,levels=c("KX004","KX003","KX007","KX008"))),scales="free",space = "free")+
  theme_bw()+
  theme(axis.text.x = element_blank(),legend.key.height = unit(1.5,"mm"),legend.key.width = unit(5,"mm"))+
  my_theme+
  labs(fill="Clone\nassignments",
       x="Clonal expansions",
       y="Clone proportions")
ggsave(filename=paste0(plots_dir,"Extended_Data_Figure_10/expansion_clone_assignment_props_plot.pdf"),expansion.assignment.proportions.plot,width=7,height=2)

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
singleton_assignment_df<-dplyr::bind_rows(Map(list=mito_data,exp_ID=names(mito_data),f=function(list,exp_ID) cbind(list$clusters,exp_ID)))%>%
  dplyr::filter(sample_id%in%unlist(singleton_samples))%>%
  dplyr::mutate(n=1)

#Not every clone has singleton colonies, and drop=FALSE keeps the full clone list
#in the legend so the colours correspond to the other panels. A level with no
#rows gets a label but no key glyph, which reads as a missing colour, so give any
#absent clone a zero-height row: the legend is then complete and the bars are
#unchanged.
absent_clones<-setdiff(as.character(0:(n_clones-1)), as.character(singleton_assignment_df$cluster_id))
if(length(absent_clones)) {
  singleton_assignment_df<-dplyr::bind_rows(singleton_assignment_df,
                                            data.frame(cluster_id=as.numeric(absent_clones),
                                                       exp_ID=names(mito_data)[1], n=0))
}

singleton.assignment.plot<-singleton_assignment_df%>%
  ggplot(aes(x=factor(exp_ID,levels=c("KX004","KX003","KX007","KX008")),y=n,fill=factor(cluster_id,levels=0:(n_clones-1))))+
  geom_bar(stat="identity",col="black",linewidth=0.1,position="stack")+
  scale_fill_manual(values = cluster_cols,drop=F)+
  theme_bw()+
  labs(fill="Clone\nassignments",
       x="Singleton colony assignments",
       y="Count")+
  #coord_flip()+
  theme(axis.text.x=element_text(angle=90),legend.position ="right",legend.direction = "vertical",legend.key.height = unit(1.5,"mm"),legend.key.width = unit(5,"mm"))+
  my_theme
ggsave(filename=paste0(plots_dir,"Extended_Data_Figure_10/singleton_clone_assignment_plot.pdf"),singleton.assignment.plot,width=3,height=2.5)


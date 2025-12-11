library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpmisc)
library(gridExtra)
library(phylosignal)
options(stringsAsFactors = F)

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

#Source functions needed for the script
my_working_directory<-"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood"
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)
genomeFile="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"
root_dir="/lustre/scratch125/casm/team268im/al28/mtDNA"
plots_dir=paste0(my_working_directory,"/plots/")

#Import the mitochondrial copy number data
mito_cn=read.csv("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/whole_genome_coverage_pileup_and_bedtools_annotated.csv",header=T)

translate=data.frame(dataset=c("KY","HL","SO","PR","LM","NW","lymph","EM"),
                     al_ref=c("lung_organoid","colon","colon_ibd","muty_mutant","endometrium","blood_MPN","immune","blood_emily"))

muty_samples=c("PD44887","PD44888","PD44889","PD44890","PD44891")

#Import metadata relating to all individuals studied
ref_df<-readxl::read_excel("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/non_blood_metadata.xlsx")

#Import colony-level data for the lymphoid dataset as there are multiple different cell types
colony_info<-read.delim("colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt")

#Define a set of 'black-listed' mutations - these are sets of artefacts that recurrently slip through the filters
#This is often due to haplotype-specific artefacts that map as SNVs
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_456_C_T","MT_567_A_C","MT_574_A_C","MT_8270_C_T","MT_16170_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C")
#exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C")

#-----------------------------------------------------------------------------------#
##----------------------IMPORT ALL DATASETS ----------------------------------------
#-----------------------------------------------------------------------------------#

all_cohorts_plus_CML<-c("KY","HL","SO","PR","LM","NW","CML","lymph","blood")
all_cohorts<-c("KY","HL","SO","PR","LM","NW","lymph","blood")
all_mito_datasets<-lapply(all_cohorts_plus_CML,function(dataset) {
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
names(all_mito_datasets)<-all_cohorts_plus_CML

all_df_tidy<-Map(dataset=all_cohorts_plus_CML,dataset_mito_data=all_mito_datasets,function(dataset,dataset_mito_data) {
  cat(dataset,sep="\n")
  
  #Generate tidy data frame of samples, mutations and vafs
  mutCN_cutoff=25 #If the mitochondrial copy number is over 25, retain mutation even if is in the "CN correlating muts" list
  vaf_cut_off<-0.03
  rho_cut_off<-0
  CN_correlating_muts<-dataset_mito_data$CN_correlating_muts
  df_tidy<-dplyr::bind_rows(Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) {
    if(is.null(list)){stop(return(NULL))}
    CN_correlating_mut_removal_mat=list$matrices$CN_correlating_mut_removal_mat
    implied_mut_CN_tidy<-list$matrices$implied_mutCN%>%
      as.data.frame()%>%
      tibble::rownames_to_column(var="mut_ref")%>%
      tidyr::gather(key="Sample",value="implied_mut_CN",-mut_ref)
    df_tidy<-(list$matrices$vaf*list$matrices$SW*CN_correlating_mut_removal_mat)%>%
      as.data.frame()%>%
      tibble::rownames_to_column(var="mut_ref")%>%
      mutate(rho_val=list$rho_vals)%>%
      tidyr::gather(key="Sample",value="vaf",-mut_ref,-rho_val)
    
    comb_tidy<-left_join(df_tidy,implied_mut_CN_tidy,by=c("Sample","mut_ref"))%>%
      dplyr::filter(!grepl("DEL|INS",mut_ref) & !mut_ref%in%exclude_muts & !mut_ref%in%list$het_oocyte_muts)%>%
      dplyr::filter(vaf>=vaf_cut_off)%>%
      mutate(exp_ID=exp_ID)
    return(comb_tidy)
  }))
  
  #Filter the CN-correlating mutations - due to mis-mapping of nuclear reads.
  #There may be some genuine mutations at these sites, in which case the implied mutation copy number will be much higher than ~2 (here an arbitrary threshold of 8 is applied)
  df_tidy<-df_tidy%>%
    left_join(mito_cn,by="Sample")%>%
    mutate(implied_mut_CN=vaf*bedtools_mtDNA_genomes)%>%
    dplyr::filter(!(mut_ref%in%dataset_mito_data$CN_correlating_muts & implied_mut_CN<=mutCN_cutoff))
  
  return(df_tidy)
})

#-----------------------------------------------------------------------------------#
##----------------------ASSESS COVERAGE & UNIFORMITY ---------------------------------------
#-----------------------------------------------------------------------------------#


name_conversion_vec=c("Normal blood","Lymphoid","MPN","Colon","IBD-affected\nColon","MUTYH-mutant\nColon","Endometrium","Bronchial\nepithelium")
names(name_conversion_vec)=c("blood","lymph","NW","HL","SO","PR","LM","KY")
tissue_cols<-c("#96e97c","#fd81c8","#145a6a","#65e6f9","#781486","#5f70cc","#aea2eb","#1d6d1f")
names(tissue_cols)<-name_conversion_vec

sample_summary_info<-Map(dataset=all_cohorts,list=all_mito_datasets[all_cohorts],function(dataset,list) {
  Map(ind=list,exp_ID=names(list),function(ind,exp_ID) {
    data.frame(dataset=dataset,exp_ID=exp_ID,SampleID=ind$tree$tip.label)
  })%>%dplyr::bind_rows()
})%>%dplyr::bind_rows()

coverage.plot<-mito_cn%>%
  filter(Tissue%in%translate$al_ref)%>%
  right_join(sample_summary_info,by=c("Sample"="SampleID"))%>%
  mutate(Tissue=factor(name_conversion_vec[dataset],levels=name_conversion_vec))%>%
  filter(!is.na(bedtools_mtDNA_coverage))%>%
  ggplot(aes(x=bedtools_mtDNA_coverage,fill=Tissue))+
  geom_histogram(col="black",linewidth=0.2)+
  scale_fill_manual(values=tissue_cols)+
  scale_x_log10(labels=scales::label_comma())+
  facet_wrap(~Tissue,scales = "free_y",nrow=2)+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = margin(1,0,1,0, "mm")))+
  labs(x="Mean mitochondrial DNA coverage",y="Count")
ggsave(filename=paste0(plots_dir,"coverage_by_tissue_plot.pdf"),coverage.plot,width=7,height=2.5)

median.coverage.plot<-mito_cn%>%
  filter(Tissue%in%translate$al_ref)%>%
  right_join(sample_summary_info,by=c("Sample"="SampleID"))%>%
  mutate(Tissue=factor(name_conversion_vec[dataset],levels=name_conversion_vec))%>%
  filter(!is.na(pileup_median_mtDNA_coverage))%>%
  ggplot(aes(x=pileup_median_mtDNA_coverage,y=Tissue,fill=Tissue))+
  #geom_histogram(col="black",linewidth=0.2)+
  ggridges::geom_density_ridges(linewidth=0.2)+
  scale_fill_manual(values=tissue_cols,guide="none")+
  scale_x_log10(labels=scales::label_comma())+
  #facet_wrap(~Tissue,scales = "free_y",nrow=2)+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = margin(1,0,1,0, "mm")))+
  labs(x="Median mitochondrial DNA coverage",y="")
ggsave(filename=paste0(plots_dir,"median.coverage.plot.pdf"),median.coverage.plot,width=3.3,height=2.5)

uniformity_plot<-mito_cn%>%
  filter(Tissue%in%translate$al_ref)%>%
  right_join(sample_summary_info,by=c("Sample"="SampleID"))%>%
  mutate(Tissue=factor(name_conversion_vec[dataset],levels=name_conversion_vec))%>%
  filter(!is.na(perc_over0.8))%>%
  ggplot(aes(x=perc_over0.8,y=Tissue,fill=Tissue))+
  #geom_histogram(col="black",linewidth=0.2)+
  ggridges::geom_density_ridges(linewidth=0.2)+
  scale_fill_manual(values=tissue_cols,guide="none")+
  scale_x_log10(limits=c(0.7,1),breaks=seq(0.7,1,0.1))+
  #facet_wrap(~Tissue,scales = "free_y",nrow=2)+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = margin(1,0,1,0, "mm")))+
  labs(x="Coverage uniformity\n(proportion of mtDNA genome with coverage >80% of mean)",y="")
ggsave(filename=paste0(plots_dir,"uniformity_plot.pdf"),uniformity_plot,width=3.3,height=2.5)

#Print coverage summary statistics
mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  group_by(exp_ID)%>%
  summarise(samples=n(),mean_mtDNA_coverage=mean(bedtools_mtDNA_coverage),Under_1000X=sum(bedtools_mtDNA_coverage<1000))

mito_cn%>%
  filter(Tissue%in%translate$al_ref)%>%
  right_join(sample_summary_info,by=c("Sample"="SampleID"))%>%
  mutate(Tissue=factor(name_conversion_vec[dataset],levels=name_conversion_vec))%>%
  filter(Tissue!="Normal blood")%>%
  #group_by(Tissue)%>%
  summarise(mean_mtDNA_coverage=mean(bedtools_mtDNA_coverage),prop_over_100=sum(bedtools_mtDNA_coverage>1000)/n())


#-----------------------------------------------------------------------------------#
##----------------------ASSESS SHARED MUTATIONS ---------------------------------------
#-----------------------------------------------------------------------------------#

## Add the shared mutations dataframe----
all_mito_datasets=Map(dataset_mito_data=all_mito_datasets[all_cohorts_plus_CML],dataset=all_cohorts_plus_CML,function(dataset_mito_data,dataset) {
  cat(dataset,sep="\n")
  
  CN_correlating_muts=dataset_mito_data[[1]]$CN_correlating_muts
  figures_dir=paste0(dataset,"/")
  
  dataset_mito_data<-Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data), function(list,exp_ID) {
    cat(exp_ID,sep="\n")
    mut_vaf_cutoff=0.03
    
    vaf.filt<-list$matrices$vaf.filt
    tree.ultra<-list$tree.ultra
    
    ##This doesn't work if the trees are too small
    if(length(tree.ultra$tip.label)<=5) {
      list$shared_muts_df<-NULL
      stop(return(list))
    }
    
    shared_muts<-rownames(vaf.filt)[rowSums(vaf.filt>mut_vaf_cutoff,na.rm = T)>1]
    
    if(exp_ID=="PD37590") {
      #Specific artifacts to this individual's unusual haplotype
      excl_muts<-c("MT_955_A_C","MT_966_A_C","MT_953_T_C","MT_961_T_C")
      shared_muts<-shared_muts[!shared_muts%in%excl_muts]
    }
    n_pos<-rowSums(vaf.filt[shared_muts,]>mut_vaf_cutoff,na.rm = T)
    print(paste("There are",length(shared_muts),"shared mutations"))
    
    if(length(shared_muts)==0) {
      list$shared_muts_df<-NULL
      stop(return(list))
    }
    
    ### Test shared muts with phylosignal----
    tree4d<-phylobase::phylo4d(drop.tip(tree.ultra,"Ancestral"),tip.data=t(vaf.filt[shared_muts,list$tree$tip.label]))
    
    ### Test for phylogenetic signal (apart from the lymphoid dataset which does not have phylogeny data)
    if(!dataset=="lymph") {res_cor=phyloSignal(tree4d)} else {res_cor=list(stat=matrix(NA,nrow=length(shared_muts),ncol=2,dimnames = list(shared_muts,c("Lambda","Cmean"))),pval=matrix(NA,nrow=length(shared_muts),ncol=2,dimnames = list(shared_muts,c("Lambda","Cmean"))))}
    
    ### Add phylosignal info onto the dataframe----
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
  saveRDS(dataset_mito_data,file=paste0("mito_mutation_data_",dataset,".RDS"))
  return(dataset_mito_data)
})

## Do the plotting ----
temp=Map(dataset_mito_data=all_mito_datasets[all_cohorts_plus_CML],dataset=all_cohorts_plus_CML,function(dataset_mito_data,dataset) {
  
  cat(dataset,sep="\n")
  CN_correlating_muts=dataset_mito_data[[1]]$CN_correlating_muts
  figures_dir=paste0(dataset,"/")
  
  ### Visualize the shared mutations, showing only VAFs that pass shearwater ----
  col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
  lowest_VAF_to_show=0.03
  tree_to_show="tree.ultra"
  temp=Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID){
    cat(exp_ID,sep="\n"); Age<-ref_df%>%filter(ID==exp_ID)%>%pull(Age)
    if(is.null(list)) {stop(return(NULL))}
    
    #Don't bother printing if the tree is <10 samples
    if(length(list[[tree_to_show]]$tip.label)<8) {stop(return(NULL))}
    
    if(is.null(list$shared_muts_df)) {
      pdf(file = paste0(figures_dir,exp_ID,"_shared_muts.pdf"),width = 7,height=4)
      plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=F,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    if(nrow(list$shared_muts_df)==0) {
      pdf(file = paste0(figures_dir,exp_ID,"_shared_muts.pdf"),width = 7,height=4)
      plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=F,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
    
    vaf.mtx<-list$matrices$vaf.filt*list$matrices$SW[rownames(list$matrices$vaf.filt),]
    
    plot_muts<-list$shared_muts_df%>%
      filter(mut%in%rownames(vaf.mtx))%>%
      pull(mut)
    
    #Remove samples if not included at all in the shearwater table - these samples are excluded due to contamination
    list[[tree_to_show]]<-keep.tip(list[[tree_to_show]],list[[tree_to_show]]$tip.label[list[[tree_to_show]]$tip.label%in%list$sample_shearwater_calls$sampleID])
    list[[tree_to_show]]$coords<-NULL
    vaf.mtx<-vaf.mtx[,list[[tree_to_show]]$tip.label]
    
    names(col_scheme)<-seq(0,1,0.01)
    hm<-matrix(0,nrow=length(plot_muts),ncol=length(list[[tree_to_show]]$tip.label),dimnames = list(plot_muts,list[[tree_to_show]]$tip.label))
    for(i in 1:length(plot_muts)) {
      mut<-plot_muts[i]
      mut_vafs<-vaf.mtx[mut,list[[tree_to_show]]$tip.label]
      mut_vafs[mut_vafs<lowest_VAF_to_show]<-0
      mut_vafs<-round(mut_vafs,digits=2)
      hm[i,]<-col_scheme[as.character(mut_vafs)]
    }
    if(length(plot_muts)>=2) {
      plot_muts.clustered<-hclust(dist(vaf.mtx[plot_muts,]))
    } else {
      plot_muts.clustered=list(order=1)
    }
    
    par(mfrow=c(1,1))
    pdf(file = paste0(figures_dir,exp_ID,"_shared_muts.pdf"),width = 7,height=4)
    list[[tree_to_show]]=plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
    add_mito_mut_heatmap(tree=list[[tree_to_show]],heatmap=hm[plot_muts.clustered$order,,drop=F],border="gray",heatmap_bar_height=0.1,cex.label = 0.25)
    dev.off()
  })
  
  ### Visualize the shared mutations again - now in 'raw' form, without shearwater filtering ----
  tree_to_show="tree"
  col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
  lowest_VAF_to_show=0.03
  temp=Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID){
    cat(exp_ID,sep="\n");Age<-ref_df%>%filter(ID==exp_ID)%>%pull(Age)
    if(is.null(list)) {stop(return(NULL))}
    
    #Don't bother printing if the tree is <8 samples
    if(length(list$tree$tip.label)<8) {stop(return(NULL))}
    
    #Plot the trees even if there are no shared muts
    if(is.null(list$shared_muts_df)) {
      pdf(file = paste0(figures_dir,exp_ID,"_raw_shared_muts.pdf"),width = 7,height=4)
      list[[tree_to_show]]=plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    if(nrow(list$shared_muts_df)==0) {
      pdf(file = paste0(figures_dir,exp_ID,"_raw_shared_muts.pdf"),width = 7,height=4)
      list[[tree_to_show]]=plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
    plot_muts<-list$shared_muts_df%>%
      #filter(Cmean_pval<0.05)%>%
      pull(mut)
    
    if(length(plot_muts)==0) {
      pdf(file = paste0(figures_dir,exp_ID,"_raw_shared_muts.pdf"),width = 7,height=4)
      list[[tree_to_show]]=plot_tree(tree = list[[tree_to_show]],cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    
    vaf.mtx<-list$matrices$vaf[rownames(list$matrices$vaf.filt),]
    
    names(col_scheme)<-seq(0,1,0.01)
    hm<-matrix(0,nrow=length(plot_muts),ncol=length(list[[tree_to_show]]$tip.label),dimnames = list(plot_muts,list[[tree_to_show]]$tip.label))
    for(i in 1:length(plot_muts)) {
      mut<-plot_muts[i]
      mut_vafs<-vaf.mtx[mut,list[[tree_to_show]]$tip.label]
      mut_vafs[mut_vafs<lowest_VAF_to_show]<-0
      mut_vafs<-round(mut_vafs,digits=2)
      hm[i,]<-col_scheme[as.character(mut_vafs)]
    }
    if(length(plot_muts)>=2) {
      plot_muts.clustered<-hclust(dist((vaf.mtx^(1/4))[plot_muts,]))
    } else {
      plot_muts.clustered=list(order=1)
    }
    
    par(mfrow=c(1,1))
    pdf(file = paste0(figures_dir,exp_ID,"_raw_shared_muts.pdf"),width = 7,height=4)
    list[[tree_to_show]]=plot_tree(tree = list[[tree_to_show]],cex.label = 0,cex.axis=0.5,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : Shared mutations (VAF>",100*lowest_VAF_to_show,"%)"))
    add_mito_mut_heatmap(tree=list[[tree_to_show]],heatmap=hm[plot_muts.clustered$order,,drop=F],border="gray",heatmap_bar_height=0.1,cex.label = 0.25)
    dev.off()
  })
  
  ### Plot all mutations that pass shearwater in at least one sample: BIG HEATMAPS ----
  #Higher granularity visualization - show down to 0.001 VAF
  lowest_VAF_to_show=0.001
  col_scheme_high<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(1001))
  names(col_scheme_high)<-seq(0,1,0.001)
  temp=Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID){
    cat(exp_ID,sep="\n"); Age<-ref_df%>%filter(ID==exp_ID)%>%pull(Age)
    if(is.null(list)) {stop(return(NULL))}
    bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
    plot_muts<-rownames(list$matrices$vaf.filt)
    
    if(length(plot_muts)==0) {
      pdf(file = paste0(figures_dir,exp_ID,"_raw_shared_muts.pdf"),width = 7,height=4)
      list$tree=plot_tree(tree = list$tree,cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : All mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    
    #Showing the raw VAF here - not filtered by SW
    vaf.mtx<-list$matrices$vaf.filt[rownames(list$matrices$vaf.filt),]
    
    hm<-matrix(0,nrow=length(plot_muts),ncol=length(list$tree$tip.label),dimnames = list(plot_muts,list$tree$tip.label))
    for(i in 1:length(plot_muts)) {
      mut<-plot_muts[i]
      mut_vafs<-vaf.mtx[mut,list$tree$tip.label]
      mut_vafs[mut_vafs<lowest_VAF_to_show]<-0
      mut_vafs<-round(mut_vafs,digits=3)
      hm[i,]<-col_scheme_high[as.character(mut_vafs)]
    }
    if(length(plot_muts)>=2) {
      plot_muts.clustered<-hclust(dist((vaf.mtx^(1/4))[plot_muts,]))
    } else {
      plot_muts.clustered=list(order=1)
    }
    
    par(mfrow=c(1,1))
    pdf(file = paste0(figures_dir,exp_ID,"_all_muts_filt.pdf"),width = 7,height=8)
    list$tree=plot_tree(tree = list$tree,cex.label = 0,plot_axis=T,vspace.reserve = 4,title=paste0(exp_ID," (",Age," years old) : All mutations (VAF>",100*lowest_VAF_to_show,"%)"))
    add_mito_mut_heatmap(tree=list$tree,heatmap=hm[plot_muts.clustered$order,,drop=F],border="gray",heatmap_bar_height=0.01,cex.label = 0.25)
    dev.off()
  })
  
  ### Plot suspected germline mutations and check them ----
  #Higher granularity visualization - show down to 0.001 VAF
  col_scheme<-c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")[2:9])(100))
  lowest_VAF_to_show=0.03
  tree_to_show="tree.ultra"
  temp=Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID){
    cat(exp_ID,sep="\n"); Age<-ref_df%>%filter(ID==exp_ID)%>%pull(Age)
    if(is.null(list)) {stop(return(NULL))}
    bb_df<-data.frame(mut=rownames(list$matrices$vaf),rho=list$rho_vals)
    plot_muts<-grep(pattern = "\\*",rownames(list$matrices$vaf),value=T)
    
    if(length(plot_muts)==0) {
      pdf(file = paste0(figures_dir,exp_ID,"_germline_muts.pdf"),width = 7,height=4)
      list$tree=plot_tree(tree = list$tree,cex.label = 0,plot_axis=T,vspace.reserve = 1.1,title=paste0(exp_ID," (",Age," years old) : germline mutations (VAF>",100*lowest_VAF_to_show,"%)"))
      dev.off()
      stop(return(NULL))
    }
    
    #Showing the raw VAF here - not filtered by SW
    vaf.mtx<-list$matrices$vaf
    
    hm<-matrix(0,nrow=length(plot_muts),ncol=length(list$tree$tip.label),dimnames = list(plot_muts,list$tree$tip.label))
    for(i in 1:length(plot_muts)) {
      mut<-plot_muts[i]
      mut_vafs<-vaf.mtx[mut,list$tree$tip.label]
      mut_vafs[mut_vafs<lowest_VAF_to_show]<-0
      mut_vafs<-round(mut_vafs,digits=3)
      hm[i,]<-col_scheme_high[as.character(mut_vafs)]
    }
    
    if(length(plot_muts)>=2) {plot_muts.clustered<-hclust(dist((vaf.mtx^(1/4))[plot_muts,]))} else {plot_muts.clustered=list(order=1)}
    
    par(mfrow=c(1,1))
    pdf(file = paste0(figures_dir,exp_ID,"_germline_muts.pdf"),width = 7,height=8)
    list$tree=plot_tree(tree = list$tree,cex.label = 0,plot_axis=T,vspace.reserve = 4,title=paste0(exp_ID," (",Age," years old) : germline mutations (VAF>",100*lowest_VAF_to_show,"%)"))
    add_mito_mut_heatmap(tree=list$tree,heatmap=hm[plot_muts.clustered$order,,drop=F],border="gray",heatmap_bar_height=0.01,cex.label = 0.25)
    dev.off()
  })
})

# --------MEASURE 'CLONAL MARKING' OF EXPANDED CLADES BY MITOCHONDRIAL MUTATIONS--------------------
dataset="CML"
dataset_mito_data<-all_mito_datasets[[dataset]]

### Basic plot to show a simplified clonal structure ----
expanded.clades.plot<-dplyr::bind_rows(Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID){
  exp_nodes<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = ifelse(dataset=="NW",50,ifelse(dataset=="CML",2,200)),min_clonal_fraction=0.01,min_samples = 2)
  if(nrow(exp_nodes)==0) {stop(return(NULL))}
  exp_nodes$exp_ID<-exp_ID
  return(exp_nodes)
}))%>%
  arrange(desc(n_samples))%>%
  ggplot(aes(x=exp_ID,y=clonal_fraction,fill=n_samples))+
  geom_bar(stat="identity",position="stack",col="black")+
  scale_fill_gradient(low="lightgrey",high="darkred")+
  labs(x="Individual",y="Clonal fraction",fill="Number of\n samples\n in clone")+
  theme_bw()+
  my_theme+
  theme(axis.text.x=element_text(angle=90))
ggsave(filename=paste0(plots_dir,dataset,"_expanded_clades_plot.pdf"),expanded.clades.plot,width=4,height=2.5)

### Get info on 'marker mutations' ----
#Slightly complex snippet of code to pull out any mutations that are:
# (1) shared by ≥2 clone samples
# (2) are not present in a large fraction of samples outside the clone (must be <20%)
expanded_clades_df<-Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID){
  cat(paste0(exp_ID,"\n"))
  
  vaf.filt<-list$matrices$vaf.filt
  
  #Review how many expanded clades have reliable mitochondrial marker mutations
  marker_mut_cutoff<-0.03
  pos_mut_cutoff<-0.03
  exp_nodes<-get_expanded_clade_nodes(list$tree,height_cut_off = ifelse(dataset=="NW",50,200),min_clonal_fraction=0.01,min_samples = 2)
  if(nrow(exp_nodes)==0|is.null(list$shared_muts_df)) {stop(return(NULL))}
  full_df<-dplyr::bind_cols(data.frame(exp_ID=rep(exp_ID,nrow(exp_nodes))),
                            exp_nodes,
                            data.frame(marker_mut_cutoff=rep(marker_mut_cutoff,nrow(exp_nodes))),
                            exp_nodes_muts<-dplyr::bind_rows(lapply(exp_nodes$nodes,function(node) {
                              cat(node,sep="\n")
                              node_samples=getTips(list$tree,node)
                              if(any(vaf.filt[,node_samples]>marker_mut_cutoff)){
                                
                                #Get mutations that are present in at least one sample in the clade above the 'marker_mut' VAF cutoff
                                node_homo_muts<-names(rowSums(vaf.filt[,node_samples,drop=F]>marker_mut_cutoff,na.rm=T)[rowSums(vaf.filt[,node_samples,drop=F]>marker_mut_cutoff)>0])
                                
                                #Test samples outside the node - exclude if mutation is also present in a significant proportion of samples outside the node in question
                                #This will primarily exclude mutations that are high level oocyte heteroplasmies that are not informative
                                outside_node_samples<-list$tree$tip.label[!list$tree$tip.label%in%node_samples]
                                pos_outside_node_samples<-rowSums(vaf.filt[node_homo_muts,outside_node_samples,drop=F]>marker_mut_cutoff)
                                prop_pos_outside_node_samples<-pos_outside_node_samples/length(outside_node_samples)
                                node_homo_muts<-node_homo_muts[prop_pos_outside_node_samples<0.2]
                                
                                if(length(node_homo_muts)==0) {
                                  stop(return(data.frame(nmuts=0,homo_muts=NA,pos_samples_per_mut=NA,max_pos_samples=0,max_pos_prop=0,mean_heteroplasmy=NA)))
                                }
                                
                                pos_samples_per_mut<-sapply(node_homo_muts,function(mut){
                                  n_samples<-sum(vaf.filt[mut,node_samples,drop=F]>pos_mut_cutoff)
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
                                                  max_pos_samples=max(pos_samples_per_mut,na.rm=T),
                                                  max_pos_prop=max(pos_samples_per_mut,na.rm=T)/length(node_samples),
                                                  mean_heteroplasmy=mean(as.numeric(vaf.filt[node_homo_muts[which.max(pos_samples_per_mut)],node_samples])),
                                                  mean_het_of_pos=mean_het_of_pos[which.max(pos_samples_per_mut)]))
                              } else {
                                return(data.frame(nmuts=0,homo_muts=NA,pos_samples_per_mut=NA,max_pos_samples=0,max_pos_prop=0,mean_heteroplasmy=NA))
                              }
                            })))
  #return(dplyr::select(full_df,-homo_muts,-pos_samples_per_mut))
  return(full_df)
})%>%dplyr::bind_rows()

node_factor_levels=expanded_clades_df%>%arrange(n_samples)%>%mutate(levels=stringr::str_c(exp_ID,nodes,sep = "_"))%>%pull(levels)
expanded.clades.marking.plot<-expanded_clades_df%>%
  mutate(n_neg_samples=n_samples-max_pos_samples)%>%
  dplyr::select(exp_ID,nodes,n_samples,max_pos_samples,n_neg_samples)%>%
  mutate(max_pos_samples=ifelse(max_pos_samples<=1,0,max_pos_samples))%>%
  mutate(n_neg_samples=n_samples-max_pos_samples)%>%
  tidyr::gather(-exp_ID,-nodes,-n_samples,key="Pos_or_neg",value="n_samples")%>%
  mutate(Pos_or_neg=ifelse(Pos_or_neg=="n_neg_samples","Absent","Present"))%>%
  mutate(Pos_or_neg=factor(Pos_or_neg,levels=c("Present","Absent")))%>%
  mutate(levels=str_c(exp_ID,nodes,sep = "_"))%>%
  ggplot(aes(x=factor(levels,levels=node_factor_levels),y=n_samples,fill=Pos_or_neg))+
  geom_bar(position="stack",stat="identity",col="black",linewidth=0.15)+
  scale_fill_brewer(palette="Set2")+
  facet_grid(~exp_ID,scales="free",space = "free")+
  #scale_y_continuous(breaks=seq(0,10,2))+
  theme_bw()+
  my_theme+
  theme(panel.spacing=unit(0.5,"mm"),axis.text.x = element_blank(),strip.text.x = element_text(angle=90))+
  labs(x="Clonal expansion",y=str_wrap("Number of samples within expansion",width=15),fill=str_wrap("Best mitochondrial marker mutation",width=10))
ggsave(filename=paste0("expanded_clades_marking_plot_",dataset,".pdf"),expanded.clades.marking.plot,width=7,height=2.3)

expanded_clades_df$n_marker<-sapply(1:nrow(expanded_clades_df),function(i) {
  pos_samples_per_mut<-as.numeric(stringr::str_split(expanded_clades_df$pos_samples_per_mut[i],pattern = ",")[[1]])
  n_perfect_markers<-sum(pos_samples_per_mut==expanded_clades_df$n_samples[i])
  return(n_perfect_markers)
})


#-----------------------------------------------------------------------------------#
##----------------------COPY NUMBER ANALYSIS ----------------------------------------
#-----------------------------------------------------------------------------------#


cn_by_tissue_plot<-mito_cn%>%
  filter(Tissue%in%translate$al_ref)%>%
  right_join(sample_summary_info,by=c("Sample"="SampleID"))%>%
  mutate(Tissue=factor(name_conversion_vec[dataset],levels=name_conversion_vec))%>%
  filter(!is.na(bedtools_mtDNA_genomes))%>%
  ggplot(aes(x=factor(exp_ID,levels=ref_df%>%arrange(Age)%>%pull(ID)%>%unique()),y=bedtools_mtDNA_genomes,col=Tissue))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1,alpha=0.5,size=0.2)+
  scale_color_manual(values=tissue_cols)+
  scale_y_log10()+
  facet_grid(~Tissue,scales="free",space="free")+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90,size = 4),legend.position="none")+
  labs(x="",y="mtDNA copy number")

ggsave(filename=paste0(plots_dir,"cn_by_tissue_plot.pdf"),cn_by_tissue_plot,width=10,height=4)

mitoCN_ridges_plot<-mito_cn%>%
  filter(Tissue%in%translate$al_ref)%>%
  right_join(sample_summary_info,by=c("Sample"="SampleID"))%>%
  mutate(Tissue=factor(name_conversion_vec[dataset],levels=name_conversion_vec))%>%
  filter(!is.na(bedtools_mtDNA_genomes))%>%
  mutate(log_mitoCN=log(bedtools_mtDNA_genomes))%>%
  ggplot(aes(y=Tissue,x=bedtools_mtDNA_genomes,fill=Tissue))+
  ggridges::geom_density_ridges(linewidth=0.3)+
  scale_fill_manual(values=tissue_cols)+
  scale_x_log10()+
  theme_classic()+
  my_theme+
  theme(legend.position="none",axis.title.y=element_blank())+
  labs(x="mtDNA copy number")

ggsave(filename=paste0(plots_dir,"mitoCN_ridges_plot.pdf"),mitoCN_ridges_plot,width=3,height=2)

median_cn_by_age_plot<-mito_cn%>%
  filter(Tissue%in%translate$al_ref)%>%
  right_join(sample_summary_info,by=c("Sample"="SampleID"))%>%
  mutate(Tissue=factor(name_conversion_vec[dataset],levels=name_conversion_vec))%>%
  filter(!is.na(bedtools_mtDNA_genomes))%>%
  group_by(Tissue,exp_ID)%>%
  dplyr::summarise(Tissue=Tissue[1],median_cn=median(bedtools_mtDNA_genomes))%>%
  left_join(ref_df,by=c("exp_ID"="ID"),relationship="many-to-many")%>%
  ggplot(aes(x=Age,y=median_cn,col=Tissue))+
  geom_point(size=0.5)+
  scale_color_manual(values=tissue_cols)+
  geom_smooth(method="lm",linewidth=0.75,alpha=0.3)+
  scale_y_log10()+
  theme_classic()+
  my_theme+
  theme(legend.position="right")+
  labs(x="Age",y="mtDNA copy number")

ggsave(filename=paste0(plots_dir,"median_cn_by_age_plot.pdf"),median_cn_by_age_plot,width=3.3,height=2)



mean_cn_plot<-mito_cn%>%
  filter(Tissue%in%translate$al_ref)%>%
  right_join(sample_summary_info,by=c("Sample"="SampleID"))%>%
  mutate(Tissue=factor(name_conversion_vec[dataset],levels=name_conversion_vec))%>%
  filter(!is.na(bedtools_mtDNA_genomes))%>%
  group_by(Tissue)%>%
  dplyr::summarise(n=n(),
                   mean_cn=mean(bedtools_mtDNA_genomes),
                   var_cn=var(bedtools_mtDNA_genomes),
                   sd_cn=sd(bedtools_mtDNA_genomes),
                   mean_logCN=mean(log(bedtools_mtDNA_genomes)),
                   var_logCN=var(log(bedtools_mtDNA_genomes)),
                   sd_logCN=sd(log(bedtools_mtDNA_genomes)))%>%
  dplyr::mutate(sem_logCN=sd_logCN/sqrt(n))%>%
  dplyr::mutate(lowerCI_mean_logCN=mean_logCN-(1.96*sem_logCN),upperCI_mean_logCN=mean_logCN+(1.96*sem_logCN))%>%
  dplyr::mutate(mean_cn2=exp(mean_logCN),lowerCI_mean_CN=exp(lowerCI_mean_logCN),upperCI_mean_CN=exp(upperCI_mean_logCN))%>%
  dplyr::select(Tissue,mean_cn,mean_cn2,lowerCI_mean_CN,upperCI_mean_CN)%>%
  ggplot(aes(y=Tissue,x=mean_cn2,xmin=lowerCI_mean_CN,xmax=upperCI_mean_CN,col=Tissue))+
  geom_point(size=0.4)+
  geom_errorbar(width=0.2)+
  scale_x_continuous(limits=c(0,2200))+
  scale_color_manual(values=tissue_cols)+
  theme_classic()+
  my_theme+
  labs(x="Mean tissue mtDNA copy number")+
  theme(legend.position="none")

ggsave(filename=paste0(plots_dir,"mean_cn_plot.pdf"),mean_cn_plot,width=3,height=2)

mito_cn_lmer_df<-mito_cn%>%
  filter(Tissue%in%translate$al_ref)%>%
  right_join(sample_summary_info,by=c("Sample"="SampleID"))%>%
  mutate(Tissue=factor(name_conversion_vec[dataset],levels=name_conversion_vec))%>%
  filter(!is.na(bedtools_mtDNA_genomes))%>%
  mutate(log_mitoCN=log(bedtools_mtDNA_genomes))%>%
  left_join(ref_df%>%dplyr::select(ID,Age,Smoking_years)%>%filter(!duplicated(.)),by=c("exp_ID"="ID"))

mito_cn.lmer_with_age<-lme4::lmer(log_mitoCN~Age+Tissue+(1|exp_ID),data=mito_cn_lmer_df)
summary(mito_cn.lmer_with_age)
confint(mito_cn.lmer_with_age)

mito_cn.lmer<-lme4::lmer(log_mitoCN~Tissue+(1|exp_ID),data=mito_cn_lmer_df)
summary(mito_cn.lmer)
lmer.CIs<-confint(mito_cn.lmer)

mito_cn.lmer_lung_smoking<-lme4::lmer(log_mitoCN~Smoking_years+(1|exp_ID),data=mito_cn_lmer_df%>%filter(dataset=="KY")%>%mutate(Smoking_years=as.numeric(Smoking_years))%>%replace_na(replace = list(Smoking_years=0)))
summary(mito_cn.lmer_lung_smoking)
lmer.CIs<-confint(mito_cn.lmer_lung_smoking)

#Visualize the model output tissue-specific coefficients
intercept_value<-summary(mito_cn.lmer)$coefficients["(Intercept)","Estimate"]
tissue_lmer_coefs_plot<-as.data.frame(summary(mito_cn.lmer)$coefficients)%>%
  tibble::rownames_to_column(var="Tissue")%>%
  left_join(as.data.frame(lmer.CIs)%>%tibble::rownames_to_column(var="Tissue"))%>%
  dplyr::select(Tissue,Estimate,lowerCI=`2.5 %`,upperCI=`97.5 %`)%>%
  mutate(Estimate=ifelse(grepl("Tissue",Tissue),Estimate+intercept_value,Estimate),
         lowerCI=ifelse(grepl("Tissue",Tissue),lowerCI+intercept_value,lowerCI),
         upperCI=ifelse(grepl("Tissue",Tissue),upperCI+intercept_value,upperCI))%>%
  filter(Tissue!="Age")%>%
  mutate(Tissue=ifelse(Tissue=="(Intercept)","Blood",gsub("Tissue","",Tissue)))%>%
  mutate(Tissue=factor(Tissue,levels=name_conversion_vec))%>%
  #mutate_at(c("Estimate","lowerCI","upperCI"),function(x) {exp(x)})%>%
  ggplot(aes(y=Tissue,x=Estimate,xmin=lowerCI,xmax=upperCI,col=Tissue))+
  geom_point(size=0.4)+
  geom_errorbar(width=0.2)+
  scale_x_continuous(limits=c(4,8))+
  #scale_x_continuous(limits=c(0,2300))+
  scale_color_manual(values=tissue_cols)+
  theme_classic()+
  my_theme+
  labs(x="Tissue-specific log(mtDNA copy number)\nfixed effect coefficients")+
  theme(legend.position="none")

ggsave(filename=paste0(plots_dir,"tissue_lmer_coefs_plot.pdf"),tissue_lmer_coefs_plot,width=3,height=2)

#Visualize the model output tissue-specific coefficients
intercept_value<-summary(mito_cn.lmer)$coefficients["(Intercept)","Estimate"]
tissue_lmer_exp_coefs_plot<-as.data.frame(summary(mito_cn.lmer)$coefficients)%>%
  tibble::rownames_to_column(var="Tissue")%>%
  left_join(as.data.frame(lmer.CIs)%>%tibble::rownames_to_column(var="Tissue"))%>%
  dplyr::select(Tissue,Estimate,lowerCI=`2.5 %`,upperCI=`97.5 %`)%>%
  mutate(Estimate=ifelse(grepl("Tissue",Tissue),Estimate+intercept_value,Estimate),
         lowerCI=ifelse(grepl("Tissue",Tissue),lowerCI+intercept_value,lowerCI),
         upperCI=ifelse(grepl("Tissue",Tissue),upperCI+intercept_value,upperCI))%>%
  filter(Tissue!="Age")%>%
  mutate(Tissue=ifelse(Tissue=="(Intercept)","Blood",gsub("Tissue","",Tissue)))%>%
  mutate(Tissue=factor(Tissue,levels=name_conversion_vec))%>%
  mutate_at(c("Estimate","lowerCI","upperCI"),function(x) {exp(x)})%>%
  ggplot(aes(y=Tissue,x=Estimate,xmin=lowerCI,xmax=upperCI,col=Tissue))+
  geom_point(size=0.4)+
  geom_errorbar(width=0.2)+
  scale_x_continuous(limits=c(0,2700))+
  scale_color_manual(values=tissue_cols)+
  theme_classic()+
  my_theme+
  labs(x="Tissue-specific mtDNA copy number\nfixed effect coefficients")+
  theme(legend.position="none")

ggsave(filename=paste0(plots_dir,"tissue_lmer_exp_coefs_plot.pdf"),tissue_lmer_exp_coefs_plot,width=3,height=2)

#-----------------------------------------------------------------------------------#
# ----------------------MUTATIONAL PROFILES----------------------
#-----------------------------------------------------------------------------------#

# read in the trinucleotide context
mtdna_trinuc_freq <- readRDS("/lustre/scratch126/casm/team154pc/ms56/mito_mutations_blood/data/mtDNA_trinuc_freqs_coding_dloop_heavy_light.Rds")

# set regions for coding and dloop regions
coding_region <- 577:16023
d_loop_region <- c(1:576,16024:16569)

trinucleotide_plot(mutations = all_df_tidy[all_cohorts]%>%
                     dplyr::bind_rows()%>%
                     filter(Sample!="global"&
                              !mut_ref%in%exclude_muts)%>%
                     separate(mut_ref,into=c("chr","pos","ref","mut"))%>%
                     dplyr::rename("donor"=exp_ID)%>%
                     mutate(pos=as.numeric(pos)),
                   analysis_region = "all_mtDNA",
                   analysis_type = "obs_exp",
                   file_name = paste0(plots_dir,"Mutational_signature_nonblood_all.pdf"))

trinucleotide_plot(mutations = all_df_tidy[all_cohorts]%>%
                     dplyr::bind_rows()%>%
                     filter(Sample!="global"&
                              !mut_ref%in%exclude_muts&
                              vaf<0.1)%>%
                     separate(mut_ref,into=c("chr","pos","ref","mut"))%>%
                     dplyr::rename("donor"=exp_ID)%>%
                     mutate(pos=as.numeric(pos)),
                   analysis_region = "all_mtDNA",
                   analysis_type = "obs_exp",
                   file_name = paste0(plots_dir,"Mutational_signature_non_blood_under0.1.pdf"))

trinucleotide_plot(mutations = df_tidy%>%
                     filter(Sample!="global"&
                              ML_Sig=="N1"&
                              !mut_ref%in%exclude_muts&
                              vaf<0.005)%>%
                     separate(mut_ref,into=c("chr","pos","ref","mut"))%>%
                     dplyr::rename("donor"=exp_ID)%>%
                     mutate(pos=as.numeric(pos)),
                   analysis_region = "all_mtDNA",
                   analysis_type = "obs_exp",
                   file_name = paste0(plots_dir,"Mutational_signature_<0.005.pdf"))


#-----------------------------------------------------------------------------------#
##----------------------Print stats for total mutation numbers----------------------
#-----------------------------------------------------------------------------------#
mutCN_cutoff=25 #If the mitochondrial copy number is over 25, retain mutation even if is in the "CN correlating muts" list
vaf_cut_off<-0.03
rho_cut_off<-0
temp=Map(dataset_mito_data=all_mito_datasets[all_cohorts],dataset=all_cohorts, function(dataset_mito_data,dataset) {
  
  #Total number of SNVs identified by shearwater - this includes the heterozygous oocyte mutations
  SW_total<-sum(unlist(lapply(dataset_mito_data,function(list) {
    if(is.null(list)){stop(return(NULL))}
    vaf.filt<-(list$matrices$vaf*list$matrices$SW)[!grepl("DEL|INS",rownames(list$matrices$vaf))&
                                                     list$rho_vals>rho_cut_off&
                                                     !rownames(list$matrices$vaf)%in%exclude_muts,list$tree$tip.label]
    total_pos=rowSums(vaf.filt>=vaf_cut_off,na.rm = T)
    return(sum(total_pos>0,na.rm=T))
  })))
  
  #Total number of SNVs identified by shearwater - once CN correlating muts excluded
  SW_total_CN_excluded<-sum(unlist(lapply(dataset_mito_data,function(list) {
    if(is.null(list)){stop(return(NULL))}
    CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%dataset_mito_data$CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
    vaf.filt<-(list$matrices$vaf*list$matrices$SW*CN_correlating_mut_removal_mat)[!grepl("DEL|INS",rownames(list$matrices$vaf))&
                                                                                    list$rho_vals>rho_cut_off&
                                                                                    !rownames(list$matrices$vaf)%in%exclude_muts,list$tree$tip.label]
    total_pos=rowSums(vaf.filt>=vaf_cut_off,na.rm = T)
    return(sum(total_pos>0,na.rm=T))
  })))
  
  cat(paste("Print stats for",dataset),sep="\n")
  cat(paste("There are",SW_total,"that pass shearwater in at least one sample"),sep="\n")
  cat(paste("There are",SW_total_CN_excluded,"that pass shearwater in at least one sample, after copy-number correlating mutations are excluded"),sep="\n")
  
})

#-----------------------------------------------------------------------------------#
##------Save 'prop of samples with mut' plots for each dataset individually----------
#-----------------------------------------------------------------------------------#

#Compile list of 'passing' samples and donors
samples_per_individual=Map(dataset_mito_data=all_mito_datasets[all_cohorts],dataset=all_cohorts, function(dataset_mito_data,dataset) {
  Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID){data.frame(exp_ID=exp_ID,Sample=list$tree$tip.label)})%>%
    dplyr::bind_rows()%>%
    dplyr::mutate(dataset=dataset)
})%>%dplyr::bind_rows()%>%group_by(dataset,exp_ID)%>%summarise(dataset=dataset[1],n=n())


#Do not do this for individuals with <8 samples
temp=Map(dataset_mito_data=all_mito_datasets,df_tidy=all_df_tidy,this_dataset=names(all_mito_datasets), function(dataset_mito_data,df_tidy,this_dataset) {
  
  min_samples_per_individual=8
  individuals_to_include<-samples_per_individual%>%filter(n>=min_samples_per_individual & dataset==this_dataset)%>%pull(exp_ID)
  exp_cols<-c("#aee39a", "#18786a", "#99def9", "#0362a0", "#e689eb", "#53348e", "#9296ee", "#9e37d0", "#32e195", "#50942f", "#a3d71e", "#7c440e", "#ecbcab", "#a53460", "#fb899b", "#ed4b04", "#ffce54", "#f82387", "#61f22d", "#d07d09","#E41A1C")
  names(exp_cols)<-c("Dummy_foetal",individuals_to_include)
  
  Individual_burden_means_df<-df_tidy%>%
    dplyr::filter(!mut_ref%in%exclude_muts)%>%
    group_by(Sample)%>%
    dplyr::summarise(exp_ID=unique(exp_ID),n_mut=sum(vaf>vaf_cut_off))%>%
    tidyr::complete(Sample,fill=list(n_mut=0))%>%
    group_by(exp_ID)%>%
    summarise(mean=mean(n_mut))
  
  samples_with_mut_df<-dplyr::bind_rows(lapply(seq(vaf_cut_off,0.99,0.01),function(cut_off) {
    cut_off_df<-df_tidy%>%
      group_by(Sample)%>%
      dplyr::summarise(exp_ID=exp_ID[1],n_mut=sum(vaf>cut_off))%>%
      tidyr::complete(Sample,fill=list(n_mut=0))%>%
      #left_join(dataset_sample_ref,by="Sample")%>%
      group_by(exp_ID)%>%
      dplyr::summarise(mean=mean(n_mut),n_with_mut=sum(n_mut>0),n_samp=n())%>%
      mutate(prop_with_mut=n_with_mut/n_samp,cut_off=cut_off)
    return(cut_off_df)
  }))
  
  prop.of.samples.with.mut<-samples_with_mut_df%>%
    filter(exp_ID%in%individuals_to_include)%>% #Only include individuals with at least 8 samples
    ggplot(aes(x=cut_off,y=prop_with_mut,col=factor(exp_ID)))+
    geom_line(alpha=0.6)+
    theme_bw()+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
    scale_color_manual(values=exp_cols)+
    guides(colour = guide_legend(ncol=1))+
    labs(x="VAF cut-off threshold",
         y=str_wrap("Proportion of samples with at least 1 mutation with VAF > cut-off",width=40),
         col="Individual")+
    my_theme+
    theme(legend.key.height = unit(2.8,"mm"),legend.title = element_text(size=7))
  
  ggsave(filename=paste0(plots_dir,"Sample_proportions_with_mut_",this_dataset,".pdf"),prop.of.samples.with.mut,width=3.5,height=2)
  
})

#-----------------------------------------------------------------------------------#
##----------------------Get sum of vaf data for all datasets ------------------------
#-----------------------------------------------------------------------------------#

#Sum of VAF mutation burden measures
mutCN_cutoff=25 #If the mitochondrial copy number is over 25, retain mutation even if is in the "CN correlating muts" list
vaf_cut_off=0.03 #only include VAFs above a cutoff
all_sum_of_vaf_df=Map(dataset_mito_data=all_mito_datasets[all_cohorts],dataset=all_cohorts, function(dataset_mito_data,dataset) {
  cat(dataset,sep="\n")
  CN_correlating_muts<-dataset_mito_data$CN_correlating_muts
  sum_of_vaf_df<-Map(list=dataset_mito_data,Exp_ID=names(dataset_mito_data),function(list,Exp_ID){
    
    #sum the vaf of mutations that:
    #(1) Pass shearwater,
    #(2) are not indels, specific artefacts or 'copy number correlating mutations' (i.e. most likely artefacts from mismapping of nuclear DNA or low-level contamination)
    #(3) are not heteroplasmic oocyte mutations - if you have these, all cells have a 'headstart' in their mtDNA mutation acquisition; therefore no fair comparison
    
    CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
    
    df<-as.data.frame(colSums((list$matrices$vaf*(list$matrices$vaf>vaf_cut_off)*list$matrices$SW*(CN_correlating_mut_removal_mat))[!grepl("DEL|INS",rownames(list$matrices$vaf))&
                                                                                                                                      !rownames(list$matrices$vaf)%in%c(exclude_muts,list$het_oocyte_muts),list$tree$tip.label]))%>%
      tibble::rownames_to_column(var="Sample")%>%
      mutate(exp_ID=Exp_ID)%>%
      dplyr::rename(sum_of_vaf=2)
    
    return(df)
  })%>%dplyr::bind_rows()%>%
    mutate(dataset=dataset,.before=1)
  
  return(sum_of_vaf_df)
})%>%dplyr::bind_rows()

#-----------------------------------------------------------------------------------#
##----------plot 'sum of vaf' data individually for all datasets --------------------
#-----------------------------------------------------------------------------------#

temp=lapply(all_cohorts,function(this_dataset) {
  cat(this_dataset,sep="\n")
  
  min_samples_per_individual=8
  sum_of_vaf_summary<-all_sum_of_vaf_df%>%
    filter(dataset==this_dataset)%>%
    group_by(exp_ID)%>%
    summarise(mean=mean(sum_of_vaf,na.rm=T),median=median(sum_of_vaf,na.rm = T))%>%
    mutate(dataset=this_dataset,.before=1)%>%
    left_join(samples_per_individual,by=c("exp_ID","dataset"))%>%
    filter(n>=min_samples_per_individual)
  
  sum_of_vaf_plot_log<-all_sum_of_vaf_df%>%
    filter(exp_ID%in%sum_of_vaf_summary$exp_ID & dataset==this_dataset)%>%
    mutate(sum_of_vaf=ifelse(sum_of_vaf==0,0.001,sum_of_vaf))%>%
    ggplot(aes(x=forcats::fct_reorder(Sample,sum_of_vaf),y=sum_of_vaf))+
    geom_point(size=0.1)+
    theme_classic()+
    my_theme+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.text.x = element_text(size=6))+
    facet_grid(cols=vars(factor(exp_ID,levels=ref_df%>%arrange(Age)%>%filter(Cohort==this_dataset&ID%in%sum_of_vaf_summary$exp_ID)%>%pull(ID))),scales = "free",space="free",drop = T)+
    scale_y_log10(breaks=c(0.005,0.01,0.03,0.06,0.12,0.25,0.5,1,2,4,8))+
    geom_hline(aes(yintercept=mean),linetype=2,data=sum_of_vaf_summary%>%mutate(mean=ifelse(mean==0,0.001,mean)),col="red",linewidth=0.5)+
    #geom_hline(aes(yintercept=median),linetype=2,data=sum_of_vaf_summary%>%mutate(median=ifelse(median==0,0.001,median)),col="red",size=0.5)+
    geom_text(aes(x=0,y=mean_pos,label=paste0("bar(x) == ",round(mean,2))),size=2,nudge_x=+5,nudge_y=+0.5,data=sum_of_vaf_summary%>%mutate(mean_pos=ifelse(mean==0,0.001,mean)),parse=T)+
    #geom_text(aes(x=0,y=median_pos,label=paste0("tilde(x) == ",round(median,3))),size=2,nudge_x=+115,nudge_y=+0.5,data=sum_of_vaf_summary%>%mutate(median_pos=ifelse(median==0,0.001,median)),parse=T)+
    labs(x="Sample",y="Mutation burden\n(sum of VAF)")+
    theme(panel.spacing.x=unit(1, "mm"),strip.text.x = element_text(angle=90))
  ggsave(filename = paste0(plots_dir,"sum_of_vaf_log_",this_dataset,".pdf"),plot=sum_of_vaf_plot_log,height=2,width=7)
  
})

#-----------------------------------------------------------------------------------#
##----------plot mean sum of vaf by age data individually for all datasets ----------
#-----------------------------------------------------------------------------------#

# Highly clonally related samples are not 'independent' measures of mtDNA mutation burden.
# Therefore, for some analyses, may choose to drop these co-correlated readings.
# For these, randomly select one of any highly clonally related clades (those with >half of shared ancestry as assessed by ultrametric tree), and drop the rest

all_samples_to_drop<-lapply(all_mito_datasets[all_cohorts],function(dataset_mito_data) {
  lapply(names(dataset_mito_data),function(this_exp_ID) {
    mut_burden<-get_mut_burden(dataset_mito_data[[this_exp_ID]]$tree.ultra)[1]
    expanded_clade_nodes<-get_expanded_clade_nodes(tree=dataset_mito_data[[this_exp_ID]]$tree.ultra,min_samples=2,height_cut_off=mut_burden/2,min_clonal_fraction = 0)
    samples_to_drop=unlist(lapply(expanded_clade_nodes$nodes,function(node) {
      clade_samples<-getTips(tree=dataset_mito_data[[this_exp_ID]]$tree.ultra,node)
      to_drop=sample(clade_samples,size=length(clade_samples)-1)
    }))
  })
})%>%unlist()

all_lms=lapply(all_cohorts,function(this_dataset) {
  
  cat(this_dataset,sep="\n")
  ref_df_cohort<-ref_df%>%filter(Cohort==this_dataset)
  
  #Do a version where mean/ median are calculated using all samples (even if highly clonally related)
  min_samples_per_individual=8
  sum_of_vaf_summary<-all_sum_of_vaf_df%>%
    filter(dataset==this_dataset)%>%
    group_by(exp_ID)%>%
    summarise(n=n(),mean=mean(sum_of_vaf,na.rm=T),median=median(sum_of_vaf,na.rm = T))%>%
    mutate(dataset=this_dataset,.before=1)%>%
    filter(n>=min_samples_per_individual)%>%
    left_join(ref_df_cohort,by=c("exp_ID"="ID"))
  
  exp_cols<-c("#aee39a", "#18786a", "#99def9", "#0362a0", "#e689eb", "#53348e", "#9296ee", "#9e37d0", "#32e195", "#50942f", "#a3d71e", "#7c440e", "#ecbcab", "#a53460", "#fb899b", "#ed4b04", "#ffce54", "#f82387", "#61f22d", "#d07d09","#E41A1C")
  names(exp_cols)<-c("Dummy_foetal",sum_of_vaf_summary$exp_ID)
  
  #If there are no very young individuals in the cohort, create a dummy dataframe to root the regression to ~0 mutations at conception
  #Based on mutation burdens of the cord blood HSPC data, assumption is that mtDNA mutation burdens likely to be fairly similar in most tissues at birth
  if(sum(sum_of_vaf_summary$Age<5)==0) {
    foetal_dummy_df<-data.frame("exp_ID"="Dummy_foetal","Cohort"=this_dataset,"mean"=0.03,"Age"=0)
  } else {
    foetal_dummy_df<-data.frame("exp_ID"=c(),"Cohort"=c(),"mean"=c(),"Age"=c())
  }
  
  mean_sum_of_vaf_by_age<-sum_of_vaf_summary%>%
    bind_rows(foetal_dummy_df)%>%
    ggplot(aes(x=Age,y=mean))+
    geom_point(aes(col=factor(exp_ID,levels = c(foetal_dummy_df$exp_ID,ref_df_cohort$ID[order(ref_df_cohort$Age)]))))+
    theme_bw()+
    scale_color_manual(values=exp_cols)+
    scale_x_continuous(limits=c(-1,82))+
    labs(col="Individual",y="Mean mutation burden (sum of VAF)")+
    guides(col=guide_legend(nrow=10))+
    stat_poly_line(formula=y~x,col="black",linewidth=0.5) +
    stat_poly_eq(use_label(c("eq", "R2")),formula=y~x,size=2) +
    #geom_smooth(aes(x=Age,y=mean),col="black",linewidth=0.5,method="lm",inherit.aes = F,fullrange=T)+
    my_theme+
    theme(legend.key.height = unit(3,"mm"),legend.title=element_text(size=7),legend.text=element_text(size=5))
  ggsave(filename=paste0(plots_dir,"mean_sum_of_vaf_by_age",this_dataset,".pdf"),plot=mean_sum_of_vaf_by_age,height=1.7,width=3.5)
  
  #Now a version where mean/ median are calculated after excluding closely clonally-related samples
  sum_of_vaf_summary<-all_sum_of_vaf_df%>%
    filter(dataset==this_dataset & !Sample%in%all_samples_to_drop)%>%
    group_by(exp_ID)%>%
    summarise(n=n(),mean=mean(sum_of_vaf,na.rm=T),median=median(sum_of_vaf,na.rm = T))%>%
    mutate(dataset=this_dataset,.before=1)%>%
    filter(n>=min_samples_per_individual)%>%
    left_join(ref_df_cohort,by=c("exp_ID"="ID"))
  
  lm.mean_mutburden_by_age=lm(mean~Age,data=sum_of_vaf_summary%>%
                                bind_rows(foetal_dummy_df))
  
  mean_sum_of_vaf_by_age<-sum_of_vaf_summary%>%
    bind_rows(foetal_dummy_df)%>%
    ggplot(aes(x=Age,y=mean))+
    geom_point(aes(col=factor(exp_ID,levels = c(foetal_dummy_df$exp_ID,ref_df_cohort$ID[order(ref_df_cohort$Age)]))))+
    theme_bw()+
    scale_color_manual(values=exp_cols)+
    scale_x_continuous(limits=c(-1,82))+
    stat_poly_line(formula=y~x,col="black",linewidth=0.5) +
    stat_poly_eq(use_label(c("eq", "R2")),formula=y~x,size=2) +
    labs(col="Individual",y="Mean mutation burden (sum of VAF)")+
    guides(col=guide_legend(title = element_blank(),ncol=1))+
    #geom_smooth(aes(x=Age,y=mean),col="black",linewidth=0.5,method="lm",inherit.aes = F,fullrange=T)+
    my_theme+
    theme(legend.spacing.x=unit(0.5,"mm"),legend.box.spacing = unit(0, "pt"),legend.key.height = unit(2.8,"mm"),legend.title=element_text(size=7),legend.text=element_text(size=5))
  ggsave(filename=paste0(plots_dir,"mean_sum_of_vaf_by_age",this_dataset,"_cocorrelated_excluded.pdf"),plot=mean_sum_of_vaf_by_age,height=2,width=2.8)
  
  return(lm.mean_mutburden_by_age)
})

names(all_lms)<-all_cohorts

#-----------------------------------------------------------------------------------#
##-------------VISUALIZE THE INDIVIDUAL LINEAR MODELS FROM ALL THE DATASETS----------
#-----------------------------------------------------------------------------------#

name_conversion_vec=c("Blood","Lymphoid","MPN","Colon","IBD-affected\nColon","MUTYH-mutant\nColon","Endometrium","Bronchial\nepithelium")
names(name_conversion_vec)=c("blood","lymph","NW","HL","SO","PR","LM","KY")
tissue_cols<-c("#96e97c","#fd81c8","#145a6a","#65e6f9","#781486","#5f70cc","#aea2eb","#1d6d1f")
names(tissue_cols)<-name_conversion_vec

all.lm.coefs<-Map(lm.mean_mutburden_by_age=all_lms,dataset=all_cohorts,function(lm.mean_mutburden_by_age,dataset) {
  coefs<-summary(lm.mean_mutburden_by_age)$coefficients[,'Estimate']
  r2<-summary(lm.mean_mutburden_by_age)$r.squared
  confints<-confint(lm.mean_mutburden_by_age)
  
  return(cbind(data.frame(dataset=dataset),term_type=c("Intercept","Age coefficient"),Estimate=coefs,data.frame(r2=r2),confints))
  
})%>%dplyr::bind_rows()

Age.coefs.plot<-all.lm.coefs%>%
  mutate(dataset=factor(name_conversion_vec[dataset],name_conversion_vec))%>%
  ggplot(aes(y=dataset,x=Estimate,xmin=`2.5 %`,xmax=`97.5 %`,col=dataset))+
  geom_point()+geom_errorbar(width=0.25)+
  scale_color_manual(values=tissue_cols)+
  theme_classic()+
  facet_wrap(~term_type,scales="free")+
  my_theme+
  theme(axis.title.y = element_blank(),legend.position="none")
ggsave(filename=paste0(plots_dir,"lm_comparison.pdf"),Age.coefs.plot,width=4.5,height=2)


#-----------------------------------------------------------------------------------#
##----------------------DO SINGLE LMER ACROSS DATA----------------------------------
#-----------------------------------------------------------------------------------#


all_sum_of_vaf_summary_df=lapply(all_cohorts,function(this_dataset) {
  
  cat(this_dataset,sep="\n")
  ref_df_cohort<-ref_df%>%filter(Cohort==this_dataset)
  
  #Do a version where mean/ median are calculated using all samples (even if highly clonally related)
  min_samples_per_individual=8
  sum_of_vaf_summary<-all_sum_of_vaf_df%>%
    filter(dataset==this_dataset & !Sample%in%all_samples_to_drop)%>%
    group_by(exp_ID)%>%
    summarise(n=n(),mean=mean(sum_of_vaf,na.rm=T),median=median(sum_of_vaf,na.rm = T))%>%
    mutate(dataset=this_dataset,.before=1)%>%
    filter(n>=min_samples_per_individual)%>%
    left_join(ref_df_cohort,by=c("exp_ID"="ID"))

  #If there are no young individuals in the cohort, create a dummy dataframe to root the regression to ~0 mutations at conception
  #Based on mutation burdens of the cord blood HSPC data, assumption is that mtDNA mutation burdens likely to be fairly similar in most tissues at birth
  if(sum(sum_of_vaf_summary$Age<5)==0 & F) {
    foetal_dummy_df<-data.frame("exp_ID"="Dummy_foetal","dataset"=this_dataset,"mean"=0.03,"Age"=0)
  } else {
    foetal_dummy_df<-data.frame("exp_ID"=c(),"dataset"=c(),"mean"=c(),"Age"=c())
  }
  
  sum_of_vaf_summary<-dplyr::bind_rows(sum_of_vaf_summary,foetal_dummy_df)
  
  return(sum_of_vaf_summary)
})%>%dplyr::bind_rows()

#lmer1<-lme4::lmer(mean~Age+Age:dataset+(1|exp_ID),data=all_sum_of_vaf_summary_df)
lmer2<-lme4::lmer(mean~Age*dataset+(1|exp_ID),data=all_sum_of_vaf_summary_df)


#Plot results of the lmer with both fixed and Age interaction term
name_conversion_vec=c("Blood","Lymphoid","MPN","Colon","IBD-affected\nColon","MUTYH-mutant\nColon","Endometrium","Bronchial\nepithelium")
names(name_conversion_vec)=c("blood","lymph","NW","HL","SO","PR","LM","KY")
tissue_cols<-c("#96e97c","#fd81c8","#145a6a","#65e6f9","#781486","#5f70cc","#aea2eb","#1d6d1f")
names(tissue_cols)<-name_conversion_vec


CI_table<-confint(lmer2)
lmer_info_for_plotting<-as.data.frame(summary(lmer2)$coefficients)%>%
  tibble::rownames_to_column(var="name")%>%
  left_join(as.data.frame(CI_table)%>%tibble::rownames_to_column(var="name"),by="name")%>%
  filter(grepl("dataset",name))%>%
  mutate(name=gsub("dataset","",name))%>%
  mutate(term_type=ifelse(grepl("Age",name),"Age interaction term","Fixed term"))%>%
  mutate(dataset=gsub("Age:","",name))%>%
  mutate(dataset=factor(name_conversion_vec[dataset],levels=name_conversion_vec))

p1.1<-lmer_info_for_plotting%>%
  ggplot(aes(y=dataset,x=Estimate,xmin=`2.5 %`,xmax=`97.5 %`,col=dataset))+
  geom_point()+
  scale_color_manual(values=tissue_cols)+
  geom_vline(xintercept=0,linetype=2)+
  facet_wrap(~term_type,scales = "free")+
  geom_errorbar(width=0.4)+
  theme_classic()+
  my_theme+
  theme(axis.title.y=element_blank(),legend.position = "none")

ggsave(filename=paste0(plots_dir,"lmer_terms.pdf"),p1.1,width=4.5,height=2)

coef_rename=c("lm: Intercept","lm: Age coefficient","lmer: Fixed term","lmer: Age interaction term")
names(coef_rename)<-c("Intercept","Age coefficient","Fixed term","Age interaction term")
regression_comb_plot<-all.lm.coefs%>%
  mutate(dataset=factor(name_conversion_vec[dataset],name_conversion_vec))%>%
  dplyr::bind_rows(lmer_info_for_plotting)%>%
  mutate(term_type=factor(coef_rename[term_type],levels=coef_rename))%>%
  ggplot(aes(y=dataset,x=Estimate,xmin=`2.5 %`,xmax=`97.5 %`,col=dataset))+
  geom_point()+
  scale_color_manual(values=tissue_cols)+
  geom_vline(xintercept=0,linetype=2)+
  facet_wrap(~term_type,scales = "free_x",nrow=1)+
  geom_errorbar(width=0.4)+
  theme_classic()+
  my_theme+
  theme(axis.title.y=element_blank(),legend.position = "none")

ggsave(filename=paste0(plots_dir,"regression_comb_plot.pdf"),regression_comb_plot,width=7,height=2)

#-----------------------------------------------------------------------------------#
##----------------------SELECTION ANALYSIS ------------------------
#-----------------------------------------------------------------------------------#

library(dndscv)
mtref_rda_path=ifelse(Sys.info()['sysname']=="Darwin",paste0(root_dir,"data/mtref.rda"),"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/mtref.rda")
input.dir <- paste0(root_dir,"dnds_tables/")

dnds_theme<-theme(panel.border = element_rect(color = "black",
                                              fill = NA,
                                              linewidth = 0.75),
                  strip.text = element_text(face="plain", size=6, colour = "black",),
                  strip.background = element_rect(fill="white", colour="black", linewidth =1),
                  axis.text.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = 0.5, face = "plain"),
                  axis.text.y = element_text(color = "black", size = 6, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
                  axis.title.x = element_text(color = "black", size = 8, angle = 0, hjust = .5, vjust = 0, face = "plain"),
                  axis.title.y = element_text(color = "black", size = 8, angle = 90, hjust = .5, vjust = .5, face = "plain", 
                                              margin = margin(t = 0, r = 10, b = 0, l = 0)), 
                  plot.title = element_text(color = "black", size = 10, hjust = .5, face = "plain"), 
                  legend.position = "right")

#exclude ND6 from analysis since its on the other strand
all_mtDNA_genes <- c("MT-CYB", "MT-ND5", "MT-ND2", "MT-ND4", "MT-ND1", "MT-CO3", "MT-ATP6","MT-ND3", "MT-ATP8", "MT-ND4L", "MT-CO2", "MT-CO1","MT-ND5","MT-ND6")
target_genes <- c("MT-CYB", "MT-ND5", "MT-ND2", "MT-ND4", "MT-ND1", "MT-CO3", "MT-ATP6","MT-ND3", "MT-ATP8", "MT-ND4L", "MT-CO2", "MT-CO1")
#target_genes <- c("MT-CYB", "MT-ND5", "MT-ND2", "MT-ND4", "MT-ND1", "MT-CO3", "MT-ATP6","MT-ND3", "MT-ATP8", "MT-ND4L", "MT-CO2", "MT-CO1","MT-ND5")

df_tidy_annotated<-Map(df_tidy=all_df_tidy,this_tissue=names(all_df_tidy),function(df_tidy,this_tissue) {
  cat(this_tissue,sep="\n")
  
  # read in the mtdna variant file and remove patient id
  mtdna.variant.data <- df_tidy%>%
    separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select("sampleID"=Sample,chr,pos,ref,mut,vaf)%>% # next, rename/reorder the columns to make them compatible with dnds input format
    dplyr::filter(!duplicated(.))%>% # remove duplicated variants
    arrange(sampleID,pos)
  
  # run dnds with the mtDNA variants to annotate them (the selection analysis isn't actually used here)
  mtdna.dndsout <- dndscv(mtdna.variant.data, gene_list=all_mtDNA_genes, 
                          refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
  
  # get the results with the annotated variants
  annotated.mtdna.variants <- left_join(mtdna.variant.data,mtdna.dndsout$annotmuts,by=c("sampleID","chr","pos","ref","mut"))%>%
    tidyr::replace_na(replace=list(impact="Non-Coding"))%>%
    mutate(tissue=this_tissue,.before=1)
  return(annotated.mtdna.variants)
})%>%dplyr::bind_rows()

complete.annotated.mutation.table<-df_tidy_annotated%>%
  left_join(all_df_tidy%>%dplyr::bind_rows()%>%dplyr::select(Sample,exp_ID),by=c("sampleID"="Sample"),relationship="many-to-many")%>%
  dplyr::rename("patientID"=exp_ID)%>%
  filter(!duplicated(.))%>%
  arrange(patientID,pos)%>% # order the dataframe
  mutate(impact=ifelse(impact%in%c("Stop_loss","Nonsense"),"Truncating",impact))%>% #rename stop loss and nonsense mutations
  tidyr::unite(col="mut_ref",chr,pos,ref,mut,sep="_",remove=F)

my_comparisons <- list(c("Missense", "Synonymous"), c("Synonymous", "Truncating"))


### Generate global dnds ---------
dndscv_by_tissue<-Map(df_tidy=all_df_tidy[all_cohorts],this_tissue=all_cohorts,function(df_tidy,this_tissue) {
  
  cat(this_tissue,sep="\n")
  tissue_info<-df_tidy%>%
    separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>%
    dplyr::filter(!duplicated(.))%>%
    arrange(sampleID,pos)
  
  mtdna.dndsout <- dndscv(tissue_info, gene_list=target_genes, 
                          refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
  return(mtdna.dndsout)
})

# processing the output to get a combined 'globaldnds' dataframe by tissue type
rename_vec=c("Missense","Nonsense","Truncating","Overall")
names(rename_vec)=c("wmis","wnon","wtru","wall")

globaldnds_res<-Map(dndsout=dndscv_by_tissue,this_tissue=all_cohorts,function(dndsout,this_tissue) {
  dndsout$globaldnds%>%
    filter(!is.na(name) & complete.cases(.))%>%
    mutate(name=rename_vec[name])%>%
    mutate(tissue=this_tissue)
})%>%dplyr::bind_rows()

nb.cols <- length(name_conversion_vec)
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(3)

max_dnds_value<-2
stats_to_include=c("Missense","Truncating")
global_dnds_nonblood<-globaldnds_res%>%
  filter(name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         tissue=factor(name_conversion_vec[tissue],levels=name_conversion_vec))%>%
  ggplot(aes(x=tissue, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Tissue") +
  theme_classic() + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6))+
  scale_y_continuous(limits=c(0,max_dnds_value))+
  scale_color_manual(values = mycolors[1:2], name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', linewidth = 0.5) +
  my_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle=90),
        legend.position="none")

ggsave(filename = paste0(plots_dir,"global_dnds_nonblood.pdf"),global_dnds_nonblood,width=2,height=2)

#-----------------------------------------------------------------------------------#
##----------------------Get gene-level dNdS info by tissue---------------
#-----------------------------------------------------------------------------------#
selcv_res_by_tissue<-Map(dndsout=dndscv_by_tissue,this_tissue=all_cohorts,function(dndsout,this_tissue) {
  temp<-dndsout$sel_cv%>%
    dplyr::select(gene_name,wmis_cv,wnon_cv,qmis_cv,qtrunc_cv)
  
  dplyr::bind_rows(temp%>%dplyr::select(gene_name,"dNdS"=wmis_cv,"qval"=qmis_cv)%>%mutate(type="Missense"),
                   temp%>%dplyr::select(gene_name,"dNdS"=wnon_cv,"qval"=qtrunc_cv)%>%mutate(type="Nonsense"))%>%
    mutate(tissue=this_tissue)
})%>%dplyr::bind_rows()

gene_order=c(paste0("MT-ND",1:6),"MT-ND4L","MT-CYB",paste0("MT-CO",1:3),"MT-ATP6","MT-ATP8")
gene_order<-gene_order[gene_order%in%target_genes]

dnds_heatmap_by_gene_by_tissue<-selcv_res_by_tissue%>%
  mutate(dnds_if_signif=ifelse(qval<0.1,paste0(sprintf(dNdS, fmt = '%#.1f'),"*"),ifelse(dNdS>2,sprintf(dNdS, fmt = '%#.1f'),"")),
         dNdS=ifelse(dNdS>2,2,dNdS),
         gene_name=factor(gene_name,levels=rev(gene_order)),
         tissue=factor(name_conversion_vec[tissue],levels=name_conversion_vec))%>%
  ggplot(aes(y=gene_name,x=tissue,fill=dNdS,label=dnds_if_signif))+
  geom_tile()+
  facet_grid(~type)+
  geom_text(col="white",size=1.6)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(n=11,name = "RdYlGn"))+
  labs(x="Variant allele fraction",y="")

ggsave(filename = paste0(plots_dir,"dnds_heatmap_by_gene_by_tissue.pdf"),dnds_heatmap_by_gene_by_tissue,width=4,height=2.5)


#-----------------------------------------------------------------------------------#
##----------------------Enhance the tidy dataframe with 'VAF bin' info---------------
#-----------------------------------------------------------------------------------#
#Define the VAF 'bins' that will be used for comparing VAF distributions
#Can set this on a log scale, or as decile bins
boundaries<-c(0,0.01,0.1,0.2,0.5,1)
new_VAF_groups=paste0(100*head(boundaries,-1),"-",100*tail(boundaries,-1),"%")
VAF_groups=data.frame(labels=new_VAF_groups,LL=head(boundaries,-1),UL=tail(boundaries,-1))

#Assign each observed VAF into its 'VAF group' bin
all_df_tidy<-lapply(all_df_tidy,function(df_tidy) {
  df_tidy$VAF_group=sapply(df_tidy$vaf,function(vaf) {
    if(vaf==0) {
      return("Absent")
    } else {
      return(VAF_groups$labels[VAF_groups$LL<vaf & VAF_groups$UL>=vaf])
    }
  })
  return(df_tidy)
})

#-----------------------------------------------------------------------------------#
##----------------------Run global dNdS by VAF group---------------
#-----------------------------------------------------------------------------------#
dndscv_by_tissue_and_vaf<-Map(df_tidy=all_df_tidy[all_cohorts],this_tissue=all_cohorts,function(df_tidy,this_tissue) {
  cat(this_tissue,sep = "\n")
  dndscv_by_vaf_tissue<-lapply(new_VAF_groups[2:5],function(VAF_bin) {
    cat(VAF_bin,sep="\n")
    tissue_vaf_info<-df_tidy%>%
      dplyr::filter(VAF_group==VAF_bin)%>%
      separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
      mutate(pos=as.numeric(pos))%>%
      dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>% 
      dplyr::filter(!duplicated(.))%>% #only count each mutation once per individual in each VAF group
      arrange(sampleID,pos)
    
    mtdna.dndsout <- dndscv(tissue_vaf_info, gene_list=target_genes, 
                            refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
    return(mtdna.dndsout)
  })
  return(dndscv_by_vaf_tissue)
})

globaldnds_res_by_tissue_and_vaf<-Map(list1=dndscv_by_tissue_and_vaf,this_tissue=all_cohorts,function(list1,this_tissue) {
  globaldnds_res_by_vaf<-Map(dndsout=list1,VAF_bin=new_VAF_groups[2:5],function(dndsout,VAF_bin) {
    dndsout$globaldnds%>%
      filter(!is.na(name) & complete.cases(.))%>%
      mutate(name=rename_vec[name])%>%
      mutate(VAF_group=VAF_bin)
  })%>%dplyr::bind_rows()%>%
    mutate(tissue=this_tissue)
})%>%dplyr::bind_rows()

nb.cols <- length(name_conversion_vec)
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(3)

max_dnds_value<-5
stats_to_include=c("Missense","Truncating")
global_dnds_by_tissue_and_vaf_plot<-globaldnds_res_by_tissue_and_vaf%>%
  filter(VAF_group%in%VAF_groups$labels[2:5] & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         VAF_group=factor(VAF_group,levels=VAF_groups$labels[2:5]),
         tissue=factor(name_conversion_vec[tissue],levels=name_conversion_vec))%>%
  ggplot(aes(x=VAF_group, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_hline(yintercept=1, linetype='dashed', col = 'black', linewidth = 0.5) +
  geom_linerange(position= position_dodge2(width=0.75), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mycolors[1:2], name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  theme_classic()+
  my_theme+
  facet_grid(cols=vars(tissue))+
  theme(legend.position="none",
        axis.text.x=element_text(angle=90),
        axis.title.x=element_blank(),
        strip.text.y = element_text(angle = 0))

ggsave(filename = paste0(plots_dir,"global_dnds_by_tissue_and_vaf_plot.pdf"),global_dnds_by_tissue_and_vaf_plot,width=5,height=2)

max_dnds_value<-5
global_dnds_by_vaf_MPN_and_lymph_onlyplot<-globaldnds_res_by_tissue_and_vaf%>%
  filter(VAF_group%in%VAF_groups$labels[2:5] & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         VAF_group=factor(VAF_group,levels=VAF_groups$labels[2:5]),
         tissue=factor(name_conversion_vec[tissue],levels=name_conversion_vec))%>%
  filter(tissue%in%c("Lymphoid","MPN"))%>%
  ggplot(aes(x=VAF_group, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_hline(yintercept=1, linetype='dashed', col = 'black', linewidth = 0.5) +
  geom_linerange(position= position_dodge2(width=0.75), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mycolors[1:2], name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  theme_classic()+
  my_theme+
  facet_grid(cols=vars(tissue))+
  theme(legend.position="none",
        axis.text.x=element_text(angle=90),
        axis.title.x=element_blank(),
        strip.text.y = element_text(angle = 0))

ggsave(filename = paste0(plots_dir,"global_dnds_by_vaf_MPN_and_lymph_onlyplot.pdf"),global_dnds_by_vaf_MPN_and_lymph_onlyplot,width=2,height=2)

max_dnds_value<-2.5
global_dnds_by_vaf_normal_nonblood_onlyplot<-globaldnds_res_by_tissue_and_vaf%>%
  filter(VAF_group%in%VAF_groups$labels[2:5] & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         VAF_group=factor(VAF_group,levels=VAF_groups$labels[2:5]),
         tissue=factor(name_conversion_vec[tissue],levels=name_conversion_vec))%>%
  filter(tissue%in%c("Colon","Endometrium","Bronchial\nepithelium"))%>%
  ggplot(aes(x=VAF_group, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_hline(yintercept=1, linetype='dashed', col = 'black', linewidth = 0.5) +
  geom_linerange(position= position_dodge2(width=0.75), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mycolors[1:2], name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  theme_classic()+
  my_theme+
  facet_grid(cols=vars(tissue))+
  theme(legend.position="none",
        axis.text.x=element_text(angle=90),
        axis.title.x=element_blank(),
        strip.text.y = element_text(angle = 0))

ggsave(filename = paste0(plots_dir,"global_dnds_by_vaf_normal_nonblood_onlyplot.pdf"),global_dnds_by_vaf_normal_nonblood_onlyplot,width=2.5,height=2)


gene_order=c(paste0("MT-ND",1:6),"MT-ND4L","MT-CYB",paste0("MT-CO",1:3),"MT-ATP6","MT-ATP8")
gene_order<-gene_order[gene_order%in%target_genes]

selcv_res_by_tissue_and_vaf<-Map(list1=dndscv_by_tissue_and_vaf,this_tissue=names(dndscv_by_tissue_and_vaf),function(list1,this_tissue) {
  
  Map(dndsout=list1,VAF_bin=new_VAF_groups[2:5],function(dndsout,VAF_bin) {
    temp<-dndsout$sel_cv%>%
      dplyr::select(gene_name,wmis_cv,wnon_cv,qmis_cv,qtrunc_cv)
    
    dplyr::bind_rows(temp%>%dplyr::select(gene_name,"dNdS"=wmis_cv,"qval"=qmis_cv)%>%mutate(type="Missense"),
                     temp%>%dplyr::select(gene_name,"dNdS"=wnon_cv,"qval"=qtrunc_cv)%>%mutate(type="Nonsense"))%>%
      mutate(VAF_group=VAF_bin)
  })%>%dplyr::bind_rows()%>%
    mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels[2:5]),tissue=this_tissue)
  
})%>%dplyr::bind_rows()

dnds_heatmap_by_tissue_and_vaf_by_gene<-selcv_res_by_tissue_and_vaf%>%
  mutate(dnds_if_signif=ifelse(qval<0.1,paste0(sprintf(dNdS, fmt = '%#.1f'),"*"),ifelse(dNdS>2,sprintf(dNdS, fmt = '%#.1f'),"")),
         dNdS=ifelse(dNdS>2,2,dNdS),
         gene_name=factor(gene_name,levels=rev(gene_order)),
         tissue=factor(name_conversion_vec[tissue],levels=name_conversion_vec))%>%
  filter(VAF_group%in%VAF_groups$labels[2:5] & tissue!="Blood")%>%
  ggplot(aes(y=gene_name,x=VAF_group,fill=dNdS,label=dnds_if_signif))+
  geom_tile()+
  facet_grid(tissue~type)+
  geom_text(col="white",size=1.5)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(n=11,name = "RdYlGn"))+
  labs(x="Variant allele fraction",y="")

ggsave(filename = paste0(plots_dir,"dnds_heatmap_by_tissue_and_vaf_by_gene.pdf"),dnds_heatmap_by_tissue_and_vaf_by_gene,width=4,height=6.5)

#Now including lymphoid and blood only
dnds_heatmap_by_tissue_and_vaf_by_gene_blood_and_lymph_only<-selcv_res_by_tissue_and_vaf%>%
  mutate(dnds_if_signif=ifelse(qval<0.1,paste0(sprintf(dNdS, fmt = '%#.1f'),"*"),ifelse(dNdS>2,sprintf(dNdS, fmt = '%#.1f'),"")),
         dNdS=ifelse(dNdS>2,2,dNdS),
         gene_name=factor(gene_name,levels=rev(gene_order)),
         tissue=factor(name_conversion_vec[tissue],levels=name_conversion_vec))%>%
  filter(VAF_group%in%VAF_groups$labels[2:5] & tissue%in%c("Lymphoid","Blood"))%>%
  ggplot(aes(y=gene_name,x=VAF_group,fill=dNdS,label=dnds_if_signif))+
  geom_tile()+
  facet_grid(tissue~type)+
  geom_text(col="white",size=1.5)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(n=11,name = "RdYlGn"))+
  labs(x="Variant allele fraction",y="")

ggsave(filename = paste0(plots_dir,"dnds_heatmap_by_tissue_and_vaf_by_gene_blood_and_lymph_only.pdf"),dnds_heatmap_by_tissue_and_vaf_by_gene_blood_and_lymph_only,width=4,height=4)

#-----------------------------------------------------------------------------------#
# -------------MPN analysis, test dNdS by driver clone/ wild-type ---------------
#-----------------------------------------------------------------------------------#

#1. Add basic nDNA driver mutation info to the object
all_mito_datasets$NW<-Map(list=all_mito_datasets$NW,exp_ID=names(all_mito_datasets$NW),function(list,exp_ID) {
  #annot_muts_file<-paste0(root_dir,"/data/nonblood/annot_files_filtered/annotated_muts_filt_",exp_ID,".Rds")
  annot_muts_file<-paste0("annotated_muts_filt_",exp_ID,".Rds")
  if(file.exists(annot_muts_file)) {
    cat(paste("Importing annotated mutations dataframe for",exp_ID),sep="\n")
    list$nDNA_mats<-readRDS(annot_muts_file)
    
    if(nrow(list$nDNA_mats$mat)==0) {
      stop(return(list))
    }
    
    tree_samples<-list$tree.ultra$tip.label[list$tree.ultra$tip.label!="Ancestral"]
    
    if(any(!tree_samples%in%colnames(list$nDNA_mats$NV))) {
      cat("Dropping samples not found in the count matrix.")
      drop_tips<-tree_samples[!tree_samples%in%colnames(list$nDNA_mats$NV)]
      list$tree.ultra<-drop.tip(list$tree.ultra,tip = drop_tips)
      list$tree<-drop.tip(list$tree,tip = drop_tips)
      tree_samples<-tree_samples[tree_samples%in%colnames(list$nDNA_mats$NV)]
    }
    
    # mtr<-as.matrix(cbind(list$nDNA_mats$NV[,tree_samples,drop=F],matrix(0,ncol = 1,nrow = nrow(list$nDNA_mats$NV),dimnames = list(rownames(list$nDNA_mats$NV),"Ancestral"))))
    # dep<-as.matrix(cbind(list$nDNA_mats$NR[,tree_samples,drop=F],matrix(10,ncol = 1,nrow = nrow(list$nDNA_mats$NR),dimnames = list(rownames(list$nDNA_mats$NV),"Ancestral"))))
    
    mtr<-as.matrix(list$nDNA_mats$NV[,tree_samples,drop=F])
    dep<-as.matrix(list$nDNA_mats$NR[,tree_samples,drop=F])
    
    
    res<-treemut::assign_to_tree(tree=list$tree.ultra,mtr=mtr,dep=dep)
    list$nDNA_mats$mat$node<-res$tree$edge[res$summary$edge_ml,2]
    return(list)
  } else {
    return(list)
  }
})

#Generate a similar 'sum of vaf' dataframe, but now incorporating driver information
mut_vs_wt_list<-Map(list=all_mito_datasets$NW,exp_ID=names(all_mito_datasets$NW),function(list,exp_ID) {
  cat(exp_ID,sep="\n")
  #Using the 'get_expanded_clade_nodes' function to effectively cut across the tree at 100 mutations of molecular time
  #Some of the clades will be singletons
  ec_df<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = 100,min_clonal_fraction = 0,min_samples=1)%>%
    left_join(list$nDNA_mats$mat%>%dplyr::select(node,mut_ref,GENE,CDS),by=c("nodes"="node"))%>%
    mutate(exp_ID=exp_ID,.before=1)
  
  # genes_of_interest=c("JAK2","TET2","DNMT3A","ASXL1")
  # mut_nodes<-ec_df%>%filter(GENE%in%genes_of_interest)%>%pull(nodes)%>%unique()
  # wt_nodes<-ec_df$nodes[!ec_df$nodes%in%mut_nodes]
  
  mut_nodes<-ec_df%>%filter(n_samples>1)%>%pull(nodes)
  wt_nodes<-ec_df%>%filter(n_samples==1)%>%pull(nodes)
  
  mut_samples<-unlist(lapply(mut_nodes,function(node) getTips(list$tree.ultra,node=node)))
  wt_samples<-unlist(lapply(wt_nodes,function(node) getTips(list$tree.ultra,node=node)))
  
  return(list(mut=mut_samples,wt=wt_samples))
})

all_mut<-unlist(lapply(mut_vs_wt_list,function(x) x$mut))
all_wt<-unlist(lapply(mut_vs_wt_list,function(x) x$wt))

dndscv_by_vaf_and_mut_status<-lapply(list(mut=all_mut,wt=all_wt),function(sample_set) {
  dndscv_by_vaf<-lapply(new_VAF_groups[2:5],function(VAF_bin) {
    cat(VAF_bin,sep="\n")
    tissue_vaf_info<-all_df_tidy$NW%>%
      dplyr::filter(Sample%in%sample_set & VAF_group==VAF_bin)%>%
      separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
      mutate(pos=as.numeric(pos))%>%
      dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>% 
      dplyr::filter(!duplicated(.))%>% #only count each mutation once per individual in each VAF group
      arrange(sampleID,pos)
    
    mtdna.dndsout <- dndscv(tissue_vaf_info, gene_list=target_genes, 
                            refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
    return(mtdna.dndsout)
  })
  names(dndscv_by_vaf)<-new_VAF_groups[2:5]
  return(dndscv_by_vaf)
})

global_dnds_by_mut_status<-Map(list1=dndscv_by_vaf_and_mut_status,status=c("mut","wt"),function(list1,status) {
  
  Map(dndsout=list1,VAF_bin=names(list1),function(dndsout,VAF_bin) {
    dndsout$globaldnds%>%
      filter(!is.na(name) & complete.cases(.))%>%
      mutate(name=rename_vec[name])%>%
      mutate(VAF_group=VAF_bin)
  })%>%dplyr::bind_rows()%>%
    mutate(mut_status=status)
  
})%>%dplyr::bind_rows()

max_dnds_value<-4
stats_to_include=c("Missense","Truncating")
rename_mut_status_vec=c("MPN clonal expansion","Singleton")
names(rename_mut_status_vec)<-c("mut","wt")

nb.cols <- 2
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)

MPN_global_dnds_by_mut_status_plot<-global_dnds_by_mut_status%>%
  filter(VAF_group%in%VAF_groups$labels[2:6] & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         VAF_group=factor(VAF_group,levels=VAF_groups$labels[2:6]),
         mut_status=rename_mut_status_vec[mut_status])%>%
  ggplot(aes(x=VAF_group, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), size = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mycolors[1:2], name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', linewidth = 0.5) +
  theme_classic()+
  my_theme+
  facet_wrap(~mut_status,nrow=1)+
  theme(legend.position="top",axis.text.x=element_text(angle=90))

ggsave(filename = paste0(plots_dir,"MPN_global_dnds_by_mut_status_plot.pdf"),MPN_global_dnds_by_mut_status_plot,width=3.5,height=2.5)

# -------------MPN analysis, test gene-level dNdS by driver clone/ wild-type and VAF ---------------

selcv_res_by_mut_status<-Map(list1=dndscv_by_vaf_and_mut_status,status=c("mut","wt"),function(list1,status) {
  
  Map(dndsout=list1,VAF_bin=names(list1),function(dndsout,VAF_bin) {
    temp<-dndsout$sel_cv%>%
      dplyr::select(gene_name,wmis_cv,wnon_cv,qmis_cv,qtrunc_cv)
    
    dplyr::bind_rows(temp%>%dplyr::select(gene_name,"dNdS"=wmis_cv,"qval"=qmis_cv)%>%mutate(type="Missense"),
                     temp%>%dplyr::select(gene_name,"dNdS"=wnon_cv,"qval"=qtrunc_cv)%>%mutate(type="Nonsense"))%>%
      mutate(VAF_group=VAF_bin)
  })%>%dplyr::bind_rows()%>%
    mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels[2:5]),mut_status=status)
  
})%>%dplyr::bind_rows()

nb.cols <- 2
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)
gene_order=c(paste0("MT-ND",1:6),"MT-ND4L","MT-CYB",paste0("MT-CO",1:3),"MT-ATP6","MT-ATP8")
gene_order<-gene_order[gene_order%in%target_genes]

MPN_dnds_heatmap_by_mut_status_by_gene<-selcv_res_by_mut_status%>%
  mutate(dnds_if_signif=ifelse(qval<0.1,paste0(sprintf(dNdS, fmt = '%#.1f'),"*"),ifelse(dNdS>2,sprintf(dNdS, fmt = '%#.1f'),"")),
         dNdS=ifelse(dNdS>2,2,dNdS),
         gene_name=factor(gene_name,levels=rev(gene_order)),
         mut_status=rename_mut_status_vec[mut_status])%>%
  filter(VAF_group%in%VAF_groups$labels[2:5])%>%
  ggplot(aes(y=gene_name,x=VAF_group,fill=dNdS,label=dnds_if_signif))+
  geom_tile()+
  facet_grid(mut_status~type)+
  geom_text(col="white",size=1.5)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(n=11,name = "RdYlGn"))+
  labs(x="Variant allele fraction",y="")

ggsave(filename = paste0(plots_dir,"MPN_dnds_heatmap_by_mut_status_by_gene.pdf"),MPN_dnds_heatmap_by_mut_status_by_gene,width=3.5,height=3.5)


#-----------------------------------------------------------------------------------#
# ----------------------HETEROPLASMIC OOCYTE MUT ANALYSIS------------------------------
#-----------------------------------------------------------------------------------#

## Identify the heteroplasmic oocyte mutations -----
# This is based on the principle that mutations present across multiple unrelated clones (till you go back to development)

all_het_oocyte_mut_df=Map(dataset_mito_data=all_mito_datasets[all_cohorts],dataset=all_cohorts, function(dataset_mito_data,dataset) {
  cat(dataset,sep="\n")
  het_oocyte_mut_df<-Map(list=dataset_mito_data,Exp_ID=names(dataset_mito_data),function(list,Exp_ID){
    cat(Exp_ID,sep="\n")
    
    #Foetal trees are analysed separately.
    #Only analyse if ≥8 samples
    if(length(list$tree$tip.label)<8){stop(return(NULL))}
    
    if(Exp_ID%in%c("8 pcw","18 pcw","CB001","CB002")) {
      threshold_vaf<-0.025; threshold_molecular_time<-5; node_height_for_independence<-20
    } else {
      threshold_vaf<-0.25; threshold_molecular_time<-30; node_height_for_independence<-100
    }
    
    het_oocyte_muts<-detect_het_oocyte_mutation(matrices=list$matrices,tree=list$tree,vaf_cutoff=threshold_vaf)
    if(length(het_oocyte_muts)==0){stop(return(NULL))}
    
    het_oocyte_muts<-het_oocyte_muts[!het_oocyte_muts%in%exclude_muts & !grepl("DEL|INS",het_oocyte_muts)]
    if(length(het_oocyte_muts)==0){stop(return(NULL))}
    
    ho_vaf_mat<-list$matrices$vaf[het_oocyte_muts,list$tree$tip.label,drop=F]

    MRCA_mol_times<-sapply(1:nrow(ho_vaf_mat),function(i) {
      pos_samples<-colnames(ho_vaf_mat)[ho_vaf_mat[i,]>threshold_vaf]
      MRCA_node<-find_latest_acquisition_node(tree = list$tree,pos_samples = pos_samples)
      nodeheight(tree=list$tree,node=MRCA_node)
    })
    
    high_vaf_het_oocyte_muts<-het_oocyte_muts[MRCA_mol_times<threshold_molecular_time]
    ho_vaf_mat<-ho_vaf_mat[MRCA_mol_times<threshold_molecular_time,]
    
    if(length(high_vaf_het_oocyte_muts)==0){stop(return(NULL))}
    
    #Most likely oocyte VAF for these mutations
    if(dataset=="lymph"){
      samples_to_include<-list$tree$tip.label
    } else {
      post_dev_nodes<-get_expanded_clade_nodes(tree=list$tree,height_cut_off = node_height_for_independence,min_samples = 1,min_clonal_fraction = 0)
      samples_to_include<-unlist(lapply(post_dev_nodes$nodes,function(node) sample(size=1,x=getTips(list$tree,node))))
    }
    
    ml_vaf<-sapply(high_vaf_het_oocyte_muts,function(mut_ref) {
      mean(unlist(ho_vaf_mat[mut_ref,samples_to_include,drop=T]))
    })
    
    df<-data.frame(exp_ID=Exp_ID,
                   vaf_threshold_used_to_define_positive=threshold_vaf,
                   molecular_time_used_to_define_embryonic=threshold_molecular_time,
                   mut_ref=high_vaf_het_oocyte_muts,
                   n_pos=sapply(1:nrow(ho_vaf_mat),function(i) {length(colnames(ho_vaf_mat)[ho_vaf_mat[i,]>threshold_vaf])}),
                   MRCA_molecular_time=MRCA_mol_times[MRCA_mol_times<threshold_molecular_time],
                   ml_vaf=ml_vaf)
     
    return(df)
  })%>%dplyr::bind_rows()%>%
    mutate(dataset=dataset,.before=1)
  
  return(het_oocyte_mut_df)
})%>%dplyr::bind_rows()

individuals_assessed<-lapply(all_mito_datasets[all_cohorts],function(list) {
  dataset_out<-Map(list2=list,exp_ID=names(list),function(list2,exp_ID){
    if(length(list2$tree$tip.label)<8){return(NULL)} else {return(data.frame(exp_ID=exp_ID,n_samp=length(list2$tree$tip.label)))}
  })%>%dplyr::bind_rows()
  return(dataset_out)
})%>%dplyr::bind_rows()

all_het_oocyte_mut_df<-all_het_oocyte_mut_df%>%filter(!exp_ID%in%c("8 pcw", "18 pcw","CB001","CB002"))
individuals_assessed<-individuals_assessed%>%filter(!exp_ID%in%c("8 pcw", "18 pcw","CB001","CB002"))

n_samples_with_ho_mut_over_threshold<-sapply(seq(0.005,0.9,0.005),function(threshold){all_het_oocyte_mut_df%>%filter(ml_vaf>threshold)%>%distinct(exp_ID)%>%nrow()})

prop_of_samples_with_ho_mut<-data.frame(threshold=seq(0.005,0.9,0.005),n=n_samples_with_ho_mut_over_threshold)%>%
  ggplot(aes(x=threshold,y=n/nrow(individuals_assessed)))+
  geom_line()+
  theme_classic()+
  scale_y_continuous(limits=c(0,0.35))+
  scale_x_log10(limits=c(0.005,1),breaks=c(0.005,0.01,0.02,0.05,0.1,0.2,0.5))+
  my_theme+
  labs(x="Heteroplasmy level",y="Proportion of individuals with\n at least one heteroplasmic oocyte\nmutation above threshold")

ggsave(filename = paste0(plots_dir,"prop_of_samples_with_ho_mut.pdf"),prop_of_samples_with_ho_mut,width=2.5,height=1.8)

library(dndscv)
library(mitovizR)

all_mtDNA_genes <- c("MT-CYB", "MT-ND5", "MT-ND2", "MT-ND4", "MT-ND1", "MT-CO3", "MT-ATP6","MT-ND3", "MT-ATP8", "MT-ND4L", "MT-CO2", "MT-CO1","MT-ND5","MT-ND6")
mtref_rda_path=ifelse(Sys.info()['sysname']=="Darwin",paste0(root_dir,"data/mtref.rda"),"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/mtref.rda")

# read in the mtdna variant file and remove patient id
mtdna.variant.data <- all_het_oocyte_mut_df%>%
  separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
  mutate(pos=as.numeric(pos))%>%
  dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>% # next, rename/reorder the columns to make them compatible with dnds input format
  dplyr::filter(!duplicated(.))%>% # remove duplicated variants
  arrange(sampleID,pos)

mtdna.dndsout <- dndscv(mtdna.variant.data, gene_list=all_mtDNA_genes, 
                        refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)

all_het_oocyte_mut_df_annotated<-all_het_oocyte_mut_df%>%
  tidyr::separate(mut_ref,into=c("chr","pos","ref","mut"),sep = "_")%>%
  mutate(pos=as.numeric(pos))%>%
  left_join(mtdna.dndsout$annotmuts,by=c("chr","pos","ref","mut"),relationship = "many-to-many")%>%
  replace_na(replace=list(impact="Non-coding"))

#Create a coordinate reference for the different regions----
mito_domain_convert<-c("tRNA","rRNA","Coding","D-loop")
names(mito_domain_convert)<-c("trna","rrna","cds","reg")
mito_coords_reference_df<-mitovizR:::mito_df()%>%
  dplyr::select(type,ymin,ymax)%>%
  dplyr::mutate(Mutation_type=mito_domain_convert[type])

all_het_oocyte_mut_df_annotated$Mutation_type<-sapply(all_het_oocyte_mut_df_annotated$pos,
                                                      function(pos) {mito_coords_reference_df$Mutation_type[mito_coords_reference_df$ymin<pos & mito_coords_reference_df$ymax>=pos]})

het_oocyte_mut_type_summary<-all_het_oocyte_mut_df_annotated%>%
  dplyr::filter(!grepl("\\*",mut))%>%
  dplyr::mutate(final_type=ifelse(is.na(impact)|impact=="Non-coding",Mutation_type,impact),
                cat=ifelse(ml_vaf>=0.01,"Heteroplasmic\noocyte\n(VAF≥1%)","Heteroplasmic\noocyte\n(VAF<1%)"))%>%
  group_by(cat,final_type)%>%
  dplyr::summarise(n=n())%>%
  dplyr::mutate(prop=n/sum(n))


#Plot barplot showing number of het oocyte muts per individual----
n_of_het_oocyte_muts_plot<-all_het_oocyte_mut_df_annotated%>%
  dplyr::filter(!grepl("\\*",mut) & ml_vaf>0.01)%>%
  group_by(exp_ID)%>%
  dplyr::summarise(nmut=n())%>%
  right_join(individuals_assessed)%>%
  replace_na(list(nmut=0))%>%
  dplyr::group_by(nmut)%>%
  dplyr::summarise(n_samples=n())%>%
  ggplot(aes(x=nmut,y=n_samples))+
  geom_bar(stat="identity",col="black",linewidth=0.2,fill="lightblue")+
  theme_classic()+
  my_theme+
  labs(x="Number of heteroplasmic oocyte mutations\nwith VAF inferred >1%",y="Number of individuals")

ggsave(filename = paste0(plots_dir,"n_of_het_oocyte_muts_plot.pdf"),n_of_het_oocyte_muts_plot,width=2,height=1.8)


#Compare mutation distribution with full somatic mutation set----
complete.annotated.mutation.table$Mutation_type<-sapply(complete.annotated.mutation.table$pos,
       function(pos) {mito_coords_reference_df$Mutation_type[mito_coords_reference_df$ymin<pos & mito_coords_reference_df$ymax>=pos]})

Mut_type_levels=c("D-loop","Synonymous","tRNA","rRNA","Missense","Stop_loss","Truncating","Inconsistent\nannotation")
Mut_cat_comparison<-het_oocyte_mut_type_summary%>%
  dplyr::bind_rows(complete.annotated.mutation.table%>%
                     dplyr::mutate(final_type=ifelse(impact=="Non-Coding"&Mutation_type=="Coding","Inconsistent\nannotation",ifelse(is.na(impact)|impact=="Non-Coding",Mutation_type,impact)))%>%
                     group_by(final_type)%>%
                     dplyr::summarise(n=n())%>%
                     dplyr::mutate(prop=n/sum(n),cat="Somatic"))%>%
  mutate(final_type=factor(final_type,levels=Mut_type_levels))%>%
  ggplot(aes(x=cat,y=prop,fill=forcats::fct_rev(final_type)))+
  geom_bar(stat="identity",position = "fill",col="black",linewidth=0.2)+
  theme_classic()+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  scale_fill_brewer(palette = "Set2",direction = -1)+
  my_theme+
  labs(fill="Mutation\ncategory",x="",y="Proportion")

het_oocyte_mut_type_summary%>%
  dplyr::bind_rows(complete.annotated.mutation.table%>%
                     dplyr::mutate(final_type=ifelse(impact=="Non-Coding"&Mutation_type=="Coding","Inconsistent\nannotation",ifelse(is.na(impact)|impact=="Non-Coding",Mutation_type,impact)))%>%
                     group_by(final_type)%>%
                     dplyr::summarise(n=n())%>%
                     dplyr::mutate(prop=n/sum(n),cat="Somatic"))%>%
  group_by(cat)%>%dplyr::summarise(n=sum(n))

ggsave(filename = paste0(plots_dir,"Mut_cat_comparison.pdf"),Mut_cat_comparison,width=3,height=2.5)


# Save plot of het oocyte mutations by genome position----
mtDNA_mut_pos<-plot_df(all_het_oocyte_mut_df_annotated%>%
                         dplyr::filter(!grepl("\\*",mut))%>%
          dplyr::mutate(SAMPLE="all")%>%
          dplyr::select(SAMPLE,pos,ref,mut,"HF"=ml_vaf),
          pos_col = "pos", ref_col = "ref", alt_col = "mut")

ggsave(filename = paste0(plots_dir,"mtDNA_mut_pos.pdf"),mtDNA_mut_pos,width=5,height=5)

#-----------------------------------------------------------------------------------#
## Do example clone inferences on an individuals with heteroplasmic oocyte muts-----
#-----------------------------------------------------------------------------------#

library(Seurat)
library(ComplexHeatmap)
library(BuenColors)

#Define function to get clusters
seuratSNN <- function(matSVD, resolution = 1, k.param = 10){ 
  set.seed(1)
  rownames(matSVD) <- make.unique(rownames(matSVD))
  obj <- FindNeighbors(matSVD, k.param = k.param, annoy.metric = "cosine")
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}
root_dir="/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood"
output.dir <- paste0(root_dir,"/mito_mut_clones")
dir.create(output.dir,showWarnings = F)

# store the information for the heatmap here
df.list <- list()
matrix.list <- list()

# iterate over multiple patients here (if there are multiple)
patient.info<-data.frame(dataset=c("NW","PR","PR","LM"),patient.id=c("PD5182","PD44887","PD44890","PD41857"))
for (i in 1:nrow(patient.info)){
  
  # read in one of the mgatk SE objects
  patient.tmp <- patient.info$patient.id[i]
  patient.dataset <- patient.info$dataset[i]
  cat(paste0("\nProcessing patient...", patient.tmp, "\n"))
  this.patient.info<-all_mito_datasets[[patient.dataset]][[patient.tmp]]
  vaf.mtx <- this.patient.info$matrices$vaf.filt
  
  #Remove mutations present in all samples & keep only those present at global VAF >0.5%
  if(is.null(this.patient.info$matrices$vaf$global)){
    this.patient.info$matrices$vaf$global<-rowSums(this.patient.info$matrices$NV)/rowSums(this.patient.info$matrices$NR)
  }
  high_global_vaf_muts<-rownames(this.patient.info$matrices$vaf)[this.patient.info$matrices$vaf$global>0.005]
  vaf.mtx<-vaf.mtx[rownames(vaf.mtx)%in%high_global_vaf_muts &
                     !grepl("INS|DEL",rownames(vaf.mtx)),]
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
  ## Define clones based on the mitochondrial mutations -------------------------------
  #-----------------------------------------------------------------------------------#
  
  # get clusters with cluster resolution 10 and knn 50
  # clustering fails for PD41857 with these parameters, so empirically altered to optimize
  #if(patient.tmp=="PD41857"){
    clusters <- seuratSNN((t(vaf.filtered.mtx)),resolution= 5,k.param=3)
  # } else {
  #   clusters <- seuratSNN((t(vaf.filtered.mtx)),resolution= 10,k.param=5)
  # }
  
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
  # vec_go[min_mut_cluster]<-"lightgrey"
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
  pdf(paste0(plots_dir, patient.tmp, "_mito_mutation_heatmap.pdf"), width=15, height=8)
  print(hm) 
  dev.off()
  
  #-----------------------------------------------------------------------------------#
  ## Create "trees" aka dendrograms ---------------------------------------------------
  #-----------------------------------------------------------------------------------#
  
  # Get group means 
  matty <- sapply(names(vec_go), function(cluster){
    cells <- df %>% dplyr::filter(cluster_id == cluster) %>% pull(sample_id) %>% as.character()
    Matrix::rowMeans(sqrt(afp[,cells]))
  })
  
  if(length(groups) > 2){
    
    # Do cosine distance; note that we used sqrt transformation already when creating the pseudo bulk-cell matrix
    mito.hc <- hclust(dist(lsa::cosine((matty))))
    
    pdf(paste0(plots_dir, patient.tmp, "_hierarchical_tree.pdf"), width = 5, height = 5)
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

for(j in 1:nrow(patient.info)) {
  
  patient.tmp <- patient.info$patient.id[j]
  patient.dataset <- patient.info$dataset[j]
  
  this.patient.info<-all_mito_datasets[[patient.dataset]][[patient.tmp]]
  exp_clones<-read.delim(paste0(root_dir,"/mito_mut_clones/",patient.tmp,"_mtdna_clone_assignment.txt"))
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
  
  pdf(file=paste0(plots_dir,"Inferred_clones_",patient.tmp,".pdf"),width = 7,height=3)
  par(mfrow=c(1,1))
  this.patient.info$tree=plot_tree(tree = this.patient.info$tree,cex.label = 0)
  clone_hm<-matrix(NA,nrow=1,ncol=length(this.patient.info$tree$tip.label),dimnames = list("Clones",this.patient.info$tree$tip.label))
  for(i in 1:nrow(exp_clones)){
    if(exp_clones$sample_id[i]%in%this.patient.info$tree$tip.label){
      clone_hm[1,exp_clones$sample_id[i]]<-exp_cluster_cols[as.character(exp_clones$cluster_id[i])]
    }
  }
  
  add_heatmap(tree=this.patient.info$tree,heatmap=clone_hm,border="gray",cex.label = 1)
  legend("topleft", inset=c(.01,.01), title="Clone no.",
         names(exp_cluster_cols), fill=exp_cluster_cols, horiz=F, cex=0.7,ncol=1)
  dev.off()
  this.patient.info$clusters<-exp_clones
  this.patient.info$cluster_cols<-exp_cluster_cols
  this.patient.info$cluster_hm<-clone_hm
}
  






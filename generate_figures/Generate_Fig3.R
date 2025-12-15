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

##----------------------LOOK FOR LOW LEVEL HETEROPLASMIC MUTATIONS-------------------

#Analysis only in the foetal/ cord blood phylogenies
young_IDs=c("8pcw","18pcw","CB001","CB002")

#Add a "stats_df" object that gives the phylogenetic signal/ beta-binomial values for each shared mutation
VAF_measure="vaf"
mito_data_young<-Map(list=mito_data[young_IDs],exp_ID=young_IDs,function(list,exp_ID) {
  print(exp_ID)
  
  list$matrices[[VAF_measure]]<-apply(list$matrices[[VAF_measure]],2,as.numeric)
  
  #If phylogenetic signal data is already embedded, use this. Otherwise just add for shared mutations.
  if(is.null(list$phylosignal)) {
    list$stats_df<-data.frame(exp_ID=exp_ID,
                                   rho_val=list$rho_vals,
                                   n_pos=rowSums(list$matrices$SW),
                                   mean_vaf=rowMeans(list$matrices[[VAF_measure]]),
                                   max_vaf=apply(list$matrices[[VAF_measure]],1,max))%>%
      tibble::rownames_to_column(var="mut_ref")
    
    print(nrow(list$stats_df))
    
    tree.noancestral<-drop.tip(list$tree.ultra,"Ancestral")
    tree.4d<-phylobase::phylo4d(tree.noancestral,tip.data=t(list$matrices[[VAF_measure]][list$stats_df$mut_ref,tree.noancestral$tip.label]))
    phylosignal<-phyloSignal(tree.4d)
    
    list$stats_df$Cmean_pvalue<-phylosignal$pvalue$Cmean
  } else {
    
    list$stats_df<-data.frame(exp_ID=exp_ID,
                                   Cmean_pvalue=list$phylosignal$pvalue$Cmean,
                                   rho_val=list$rho_vals,
                                   n_pos=rowSums(list$matrices$SW),
                                   mean_vaf=rowMeans(list$matrices[[VAF_measure]]),
                                   max_vaf=apply(list$matrices[[VAF_measure]],1,max))%>%
      tibble::rownames_to_column(var="mut_ref")%>%
      dplyr::filter(n_pos>1)
  }
  return(list)
})

#-----------------------------------------------------------------------------------#
### Generate FIG. 3A-D ---------
#-----------------------------------------------------------------------------------#

#Now pull out the shared mutations for which the MRCA is the tree root and visualize
resave_plots=T
het_oocyte_muts<-lapply(young_IDs,function(Exp_ID) {
  low_level_het_muts<-mito_data_young[[Exp_ID]]$stats_df%>%
    dplyr::bind_rows()%>%
    dplyr::filter(n_pos>=2 & !grepl("INS|DEL",mut_ref))%>%
    mutate(Cmean_qvalue=p.adjust(Cmean_pvalue,method = "BH"))%>%
    dplyr::filter(!mut_ref%in%exclude_muts)%>%
    dplyr::filter(max_vaf>0.01 & Cmean_pvalue<0.05)%>%
    dplyr::filter(!(mean_vaf>0.003&Cmean_pvalue>=0.005&rho_val<=5e-3))%>% #Final filter to take out a few remaining artefacts that are present at fairly high global vaf but are not phylocorrelated and have low dispersion
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
  
  #Define the root of the tree - though as the 18pcw tree is so assymmetric, accept the MRCA of the major clade
  tree_root=ifelse(Exp_ID=="18pcw",236,getRoot(tree.no.ancestral))
  MRCA_is_root_muts<-low_level_het_muts[mut_MRCA==tree_root]
  
  if(resave_plots) {
    pdf(file = paste0(plots_dir,"Figure_03/Low_level_hetmut_",Exp_ID,"_phylos.pdf"),width = 18,height=5)
    par(mfrow=c(2,6))
    temp=lapply(MRCA_is_root_muts,function(mut){
      plot_tree(mito_data[[Exp_ID]]$tree.ultra,cex.label = 0,bars = mito_data[[Exp_ID]]$matrices$vaf[mut,],plot_axis = T,title = mut)
      text(x = 0,
           y=-0.05*par()[['yaxp']][2],
           cex = 0.75,
           font=3,
           col="#00000095",
           paste0("Max VAF: ",round(max(mito_data[[Exp_ID]]$matrices$vaf[mut,]),digits = 3),"; Mean VAF: ",round(mito_data_young[[Exp_ID]]$stats_df$mean_vaf[mito_data_young[[Exp_ID]]$stats_df$mut==mut],digits=3),"; Cmean pval = ",mito_data_young[[Exp_ID]]$stats_df$Cmean_pval[mito_data_young[[Exp_ID]]$stats_df$mut==mut]),pos = 4)
    })
    dev.off()
  }
  
  return(MRCA_is_root_muts)
})

names(het_oocyte_muts)<-young_IDs
saveRDS(het_oocyte_muts,paste0(root_dir,"/data/het_oocyte_muts.Rds"))

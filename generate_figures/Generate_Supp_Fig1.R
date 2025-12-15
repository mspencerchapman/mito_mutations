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

#-----------------------------------------------------------------------------------#
#----------------MITOCHONDRIAL COVERAGE STATISTICS----------------
#-----------------------------------------------------------------------------------#

#Generate the coverage histograms
coverage.plot<-mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  ggplot(aes(x=bedtools_mtDNA_coverage))+
  geom_histogram(fill="lightblue",col="black",linewidth=0.2)+
  scale_x_log10(labels=scales::label_comma())+
  facet_wrap(~exp_ID,scales = "free_y")+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = margin(1,0,1,0, "mm")))+
  labs(x="Mean mitochondrial DNA coverage",y="Count")
ggsave(filename=paste0(plots_dir,"Supp_Figure_01/Coverage_mean_plot.pdf"),coverage.plot,width=4,height=3)

median.coverage.plot<-mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  ggplot(aes(x=pilup_median_mtDNA_coverage))+
  geom_histogram(fill="lightblue",col="black",linewidth=0.2)+
  scale_x_log10(labels=scales::label_comma())+
  facet_wrap(~exp_ID,scales = "free_y",nrow=2)+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = margin(1,0,1,0, "mm")))+
  labs(x="Median mitochondrial DNA coverage",y="Count")
ggsave(filename=paste0(plots_dir,"Supp_Figure_01/Coverage_median_plot.pdf"),median.coverage.plot,width=7,height=2)

### Generate SUPPLEMENTARY FIG. 1A ---------
median.coverage.ridges.plot<-mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  ggplot(aes(x=pilup_median_mtDNA_coverage,y=exp_ID))+
  ggridges::geom_density_ridges(linewidth=0.2,fill="lightblue")+
  scale_x_log10(labels=scales::label_comma())+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = margin(1,0,1,0, "mm")))+
  labs(x="Median mitochondrial DNA coverage",y="")
ggsave(filename=paste0(plots_dir,"Supp_Figure_01/median.coverage.ridges.plot.pdf"),median.coverage.ridges.plot,width=3.3,height=2.5)

### Generate SUPPLEMENTARY FIG. 1B ---------
uniformity_ridges.plot<-mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  ggplot(aes(x=perc_over0.8,y=exp_ID))+
  ggridges::geom_density_ridges(linewidth=0.2,fill="lightblue")+
  scale_x_log10()+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = margin(1,0,1,0, "mm")))+
  labs(x="Coverage uniformity\n(proportion of mtDNA genome with coverage >80% of mean)",y="")
ggsave(filename=paste0(plots_dir,"Supp_Figure_01/Uniformity_perc80_ridges_plot.pdf"),uniformity_ridges.plot,width=3.3,height=2.5)

#Print coverage summary statistics
mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  group_by(exp_ID)%>%
  summarise(samples=n(),mean_mtDNA_coverage=mean(bedtools_mtDNA_coverage),Under_1000X=sum(bedtools_mtDNA_coverage<1000))

mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  summarise(mean_mtDNA_coverage=mean(bedtools_mtDNA_coverage))
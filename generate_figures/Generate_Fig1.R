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
#----------------MITOCHONDRIAL COPY NUMBER STATISTICS----------------
#-----------------------------------------------------------------------------------#

#Vector to convert the foetal cell type codes into basic 'HSC'/ 'Progenitor' types
convert_vec=c("HSC","HSC","HPC","HPC","HPC","HPC")
names(convert_vec)=c("H","HSC","C","M","HSPC","Progenitor")


Individual_means_df<-mito_cn%>%
  mutate(Cell_type=convert_vec[Cell_type])%>%
  filter(Sample%in%unlist(lapply(mito_data,function(list) list$tree$tip.label))&
           !is.na(bedtools_mtDNA_genomes)&
           #Cell_type=="HSC"&
           !is.na(exp_ID))%>%
  group_by(exp_ID)%>% 
  summarise(mean=mean(bedtools_mtDNA_genomes,na.rm=T),median=median(bedtools_mtDNA_genomes,na.rm=T),sd=sd(bedtools_mtDNA_genomes,na.rm=T),IQR_lower=quantile(bedtools_mtDNA_genomes,0.25),IQR_upper=quantile(bedtools_mtDNA_genomes,0.75))
Individual_means_df$exp_ID<-factor(Individual_means_df$exp_ID,levels=ref_df$Sample[order(ref_df$Age)])

studies_to_include=ref_df$Canapps.project

#-----------------------------------------------------------------------------------#
### Generate FIG. 1A ---------
#-----------------------------------------------------------------------------------#

#Colour points by the colony phenotype
phenotype_cols=c("red",RColorBrewer::brewer.pal(n=length(pheno_levels)-1,name = "Set2"),"darkgray") #Make the erythroid red, and 'unknown' dark grey
names(phenotype_cols)<-c(pheno_levels,"Unknown")
mtDNA.copy.number.logscale.by.phenotype<-mito_cn%>%
  mutate(Cell_type=convert_vec[Cell_type])%>%
  filter(Sample%in%unlist(lapply(mito_data,function(list) list$tree$tip.label))&
           !is.na(bedtools_mtDNA_genomes)&
           !is.na(exp_ID))%>%
  replace_na(replace=list(Phenotype="Unknown"))%>%
  ggplot(aes(x=factor(Sample,levels=mito_cn%>%dplyr::filter(Study%in%studies_to_include)%>%arrange(bedtools_mtDNA_genomes)%>%pull(Sample)%>%unique()),
             y=bedtools_mtDNA_genomes,
             col=Phenotype))+
  geom_point(size = 0.3, alpha=0.7,stroke = 0, shape = 16)+
  scale_y_log10(limits=c(8,20000),breaks=c(10,100,1000,10000))+
  scale_color_manual(values=phenotype_cols)+
  theme_classic()+
  facet_grid(cols=vars(factor(exp_ID,levels = ref_df$Sample[order(ref_df$Age)])),scales="free",space = "free")+
  geom_hline(aes(yintercept=median),linetype=2,data=Individual_means_df,col="red",linewidth=0.5)+
  geom_text(aes(x=100,y=mean,label=paste0("tilde(x) == ",round(median))),nudge_y=+0.1,size=2,data=Individual_means_df,parse=T,inherit.aes = F)+
  labs(x="Individual",y="Mitochondrial copy number",col="Colony\nphenotype")+
  my_theme+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size=8),
        legend.key.height = unit(2,"mm"),
        strip.text.x = element_text(size=6),
        panel.spacing.x=unit(1, "mm"))+
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave(filename = paste0(plots_dir,"Figure_01/Fig1a.mtDNA.copy.number.logscale.by.phenotype.pdf"),plot=mtDNA.copy.number.logscale.by.phenotype,height=2,width=7)

mtDNA_by_pheno_and_ID_summary<-mito_cn%>%
  filter(!is.na(Phenotype) & Cell_type=="HSC")%>%
  group_by(exp_ID,Phenotype)%>%
  summarise(n=n(),
            mean_mtDNA_genomes=mean(bedtools_mtDNA_genomes),
            median_mtDNA_genomes=median(bedtools_mtDNA_genomes))

mtDNA_by_pheno_summary<-mito_cn%>%
  filter(!is.na(Phenotype) & Cell_type=="HSC")%>%
  group_by(Phenotype)%>%
  summarise(n=n(),
            mean_mtDNA_genomes=mean(bedtools_mtDNA_genomes),
            median_mtDNA_genomes=median(bedtools_mtDNA_genomes))

#Review the average copy number of colonies with different phenotypes relative to the average copy number of erythroid colonies
pheno_levels=c("Ery","EryMy","Gran","MyGran","MyMono","Mono","NKMy")
mtDNA_by_pheno_and_ID_summary%>%
  dplyr::select(-n)%>%
  pivot_wider(id_cols=c("exp_ID"),names_from="Phenotype",values_from="median_mtDNA_genomes")%>%
  mutate(across(all_of(pheno_levels),.fns = function(x) x/Ery))

#-----------------------------------------------------------------------------------#
### Generate FIG. 1B ---------
#-----------------------------------------------------------------------------------#

mito_cn_by_pheno<-mito_cn%>%
  filter(!is.na(Phenotype) & Cell_type=="HSC")%>%
  dplyr::select(exp_ID,Sample,bedtools_mtDNA_genomes,Phenotype)%>%
  ggplot(aes(x=factor(Phenotype,levels=pheno_levels),y=bedtools_mtDNA_genomes))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(col=factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])),width=0.2,height=0,alpha=0.3,size=0.025)+
  scale_y_log10(limits=c(50,3000))+
  theme_bw()+
  scale_color_manual(values = Individual_cols[-c(1,2)])+
  guides(col=guide_legend(override.aes = list(alpha=1,size=0.5)))+
  labs(x="Colony cell composition",y="mtDNA copies per cell",col="Individual")+
  geom_text(data=mtDNA_by_pheno_summary,aes(label=paste0("µ = ",round(mean_mtDNA_genomes)),x=Phenotype,y=75),size=1.6,inherit.aes = F)+
  my_theme+
  theme(#axis.text.x = element_text(angle=90),
    strip.text.x = element_text(size=7,margin = unit(c(1,0,1,0),"mm")),
    legend.key.size = unit(0,"mm"),
    legend.position="right")
ggsave(filename=paste0(plots_dir,"Figure_01/Fig1b.mito_cn_by_pheno.pdf"),mito_cn_by_pheno,width=4,height=2)

#-----------------------------------------------------------------------------------#
# ----------------------ANALYSIS OF MITOCHONDRIAL MUTATION BURDENS------------------
#-----------------------------------------------------------------------------------#

## Generate tidy data frame of samples, mutations and vafs-------------------------------
# This combines various pieces of information to create a filtered mutations set
# 1. If the mutation is in the 'copy number correlated' list, and the inferred mtDNA copy number is <25 -> remove
# 2. Only keep if it has passed shearwater
# 3. Only keep if the most likely mutational signature is the 'N1' signature
# 4. Remove insertions and deletions

mutCN_cutoff=25 #If the mitochondrial copy number is over 25, retain mutation even if is in the "CN correlating muts" list
df_tidy<-dplyr::bind_rows(Map(list=mito_data,exp_ID=names(mito_data),function(list,exp_ID) {
  if(is.null(list)){stop(return(NULL))}
  CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
  implied_mut_CN_tidy<-list$matrices$implied_mutCN%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var="mut_ref")%>%
    dplyr::select(-global)%>%
    tidyr::gather(key="Sample",value="implied_mut_CN",-mut_ref)
  
  df_tidy<-(list$matrices$vaf*list$matrices$SW*(list$matrices$ML_Sig=="N1")*CN_correlating_mut_removal_mat)%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var="mut_ref")%>%
    mutate(rho_val=list$rho_vals)%>%
    dplyr::select(-global)%>%
    tidyr::gather(key="Sample",value="vaf",-mut_ref,-rho_val)
  
  comb_tidy<-left_join(df_tidy,implied_mut_CN_tidy,by=c("Sample","mut_ref"))%>%
    dplyr::filter(!grepl("DEL|INS",mut_ref) & !mut_ref%in%exclude_muts)%>% #Exclude indels, and the 'black listed' mutations
    dplyr::filter(vaf>0)%>%
    mutate(exp_ID=exp_ID)
  return(comb_tidy)
  }))

#Filter the CN-correlating mutations - due to mis-mapping of nuclear reads.
#There may be some genuine mutations at these sites, in which case the implied mutation copy number will be much higher than ~2 (here an arbitrary threshold of 8 is applied)
df_tidy<-df_tidy%>%
  left_join(mito_cn%>%mutate(exp_ID=gsub("8pcw","8 pcw",exp_ID)))%>%
  mutate(implied_mut_CN=vaf*bedtools_mtDNA_genomes)%>%
  dplyr::filter(!(mut_ref%in%CN_correlating_muts & implied_mut_CN<=mutCN_cutoff))


## ----------------------Sum of VAF mutation burden measures----------------------

sum_of_vaf_df<-Map(list=mito_data,Exp_ID=names(mito_data),function(list,Exp_ID){
  #sum the vaf of mutations that:
  #(1) Pass shearwater,
  #(2) are most likely signature N1 (the real mutatation signature),
  #(3) are not indels, specific artefacts or 'copy number correlating mutations' (i.e. most likely artefacts from mismapping of nuclear DNA or low-level contamination)
  CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
  df<-as.data.frame(colSums((list$matrices$vaf*list$matrices$SW*(list$matrices$ML_Sig=="N1")*(CN_correlating_mut_removal_mat))[!grepl("DEL|INS",rownames(list$matrices$vaf))&
                                                                                                !rownames(list$matrices$vaf)%in%exclude_muts,],na.rm = T))%>%
    tibble::rownames_to_column(var="Sample")%>%
    mutate(exp_ID=Exp_ID)%>%
    dplyr::rename(sum_of_vaf=2)%>%
    dplyr::filter(Sample!="global")
  return(df)
})%>%
  dplyr::bind_rows()%>%
  left_join(mito_cn%>%dplyr::select(Sample,Cell_type,Phenotype),by="Sample")%>%
  dplyr::mutate(Cell_type=convert_vec[Cell_type])

sum_of_vaf_summary<-sum_of_vaf_df%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
  group_by(exp_ID)%>%
  summarise(mean=mean(sum_of_vaf,na.rm=T),median=median(sum_of_vaf,na.rm = T))


sum_of_vaf_plot_log<-sum_of_vaf_df%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
  mutate(sum_of_vaf=ifelse(sum_of_vaf==0,0.001,sum_of_vaf))%>%
  ggplot(aes(x=forcats::fct_reorder(Sample,sum_of_vaf),y=sum_of_vaf))+
  geom_point(size=0.1)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.text.x = element_text(size=6))+
  facet_grid(cols=vars(factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])),scales = "free",space="free")+
  scale_y_log10(breaks=c(0.005,0.01,0.03,0.06,0.12,0.25,0.5,1,2,4,8))+
  geom_hline(aes(yintercept=median),linetype=2,data=sum_of_vaf_summary%>%mutate(median=ifelse(median==0,0.001,median)),col="red",linewidth=0.5)+
  geom_text(aes(x=0,y=median_pos,label=paste0("tilde(x) == ",round(median,3))),size=2,nudge_x=+115,nudge_y=+0.5,data=sum_of_vaf_summary%>%mutate(median_pos=ifelse(median==0,0.001,median)),parse=T)+
  labs(x="Sample",y="Mutation burden\n(sum of VAF)")+
  theme(panel.spacing.x=unit(1, "mm"))
ggsave(filename = paste0(plots_dir,"Figure_01/sum_of_vaf_log.pdf"),plot=sum_of_vaf_plot_log,height=2,width=7)

#-----------------------------------------------------------------------------------#
## Generate FIG. 1C ---------
#-----------------------------------------------------------------------------------#

#Now coloured by cell type
sum_of_vaf_plot_log_by_celltype<-sum_of_vaf_df%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
  mutate(sum_of_vaf=ifelse(sum_of_vaf==0,0.001,sum_of_vaf))%>%
  ggplot(aes(x=forcats::fct_reorder(Sample,sum_of_vaf),y=sum_of_vaf))+
  geom_point(aes(col=Cell_type),size=0.1,alpha=0.6)+
  scale_color_manual(values=c("#1a80bb", "#ea801c"))+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.text.x = element_text(size=6))+
  facet_grid(cols=vars(factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])),scales = "free",space="free")+
  scale_y_log10(breaks=c(0.005,0.01,0.03,0.06,0.12,0.25,0.5,1,2,4,8))+
  geom_hline(aes(yintercept=median),linetype=2,data=sum_of_vaf_summary%>%mutate(median=ifelse(median==0,0.001,median)),col="red",size=0.5)+
  geom_text(aes(x=0,y=median_pos,label=paste0("tilde(x) == ",round(median,3))),size=2,nudge_x=+115,nudge_y=+0.5,data=sum_of_vaf_summary%>%mutate(median_pos=ifelse(median==0,0.001,median)),parse=T)+
  labs(x="Sample",y="Mutation burden\n(sum of VAF)",col="Cell type")+
  guides(color=guide_legend(override.aes = list(size=0.4,alpha=1)))+
  theme(panel.spacing.x=unit(1, "mm"),legend.position = "bottom")
ggsave(filename = paste0(plots_dir,"Figure_01/Fig1c.sum_of_vaf_log_by_celltype.pdf"),plot=sum_of_vaf_plot_log_by_celltype,height=2.2,width=7)

#Comparison of mutaiton burden by cell type
individuals_with_both=c("18pcw","CB002","SX001","AX001","KX004")
sum_of_vaf_by_celltype<-sum_of_vaf_df%>%
  filter(exp_ID%in%individuals_with_both)%>%
  ggplot(aes(x=Cell_type,y=sum_of_vaf))+
  geom_violin(aes(fill=Cell_type),alpha=0.2)+
  scale_color_manual(values=c("#1a80bb", "#ea801c"))+
  scale_fill_manual(values=c("#1a80bb", "#ea801c"))+
  stat_compare_means(method = "wilcox.test",label = "p.format",fontface="italic",vjust = +1,size=2.6)+
  geom_jitter(aes(col=Cell_type),width=0.2,alpha=0.5)+
  facet_grid(cols=vars(factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])),scales = "free",space="free")+
  theme_classic()+
  my_theme+
  labs(y="Mutation burden\n(sum of VAF)")+
  theme(axis.title.x=element_blank(),legend.position = "none")
ggsave(filename = paste0(rebuttal_figs_dir,"sum_of_vaf_by_celltype_comparison.pdf"),plot=sum_of_vaf_by_celltype,height=2.2,width=7)

#-----------------------------------------------------------------------------------#
## Generate FIG. 1D ---------
#-----------------------------------------------------------------------------------#

mean_sum_of_vaf_by_age<-sum_of_vaf_summary%>%
  left_join(ref_df,by=c("exp_ID"="Sample"))%>%
  ggplot(aes(x=Age,y=mean,col=factor(exp_ID,levels = ref_df$Sample[order(ref_df$Age)])))+
  geom_point()+
  theme_bw()+
  scale_color_brewer(palette="Paired")+
  labs(col="Individual",y="Mean mutation burden (sum of VAF)")+
  geom_smooth(aes(x=Age,y=mean),col="black",linewidth=0.5,method="lm",inherit.aes = F)+
  my_theme+
  theme(legend.key.height = unit(3,"mm"),legend.title=element_text(size=7),legend.text=element_text(size=5))
ggsave(filename=paste0(plots_dir,"Figure_01/Fig1d.mean_sum_of_vaf_by_age.pdf"),plot=mean_sum_of_vaf_by_age,height=2,width=3)

lm.mean_mutburden_by_age=lm(mean~Age,data=sum_of_vaf_summary%>%
                              left_join(ref_df,by=c("exp_ID"="Sample")))
summary(lm.mean_mutburden_by_age)$coefficients
summary(lm.mean_mutburden_by_age)$r.squared
confint(lm.mean_mutburden_by_age)


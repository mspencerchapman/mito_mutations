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

#geom_jitter() displaces points at random, so fix the seed to make the saved
#panels reproducible between runs.
set.seed(42)

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"/data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
#plots_dir, rebuttal_figs_dir and my_theme all come from config.R

#Create the figure output directories if they do not already exist
for(d in c("Extended_Data_Figure_05","Figure_02")) dir.create(paste0(plots_dir,d),showWarnings=FALSE,recursive=TRUE)

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

#Vector to convert the foetal cell type codes into basic 'HSC'/ 'Progenitor' types
convert_vec=c("HSC","HSC","HPC","HPC","HPC","HPC")
names(convert_vec)=c("H","HSC","C","M","HSPC","Progenitor")

#-----------------------------------------------------------------------------------#
## ----------------------ANALYSIS OF MITOCHONDRIAL MUTATION BURDENS----------------------
#-----------------------------------------------------------------------------------#

## Generate tidy data frame of samples, mutations and vafs-------------------------------

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

#-----------------------------------------------------------------------------------#
# ----------------------Sum of VAF mutation burden measures----------------------
#-----------------------------------------------------------------------------------#

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

#-----------------------------------------------------------------------------------#
# ------------------------------SELECTION ANALYSIS----------------------------------
#-----------------------------------------------------------------------------------#

mtref_rda_path=ifelse(Sys.info()['sysname']=="Darwin",paste0(root_dir,"/data/mtref.rda"),"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/mtref.rda")
input.dir <- paste0(root_dir,"/dnds_tables/")

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
# (therefore the mutational signature normalization will not work)
#Maps the dndscv globaldnds statistic names onto readable labels
rename_vec=c("Missense","Nonsense","Truncating","Overall")
names(rename_vec)=c("wmis","wnon","wtru","wall")

#Colours for mutation consequence classes, shared by the panels below
mut_type_cols<-c("#00468B","#925E9F","#925E9F","#810F7C","black","#E10C04","black")
names(mut_type_cols)<-c("Missense","Nonsense","Truncating","Stop_loss","Non-Coding","Synonymous","Overall")

all_mtDNA_genes <- c("MT-CYB", "MT-ND5", "MT-ND2", "MT-ND4", "MT-ND1", "MT-CO3", "MT-ATP6","MT-ND3", "MT-ATP8", "MT-ND4L", "MT-CO2", "MT-CO1","MT-ND5","MT-ND6")
target_genes <- c("MT-CYB", "MT-ND5", "MT-ND2", "MT-ND4", "MT-ND1", "MT-CO3", "MT-ATP6","MT-ND3", "MT-ATP8", "MT-ND4L", "MT-CO2", "MT-CO1")

df_tidy$Tissue<-ifelse(grepl("pcw",df_tidy$exp_ID),"Foetal blood",ifelse(grepl("CB",df_tidy$exp_ID),"Cord blood","Adult blood"))
valid.tissues <- unique(df_tidy$Tissue)

df_tidy_annotated<-lapply(valid.tissues,function(tissue.id) {
  cat(tissue.id,sep="\n")
  
  # read in the mtdna variant file and remove patient id
  mtdna.variant.data <- df_tidy%>%
    dplyr::filter(Tissue==tissue.id)%>%
    separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select("sampleID"=Sample,chr,pos,ref,mut,vaf)%>% # next, rename/reorder the columns to make them compatible with dnds input format
    dplyr::filter(!duplicated(.))%>% # remove duplicated variants
    arrange(sampleID,pos)
  
  # run dnds with the mtDNA variants to annotate them (the selection analysis isn't actually used here)
  mtdna.dndsout <- dndscv(mtdna.variant.data, gene_list=all_mtDNA_genes, 
                          refdb = mtref_rda_path, numcode = 2, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
  # get the results with the annotated variants
  annotated.mtdna.variants <- left_join(mtdna.variant.data,mtdna.dndsout$annotmuts,by=c("sampleID","chr","pos","ref","mut"))%>%
    tidyr::replace_na(replace=list(impact="Non-Coding"))%>%
    mutate(tissue=tissue.id,.before=1)
  return(annotated.mtdna.variants)
})%>%dplyr::bind_rows()

complete.annotated.mutation.table<-df_tidy_annotated%>%
  left_join(df_tidy%>%dplyr::select(Sample,exp_ID),by=c("sampleID"="Sample"),relationship="many-to-many")%>%
  dplyr::rename("patientID"=exp_ID)%>%
  filter(!duplicated(.))%>%
  arrange(patientID,pos)%>% # order the dataframe
  tidyr::unite(col="mut_ref",chr,pos,ref,mut,sep="_",remove=F)

my_comparisons <- list(c("Missense", "Synonymous"), c("Synonymous", "Truncating"))

#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 5A ---------------
## Variant counts per gene by predicted impact, blood dataset
#-----------------------------------------------------------------------------------#

mutation_type_by_gene<-complete.annotated.mutation.table%>%
  filter(impact!="Non-Coding")%>%
  distinct(sampleID,mut_ref,gene,impact)%>%
  group_by(gene,impact)%>%
  dplyr::summarise(n=n())%>%
  ggplot(aes(y=gene,x=n,fill=impact))+
  geom_bar(position="stack",stat="identity",col="black",linewidth=0.1)+
  scale_fill_manual(values=mut_type_cols)+
  scale_x_continuous(breaks=seq(0,1400,200))+
  theme_classic()+
  my_theme+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(face="italic"))+
  labs(x="Number of variants",fill="Mutation\nconsequence")

ggsave(filename = paste0(plots_dir,"Extended_Data_Figure_05/ExtDataFig5a.mutation_type_by_gene_blood.pdf"),mutation_type_by_gene,width=3.5,height=2.5)


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 5B ---------------
## Genome-wide dN/dS by VAF level, missense and truncating, blood dataset
#-----------------------------------------------------------------------------------#

#Define the VAF 'bins' that will be used for comparing VAF distributions
#Can set this on a log scale, or as decile bins
boundaries<-c(0,0.01,0.1,0.2,0.5,1)
new_VAF_groups=paste0(100*head(boundaries,-1),"-",100*tail(boundaries,-1),"%")
VAF_groups=data.frame(labels=new_VAF_groups,LL=head(boundaries,-1),UL=tail(boundaries,-1))

#Assign each observed VAF into its 'VAF group' bin
df_tidy$VAF_group=sapply(df_tidy$vaf,function(vaf) {
  if(vaf==0) {
    return("Absent")
  } else {
    return(VAF_groups$labels[VAF_groups$LL<vaf & VAF_groups$UL>=vaf])
  }
})


dndscv_by_vaf_blood<-lapply(new_VAF_groups,function(VAF_bin) {
  cat(VAF_bin,sep="\n")
  tissue_vaf_info<-df_tidy%>%
    dplyr::filter(VAF_group==VAF_bin)%>%
    separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>% 
    dplyr::filter(!duplicated(.))%>% #only count each mutation once per individual in each VAF group
    arrange(sampleID,pos)
  
  mtdna.dndsout <- dndscv(tissue_vaf_info, gene_list=target_genes, 
                          refdb = mtref_rda_path,
                          max_coding_muts_per_sample = Inf,
                          numcode = 2,
                          max_muts_per_gene_per_sample = Inf)
  return(mtdna.dndsout)
})

globaldnds_res_by_vaf<-Map(dndsout=dndscv_by_vaf_blood,VAF_bin=new_VAF_groups,function(dndsout,VAF_bin) {
  dndsout$globaldnds%>%
    filter(!is.na(name) & complete.cases(.))%>%
    mutate(name=rename_vec[name])%>%
    mutate(VAF_group=VAF_bin)
})%>%dplyr::bind_rows()

max_dnds_value<-1.5
stats_to_include=c("Overall")
global_dnds_by_vaf_plot<-globaldnds_res_by_vaf%>%
  filter(VAF_group%in%VAF_groups$labels[-1] & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         VAF_group=factor(VAF_group,levels=VAF_groups$labels[-1]))%>%
  ggplot(aes(x=VAF_group, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Mitochondrial genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mut_type_cols, name="Mutation type") +
  scale_shape_manual(values = 17, name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', linewidth = 0.5) +
  theme_classic()+
  my_theme+
  #facet_wrap(~tissue,nrow=1)+
  theme(legend.position="none",axis.text.x=element_text(angle=90))



max_dnds_value<-2
stats_to_include=c("Truncating","Missense")
global_dnds_bytype_by_vaf_plot<-globaldnds_res_by_vaf%>%
  filter(VAF_group%in%VAF_groups$labels[-1] & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         VAF_group=factor(VAF_group,levels=VAF_groups$labels[-1]))%>%
  ggplot(aes(x=VAF_group, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.5), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.5), size = 1.5) +
  labs(y = "Mitochondrial genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mut_type_cols, name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', linewidth = 0.5) +
  theme_classic()+
  my_theme+
  facet_grid(cols=vars(name))+
  theme(legend.position="none",axis.text.x=element_text(angle=90))



max_dnds_value<-1.5
stats_to_include=c("Missense")
global_dnds_missense_by_vaf_plot<-globaldnds_res_by_vaf%>%
  filter(VAF_group%in%VAF_groups$labels[-1] & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         VAF_group=factor(VAF_group,levels=VAF_groups$labels[-1]))%>%
  ggplot(aes(x=VAF_group, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Mitochondrial genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mut_type_cols, name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', linewidth = 0.5) +
  theme_classic()+
  my_theme+
  #facet_wrap(~tissue,nrow=1)+
  theme(legend.position="none",axis.text.x=element_text(angle=90))

ggsave(filename = paste0(plots_dir,"Extended_Data_Figure_05/ExtDataFig5b.global_dnds_bytype_by_vaf.pdf"),global_dnds_bytype_by_vaf_plot,width=2.5,height=2.5)


#-----------------------------------------------------------------------------------#
### CROSS TISSUE ANALYSIS----
#-----------------------------------------------------------------------------------#

nonblood_ref_file=paste0(root_dir,"/data/metadata/non_blood_metadata.xlsx")
colony_info_file=paste0(root_dir,"/data/metadata/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt")

#Read in the mitochondrial copy number data
mito_cn=read.csv(paste0(root_dir,"/data/whole_genome_coverage_pileup_and_bedtools_annotated.csv"),header=T)%>%
  dplyr::rename("pileup_median_mtDNA_coverage"=pilup_median_mtDNA_coverage) #This column has a typo

#Import individual level metadata for the adult and foetal blood samples
ref_df<-readxl::read_excel(nonblood_ref_file) #Import metadata relating to all individuals studied
colony_info<-read.delim(colony_info_file) #Import colony-level data for the lymphoid dataset as there are multiple different cell types

#Define the dataframe for converting IDs between my IDs, and Andrew's for the different datasets
translate=data.frame(dataset=c("KY","HL","SO","PR","LM","NW","lymph","EM"),
                     al_ref=c("lung_organoid","colon","colon_ibd","muty_mutant","endometrium","blood_MPN","immune","blood_emily"))

muty_samples=c("PD44887","PD44888","PD44889","PD44890","PD44891")

#Define a set of 'black-listed' mutations - these are sets of artefacts that recurrently slip through the filters
#This is often due to haplotype-specific artefacts that map as SNVs
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_456_C_T","MT_567_A_C","MT_574_A_C","MT_8270_C_T","MT_16170_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C")

#-----------------------------------------------------------------------------------#
# IMPORT ALL CROSS TISSUE DATASETS ----------------------------------------
#-----------------------------------------------------------------------------------#

all_cohorts_plus_CML<-c("KY","HL","SO","PR","LM","NW","CML","lymph","blood")
all_cohorts<-c("KY","HL","SO","PR","LM","NW","lymph","blood")

name_conversion_vec=c("Normal blood","Lymphoid","MPN","Colon","IBD-affected\nColon","MUTYH-mutant\nColon","Endometrium","Bronchial\nepithelium")
names(name_conversion_vec)=c("blood","lymph","NW","HL","SO","PR","LM","KY")
tissue_cols<-c("#96e97c","#fd81c8","#145a6a","#65e6f9","#781486","#5f70cc","#aea2eb","#1d6d1f")
names(tissue_cols)<-name_conversion_vec

all_mito_datasets<-lapply(all_cohorts_plus_CML,function(dataset) {
  dataset_mito_data<-readRDS(paste0(root_dir,"/data/nonblood/mito_mutation_data_",dataset,".RDS"))
  CN_correlating_muts<-readRDS(paste0(root_dir,"/data/nonblood/CN_correlation_",dataset,".RDS"))
  
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
## Annotate the cross-tissue variants with dndscv. NB this redefines
## complete.annotated.mutation.table for the cross-tissue panels below; the blood
## panels above are already written.
#-----------------------------------------------------------------------------------#

library(dndscv)
mtref_rda_path=ifelse(Sys.info()['sysname']=="Darwin",paste0(root_dir,"/data/mtref.rda"),"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/mtref.rda")
input.dir <- paste0(root_dir,"/dnds_tables/")

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

#Define a gene order by which respiratory chain complex they are part of
gene_order=c(paste0("MT-ND",1:6),"MT-ND4L","MT-CYB",paste0("MT-CO",1:3),"MT-ATP6","MT-ATP8")
gene_order<-gene_order[gene_order%in%target_genes]

# Vector to rename the dnds types in more understandable format
rename_vec=c("Missense","Nonsense","Truncating","Overall")
names(rename_vec)=c("wmis","wnon","wtru","wall")

df_tidy_annotated<-Map(df_tidy=all_df_tidy,this_tissue=names(all_df_tidy),function(df_tidy,this_tissue) {
  
  #cat(this_tissue,sep="\n")
  
  # read in the mtdna variant file and remove patient id
  mtdna.variant.data <- df_tidy%>%
    separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select("sampleID"=Sample,chr,pos,ref,mut,vaf)%>% # next, rename/reorder the columns to make them compatible with dnds input format
    dplyr::filter(!duplicated(.))%>% # remove duplicated variants
    arrange(sampleID,pos)
  
  # run dnds with the mtDNA variants to annotate them (the selection analysis isn't actually used here)
  mtdna.dndsout <- suppressMessages(dndscv(mtdna.variant.data,
                          gene_list=all_mtDNA_genes, 
                          refdb = mtref_rda_path,
                          numcode = 2,
                          max_coding_muts_per_sample = Inf,
                          max_muts_per_gene_per_sample = Inf))
  
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
  tidyr::unite(col="mut_ref",chr,pos,ref,mut,sep="_",remove=F)


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 5C ---------------
## Variant counts in protein-coding genes by impact, cross-tissue datasets
#-----------------------------------------------------------------------------------#

mutation_type_by_tissue<-complete.annotated.mutation.table%>%
  filter(impact!="Non-Coding" & !tissue%in%c("blood","CML"))%>%
  distinct(tissue,sampleID,mut_ref,gene,impact)%>%
  group_by(tissue,impact)%>%
  dplyr::summarise(n=n())%>%
  dplyr::mutate(tissue=factor(name_conversion_vec[tissue],levels=name_conversion_vec))%>%
  ggplot(aes(y=tissue,x=n,fill=impact))+
  geom_bar(position="stack",stat="identity",col="black",linewidth=0.1)+
  scale_fill_manual(values=mut_type_cols)+
  scale_x_continuous(breaks=seq(0,1400,200))+
  theme_classic()+
  my_theme+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(face="italic"))+
  labs(x="Number of coding variants",fill="Mutation\nconsequence")

ggsave(filename = paste0(plots_dir,"Extended_Data_Figure_05/ExtDataFig5c.mutation_type_by_tissue.pdf"),mutation_type_by_tissue,width=3.5,height=2)


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 5D ---------------
## Genome-wide dN/dS by tissue
#-----------------------------------------------------------------------------------#

dndscv_by_tissue<-Map(df_tidy=all_df_tidy[all_cohorts],this_tissue=all_cohorts,function(df_tidy,this_tissue) {
  
  #cat(this_tissue,sep="\n")
  tissue_info<-df_tidy%>%
    separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>%
    dplyr::filter(!duplicated(.))%>%
    arrange(sampleID,pos)
  
  mtdna.dndsout <- dndscv(tissue_info,
                          gene_list=target_genes,
                          refdb = mtref_rda_path,
                          numcode=2,
                          max_coding_muts_per_sample = Inf,
                          max_muts_per_gene_per_sample = Inf)
  return(mtdna.dndsout)
})

globaldnds_res<-Map(dndsout=dndscv_by_tissue,this_tissue=all_cohorts,function(dndsout,this_tissue) {
  dndsout$globaldnds%>%
    filter(!is.na(name) & complete.cases(.))%>%
    mutate(name=rename_vec[name])%>%
    mutate(tissue=this_tissue)
})%>%dplyr::bind_rows()

nb.cols <- length(name_conversion_vec)
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)

max_dnds_value<-2
stats_to_include=c("Overall")
global_dnds_nonblood<-globaldnds_res%>%
  filter(name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         tissue=factor(name_conversion_vec[tissue],levels=name_conversion_vec))%>%
  ggplot(aes(x=tissue, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Mitochondrial genome-wide dN/dS") +
  theme_classic() + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6))+
  scale_y_continuous(limits=c(0,max_dnds_value))+
  scale_color_manual(values = mut_type_cols, name="Mutation type") +
  scale_shape_manual(values = 17, name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', linewidth = 0.5) +
  my_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle=90),
        legend.position="none")

ggsave(filename = paste0(plots_dir,"Extended_Data_Figure_05/ExtDataFig5d.global_dnds_by_tissue.pdf"),global_dnds_nonblood,width=2,height=2)


#-----------------------------------------------------------------------------------#
## Assign the cross-tissue variants to VAF bins, used by panel e.
#-----------------------------------------------------------------------------------#

#Define the VAF 'bins' that will be used for comparing VAF distributions
#Can set this on a log scale, or as decile bins, or dichotomous (>/< 20%)
#Only mutations >3% are actually included
boundaries<-c(0,0.03,0.1,0.2,0.5,1)
new_VAF_groups=paste0(round(100*head(boundaries,-1)),"-",round(100*tail(boundaries,-1)),"%")
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
##--------Generate EXTENDED DATA FIG. 5E ---------------
## Genome-wide dN/dS by tissue and VAF
#-----------------------------------------------------------------------------------#

dndscv_by_tissue_and_vaf<-Map(df_tidy=all_df_tidy[all_cohorts],this_tissue=all_cohorts,function(df_tidy,this_tissue) {
  #cat(this_tissue,sep = "\n")
  dndscv_by_vaf_tissue<-lapply(new_VAF_groups[2:length(new_VAF_groups)],function(VAF_bin) {
    #cat(VAF_bin,sep="\n")
    tissue_vaf_info<-df_tidy%>%
      dplyr::filter(VAF_group==VAF_bin)%>%
      separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
      mutate(pos=as.numeric(pos))%>%
      dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>% 
      dplyr::filter(!duplicated(.))%>% #only count each mutation once per individual in each VAF group
      arrange(sampleID,pos)
    
    mtdna.dndsout <- dndscv(tissue_vaf_info,
                            gene_list=target_genes, 
                            refdb = mtref_rda_path,
                            numcode = 2,
                            max_coding_muts_per_sample = Inf,
                            max_muts_per_gene_per_sample = Inf)
    return(mtdna.dndsout)
  })
  return(dndscv_by_vaf_tissue)
})

globaldnds_res_by_tissue_and_vaf<-Map(list1=dndscv_by_tissue_and_vaf,this_tissue=all_cohorts,function(list1,this_tissue) {
  globaldnds_res_by_vaf<-Map(dndsout=list1,VAF_bin=new_VAF_groups[2:length(new_VAF_groups)],function(dndsout,VAF_bin) {
    dndsout$globaldnds%>%
      filter(!is.na(name) & complete.cases(.))%>%
      mutate(name=rename_vec[name])%>%
      mutate(VAF_group=VAF_bin)
  })%>%dplyr::bind_rows()%>%
    mutate(tissue=this_tissue)
})%>%dplyr::bind_rows()

nb.cols <- length(name_conversion_vec)
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(3)
#mycolors<-ggsci::pal_lancet(palette = "lanonc")(9)

max_dnds_value<-5
stats_to_include=c("Overall")
global_dnds_by_tissue_and_vaf_plot<-globaldnds_res_by_tissue_and_vaf%>%
  filter(VAF_group%in%VAF_groups$labels[-1] & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         VAF_group=factor(VAF_group,levels=VAF_groups$labels[-1]),
         tissue=factor(name_conversion_vec[tissue],levels=name_conversion_vec))%>%
  ggplot(aes(x=VAF_group, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_hline(yintercept=1, linetype='dashed', col = 'black', linewidth = 0.5) +
  geom_linerange(position= position_dodge2(width=0.75), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mut_type_cols, name="Mutation type") +
  scale_shape_manual(values = 17, name = "Mutation type") + 
  theme_classic()+
  my_theme+
  facet_grid(cols=vars(tissue))+
  theme(legend.position="none",
        axis.text.x=element_text(angle=90),
        axis.title.x=element_blank(),
        strip.text.y = element_text(angle = 0))

ggsave(filename = paste0(plots_dir,"Extended_Data_Figure_05/ExtDataFig5e.global_dnds_by_tissue_and_vaf.pdf"),global_dnds_by_tissue_and_vaf_plot,width=7,height=2)


#-----------------------------------------------------------------------------------#
## Combine the epithelial tissues (colon, endometrium, bronchial epithelium) for
## panels f and g.
#-----------------------------------------------------------------------------------#

epithelial_datasets=c("HL","KY","LM","PR","SO")
combined_epithelial_for_dnds<-complete.annotated.mutation.table%>%
  dplyr::filter(tissue%in%epithelial_datasets)%>%
  separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
  mutate(pos=as.numeric(pos))%>%
  dplyr::distinct(sampleID,chr,pos,ref,mut,.keep_all = F)%>%
  arrange(sampleID,pos)


mtdna.dndsout.combined.epithelial <- dndscv(combined_epithelial_for_dnds,
                        gene_list=target_genes,
                        refdb = mtref_rda_path,
                        numcode=2,
                        max_coding_muts_per_sample = Inf,
                        max_muts_per_gene_per_sample = Inf)

mtdna.dndsout.combined.epithelial$globaldnds%>%
  filter(!is.na(name) & complete.cases(.))%>%
  mutate(name=rename_vec[name])%>%
  filter(name%in%c("Overall","Nonsense","Missense"))%>%
  ggplot(aes(x=name,y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Tissue") +
  theme_classic() + 
  scale_color_manual(values = mut_type_cols, name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', linewidth = 0.5) +
  my_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle=90),
        legend.position="none")


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 5F ---------------
## Empirical cumulative distribution of VAF by consequence, combined epithelial tissues
#-----------------------------------------------------------------------------------#

#Gives some sense of differing heteroplasmy distributions
#But bear in mind there are multiple non-independent measures for some variants
#Would need to account for this e.g. by taking single measure of e.g. max VAF of variant per individual
syn_vals <- complete.annotated.mutation.table %>%
  dplyr::filter(tissue%in%epithelial_datasets)%>%
  dplyr::filter(impact == "Synonymous"&vaf>0.01) %>% pull(vaf)

mis_vals <- complete.annotated.mutation.table %>%
  dplyr::filter(tissue%in%epithelial_datasets)%>%
  dplyr::filter(impact == "Missense"&vaf>0.01) %>% pull(vaf)

nons_vals <- complete.annotated.mutation.table %>%
  dplyr::filter(tissue%in%epithelial_datasets)%>%
  dplyr::filter(impact == "Nonsense"&vaf>0.01) %>% pull(vaf)

ks_syn_mis <- ks.test(syn_vals, mis_vals)
ks_syn_non <- ks.test(syn_vals, nons_vals)

find_ks_coords <- function(vals1, vals2) {
  ecdf1 <- ecdf(vals1)
  ecdf2 <- ecdf(vals2)
  all_y <- sort(c(vals1, vals2))
  diffs <- abs(ecdf1(all_y) - ecdf2(all_y))
  max_y <- all_y[which.max(diffs)]
  data.frame(
    y    = max_y,
    xmin = min(ecdf1(max_y), ecdf2(max_y)),
    xmax = max(ecdf1(max_y), ecdf2(max_y))
  )
}


fmt_p <- function(p) {
  ifelse(p < 0.001, "p < 0.001",
         ifelse(p < 0.01,  paste0("p = ", round(p, 3)),
                paste0("p = ", round(p, 2))))
}

coords_non <- find_ks_coords(syn_vals, nons_vals)
coords_mis <- find_ks_coords(syn_vals, mis_vals)


mutation_type_ecdf_plot<-complete.annotated.mutation.table %>%
  filter(impact != "Non-Coding" & impact != "Stop_loss") %>%
  ggplot(aes(y = vaf, col = impact)) +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  # Syn vs Missense segment
  annotate("segment",
           x = coords_mis$xmin, xend = coords_mis$xmax,
           y = coords_mis$y,    yend = coords_mis$y,
           colour = "black", linewidth = 0.5, linetype = "dashed") +
  annotate("text",
           x = mean(c(coords_mis$xmin, coords_mis$xmax)),
           y = coords_mis$y,
           vjust = -0.5, size = 2,
           label = paste0("Syn vs Mis\n", fmt_p(ks_syn_mis$p.value))) +
  # Syn vs Nonsense segment
  annotate("segment",
           x = coords_non$xmin, xend = coords_non$xmax,
           y = coords_non$y,    yend = coords_non$y,
           colour = "black", linewidth = 0.5, linetype = "dashed") +
  annotate("text",
           x = mean(c(coords_non$xmin, coords_non$xmax)),
           y = coords_non$y,
           vjust = -0.5, size = 2,
           label = paste0("Syn vs Nons\n", fmt_p(ks_syn_non$p.value)))+
  scale_y_continuous(limits = c(0.01, 1),
                     labels = scales::percent_format(accuracy = 1),
                     name = "Variant Allele Fraction (VAF)") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     name = "Cumulative fraction of variants") +
  theme_classic(base_size = 13) +
  scale_color_manual(values = mut_type_cols) +
  my_theme +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.3, 0.8)) +
  labs(col = "Mutation type")

ggsave(filename = paste0(plots_dir,"Extended_Data_Figure_05/ExtDataFig5f.mutation_type_ecdf_epithelial.pdf"),mutation_type_ecdf_plot,width=2.5,height=2)


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 5G ---------------
## Genome-wide dN/dS above and below 20% heteroplasmy, combined epithelial tissues
#-----------------------------------------------------------------------------------#

cutoffs=list(low_vaf=c(0.03,0.2),
             high_vaf=c(0.2,1.1))
cutoff_names=c("<20%",">20%")

# cutoffs=list(vaf1=c(0.03,0.1),
#              vaf2=c(0.1,0.2),
#              vaf3=c(0.2,0.5),
#              vaf4=c(0.5,1.011))
# cutoff_names=c("<10%","10-20%","20-50%",">50%")


dndsout.epithelial<-lapply(cutoffs,function(cutoffs) {
  combined_nonblood_for_dnds<-complete.annotated.mutation.table%>%
    dplyr::filter(tissue%in%epithelial_datasets & vaf>=cutoffs[1] & vaf<cutoffs[2])%>%
    separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select(sampleID,chr,pos,ref,mut)%>%
    dplyr::filter(!duplicated(.))%>%
    arrange(sampleID,pos)
  
  mtdna.dndsout.combined.epithelial <- dndscv(combined_nonblood_for_dnds,
                                            gene_list=target_genes,
                                            refdb = mtref_rda_path,
                                            numcode=2,
                                            max_coding_muts_per_sample = Inf,
                                            max_muts_per_gene_per_sample = Inf)
  return(mtdna.dndsout.combined.epithelial)
})

globaldnds_by_vaf_epithelial<-Map(dndsout=dndsout.epithelial,vaf_group=cutoff_names,function(dndsout,vaf_group) {
  dndsout$globaldnds%>%
    filter(!is.na(name) & complete.cases(.))%>%
    mutate(name=rename_vec[name])%>%
    mutate(VAF=vaf_group)
})%>%dplyr::bind_rows()

dnds_global_by_dichotomous_vaf_epithelial_plot<-globaldnds_by_vaf_epithelial%>%
  filter(!is.na(name) & complete.cases(.))%>%
  filter(name%in%c("Overall","Nonsense","Missense"))%>%
  mutate(VAF=factor(VAF,levels = cutoff_names))%>%
  ggplot(aes(x=VAF,y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), linewidth = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Mitochondrial genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  facet_grid(~name)+
  scale_color_manual(values = mut_type_cols, name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', linewidth = 0.5) +
  my_theme+
  theme(axis.text.x = element_text(angle=0),
        legend.position="none")

ggsave(filename = paste0(plots_dir,"Extended_Data_Figure_05/ExtDataFig5g.dnds_global_by_dichotomous_vaf_epithelial.pdf"),dnds_global_by_dichotomous_vaf_epithelial_plot,width=3,height=2)

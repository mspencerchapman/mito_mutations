library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(phylosignal)
options(stringsAsFactors = F)

my_working_directory<-"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood"
R_function_files = list.files("/lustre/scratch126/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)
genomeFile="/nfs/cancer_ref02/human/GRCh37d5/genome.fa"
root_dir="/lustre/scratch125/casm/team268im/al28/mtDNA"

#Import the mitochondrial copy number data
mito_cn=read.csv("/lustre/scratch125/casm/team268im/al28/mtDNA/whole_genome_coverage_pileup_and_bedtools_annotated.csv",header=T)

#Import sample metadata
ref_df<-readxl::read_excel("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/non_blood_metadata.xlsx")

translate=data.frame(dataset=c("KY","HL","SO","PR","LM","NW","lymph"),
                     al_ref=c("lung_organoid","colon","colon_ibd","muty_mutant","endometrium","blood_MPN","immune"))

muty_samples=c("PD44887","PD44888","PD44889","PD44890","PD44891")


##CHOOSE DATASET AND START ANALYSIS
dataset="LM" #Can be any of c("KY","HL","SO","PR","LM","NW","lymph")
CN_correlating_muts=readRDS(paste0("CN_correlation_",dataset,".RDS"))
dataset_mito_data=readRDS(paste0("mito_mutation_data_",dataset,".RDS"))
figures_dir=paste0(dataset,"/")

rho_cut_off=0
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C") #These are the ones to exclude for the endometrial analysis

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

#In the lymph dataset, the 'trees' are only dummy phylogenies and don't reflect the true clonal relationships
# Therefore, only the shared mutations in the youngest individuals (TX001/ TX002) likely reflect true heteroplasmic oocyte mutations (as opposed to shared somatically-acuiqred variants)
if(dataset=="lymph") {
  dataset_mito_data<-Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) {
    if(!exp_ID=="TX001") {list$het_oocyte_muts<-c()}
    return(list)
  })
}

####
#Generate tidy data frame of samples, mutations and vafs
mutCN_cutoff=25 #If the mitochondrial copy number is over 25, retain mutation even if is in the "CN correlating muts" list
vaf_cut_off<-0.03
rho_cut_off<-0
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
  left_join(mito_cn)%>%
  mutate(implied_mut_CN=vaf*bedtools_mtDNA_genomes)%>%
  dplyr::filter(!(mut_ref%in%CN_correlating_muts & implied_mut_CN<=mutCN_cutoff))

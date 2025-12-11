library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpmisc)
library(gridExtra)
library(phylosignal)
library(treemut)
options(stringsAsFactors = F)

#Set these file paths before running the script
genomeFile="~/Documents/Reference_files/hs37d5.fa"
root_dir="~/R_work/mito_mutations_blood"
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"/data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
nonblood_ref_file=paste0(root_dir,"/data/metadata/non_blood_metadata.xlsx")
colony_info_file=paste0(root_dir,"/data/metadata/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt")
figures_dir=paste0(root_dir,"/figures/")
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
##----------------------CHOOSE DATASET AND START ANALYSIS---------------------
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
##----------------------Enhance the tidy dataframe with 'VAF bin' info---------------
#-----------------------------------------------------------------------------------#
#Define the VAF 'bins' that will be used for comparing VAF distributions
#Can set this on a log scale, or as decile bins
log_scale=F
if(log_scale) {
  new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
  VAF_groups=data.frame(labels=new_VAF_groups,LL=c(0,2^(-10:-1)),UL=2^(-10:0))
} else {
  new_VAF_groups=paste0(100*seq(0,0.95,0.05),"-",100*seq(0.05,1,0.05),"%")
  VAF_groups=data.frame(labels=new_VAF_groups,LL=seq(0,0.95,0.05),UL=seq(0.05,1,0.05))
}

all_df_tidy<-lapply(all_df_tidy,function(df_tidy) {
  
  #Assign each observed VAF into its 'VAF group' bin
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
##----------------------Generate VAF distribution matrices---------------
#-----------------------------------------------------------------------------------#
# For each data set, generate a matrix where each row is an individual, and each column is a VAF bin
# The matrix is filled with the average number of mutations in that bin per sample
# If the number of samples is too low, this values will have a large degree of sampling error, therefore limit to individuals with >15 samples

minimum_samples_to_include=15

#Iterate across samples/ VAF groups to get numbers of mutations per sample/ VAF group
VAF_distribution_mats<-Map(dataset_mito_data=all_mito_datasets[all_cohorts],df_tidy=all_df_tidy[all_cohorts],function(dataset_mito_data,df_tidy) {
  all_samples=dplyr::bind_rows(Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) data.frame(exp_ID=exp_ID,Sample=list$tree$tip.label)))
  samples_vaf_groups_df<-right_join(all_samples,tidyr::expand_grid(Sample=all_samples$Sample,VAF_group=new_VAF_groups),by="Sample")
  nmut_by_vaf_by_sample<-left_join(samples_vaf_groups_df,df_tidy%>%
                                     group_by(Sample,VAF_group)%>%
                                     dplyr::summarise(n=n()),by=c("Sample","VAF_group"))%>%
    replace_na(replace = list(n=0))
  
    donors_to_include=names(dataset_mito_data)[sapply(dataset_mito_data,function(list) length(list$tree$tip.label)>=minimum_samples_to_include)]
  
  tissue_mat<-lapply(donors_to_include,function(this_exp_ID) {
    
    #Highly clonally related samples are not 'independent' measures of mtDNA mutation burden. Therefore, randomly select one of
    #any highly clonally related clades (those with >half of mutations related)
    mut_burden<-get_mut_burden(dataset_mito_data[[this_exp_ID]]$tree.ultra)[1]
    expanded_clade_nodes<-get_expanded_clade_nodes(tree=dataset_mito_data[[this_exp_ID]]$tree.ultra,min_samples=2,height_cut_off=mut_burden/2)
    samples_to_drop=unlist(lapply(expanded_clade_nodes$nodes,function(node) {
      clade_samples<-getTips(tree=dataset_mito_data[[this_exp_ID]]$tree.ultra,node)
      to_drop=sample(clade_samples,size=length(clade_samples)-1)
    }))
    
    donor_df<-nmut_by_vaf_by_sample%>%filter(exp_ID==this_exp_ID & !Sample%in%samples_to_drop)
    
    if(nrow(donor_df)>0) {
      donor_df%>%
        pivot_wider(names_from = Sample,values_from = n)%>%
        dplyr::select(-exp_ID)%>%
        tibble::column_to_rownames(var="VAF_group")%>%
        as.matrix()%>%rowMeans()%>%
        as.data.frame()%>%
        rlang::set_names(this_exp_ID)
    } else {
      NULL
    }
    
  })%>%bind_cols()%>%
    as.matrix()%>%t()
  
  return(tissue_mat)
})

#-----------------------------------------------------------------------------------#
##-----------------Generate plots of these VAF distributions by individual---------
#-----------------------------------------------------------------------------------#

# Set up vector to rename the datasets based on their actual tissue type
name_conversion_vec=c("MPN","Colon","IBD Colon","Lung","Lymphoid","MUTYH-mutant\nColon","Endometrium","Normal blood")
names(name_conversion_vec)=c("NW","HL","SO","KY","lymph","PR","LM","blood")
tissue_cols<-c("#96e97c","#145a6a","#65e6f9","#781486","#fd81c8","#5f70cc","#aea2eb","#1d6d1f")
names(tissue_cols)<-name_conversion_vec

# Now generate the plots one dataset at a time
temp=Map(mat=VAF_distribution_mats,dataset=all_cohorts,function(mat,dataset) {
  cat(dataset,sep="\n")
  figures_dir=paste0(dataset,"/")
  exp_ID_by_age<-ref_df%>%filter(ID%in%rownames(mat) & Cohort==dataset)%>%
    arrange(Age)%>%pull(ID)
  if(length(exp_ID_by_age)==0) {stop(return(NULL))}
  
  VAF_dist_plot<-as.data.frame(mat)%>%t()%>%
    as.data.frame()%>%tibble::rownames_to_column(var="VAF_group")%>%
    tidyr::gather(-VAF_group,key="exp_ID",value="n_per_sample")%>%
    ggplot(aes(x=factor(VAF_group,levels=new_VAF_groups),y=n_per_sample))+
    geom_bar(stat="identity",fill=tissue_cols[name_conversion_vec[dataset]],col="black",linewidth=0.1)+
    facet_grid(rows=vars(factor(exp_ID,levels=exp_ID_by_age)))+
    theme_classic()+
    my_theme+
    theme(axis.text.x=element_text(angle=90),strip.text.y=element_text(angle=0))+
    labs(x="Heteroplasmy level",y="mtDNA mutations per cell")
  
  ggsave(filename=paste0(plots_dir,"VAF_dist_plot2_",dataset,".pdf"),VAF_dist_plot,width=2.5,height=(0.6*nrow(mat)))
  
})

# Now a single plot across all datasets
comb_VAF_distribution=Map(mat=VAF_distribution_mats,data_set=all_cohorts,function(mat,data_set) {
  cat(data_set,sep="\n")
  exp_ID_by_age<-ref_df%>%filter(ID%in%rownames(mat) & Cohort==data_set)%>%
    arrange(Age)%>%pull(ID)
  if(length(exp_ID_by_age)==0) {stop(return(NULL))}
  
  as.data.frame(mat)%>%t()%>%
    as.data.frame()%>%tibble::rownames_to_column(var="VAF_group")%>%
    tidyr::gather(-VAF_group,key="exp_ID",value="n_per_sample")%>%
    dplyr::mutate(dataset=data_set)
})%>%dplyr::bind_rows()

cohort_exp_ID_by_age<-ref_df%>%
  mutate(cohort_exp_ID=paste(Cohort,ID,sep="_"))%>%
  filter(ID%in%comb_VAF_distribution$exp_ID)%>%arrange(Cohort,Age)%>%pull(cohort_exp_ID)
combined_VAF_dist_plot<-comb_VAF_distribution%>%
  mutate(cohort_exp_ID=factor(paste(dataset,exp_ID,sep="_"),levels=cohort_exp_ID_by_age))%>%
  mutate(dataset=name_conversion_vec[dataset])%>%
  ggplot(aes(fill=dataset,x=factor(VAF_group,levels=new_VAF_groups),y=n_per_sample))+
  geom_bar(stat="identity",col="black",linewidth=0.1)+
  facet_wrap(~cohort_exp_ID,ncol=4,
             labeller=labeller(cohort_exp_ID=function(x) {stringr::str_split(x,pattern="_",simplify=T)[,2]}))+
  scale_fill_manual(values=tissue_cols)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),
        strip.text=element_text(angle=0,size = 7, margin = margin(1,0.5,1,0.5)))+
  labs(x="Heteroplasmy level",y="mtDNA mutations per cell")

ggsave(filename=paste0(plots_dir,"combined_VAF_distirbution_plot.pdf"),combined_VAF_dist_plot,width=7,height=12)

#Now do a comparison of normal blood vs lymphoid in the individuals with both
both_dataset_IDs<-intersect(rownames(VAF_distribution_mats$lymph),rownames(VAF_distribution_mats$blood))
both_dataset_IDs<-ref_df%>%dplyr::filter(ID%in%both_dataset_IDs&Cohort=="blood")%>%arrange(Age)%>%pull(ID)

blood_lymph_common_data<-lapply(c("lymph","blood"),function(tissue) {
  as.data.frame(VAF_distribution_mats[[tissue]])%>%
    tibble::rownames_to_column(var="exp_ID")%>%
    tidyr::gather(-exp_ID,key="VAF_group",value=n_per_sample)%>%
    filter(exp_ID%in%both_dataset_IDs)%>%
    mutate(dataset=tissue,.before=1)
})%>%bind_rows()

blood_vs_lymph_VAF_distribution<-blood_lymph_common_data%>%
  mutate(dataset=factor(name_conversion_vec[dataset],levels=c("Normal blood","Lymphoid")),exp_ID=factor(exp_ID,levels=both_dataset_IDs))%>%
  ggplot(aes(x=factor(VAF_group,levels=new_VAF_groups),y=n_per_sample,fill=dataset))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=tissue_cols)+
  facet_grid(exp_ID~dataset)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),strip.text.y=element_text(angle=0),legend.position="none")+
  labs(x="Heteroplasmy level",y="mtDNA mutations per cell")
ggsave(filename=paste0(plots_dir,"blood_vs_lymph_VAF_dist.pdf"),blood_vs_lymph_VAF_distribution,width=3.3,height=2.5)

#Kolmogorov Smirnov test to compare the distribution of mutations in the lymphocytes vs HSCs for each of these individuals
for(i in 1:length(both_dataset_IDs)) {
  this_ID<-both_dataset_IDs[i]
  
  blood_vafs<-all_df_tidy[["blood"]]%>%filter(exp_ID==this_ID)%>%pull(vaf)
  lymph_vafs<-all_df_tidy[["lymph"]]%>%filter(exp_ID==this_ID)%>%pull(vaf)
  ks.res<-ks.test(x=blood_vafs,y=lymph_vafs)
  print(ks.res)
}

#-----------------------------------------------------------------------------------#
##----------------------Generate tidy VAF distribution dataframe---------------
#-----------------------------------------------------------------------------------#
#Generate a single tidy dataframe with numbers of mutations in each VAF category per sample, with each VAF category as a separate row
#Because this doesn't average across samples, can have a more liberal inclusion of individuals with fewer samples
minimum_samples_to_include=8

VAF_distribution_by_sample_df<-Map(dataset=all_cohorts,dataset_mito_data=all_mito_datasets[all_cohorts],df_tidy=all_df_tidy[all_cohorts],function(dataset,dataset_mito_data,df_tidy) {
  all_samples=dplyr::bind_rows(Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) data.frame(exp_ID=exp_ID,Sample=list$tree$tip.label)))
  samples_vaf_groups_df<-right_join(all_samples,tidyr::expand_grid(Sample=all_samples$Sample,VAF_group=new_VAF_groups),by="Sample")
  nmut_by_vaf_by_sample<-left_join(samples_vaf_groups_df,df_tidy%>%
                                     group_by(Sample,VAF_group)%>%
                                     dplyr::summarise(n=n()),by=c("Sample","VAF_group"))%>%
    replace_na(replace = list(n=0))
  
  
  donors_to_include=names(dataset_mito_data)[sapply(dataset_mito_data,function(list) length(list$tree$tip.label)>=minimum_samples_to_include)]
  
  tissue_df<-lapply(donors_to_include,function(this_exp_ID) {
    
    #Highly clonally related samples are not 'independent' measures of mtDNA mutation burden. Therefore, randomly select one of
    #any highly clonally related clades (those with >half of mutations related)
    mut_burden<-get_mut_burden(dataset_mito_data[[this_exp_ID]]$tree.ultra)[1]
    expanded_clade_nodes<-get_expanded_clade_nodes(tree=dataset_mito_data[[this_exp_ID]]$tree.ultra,min_samples=2,height_cut_off=mut_burden/2)
    samples_to_drop=unlist(lapply(expanded_clade_nodes$nodes,function(node) {
      clade_samples<-getTips(tree=dataset_mito_data[[this_exp_ID]]$tree.ultra,node)
      to_drop=sample(clade_samples,size=length(clade_samples)-1)
    }))
    
    donor_df<-nmut_by_vaf_by_sample%>%filter(exp_ID==this_exp_ID & !Sample%in%samples_to_drop)
    return(donor_df)
  })%>%dplyr::bind_rows()%>%
    mutate(dataset=dataset,.before=1)
  return(tissue_df)
})%>%dplyr::bind_rows()

#------------------------------------------------------------------------------------------------#
##-----Analyse drift using the number of homoplasmic/ near-homoplasmic mutations per sample------
#------------------------------------------------------------------------------------------------#

# Select the VAF bins >85%, to include only the number of 'homoplasmic' or 'near homoplasmic' mutations
n_near_homoplasmic_per_sample<-VAF_distribution_by_sample_df%>%
  filter(VAF_group%in%new_VAF_groups[18:20])%>%
  group_by(dataset,exp_ID,Sample)%>%
  dplyr::summarise(dataset=dataset[1],exp_ID=exp_ID[1],Sample=Sample[1],n_homo=sum(n))%>%
  left_join(ref_df%>%dplyr::select(ID,Age),relationship = "many-to-many",by=c("exp_ID"="ID"))

Sample_order<-ref_df%>%mutate(Cohort=name_conversion_vec[Cohort])%>%arrange(Cohort,Age)%>%unite(col = "Cohort_ID",Cohort,ID,sep="_")%>%pull(Cohort_ID)

# Plot a basic summary of the number of homoplasmic mutations per sample by individual in a line plot  
homo_muts_per_sample<-n_near_homoplasmic_per_sample%>%
  mutate(dataset=name_conversion_vec[dataset])%>%
  unite(col = "Cohort_ID",dataset,exp_ID,sep="_",remove = F)%>%
  ggplot(aes(x=forcats::fct_reorder(factor(Sample),n_homo),y=n_homo,col=factor(dataset)))+
  #geom_point(size=0.3)+
  geom_line(aes(group=Cohort_ID))+
  facet_grid(~factor(Cohort_ID,levels=Sample_order),scales = "free",space = "free")+
  scale_color_manual(values = tissue_cols)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_blank(),
        strip.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(x="Individual sample",col="",y="Burden of homoplasmic mutations")

ggsave(plot = homo_muts_per_sample,filename = paste0(plots_dir,"homo_muts_per_sample.pdf"),width=15,height=3)

#------------------------------------------------------------------------------------------------#
##----------------------------------- LINEAR REGRESSION MODEL-----------------------------------
#------------------------------------------------------------------------------------------------#
#Do linear model just including the mean values from the adult blood and lymphoid data
#Because these have large numbers of samples, the mean values are fairly accurate and the linear model is a reasonable way
# to assess the increased number of homoplasmic muts in lymphoid cells compare to HSCs of the same age

studies_to_include=c("blood","lymph")
blood_vs_lymph_homo_muts_dat<-n_near_homoplasmic_per_sample%>%
  group_by(dataset,exp_ID)%>%
  dplyr::summarise(n=n(),mean=mean(n_homo),var=var(n_homo),sd=sd(n_homo))%>%
  mutate(sem=sd/sqrt(n))%>%
  left_join(ref_df,by=c("dataset"="Cohort","exp_ID"="ID"))%>%
  filter(dataset%in%studies_to_include & Age>20)

lmer.blood_vs_lymph<-lme4::lmer(mean~Age+dataset + (1|exp_ID),data = blood_vs_lymph_homo_muts_dat)
coefs<-summary(lmer.blood_vs_lymph)$coefficients
coef_df<-data.frame(dataset=c("Normal blood","Lymphoid"),slope=rep(coefs['Age','Estimate'],2),intercept=c(coefs['(Intercept)','Estimate'],coefs['(Intercept)','Estimate']+coefs['datasetlymph','Estimate']))

HSC_vs_lymph_mean_homo_burden<-blood_vs_lymph_homo_muts_dat%>%
  mutate(dataset=name_conversion_vec[dataset])%>%
  ggplot(aes(x=Age,y=mean,ymin=(mean-1.96*sem),ymax=(mean+1.96*sem),col=dataset))+
  geom_point(size=0.5)+
  geom_errorbar(width=2,alpha=0.5)+
  geom_abline(aes(slope=slope,intercept=intercept,col=dataset),data = coef_df,show.legend = F)+
  scale_color_manual(values=tissue_cols)+
  scale_x_continuous(limits=c(20,90))+
  labs(x="Age",y="Mean homoplasmic mutations per cell",col="")+
  theme_classic()+
  my_theme

ggsave(plot = HSC_vs_lymph_mean_homo_burden,filename = paste0(plots_dir,"HSC_vs_lymph_mean_homo_burden.pdf"),width=3,height=2)


#------------------------------------------------------------------------------------------------#
##----------------------------------- POISSON REGRESSION MODEL-----------------------------------
#------------------------------------------------------------------------------------------------#
# Perform a poisson regression model for this data to get a sense of different rates of drift in different tissues
# Uses a GLM mixed-effects model
# Fixed effects are 'Age' and the 'Tissue' (defined by the dataset)
# Random effects are the individual ID and the Age
# Include only the adult data as this is the more linear region of the data

#Check mean vs variance: generally very similar, therefore poisson regression is reasonable
Mean_vs_Variance_plot<-n_near_homoplasmic_per_sample%>%
  group_by(exp_ID)%>%
  dplyr::summarise(dataset=dataset[1],n=n(),Mean=mean(n_homo),Variance=var(n_homo))%>%
  filter(n>10)%>%
  mutate(dataset=name_conversion_vec[dataset])%>%
  ggplot(aes(x=Mean,y=Variance,col=dataset))+
  ggrepel::geom_label_repel(aes(label = exp_ID),box.padding   = 0.35, point.padding = 0.5, size=2, segment.color = 'grey50',show.legend=F)+
  geom_point(size=0.5)+
  geom_abline(linetype=2)+
  scale_color_manual(values=tissue_cols)+
  scale_x_continuous(limits=c(0,1))+
  scale_y_continuous(limits=c(0,1))+
  theme_classic()+
  my_theme+
  labs(col="")
ggsave(plot = Mean_vs_Variance_plot,filename = paste0(plots_dir,"homo_muts_mean_vs_variance.pdf"),width=3,height=2)

minimum_age_for_regression=20
n_near_homoplasmic_per_sample_FILT<-n_near_homoplasmic_per_sample%>%filter(Age>=minimum_age_for_regression)%>%mutate(logAge=log(Age+0.1))
glmer.res<-lme4::glmer(n_homo~Age+dataset+(1+Age|exp_ID),family=poisson(link="log"),data=n_near_homoplasmic_per_sample_FILT)
summary(glmer.res)

#glmer.with.interaction.res<-lme4::glmer(n_homo~Age+Age*dataset+dataset+(1|Sample),family=poisson(link="log"),data=n_near_homoplasmic_per_sample)

#Now generate the regression curves for plotting by using the 'predict' function on test data
ages_to_cover<-seq(minimum_age_for_regression,80,2)
tissues_to_cover<-c("blood","lymph","KY","LM")
test_data<-data.frame(dataset=rep(tissues_to_cover,each=length(ages_to_cover)),exp_ID=ids::ids(length(ages_to_cover)*length(tissues_to_cover)),Age=rep(ages_to_cover,times=length(tissues_to_cover)))%>%
  mutate(logAge=log(Age))
test_data$n_homo<-exp(predict(glmer.res,newdata=test_data,allow.new.levels =T))

#Visualize the regression model on top of the data for the selected tissues
#Functions to calculate the upper/ lower CI of the mean from count data
get_upper_CI<-function (X, conf.level=0.95) {alpha = 1 - conf.level; upper <- 0.5 * qchisq(1-alpha/2, 2*X+2)}
get_lower_CI<-function (X, conf.level=0.95) {alpha = 1 - conf.level; lower <- 0.5 * qchisq(alpha/2, 2*X)}

data_vs_regression<-n_near_homoplasmic_per_sample%>%
  group_by(dataset,exp_ID)%>%
  dplyr::summarise(n_samp=n(),mean=mean(n_homo),total_homo=sum(n_homo))%>%
  dplyr::mutate(CI_low=get_lower_CI(total_homo)/n_samp,CI_high=get_upper_CI(total_homo)/n_samp)%>%
    left_join(ref_df,by=c("dataset"="Cohort","exp_ID"="ID"))%>%
  filter(Age>=minimum_age_for_regression)%>%
  filter(dataset%in%tissues_to_cover)%>%
  mutate(dataset=name_conversion_vec[dataset])%>%
  ggplot(aes(x=Age,y=mean,ymin=CI_low,ymax=CI_high,col=dataset))+
  geom_point(size=0.5)+
  geom_line(aes(x=Age,y=mean,col = dataset),data = test_data%>%dplyr::rename("mean"=n_homo)%>%mutate(dataset=name_conversion_vec[dataset]),inherit.aes = F)+
  geom_errorbar(width=2,alpha=0.5)+
  scale_color_manual(values = tissue_cols)+
  labs(x="Age",y="Mean homoplasmic mutations per cell",col="")+
  theme_classic()+
  my_theme
ggsave(plot = data_vs_regression,filename = paste0(plots_dir,"data_vs_regression.pdf"),width=3,height=2)

# Visualize the tissue-specific coefficients (relative to blood)
homo_muts_glm_coefficients<-summary(glmer.res)$coefficients%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="Coefficient")%>%
  filter(grepl("dataset",Coefficient))%>%
  mutate(Coefficient=name_conversion_vec[gsub("dataset","",Coefficient)])%>%
  mutate(Coefficient=factor(Coefficient))%>%
  ggplot(aes(y=forcats::fct_reorder(Coefficient,Estimate),
             col=Coefficient,
             x=Estimate,xmin=(Estimate-1.96*`Std. Error`),
             xmax=(Estimate+1.96*`Std. Error`)))+
  geom_point(size=0.5)+
  scale_color_manual(values = tissue_cols)+
  geom_errorbar(width = 0.3)+
  geom_vline(xintercept = 0,linetype=2)+
  theme_classic()+
  my_theme+
  theme(axis.title.y=element_blank(),legend.position = "none")+
  labs(x="Tissue drift coefficient\n(relative to normal blood)")
ggsave(plot = homo_muts_glm_coefficients,filename = paste0(plots_dir,"homo_muts_glm_coefficients.pdf"),width=2,height=2.2)

#------------------------------------------------------------------------------------------------#
##-------------- Cell type-specific analysis for the lymphoid cells------------------------------
#------------------------------------------------------------------------------------------------#

#Filter out only the blood & lymphoid cell types & enhance with the cell type information for the lymphoid cells
minimum_age_for_regression=20
n_homoplasmic_lymph_type<-n_near_homoplasmic_per_sample_FILT%>%
  filter(dataset%in%c("blood","lymph"))%>%
  left_join(colony_info%>%dplyr::select("Sample"=colony,"exp_ID"=Donor,CellType,Cell.type2),by=c("Sample","exp_ID"))%>%
  mutate(Cell.type2=ifelse(dataset=="blood","HSPC",Cell.type2))

Cell.types<-c("HSPC" ,"Naive B","Naive T","Memory B","Memory T","Treg")
cell_type_cols<-c("#1d6d1f",colorRampPalette(c("#fd81c8","#5f70cc","#aea2eb"))(5))
names(cell_type_cols)<-Cell.types

#Now do the glmer poisson regression model using the cell type, not just lymphoid vs blood
glmer.res_blood_vs_lymph_subdivided<-lme4::glmer(n_homo~Age+Cell.type2+(1+Age|exp_ID),family=poisson(link="log"),data=n_homoplasmic_lymph_type)
summary(glmer.res_blood_vs_lymph_subdivided)

#Now generate the regression curves for plotting by using the 'predict' function on test data
ages_to_cover<-seq(minimum_age_for_regression,80,2)
tissues_to_cover<-unique(n_homoplasmic_lymph_type$Cell.type2)
test_data<-data.frame(Cell.type2=rep(tissues_to_cover,each=length(ages_to_cover)),exp_ID=ids::ids(length(ages_to_cover)*length(tissues_to_cover)),Age=rep(ages_to_cover,times=length(tissues_to_cover)))%>%
  mutate(logAge=log(Age))
test_data$n_homo<-exp(predict(glmer.res_blood_vs_lymph_subdivided,newdata=test_data,allow.new.levels =T))

#Visualize the regression model on top of the data for the selected tissues
#Functions to calculate the upper/ lower CI of the mean from count data
get_upper_CI<-function (X, conf.level=0.95) {alpha = 1 - conf.level; upper <- 0.5 * qchisq(1-alpha/2, 2*X+2)}
get_lower_CI<-function (X, conf.level=0.95) {alpha = 1 - conf.level; lower <- 0.5 * qchisq(alpha/2, 2*X)}

data_vs_regression_lymph_subdivided<-n_homoplasmic_lymph_type%>%
  group_by(Age,Cell.type2,exp_ID)%>%
  dplyr::summarise(Age=Age[1],n_samp=n(),mean=mean(n_homo),total_homo=sum(n_homo))%>%
  dplyr::mutate(CI_low=get_lower_CI(total_homo)/n_samp,CI_high=get_upper_CI(total_homo)/n_samp)%>%
  filter(Age>=minimum_age_for_regression)%>%
  ggplot(aes(x=Age,y=mean,ymin=CI_low,ymax=CI_high,col=Cell.type2))+
  geom_point(size=0.5)+
  geom_line(aes(x=Age,y=n_homo,col = Cell.type2),data = test_data,inherit.aes = F)+
  geom_errorbar(width=2,alpha=0.5)+
  scale_color_manual(values = cell_type_cols)+
  labs(x="Age",y="Mean homoplasmic mutations per cell",col="")+
  theme_classic()+
  my_theme
ggsave(plot = data_vs_regression_lymph_subdivided,filename = paste0(plots_dir,"data_vs_regression_lymph_subdivided.pdf"),width=3,height=2)


# Visualize the lymphoid celltype-specific coefficients (relative to blood)
blood_vs_lymph_subdivided_homo_muts_glm_coefficients<-summary(glmer.res_blood_vs_lymph_subdivided)$coefficients%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="Coefficient")%>%
  filter(grepl("Cell.type2",Coefficient))%>%
  mutate(Coefficient=gsub("Cell.type2","",Coefficient))%>%
  mutate(Coefficient=factor(Coefficient))%>%
  ggplot(aes(y=forcats::fct_reorder(Coefficient,Estimate),
             col=Coefficient,
             x=Estimate,xmin=(Estimate-1.96*`Std. Error`),
             xmax=(Estimate+1.96*`Std. Error`)))+
  geom_point(size=0.5)+
  scale_color_manual(values = cell_type_cols)+
  geom_errorbar(width = 0.3)+
  geom_vline(xintercept = 0,linetype=2)+
  theme_classic()+
  my_theme+
  theme(axis.title.y=element_blank(),legend.position = "none")+
  labs(x="Tissue drift coefficient\n(relative to normal HSPCs)")
ggsave(plot = blood_vs_lymph_subdivided_homo_muts_glm_coefficients,filename = paste0(plots_dir,"blood_vs_lymph_subdivided_homo_muts_glm_coefficients.pdf"),width=2,height=2.2)


#------------------------------------------------------------------------------------------------#
##-------------- Driver mutation-specific analysis for the HSPCs cells------------------------------
#------------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------------#
##--------OLD ANALYSIS - ?CAN DELETE, FROM WHEN ATTEMPTING ABC FRAMEWORK------------------------
#------------------------------------------------------------------------------------------------#

#Exclude the lowest VAF categories from abc that may be unreliably captured by sequencing
new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
VAF_categories_to_include_in_ABC=new_VAF_groups[7:11]

exp_ID_by_age<-ref_df%>%filter(ID%in%rownames(tissue_df) & Cohort==dataset)%>%
  arrange(Age)%>%pull(ID)

VAF_dist_plot<-as.data.frame(tissue_df[,VAF_categories_to_include_in_ABC])%>%t()%>%
  as.data.frame()%>%tibble::rownames_to_column(var="VAF_group")%>%
  tidyr::gather(-VAF_group,key="exp_ID",value="n_per_sample")%>%
  ggplot(aes(x=factor(VAF_group,levels=new_VAF_groups),y=n_per_sample))+
  geom_bar(stat="identity")+
  facet_grid(rows=vars(factor(exp_ID,levels=exp_ID_by_age)))+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),strip.text.y=element_text(angle=0))+
  labs(x="Heteroplasmy level",y="mtDNA mutations per cell")

ggsave(filename=paste0(figures_dir,"VAF_dist_plot_",dataset,".pdf"),VAF_dist_plot,width=2.5,height=10)


##----------------------PERFORM THE ABC----------------------
#This section uses summary statistics collected from many simulations and compares the results to those from the
#data stastics for each individual (as collected above).
#The actual simulations are run in the separate script "Mito_VAF_distribution_simulations.R"


library(abc)
root_dir="/lustre/scratch126/casm/team154pc/ms56/mito_mutations_blood/"
all_sumstats=readRDS(file=paste0(root_dir,"data/Drift_ABC/VAF_distribution_ABC_simulation_sumstats_combined.Rds"))
#all_sumstats=readRDS(file = "/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/VAF_distribution_ABC_simulation_sumstats_COMBINED.Rds")

#sim_framework="varying_mtDNA_CN"
sim_framework="varying_mut_rate"

if(sim_framework=="varying_mtDNA_CN") {
  all_sumstats=readRDS(file = "/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/VAF_distribution_v2_COMBINED.Rds")
} else {
  all_sumstats=readRDS(file = "/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/VAF_distribution_ABC_simulation_sumstats_COMBINED.Rds")
  all_sumstats<-all_sumstats%>%mutate(mito_copy_number=600,log10_mito_copy_number=2.778,.before=1)
}

all_sumstats<-all_sumstats%>%
  replace(is.na(.),0)%>%
  filter(muts_per_mitochondria_per_generation>0.0004 & muts_per_mitochondria_per_generation<0.0006)

normalize=F
param_select_vec=c(3,4)
abc_res_list<-lapply(1:nrow(tissue_df),function(i) {
  cat(i,sep="\n")
  Exp_ID=rownames(tissue_df)[i]
  
  total_samples<-length(dataset_mito_data[[Exp_ID]]$tree$tip.label)
  mut_burden<-get_mut_burden(dataset_mito_data[[Exp_ID]]$tree.ultra)[1]
  expanded_clade_nodes<-get_expanded_clade_nodes(tree=dataset_mito_data[[Exp_ID]]$tree.ultra,min_samples=2,height_cut_off=mut_burden/2)
  if(nrow(expanded_clade_nodes)>0) {
    samples_to_drop=unlist(lapply(expanded_clade_nodes$nodes,function(node) {
      clade_samples<-getTips(tree=dataset_mito_data[[Exp_ID]]$tree.ultra,node)
      to_drop=sample(clade_samples,size=length(clade_samples)-1)
    }))
    n_samples<-total_samples-length(samples_to_drop)
  } else {
    n_samples<-total_samples
  }
  
  
  #Due to low sample numbers per individual, the counts per cell numbers will have more statistical variation.
  #Need to reflect this in the simulations
  all_sumstats_mod<-lapply(1:nrow(all_sumstats),function(j) {
    num_vec<-as.numeric(all_sumstats[j,VAF_categories_to_include_in_ABC])
    new_vec<-sapply(num_vec,function(p) {sum(rpois(n=n_samples,lambda = p))/n_samples})
    names(new_vec)<-VAF_categories_to_include_in_ABC
    return(new_vec)
  })%>%bind_rows()
  
  data_target<-as.numeric(tissue_df[i,VAF_categories_to_include_in_ABC])
  
  
  sumstats<-as.matrix(all_sumstats_mod[,VAF_categories_to_include_in_ABC])
  
  if(normalize) {
    data_target=data_target/sum(data_target)
    sumstats<-t(apply(sumstats,1,function(x) x/sum(x)))
  }
  
  res<-abc(target = data_target,
           param = all_sumstats[,param_select_vec],
           sumstat = sumstats,
           tol = 0.05,
           transf = c("none","log"),
           method = "neuralnet")
  
  return(res)
})

#Now extract either adjusted (i.e. with neural network regression) or unadjusted (i.e. with 'rejection' method) posterior values
method="rejection"
if(sim_framework=="varying_mtDNA_CN") {
  parameter_values<-c("log10_mito_copy_number","total_generations","drift_param")
} else {
  parameter_values<-c("muts_per_mitochondria_per_generation","total_generations")
}

abc_stats_res_df<-Map(res=abc_res_list,Exp_ID=rownames(tissue_df),function(res,Exp_ID) {
  if(method=="nn") {res_df<-as.data.frame(res$adj.values)} else if(method=="rejection") {res_df<-as.data.frame(res$unadj.values)}
  
  if(sim_framework=="varying_mtDNA_CN"){
    res_df$drift_param<-apply(res_df,1,function(x) x[2]/x[1])
  }
  
  stats=apply(res_df,2,quantile,c(0.025,0.5,0.975))
  stats=as.data.frame(stats)%>%tibble::rownames_to_column(var="quantile")%>%pivot_wider(names_from = "quantile",values_from=all_of(parameter_values))
  return(cbind(data.frame(exp_ID=Exp_ID),stats))
})%>%dplyr::bind_rows()%>%
  left_join(ref_df,by=c("exp_ID"="ID"))


abc_stats_res_df%>%
  ggplot(aes(x=Age,y=`drift_param_50%`,ymin=`drift_param_2.5%`,ymax=`drift_param_97.5%`))+
  geom_point()+
  geom_errorbar()+
  scale_x_continuous(limits=c(0,80))+
  scale_y_continuous(limits=c(0,2000))+
  theme_classic()+
  my_theme

abc_stats_res_df%>%
  ggplot(aes(x=Age,y=`total_generations_50%`,ymin=`total_generations_2.5%`,ymax=`total_generations_97.5%`))+
  geom_point()+
  geom_errorbar()+
  scale_x_continuous(limits=c(0,80))+
  scale_y_continuous(limits=c(0,3000))+
  theme_classic()+
  my_theme

method="rejection"
Map(res=abc_res_list,Exp_ID=rownames(tissue_df),function(res,Exp_ID) {
  if(method=="nn") {res_df<-as.data.frame(res$adj.values)} else if(method=="rejection") {res_df<-as.data.frame(res$unadj.values)}
  res_df%>%mutate(exp_ID=Exp_ID,.before=1)
})%>%dplyr::bind_rows()%>%
  ggplot(aes(x=log10_mito_copy_number,y=total_generations,col=factor(exp_ID)))+
  geom_point(size=0.5)+
  theme_classic()+
  my_theme+
  facet_wrap(~exp_ID)

method="nn"
Map(res=abc_res_list,Exp_ID=rownames(tissue_df),function(res,Exp_ID) {
  if(method=="nn") {res_df<-as.data.frame(res$adj.values)} else if(method=="rejection") {res_df<-as.data.frame(res$unadj.values)}
  res_df%>%mutate(exp_ID=Exp_ID,.before=1)
})%>%dplyr::bind_rows()%>%
  ggplot(aes(x=muts_per_mitochondria_per_generation,y=total_generations,col=factor(exp_ID)))+
  geom_point(size=0.5)+
  theme_classic()+
  scale_y_continuous(limits=c(0,3000))+
  my_theme+
  theme(axis.text.x=element_text(angle=90))+
  labs(x="Mutations per mtDNA molecular per generation",y="Total generations")+
  facet_wrap(~factor(exp_ID,levels=exp_ID_by_age))

abc_res_plot_mut_rate<-abc_stats_res_df%>%
  filter(Age>=0)%>%
  ggplot(aes(x=forcats::fct_reorder(exp_ID,Age),y=`muts_per_mitochondria_per_generation_50%`,
             ymin=`muts_per_mitochondria_per_generation_2.5%`,
             ymax=`muts_per_mitochondria_per_generation_97.5%`))+
  geom_smooth(method="lm",col="black",linewidth=0.5)+
  geom_point(alpha=0.75,size=0.5)+
  geom_errorbar(width=0.3,alpha=0.5)+
  theme_classic()+
  #scale_y_continuous(limits=c(1e-4,8e-4),breaks=seq(2e-4,1e-3,2e-4))+
  labs(x="Individual",y="Mutations per mtDNA molecule \nper generation")+
  my_theme+
  coord_flip()

ggsave(filename=paste0(figures_dir,"abc_res_plot_mut_rate_",dataset,"_",method,".pdf"),abc_res_plot_mut_rate,width=3.5,height=2)

##Now do linear mixed effects modelling using all values of the posterior distribution and using individual as a random effect
abc_res_df<-Map(res=abc_res_list,Exp_ID=rownames(tissue_df),function(res,Exp_ID) {
  if(method=="nn") {res_df<-as.data.frame(res$adj.values)} else if(method=="rejection") {res_df<-as.data.frame(res$unadj.values)}
  res_df%>%mutate(exp_ID=Exp_ID,.before=1)
})%>%dplyr::bind_rows()%>%
  left_join(ref_df%>%dplyr::select(ID,Age),by=c("exp_ID"="ID"))


gens_by_age_plot<-abc_res_df%>%
  mutate(exp_ID=factor(exp_ID,levels=ref_df%>%filter(Cohort==dataset & ID%in% donors_to_include)%>%arrange(Age)%>%pull(ID)))%>%
  ggplot(aes(x=Age,y=total_generations,col=exp_ID))+
  geom_jitter(size=0.1,height=0,width=0.2,alpha=0.3)+
  #geom_point(size=0.1)+
  theme_classic()+
  my_theme+
  scale_x_continuous(limits = c(0,100))+
  scale_y_continuous(limits=c(0,4000))+
  labs(col="Sample",x="Age",y="Total generations")

ggsave(filename=paste0(figures_dir,"gens_by_age_plot_",dataset,"_",method,".pdf"),gens_by_age_plot,width=3.5,height=2)

mut_rate_vs_gens_plot<-abc_res_df%>%
  mutate(exp_ID=factor(exp_ID,levels=ref_df%>%filter(Cohort==dataset & ID%in% donors_to_include)%>%arrange(Age)%>%pull(ID)))%>%
  ggplot(aes(x=muts_per_mitochondria_per_generation,y=total_generations,col=exp_ID))+
  #geom_jitter(size=0.1,height=0,width=2)+
  geom_point(size=1,alpha=0.3,stroke=NA)+
  theme_classic()+
  my_theme+
  facet_wrap(~exp_ID)+
  theme(axis.text.x=element_text(angle=90))+
  labs(col="Sample",x="Mutations per mtDNA molecule/ generation",y="Total generations")

ggsave(filename=paste0(figures_dir,"mut_rate_vs_gens_plot_",dataset,"_",method,".pdf"),mut_rate_vs_gens_plot,width=3.5,height=4)

#Perform the LMER
lme.gens_by_age<-lme4::lmer(total_generations~Age+(1|exp_ID),data=abc_res_df)
summary(lme.gens_by_age)
confint(lme.gens_by_age)
mtDNA_CN*365/lme.gens_by_age@beta[2]
mtDNA_CN*365/confint(lme.gens_by_age)["Age",]
















#Group these to get a profile for an individual/ tissue
mut_numbers=data.frame(SampleID=rownames(mutation_profiles_mat),nmuts=rowSums(mutation_profiles_mat))

#Fill in the rare bins with too low mutation numbers - assume all are N1, as these are high VAF mutations
excluded_cats<-rownames(mutation_profiles_mat)[!rownames(mutation_profiles_mat)%in%colnames(exposures)]
all_N1=c(0,1,rep(0,nrow(exposures)-2))
excluded_mat<-matrix(rep(all_N1,times=length(excluded_cats)),nrow=length(all_N1),dimnames = list(rownames(exposures),excluded_cats))

real_muts_dist_df<-t(cbind(exposures,excluded_mat))%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="SampleID")%>%
  left_join(mut_numbers)%>%
  mutate(exp_ID=stringr::str_split(SampleID,pattern="_",simplify=T)[,1])%>%
  mutate(VAF_range=sapply(stringr::str_split(SampleID,pattern="_"),function(vec) paste(vec[2:3],collapse="_")))%>%
  mutate(VAF_range=factor(new_VAF_groups[VAF_range],levels=new_VAF_groups))%>%
  dplyr::select(-SampleID)%>%
  gather(-VAF_range,-exp_ID,-nmuts,key="Signature",value="Exposure")%>%
  mutate(abs_muts=Exposure*nmuts)%>%
  mutate(Signature=factor(Signature,levels=c(rownames(exposures)[nrow(exposures):3],"N0","N1")))%>%
  dplyr::filter(Signature=="N1")%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))

real_mut_distribution_plot<-real_muts_dist_df%>%
  ggplot(aes(x=VAF_range,y=abs_muts))+
  geom_bar(stat="identity",fill="#FDBF6F",col="black",size=0.25)+
  geom_line(aes(group=exp_ID))+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(size=5,angle=90),strip.text.x = element_text(size=6))+
  #facet_grid(rows=vars(factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])),scales="free")+
  facet_wrap(~factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)]),scales="fixed",dir="v",ncol=6)+
  labs(x="VAF range",y="Number of N1 mutations detected")

#Normalize mutation numbers by the number of samples included
nsamp_df<-Map(list=mito_data,Exp_ID=names(mito_data),function(list,Exp_ID) {
  data.frame(exp_ID=Exp_ID,n_samp=length(colnames(list$matrices$SW)))
})%>%dplyr::bind_rows()%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))

real_mut_distribution_plot_normalized<-left_join(real_muts_dist_df,nsamp_df)%>%
  mutate(muts_per_samp=abs_muts/n_samp)%>%
  ggplot(aes(x=VAF_range,y=muts_per_samp))+
  geom_bar(stat="identity",fill="#FDBF6F",col="black",size=0.25)+
  geom_line(aes(group=exp_ID))+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(size=5,angle=90),strip.text.x = element_text(size=6))+
  #facet_grid(rows=vars(factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])),scales="free")+
  facet_wrap(~factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)]),scales="fixed",dir="v",ncol=6)+
  labs(x="VAF range",y="Number of N1 mutations detected per sample")
ggsave(filename=paste0(figures_dir,"Figure_04/real_mut_distribution_plot.pdf"),real_mut_distribution_plot_normalized,width=7,height=2.5)

sumstats.data<-left_join(real_muts_dist_df,nsamp_df)%>%
  mutate(muts_per_samp=abs_muts/n_samp)%>%
  dplyr::select(exp_ID,VAF_range,muts_per_samp)%>%
  pivot_wider(names_from="VAF_range",values_from="muts_per_samp")%>%
  replace(is.na(.), 0)

##----------------------VISUALIZE WRIGHT-FISHER MODELS OF DRIFT----------------------

#This code is to plot the shifting distribution of VAFs over time - single run of simulation but storing information after each generation
muts_per_mitochondria_per_generation=5e-4
print(muts_per_mitochondria_per_generation)
mito_copy_number=600 #This is the
number_of_cells_in_simulation=1000 #Need enough 'cells' to have adequate mutation numbers to define the distribution
mutation_introductions_per_generation=muts_per_mitochondria_per_generation*number_of_cells_in_simulation*mito_copy_number

#Define the VAF 'bins' that will be used for comparing VAF distributions
new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
VAF_groups=data.frame(
  labels=new_VAF_groups,
  lower_limit=c(0,2^(-10:-1)),
  upper_limit=2^(-10:0)
)

#Record the distribution of VAFs every 10 generations, though the generation of mutation acquisition is recorded exactly
gens_to_include=seq(10,1500,10)
gens_record=vector(mode="list",length = length(gens_to_include))
vaf_by_gen_list=vector(mode="list")
last_gen=0
curr_vafs=c()
for(ngen in gens_to_include) {
  print(ngen)
  if(length(vaf_by_gen_list)>0) {
    drifted_old_muts=lapply(vaf_by_gen_list,function(gen_vafs) {
      new_gen_vafs<-sapply(gen_vafs,function(vaf) fisher_wright_drift(vaf,population_size = mito_copy_number,1,ngen-last_gen))
      return(new_gen_vafs[new_gen_vafs>0])
    })
  } else {
    drifted_old_muts<-NULL
  }
  
  new_muts_gen_list=vector(mode="list",length=ngen-last_gen)
  names(new_muts_gen_list)<-(last_gen+1):ngen
  for(j in 1:(ngen-last_gen)){
    gen_mut_vafs<-vector(length=mutation_introductions_per_generation)
    for(i in 1:mutation_introductions_per_generation) {
      this_mut_final_vaf<- fisher_wright_drift(1/mito_copy_number,population_size = mito_copy_number,1,j)
      gen_mut_vafs[i]<-this_mut_final_vaf
    }
    new_muts_gen_list[[j]]<-gen_mut_vafs[gen_mut_vafs>0]
  }
  vaf_by_gen_list<-c(drifted_old_muts,new_muts_gen_list)
  gen_summary=Map(vafs=vaf_by_gen_list,gen=names(vaf_by_gen_list),function(vafs,gen) if(length(vafs)>0){data.frame(gen=gen,vaf=vafs)}else {NULL})%>%dplyr::bind_rows()%>%mutate(total_gens=ngen)
  gens_record[[which(gens_to_include==ngen)]]<-gen_summary
  last_gen<-ngen
}

#Here, 7.5 is the 'haploid sequencing coverage' in this simulated experiment
#This approximates the ~15X diploid coverage of most of the experiments
haploid_coverage=7.5
gens_record<-dplyr::bind_rows(gens_record)
gens_record$observed_vaf=sapply(gens_record$vaf,function(vaf) rbinom(n=1,size=(haploid_coverage*mito_copy_number),prob=vaf)/(haploid_coverage*mito_copy_number))

gens_record$VAF_group=sapply(1:nrow(gens_record),function(i) {
  if(gens_record$observed_vaf[i]==0) {
    return("Absent")
  } else {
    return(VAF_groups$labels[VAF_groups$lower_limit<gens_record$observed_vaf[i] & VAF_groups$upper_limit>=gens_record$observed_vaf[i]])
  }
})

#Visualize the distribution at a selection of "generations" to get a sense of the evolution of the distribution through time
modelled_VAF_dist<-gens_record%>%
  filter(total_gens%in%c(100,200,400,800,1200))%>%
  mutate(total_gens=factor(paste(total_gens,"generations"),levels=paste(c(100,200,400,800,1200),"generations")))%>%
  mutate(VAF_group=factor(VAF_group,levels=new_VAF_groups),gen=factor(gen,levels=1:1200))%>%
  ggplot(aes(x=VAF_group,fill=gen))+
  geom_bar(stat="count")+
  facet_wrap(~total_gens,ncol=5)+
  theme_bw()+
  theme(legend.position = "none",axis.text.x = element_text(size=5,angle=90),strip.text.x = element_text(size=6))+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8,"Spectral"))(1200))+
  labs(x="VAF group",y="Count")+
  my_theme

ggsave(filename=paste0(figures_dir,"Figure_04/Modelled_VAF_dist.pdf"),width = 6,height=2)

#One notable feature is that these simulations suggest that high VAF mutations were acquired early in life
#whereas low VAF mutations were invariably acquired very recently.
acquisition_time_by_VAF<-gens_record%>%
  mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels))%>%
  filter(total_gens==1250)%>%
  ggplot(aes(x=as.numeric(gen)))+
  geom_density(fill="light blue")+
  facet_wrap(~VAF_group,scales="free_y",ncol=6,dir = "v")+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),strip.text.x = element_text(size=7,margin = unit(c(1,0,1,0),"mm")))+
  labs(x="Time of mutation acquisition\n(WF Generations)",y="Density")
ggsave(filename = paste0(figures_dir,"Supp_Figure_04/Acquisition_time_by_VAF.pdf"),acquisition_time_by_VAF,width=7,height=2.5)

## Can also be informative to visualize the evolving distribution over time as an animation
gens_record_summary<-gens_record%>%
  group_by(total_gens,gen,VAF_group)%>%
  summarise(n=n())
gens_record_summary$VAF_group=factor(gens_record_summary$VAF_group,levels=VAF_groups$labels)
gens_record_summary$gen=as.numeric(gens_record_summary$gen)

library(gganimate)
drift_with_age_animation<-gens_record_summary%>%
  filter(!is.na(VAF_group))%>%
  arrange(gen)%>%
  ggplot(aes(x=VAF_group,y=n,fill=gen))+
  geom_bar(stat="identity")+
  transition_manual(frames=total_gens)+
  theme_bw()+
  theme(title = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle=90,size=15),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(8,"Spectral"))+
  labs(title='{current_frame} generations',x="Observed VAF",y="Count")
anim_save(filename=paste0(plots_dir,"driftwithage.gif"),animation = drift_with_age_animation)

#Plot the colour scale for these figures
max_gen=1500
lut=colorRampPalette(RColorBrewer::brewer.pal(8,"Spectral"))(max_gen)
scale = (length(lut)-1)
pdf(paste0(plots_dir,"mito_muts_by_gen_scale.pdf"),height=6,width = 3)
plot(c(0,10), c(0,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
axis(side=2,pos=1,las=1,at=seq(0,max_gen,200)/scale,labels=seq(0,max_gen,200))
for (i in 1:(length(lut)-1)) {
  y = (i-1)/scale
  rect(1,y,2,y+1/scale, col=lut[i], border=NA)
}
dev.off()

##----------------------PERFORM THE ABC----------------------
#This section uses summary statistics collected from many simulations and compares the results to those from the
#data stastics for each individual (as collected above).
#The actual simulations are run in the separate script "Mito_VAF_distribution_simulations.R"

library(abc)
all_sumstats=readRDS(file=paste0(root_dir,"data/Drift_ABC/VAF_distribution_ABC_simulation_sumstats_combined.Rds"))

new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")

#Exclude the lowest VAF categories from abc that may be unreliably captured by sequencing
VAF_categories_to_include_in_ABC=new_VAF_groups[3:11]

mtDNA_CN=600
abc_res_df<-lapply(1:nrow(sumstats.data),function(i) {
  Exp_ID=sumstats.data[i,1]
  res<-abc(target = as.numeric(sumstats.data[i,VAF_categories_to_include_in_ABC]),
           param = all_sumstats[,1:2],
           sumstat = all_sumstats[,VAF_categories_to_include_in_ABC],
           tol = 0.02,
           transf = c("log","log"),
           method = "neuralnet")
  stats=apply(res$adj.values,2,quantile,c(0.025,0.5,0.975))
  stats=as.data.frame(stats)%>%tibble::rownames_to_column(var="quantile")%>%pivot_wider(names_from = "quantile",values_from=c("muts_per_mitochondria_per_generation","total_generations"))
  return(cbind(data.frame(exp_ID=Exp_ID),stats))
})%>%dplyr::bind_rows()%>%
  left_join(ref_df,by=c("exp_ID"="Sample"))

abc_res_plot_mut_rate<-abc_res_df%>%
  filter(Age>=0)%>%
  ggplot(aes(x=forcats::fct_reorder(exp_ID,Age),y=`muts_per_mitochondria_per_generation_50%`,
             ymin=`muts_per_mitochondria_per_generation_2.5%`,
             ymax=`muts_per_mitochondria_per_generation_97.5%`))+
  geom_smooth(method="lm",col="black",size=0.5)+
  geom_point(alpha=0.75,size=0.5)+
  geom_errorbar(width=0.3,alpha=0.5)+
  theme_classic()+
  scale_y_continuous(limits=c(1e-4,8e-4),breaks=seq(2e-4,1e-3,2e-4))+
  labs(x="Individual",y="Mutations per mitochondria\nper generation")+
  my_theme+
  coord_flip()

ggsave(filename=paste0(figures_dir,"Figure_04/abc_res_plot_mut_rate.pdf"),abc_res_plot_mut_rate,width=1.75,height=2)

##Now do linear mixed effects modelling using all values of the posterior distribution and using individual as a random effect
abc_res_df<-lapply(1:nrow(sumstats.data),function(i) {
  Exp_ID=sumstats.data[i,1]
  res<-abc(target = as.numeric(sumstats.data[i,VAF_categories_to_include_in_ABC]),
           param = all_sumstats[,1:2],
           sumstat = all_sumstats[,VAF_categories_to_include_in_ABC],
           tol = 0.05,
           transf = c("log","log"),
           method = "neuralnet")
  
  return(as.data.frame(res$adj.values)%>%mutate(exp_ID=Exp_ID$exp_ID))
})%>%dplyr::bind_rows()%>%
  left_join(ref_df%>%dplyr::select(Sample,Age),by=c("exp_ID"="Sample"))%>%
  dplyr::filter(!exp_ID%in%c("8pcw","18pcw"))

#Perform the LMER
lme.gens_by_age<-lme4::lmer(total_generations~Age+(1|exp_ID),data=abc_res_df)
summary(lme.gens_by_age)
mtDNA_CN*365/lme.gens_by_age@beta[2]
mtDNA_CN*365/confint(lme.gens_by_age)["Age",]

#Plot these results
abc_res_plot_total_generations<-abc_res_df%>%
  dplyr::filter(!exp_ID%in%c("8pcw","18pcw"))%>%
  mutate(exp_ID=factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)]))%>%
  ggplot(aes(x=Age,y=total_generations))+
  geom_point(aes(col=exp_ID),alpha=0.05,size=0.25)+
  scale_y_continuous(limits=c(0,1700))+
  geom_abline(slope=lme.gens_by_age@beta[2],intercept = lme.gens_by_age@beta[1],linetype=1)+
  scale_color_manual(values=Individual_cols[-c(1:2)])+
  theme_classic()+
  guides(colour=guide_legend(override.aes = list(alpha=1)))+
  labs(x="Age",y="Total WF generations\n(Posterior distribution from ABC)",col="")+
  my_theme+
  theme(legend.box.spacing = unit(0,"mm"),
        legend.key.size = unit(0.5,"mm"),
        legend.spacing = unit(0,"mm"),
        legend.box.margin = margin(c(0,0,0,0)),
        legend.text = element_text(margin = margin(t=0)))

ggsave(filename=paste0(figures_dir,"Figure_04/abc_res_plot_total_generations.pdf"),abc_res_plot_total_generations,width=2,height=2)









ggsave(filename=paste0("blood/VAF_dist_plot2_lymph_IDs_only.pdf"),VAF_dist_plot,width=2.5,height=(0.8*nrow(tissue_df)))
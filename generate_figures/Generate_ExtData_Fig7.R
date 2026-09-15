#-----------------------------------------------------------------------------------#
# Generate_ExtData_Fig7.R
#
# Extended Data Fig. 7 - Mitochondrial mutation drift across tissues.
#
#   a  Wright-Fisher drift: time of acquisition by VAF range
#      (Mitochondrial_drift_analysis.Rmd, also written by Generate_Fig4.R)
#   b  heteroplasmy distribution by driver status, bootstrap CIs
#   c  heteroplasmy distribution by clonal expansion status, bootstrap CIs
#      (b and c follow full_analysis_scripts/archive/Mitochondrial_mut_analysis_Apr2026.R)
#   d  mutations per cell by heteroplasmy, HSPC vs lymphoid
#   e  mean homoplasmic burden against age, by blood compartment
#   f  Poisson mixed model across tissues
#   g  as f, lymphoid split by subset
#   h  tissue-specific model parameters
#   i  lymphoid subset model parameters
#      (d-i follow Nonblood_mtDNA_drift_analysis.Rmd)
#-----------------------------------------------------------------------------------#

options(stringsAsFactors = F)
source(here::here("config.R")) #sets root_dir, genomeFile, project paths and my_theme
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))
#plots_dir, rebuttal_figs_dir and my_theme all come from config.R

load_packages(cran = c("dplyr","tidyr","ggplot2","stringr","readxl","forcats",
                       "ape","phangorn","lme4","lmerTest","ggridges","RColorBrewer","patchwork"))

ed7_dir <- paste0(plots_dir,"Extended_Data_Figure_07/")
dir.create(ed7_dir, showWarnings = FALSE, recursive = TRUE)

#The bootstraps in panels b and c resample with replacement; fix the seed so the
#intervals are reproducible.
set.seed(42)



#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 7A ---------------
## Wright-Fisher drift: time of mutation acquisition by VAF range
##
## The simulation is the one used for Fig. 4b; it is repeated here so that this
## script produces its own panel rather than depending on Generate_Fig4.R.
#-----------------------------------------------------------------------------------#


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

#The Wright-Fisher simulation below is stochastic (fisher_wright_drift), as is
#the binomial resampling of its output, so fix the seed here - before the loop -
#to make Fig 4b and Extended Data Fig 7a exactly reproducible.
set.seed(42)

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

acquisition_time_by_VAF_ridges<-gens_record%>%
  mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels))%>%
  filter(total_gens==1250)%>%
  ggplot(aes(x=as.numeric(gen),y=VAF_group,fill = factor(stat(quantile),levels=1:4))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,quantiles = 4, quantile_lines = TRUE,scale=4) +
  scale_fill_viridis_d(name = "Quartiles")+
  scale_x_continuous(limits=c(0,1350),breaks=seq(0,1250,250))+
  theme_classic()+
  my_theme+
  labs(x="Time of mutation acquisition\n(WF Generations)",y="VAF level")

ggsave(filename = paste0(ed7_dir,"ExtDataFig7a.acquisition_time_by_VAF_ridges.pdf"),acquisition_time_by_VAF_ridges,width=3.3,height=2.5)


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 7B and 7C ---------------
## Heteroplasmy distribution split by driver status (b) and by clonal expansion
## status (c), with confidence intervals from 100 bootstrap resamples of the
## sample set. Follows archive/Mitochondrial_mut_analysis_Apr2026.R.
##
## These need the blood mutation data and the nuclear driver annotation, which is
## attached below from data/blood_adult/annot_files_filtered/.
##
## Both panels use all eight adult donors, mito_data[3:10] (the archived script
## used mito_data[5:8], an arbitrary four of them; with only four donors just two
## driver-carrying clades are found, which makes panel b far noisier than the
## submitted version).
#-----------------------------------------------------------------------------------#

mito_data<-readRDS(paste0(root_dir,"/data/mito_data.Rds"))
CN_correlating_muts<-readRDS(paste0(root_dir,"/data/CN_correlation.RDS"))
#Blood sample metadata, needed to give mito_cn its exp_ID before the join below
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
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
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C")

#Attach the nuclear driver annotation, which panel b splits on. nDNA_mats is not
#part of the deposited mito_data.Rds, so it is read from the annotation files.
mito_data<-Map(list=mito_data,exp_ID=names(mito_data),function(list,exp_ID) {
  annot_muts_file<-paste0(root_dir,"/data/blood_adult/annot_files_filtered/annotated_muts_filt_",exp_ID,".Rds")
  if(file.exists(annot_muts_file)) list$nDNA_mats<-readRDS(annot_muts_file)
  return(list)
})

## Generate tidy data frame of samples, mutations and vafs-------------------------------

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

boundaries<-signif(c(0,2^(-10:0)),2)
new_VAF_groups=paste0(100*head(boundaries,-1),"-",100*tail(boundaries,-1),"%")
VAF_groups=data.frame(labels=new_VAF_groups,LL=head(boundaries,-1),UL=tail(boundaries,-1))

## 2. Assign each observed VAF into its 'VAF group' bin ----
df_tidy$VAF_group=sapply(df_tidy$vaf,function(vaf) {
  if(vaf==0) {
    return("Absent")
  } else {
    return(VAF_groups$labels[VAF_groups$LL<vaf & VAF_groups$UL>=vaf])
  }
})

## 3. Make VAF distribution matrices for each individual: each row is a sample, and each column is a VAF group, each value is number of mutations
muts_per_VAF_per_sample_mat<-Map(exp_ID=names(mito_data)[3:10],list=mito_data[3:10],function(exp_ID,list) {
  cat(exp_ID)
  exp_ID_list<-lapply(list$tree$tip.label,function(SampleID) {
    
    mat<-dplyr::select(VAF_groups,"VAF_group"=labels)%>%
      left_join(df_tidy%>%filter(Sample==SampleID)%>%group_by(VAF_group)%>%summarise(n=n()),by="VAF_group")%>%
      replace_na(list(n=0))%>%
      tibble::column_to_rownames(var="VAF_group")%>%t()
    rownames(mat)<-SampleID
    return(mat)
  })
  exp_ID_mat<-Reduce(rbind,exp_ID_list)
  return(exp_ID_mat)
})

muts_per_VAF_per_sample_mat_combined<-Reduce(rbind,muts_per_VAF_per_sample_mat)


mut_nodes_list<-Map(list=mito_data[3:10],exp_ID=names(mito_data)[3:10],function(list,exp_ID) {
  cat(exp_ID,sep="\n")
  #Using the 'get_expanded_clade_nodes' function to effectively cut across the tree at 100 mutations of molecular time
  #Some of the clades will be singletons
  ec_df<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = 100,min_clonal_fraction = 0,min_samples=1)%>%
    left_join(list$nDNA_mats$mat%>%dplyr::select(node,mut_ref,Gene,CDS,Protein),by=c("nodes"="node"))%>%
    mutate(exp_ID=exp_ID,.before=1)

  genes_of_interest=c("TET2","DNMT3A","ASXL1")
  mut_nodes<-ec_df%>%filter(Gene%in%genes_of_interest)
  return(mut_nodes)
})%>%dplyr::bind_rows()

## 3. Make a list of samples with/ without driver mutations----
# Each cell from a clone is a non-independent 
nboot=100
all_boot_res<-lapply(1:nboot,function(i) {
  cat(i,sep="\n")
  mut_vs_wt_list<-Map(list=mito_data[3:10],exp_ID=names(mito_data)[3:10],function(list,exp_ID) {
    cat(exp_ID,sep="\n")
    #Using the 'get_expanded_clade_nodes' function to effectively cut across the tree at 100 mutations of molecular time
    #Some of the clades will be singletons
    ec_df<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = 100,min_clonal_fraction = 0,min_samples=1)%>%
      left_join(list$nDNA_mats$mat%>%dplyr::select(node,mut_ref,Gene,CDS,Protein),by=c("nodes"="node"))%>%
      mutate(exp_ID=exp_ID,.before=1)
    
    genes_of_interest=c("TET2","DNMT3A","ASXL1")
    mut_nodes<-ec_df%>%filter(Gene%in%genes_of_interest)%>%pull(nodes)%>%unique()
    wt_nodes<-ec_df$nodes[!ec_df$nodes%in%mut_nodes]
    
    # mut_nodes<-ec_df%>%filter(n_samples>1)%>%pull(nodes)
    # wt_nodes<-ec_df%>%filter(n_samples==1)%>%pull(nodes)
    
    mut_samples<-unlist(lapply(mut_nodes,function(node) sample(x=getTips(list$tree.ultra,node=node),size=1)))
    wt_samples<-unlist(lapply(wt_nodes,function(node) sample(x=getTips(list$tree.ultra,node=node),size=1)))
    
    return(list(mut=mut_samples,wt=wt_samples))
  })
  
  all_mut<-unlist(lapply(mut_vs_wt_list,function(x) x$mut))
  all_wt<-unlist(lapply(mut_vs_wt_list,function(x) x$wt))
  
  ## Make a VAF distribution matrix of the mutant vs wildtype ----
  #Iterate across samples/ VAF groups to get numbers of mutations per sample/ VAF group
  VAF_distribution_by_mut_status<-Map(sample_set=list(all_mut,all_wt),status=c("mutated","wild-type"), function(sample_set,status) {
    
    sample_set=sample(sample_set,replace=T)
    
    as.data.frame(colSums(muts_per_VAF_per_sample_mat_combined[sample_set,])/length(sample_set))%>%
      tibble::rownames_to_column(var="VAF_group")%>%
      dplyr::rename("n_per_sample"=2)%>%
      mutate(status=status)
    
  })%>%dplyr::bind_rows()

  return(VAF_distribution_by_mut_status%>%dplyr::select(n_per_sample))
})%>%bind_cols()

#Plot the results of these bootstraps
VAF_dist_mut_vs_wt<-as.matrix(all_boot_res)%>%
  apply(1,function(x) quantile(x, c(0.025,0.5,0.975)))%>%t()%>%
  as.data.frame()%>%
  bind_cols(data.frame(status=rep(c("DNMT3A/TET2/ASXL1\nmutated","Wild-type"),each=11),VAF_group=rep(VAF_groups$labels,times=2)))%>%
  mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels))%>%
  filter(VAF_group%in%VAF_groups$labels[4:11])%>%
  ggplot(aes(y=`50%`,ymin=`2.5%`,ymax=`97.5%`,x=VAF_group,fill=status))+
  geom_bar(stat="identity",position="dodge",alpha=0.5)+
  scale_y_continuous(limits=c(0,1))+
  #facet_grid(~status)+
  geom_point(position= position_dodge2(width=0.75),size=0.75)+
  geom_linerange(position= position_dodge2(width=0.75), linewidth=0.5,color="gray50") +
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),legend.position="right")+
  labs(x="VAF group",y="Average number of mutations per cell",fill="")

ggsave(filename = paste0(ed7_dir,"ExtDataFig7b.VAF_dist_mut_vs_wt.pdf"),VAF_dist_mut_vs_wt,width=3.3,height=2.5)


#Repeat analysis but dividing by singleton vs member of clonal expansion
nboot=100
all_boot_res_clone_vs_singleton<-lapply(1:nboot,function(i) {
  cat(i,sep="\n")
  mut_vs_wt_list<-Map(list=mito_data[3:10],exp_ID=names(mito_data)[3:10],function(list,exp_ID) {
    cat(exp_ID,sep="\n")
    #Using the 'get_expanded_clade_nodes' function to effectively cut across the tree at 100 mutations of molecular time
    #Some of the clades will be singletons
    ec_df<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = 100,min_clonal_fraction = 0,min_samples=1)%>%
      left_join(list$nDNA_mats$mat%>%dplyr::select(node,mut_ref,Gene,CDS,Protein),by=c("nodes"="node"))%>%
      mutate(exp_ID=exp_ID,.before=1)
    
    # genes_of_interest=c("TET2","DNMT3A","ASXL1")
    # mut_nodes<-ec_df%>%filter(Gene%in%genes_of_interest)%>%pull(nodes)%>%unique()
    # wt_nodes<-ec_df$nodes[!ec_df$nodes%in%mut_nodes]
    
    mut_nodes<-ec_df%>%filter(n_samples>1)%>%pull(nodes)
    wt_nodes<-ec_df%>%filter(n_samples==1)%>%pull(nodes)
    
    mut_samples<-unlist(lapply(mut_nodes,function(node) sample(x=getTips(list$tree.ultra,node=node),size=1)))
    wt_samples<-unlist(lapply(wt_nodes,function(node) sample(x=getTips(list$tree.ultra,node=node),size=1)))
    
    return(list(mut=mut_samples,wt=wt_samples))
  })
  
  all_mut<-unlist(lapply(mut_vs_wt_list,function(x) x$mut))
  all_wt<-unlist(lapply(mut_vs_wt_list,function(x) x$wt))
  
  ## Make a VAF distribution matrix of the mutant vs wildtype ----
  #Iterate across samples/ VAF groups to get numbers of mutations per sample/ VAF group
  VAF_distribution_by_mut_status<-Map(sample_set=list(all_mut,all_wt),status=c("mutated","wild-type"), function(sample_set,status) {
    
    sample_set=sample(sample_set,replace=T)
    
    as.data.frame(colSums(muts_per_VAF_per_sample_mat_combined[sample_set,])/length(sample_set))%>%
      tibble::rownames_to_column(var="VAF_group")%>%
      dplyr::rename("n_per_sample"=2)%>%
      mutate(status=status)
    
  })%>%dplyr::bind_rows()
  
  return(VAF_distribution_by_mut_status%>%dplyr::select(n_per_sample))
  
})%>%bind_cols()

VAF_dist_exp_vs_singleton<-as.matrix(all_boot_res_clone_vs_singleton)%>%
  apply(1,function(x) quantile(x, c(0.025,0.5,0.975)))%>%t()%>%
  as.data.frame()%>%
  bind_cols(data.frame(status=rep(c("clonal expansion","singleton"),each=11),VAF_group=rep(VAF_groups$labels,times=2)))%>%
  mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels[-1]))%>%
  filter(VAF_group%in%VAF_groups$labels[4:11])%>%
  ggplot(aes(y=`50%`,ymin=`2.5%`,ymax=`97.5%`,x=VAF_group,fill=status))+
  geom_bar(stat="identity",position="dodge",alpha=0.5)+
  scale_y_continuous(limits=c(0,1))+
  #facet_grid(~status)+
  geom_point(position= position_dodge2(width=0.75),size=0.75)+
  geom_linerange(position= position_dodge2(width=0.75), linewidth=0.5,color="gray50") +
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),legend.position="right")+
  labs(x="VAF group",y="Average number of mutations per cell",fill="")

ggsave(filename = paste0(ed7_dir,"ExtDataFig7c.VAF_dist_exp_vs_singleton.pdf"),VAF_dist_exp_vs_singleton,width=3.3,height=2.5)

#-----------------------------------------------------------------------------------#
# Cross-tissue data: import, tidy, and per-sample VAF distributions
#-----------------------------------------------------------------------------------#

#Source the necessary functions - included in the repo
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"/data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
nonblood_ref_file=paste0(root_dir,"/data/metadata/non_blood_metadata.xlsx")
colony_info_file=paste0(root_dir,"/data/metadata/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt")

# read in the trinucleotide context data
mtdna_trinuc_freq_path=paste0(root_dir,"/data/mtDNA_trinuc_freqs_coding_dloop_heavy_light.Rds")
mtdna_trinuc_freq <- readRDS(mtdna_trinuc_freq_path)

# set regions for coding and dloop regions - important for the mutation profiles
coding_region <- 577:16023
d_loop_region <- c(1:576,16024:16569)

#Set the basic plotting theme for ggplot2
#my_markdown_theme is defined in config.R


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

#Set up objects for colour schemes and renaming tissues for plots
mut_type_cols<-c("#00468B","#925E9F","#925E9F","#810F7C","black","#E10C04","black")
names(mut_type_cols)<-c("Missense","Nonsense","Truncating","Stop_loss","Non-Coding","Synonymous","Overall")

name_conversion_vec=c("Normal blood","Lymphoid","MPN","Colon","IBD-affected\nColon","MUTYH-mutant\nColon","Endometrium","Bronchial\nepithelium")
names(name_conversion_vec)=c("blood","lymph","NW","HL","SO","PR","LM","KY")
tissue_cols<-c("#96e97c","#fd81c8","#145a6a","#65e6f9","#781486","#5f70cc","#aea2eb","#1d6d1f")
names(tissue_cols)<-name_conversion_vec

all_cohorts_plus_CML<-c("KY","HL","SO","PR","LM","NW","CML","lymph","blood")
all_cohorts<-c("KY","HL","SO","PR","LM","NW","lymph","blood")

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
  
  
  #cat("Adding the filtered VAF matrix.",sep="\n")
  
  #Annotate specific mutations with their most likely signature using the 'sig_ref' dataframe
  dataset_mito_data<-Map(list=dataset_mito_data,this_exp_ID=names(dataset_mito_data),function(list,this_exp_ID) {
    
    #cat(this_exp_ID,sep="\n")
    
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


minimum_samples_to_include=15

#VAF bin definitions used by the distribution matrices below
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

#Iterate across samples/ VAF groups to get numbers of mutations per sample/ VAF group
VAF_distribution_mats<-Map(dataset_mito_data=all_mito_datasets[all_cohorts],df_tidy=all_df_tidy[all_cohorts],function(dataset_mito_data,df_tidy) {
  all_samples=dplyr::bind_rows(Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) data.frame(exp_ID=exp_ID,Sample=list$tree$tip.label)))
  samples_vaf_groups_df<-right_join(all_samples,tidyr::expand_grid(Sample=all_samples$Sample,VAF_group=new_VAF_groups),by="Sample")
  nmut_by_vaf_by_sample<-left_join(samples_vaf_groups_df,df_tidy%>%
                                     group_by(Sample,VAF_group)%>%
                                     dplyr::summarise(n=n(),.groups = "drop_last"),by=c("Sample","VAF_group"))%>%
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

minimum_samples_to_include=8
VAF_distribution_by_sample_df<-Map(dataset=all_cohorts,dataset_mito_data=all_mito_datasets[all_cohorts],df_tidy=all_df_tidy[all_cohorts],function(dataset,dataset_mito_data,df_tidy) {
  all_samples=dplyr::bind_rows(Map(list=dataset_mito_data,exp_ID=names(dataset_mito_data),function(list,exp_ID) data.frame(exp_ID=exp_ID,Sample=list$tree$tip.label)))
  samples_vaf_groups_df<-right_join(all_samples,tidyr::expand_grid(Sample=all_samples$Sample,VAF_group=new_VAF_groups),by="Sample")
  nmut_by_vaf_by_sample<-left_join(samples_vaf_groups_df,df_tidy%>%
                                     group_by(Sample,VAF_group)%>%
                                     dplyr::summarise(n=n(),.groups = "drop_last"),by=c("Sample","VAF_group"))%>%
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


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 7D ---------------
## Mutations per cell by heteroplasmy, HSPC vs lymphoid, donors with both
#-----------------------------------------------------------------------------------#

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
  my_markdown_theme+
  theme(axis.text.x=element_text(angle=90),strip.text.y=element_text(angle=0),legend.position="none")+
  labs(x="Heteroplasmy level",y="mtDNA mutations per cell")

blood_vs_lymph_VAF_distribution

ggsave(filename = paste0(ed7_dir,"ExtDataFig7d.blood_vs_lymph_VAF_dist.pdf"),blood_vs_lymph_VAF_distribution,width=3.3,height=2.5)


#Per-sample homoplasmic burden, used by panels e-i
n_near_homoplasmic_per_sample<-VAF_distribution_by_sample_df%>%
  filter(VAF_group%in%new_VAF_groups[18:20])%>%
  group_by(dataset,exp_ID,Sample)%>%
  dplyr::summarise(dataset=dataset[1],exp_ID=exp_ID[1],Sample=Sample[1],n_homo=sum(n),.groups = "drop_last")%>%
  left_join(ref_df%>%dplyr::select(ID,Age),relationship = "many-to-many",by=c("exp_ID"="ID"))

Sample_order<-ref_df%>%mutate(Cohort=name_conversion_vec[Cohort])%>%arrange(Cohort,Age)%>%unite(col = "Cohort_ID",Cohort,ID,sep="_")%>%pull(Cohort_ID)
studies_to_include=c("blood","lymph")
blood_vs_lymph_homo_muts_dat<-n_near_homoplasmic_per_sample%>%
  group_by(dataset,exp_ID)%>%
  dplyr::summarise(n=n(),mean=mean(n_homo),var=var(n_homo),sd=sd(n_homo),.groups = "drop_last")%>%
  mutate(sem=sd/sqrt(n))%>%
  left_join(ref_df,by=c("dataset"="Cohort","exp_ID"="ID"))%>%
  filter(dataset%in%studies_to_include & Age>20)

lmer.blood_vs_lymph<-lme4::lmer(mean~Age+dataset + (1|exp_ID),data = blood_vs_lymph_homo_muts_dat)
coefs<-summary(lmer.blood_vs_lymph)$coefficients
coef_df<-data.frame(dataset=c("Normal blood","Lymphoid"),slope=rep(coefs['Age','Estimate'],2),intercept=c(coefs['(Intercept)','Estimate'],coefs['(Intercept)','Estimate']+coefs['datasetlymph','Estimate']))

coef_df


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 7E ---------------
## Mean homoplasmic burden against age by blood compartment
#-----------------------------------------------------------------------------------#

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
  my_markdown_theme

HSC_vs_lymph_mean_homo_burden

ggsave(filename = paste0(ed7_dir,"ExtDataFig7e.HSC_vs_lymph_mean_homo_burden.pdf"),HSC_vs_lymph_mean_homo_burden,width=3.5,height=2.5)


#Poisson mixed model of homoplasmic burden against age and tissue
minimum_age_for_regression=20
n_near_homoplasmic_per_sample_FILT<-n_near_homoplasmic_per_sample%>%filter(Age>=minimum_age_for_regression)%>%mutate(logAge=log(Age+0.1))

n_near_homoplasmic_per_sample_FILT$Age_scaled<-scale(n_near_homoplasmic_per_sample_FILT$Age)
glmer.res<-lme4::glmer(n_homo~Age_scaled+dataset+(1+Age_scaled|exp_ID),family=poisson(link="log"),data=n_near_homoplasmic_per_sample_FILT)
summary(glmer.res)
car::vif(glmer.res)
overdispersion_ratio <- sum(residuals(glmer.res, type = "pearson")^2) / df.residual(glmer.res)
overdispersion_ratio



#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 7F ---------------
## Poisson regression across tissues
#-----------------------------------------------------------------------------------#

#Now generate the regression curves for plotting by using the 'predict' function on test data

# get the scaling parameters used when creatining Age_scaled in the model data
age_center <- attr(n_near_homoplasmic_per_sample_FILT$Age_scaled, "scaled:center")
age_scale <- attr(n_near_homoplasmic_per_sample_FILT$Age_scaled, "scaled:scale")

ages_to_cover <- seq(minimum_age_for_regression, 80, 2)
tissues_to_cover <- c("blood", "lymph", "KY", "LM")

test_data <- data.frame(
  dataset = rep(tissues_to_cover, each = length(ages_to_cover)),
  exp_ID = ids::ids(length(ages_to_cover) * length(tissues_to_cover)),
  Age = rep(ages_to_cover, times = length(tissues_to_cover))
) %>%
  mutate(
    logAge = log(Age),
    Age_scaled = (Age - age_center) / age_scale
  )

test_data$n_homo <- exp(predict(glmer.res, newdata = test_data, allow.new.levels = TRUE))

#Visualize the regression model on top of the data for the selected tissues
#Functions to calculate the upper/ lower CI of the mean from count data
get_upper_CI<-function (X, conf.level=0.95) {alpha = 1 - conf.level; upper <- 0.5 * qchisq(1-alpha/2, 2*X+2)}
get_lower_CI<-function (X, conf.level=0.95) {alpha = 1 - conf.level; lower <- 0.5 * qchisq(alpha/2, 2*X)}

data_vs_regression<-n_near_homoplasmic_per_sample%>%
  group_by(dataset,exp_ID)%>%
  dplyr::summarise(n_samp=n(),mean=mean(n_homo),total_homo=sum(n_homo),.groups = "drop_last")%>%
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
  my_markdown_theme

data_vs_regression

ggsave(filename = paste0(ed7_dir,"ExtDataFig7f.data_vs_regression.pdf"),data_vs_regression,width=3.5,height=2.5)


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 7H ---------------
## Tissue-specific parameters from the Poisson model
#-----------------------------------------------------------------------------------#

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
  my_markdown_theme+
  theme(axis.title.y=element_blank(),legend.position = "none")+
  labs(x="Tissue drift coefficient\n(relative to normal blood)")

homo_muts_glm_coefficients

ggsave(filename = paste0(ed7_dir,"ExtDataFig7h.homo_muts_glm_coefficients.pdf"),homo_muts_glm_coefficients,width=3,height=2.5)


#Refit with the lymphoid compartment split by subset, for panels g and i
minimum_age_for_regression=20
n_homoplasmic_lymph_type<-n_near_homoplasmic_per_sample_FILT%>%
  filter(dataset%in%c("blood","lymph"))%>%
  left_join(colony_info%>%dplyr::select("Sample"=colony,"exp_ID"=Donor,CellType,Cell.type2),by=c("Sample","exp_ID"))%>%
  mutate(Cell.type2=ifelse(dataset=="blood","HSPC",Cell.type2))

Cell.types<-c("HSPC" ,"Naive B","Naive T","Memory B","Memory T","Treg")
cell_type_cols<-c("#1d6d1f",colorRampPalette(c("#fd81c8","#5f70cc","#aea2eb"))(5))
names(cell_type_cols)<-Cell.types

n_homoplasmic_lymph_type$Age_scaled=scale(n_homoplasmic_lymph_type$Age)

#Now do the glmer poisson regression model using the cell type, not just lymphoid vs blood
glmer.res_blood_vs_lymph_subdivided<-lme4::glmer(n_homo~Age_scaled+Cell.type2+(1+Age_scaled|exp_ID),family=poisson(link="log"),data=n_homoplasmic_lymph_type)

summary(glmer.res_blood_vs_lymph_subdivided)


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 7G ---------------
## As f, lymphoid split by subset
#-----------------------------------------------------------------------------------#

# get the scaling parameters used when creatining Age_scaled in the model data
age_center <- attr(n_homoplasmic_lymph_type$Age_scaled, "scaled:center")
age_scale <- attr(n_homoplasmic_lymph_type$Age_scaled, "scaled:scale")

ages_to_cover <- seq(minimum_age_for_regression, 80, 2)
tissues_to_cover<-unique(n_homoplasmic_lymph_type$Cell.type2)

test_data <- data.frame(
  Cell.type2 = rep(tissues_to_cover, each = length(ages_to_cover)),
  exp_ID = ids::ids(length(ages_to_cover) * length(tissues_to_cover)),
  Age = rep(ages_to_cover, times = length(tissues_to_cover))
) %>%
  mutate(
    logAge = log(Age),
    Age_scaled = (Age - age_center) / age_scale
  )

test_data$n_homo <- exp(predict(glmer.res_blood_vs_lymph_subdivided, newdata = test_data, allow.new.levels = TRUE))

#Visualize the regression model on top of the data for the selected tissues
#Functions to calculate the upper/ lower CI of the mean from count data
get_upper_CI<-function (X, conf.level=0.95) {alpha = 1 - conf.level; upper <- 0.5 * qchisq(1-alpha/2, 2*X+2)}
get_lower_CI<-function (X, conf.level=0.95) {alpha = 1 - conf.level; lower <- 0.5 * qchisq(alpha/2, 2*X)}

data_vs_regression_lymph_subdivided<-n_homoplasmic_lymph_type%>%
  group_by(Age,Cell.type2,exp_ID)%>%
  dplyr::summarise(Age=Age[1],n_samp=n(),mean=mean(n_homo),total_homo=sum(n_homo),.groups = "drop_last")%>%
  dplyr::mutate(CI_low=get_lower_CI(total_homo)/n_samp,CI_high=get_upper_CI(total_homo)/n_samp)%>%
  filter(Age>=minimum_age_for_regression)%>%
  ggplot(aes(x=Age,y=mean,ymin=CI_low,ymax=CI_high,col=Cell.type2))+
  geom_point(size=0.5)+
  geom_line(aes(x=Age,y=n_homo,col = Cell.type2),data = test_data,inherit.aes = F)+
  geom_errorbar(width=2,alpha=0.5)+
  scale_color_manual(values = cell_type_cols)+
  labs(x="Age",y="Mean homoplasmic mutations per cell",col="")+
  theme_classic()+
  my_markdown_theme
data_vs_regression_lymph_subdivided

ggsave(filename = paste0(ed7_dir,"ExtDataFig7g.data_vs_regression_lymph_subdivided.pdf"),data_vs_regression_lymph_subdivided,width=3.5,height=2.5)


#-----------------------------------------------------------------------------------#
##--------Generate EXTENDED DATA FIG. 7I ---------------
## Lymphoid subset parameters from the Poisson model
#-----------------------------------------------------------------------------------#

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
  my_markdown_theme+
  theme(axis.title.y=element_blank(),legend.position = "none")+
  labs(x="Tissue drift coefficient\n(relative to normal HSPCs)")

blood_vs_lymph_subdivided_homo_muts_glm_coefficients

ggsave(filename = paste0(ed7_dir,"ExtDataFig7i.blood_vs_lymph_subdivided_homo_muts_glm_coefficients.pdf"),blood_vs_lymph_subdivided_homo_muts_glm_coefficients,width=3,height=2.5)

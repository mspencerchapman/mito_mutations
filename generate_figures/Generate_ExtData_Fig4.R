#-----------------------------------------------------------------------------------#
# Generate_ExtData_Fig4.R
#
# Extended Data Fig. 4 | Mutation burden analysis.
#
#   a  Proportion of samples per individual with at least one mutation above the
#      VAF cut-off on the x-axis.
#   b  Mutation burden (sum of VAF) in HSCs vs HPCs, by individual.
#   c  Burden in wild-type cells vs those with a DNMT3A/TET2/ASXL1 driver.
#   d  Burden in cells from clonal expansions vs singletons.
#   e  Rate of mutation acquisition across datasets, from a linear mixed effects
#      model with dataset-specific interaction terms.
#
# All five panels are built here.
#
# Panel code follows mtDNA_mutations_blood.Rmd (the current blood analysis).
#-----------------------------------------------------------------------------------#

cran_packages=c("dplyr","tidyr","tibble","ggplot2","stringr","RColorBrewer","ape","lme4","lmerTest","readxl")
for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

options(stringsAsFactors = F)

source(here::here("config.R")) #sets root_dir, genomeFile, project paths and my_theme
#plots_dir, rebuttal_figs_dir and my_theme all come from config.R
dir.create(paste0(plots_dir,"Extended_Data_Figure_04"),showWarnings=FALSE,recursive=TRUE)
ed4_dir=paste0(plots_dir,"Extended_Data_Figure_04/")

#-----------------------------------------------------------------------------------#
# Data
#-----------------------------------------------------------------------------------#

source(paste0(root_dir,"/data/mito_mutations_blood_functions.R")) #supplies get_expanded_clade_nodes, getTips
library(treemut)

mito_data<-readRDS(paste0(root_dir,"/data/mito_data.Rds"))
ref_df=read.csv(paste0(root_dir,"/data/Samples_metadata_ref.csv"))%>%filter(Dataset!="Lymphocyte")

#Copy number table, augmented with sample-level metadata. Cell_type and Phenotype
#are not in the coverage file itself - they come from the sample metadata and the
#phenotype summary, joined here exactly as in mtDNA_mutations_blood.Rmd.
sample_level_metadata_EM<-read.csv(paste0(root_dir,"/data/EM_sample_level_metadata.csv"))%>%
  dplyr::rename("exp_ID"="donor_id","Sample"="PDID","Cell_type"="cell_type","coverage"="mean_depth")
sample_level_metadata_foetal<-read.csv(paste0(root_dir,"/data/foetal_sample_level_metadata.csv"))%>%
  mutate(Sample=paste(Donor_ID,"hum",sep="_"))
sample_level_metadata<-dplyr::bind_rows(sample_level_metadata_EM,sample_level_metadata_foetal)
phenotype_data_EM<-read.csv(paste0(root_dir,"/data/Summary_pheno_pdid.csv"),stringsAsFactors=F)%>%
  dplyr::rename("exp_ID"=Donor_ID,"Sample"=PDID)

mito_cn=read.csv(paste0(root_dir,"/data/whole_genome_coverage_pileup_and_bedtools_annotated.csv"),header=T)%>%
  left_join(sample_level_metadata,by="Sample")%>%
  left_join(ref_df%>%dplyr::rename("exp_ID"=Sample),by=intersect(names(.),names(ref_df%>%dplyr::rename("exp_ID"=Sample))))
mito_cn<-left_join(mito_cn,phenotype_data_EM,by=intersect(names(mito_cn),names(phenotype_data_EM)))

#Mutations that recurrently slip through the filters as artefacts
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C",
               "MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C")

#Mutations whose VAF tracks mtDNA copy number - most likely mis-mapped nuclear reads
CN_correlating_muts<-readRDS(paste0(root_dir,"/data/CN_correlation.RDS"))

#If implied mutation copy number exceeds this, keep the mutation even if it is in
#the CN-correlating list: a genuine mutation at such a site gives a much higher
#implied copy number than a mis-mapping artefact does.
mutCN_cutoff=25

#Collapse the foetal cell type codes to HSC / HPC
convert_vec=c("HSC","HSC","HPC","HPC","HPC","HPC")
names(convert_vec)=c("H","HSC","C","M","HSPC","Progenitor")

#' Per-donor matrix of mutations that pass every filter
#'
#' Retains calls that (1) pass shearwater, (2) are assigned to the genuine
#' mutational signature N1, and (3) are not indels, blacklisted artefacts or
#' unrescued copy-number-correlating mutations.
filtered_vaf_matrix<-function(list) {
  CN_removal<-list$matrices$implied_mutCN>mutCN_cutoff |
    (matrix(!rownames(list$matrices$vaf)%in%CN_correlating_muts,ncol=1)%*%
       matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
  list$matrices$vaf*list$matrices$SW*(list$matrices$ML_Sig=="N1")*CN_removal
}

#-----------------------------------------------------------------------------------#
# Fig ED4a | Proportion of samples with a mutation above a VAF cut-off
#-----------------------------------------------------------------------------------#

df_tidy<-dplyr::bind_rows(Map(list=mito_data,exp_ID=names(mito_data),function(list,exp_ID) {
  if(is.null(list)) return(NULL)
  as.data.frame(filtered_vaf_matrix(list))%>%
    tibble::rownames_to_column(var="mut_ref")%>%
    dplyr::select(-any_of("global"))%>%
    tidyr::gather(key="Sample",value="vaf",-mut_ref)%>%
    dplyr::filter(!grepl("DEL|INS",mut_ref) & !mut_ref%in%exclude_muts)%>%
    dplyr::filter(vaf>0)%>%
    mutate(exp_ID=exp_ID)
}))
cat("df_tidy:",nrow(df_tidy),"mutation-sample observations\n")

#Sweep the VAF threshold and record, per individual, the proportion of samples
#carrying at least one mutation above it
samples_with_mut_df<-dplyr::bind_rows(lapply(seq(0.01,0.99,0.01),function(cut_off) {
  df_tidy%>%
    dplyr::filter(Sample!="global")%>%
    group_by(Sample)%>%
    dplyr::summarise(n_mut=sum(vaf>cut_off),.groups="drop")%>%
    tidyr::complete(Sample,fill=list(n_mut=0))%>%
    mutate(exp_ID=sapply(Sample,function(SampleID) {mito_cn$exp_ID[mito_cn$Sample==SampleID][1]}))%>%
    group_by(exp_ID)%>%
    dplyr::summarise(mean=mean(n_mut),n_with_mut=sum(n_mut>0),n_samp=n(),.groups="drop")%>%
    mutate(prop_with_mut=n_with_mut/n_samp,cut_off=cut_off)
}))

prop.of.samples.with.mut<-samples_with_mut_df%>%
  dplyr::mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
  dplyr::filter(!is.na(exp_ID))%>%
  ggplot(aes(x=cut_off,y=prop_with_mut,col=factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])))+
  geom_line(alpha=0.6)+
  theme_bw()+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  scale_color_brewer(palette="Paired")+
  labs(x="Cut-off",
       y=str_wrap("Proportion of samples with at least 1 mutation with VAF > cut-off",width=40),
       col="Individual")+
  my_theme+
  theme(legend.key.height=unit(3,"mm"),legend.title=element_text(size=7))

ggsave(filename=paste0(ed4_dir,"ExtDataFig4a.Sample_proportions_with_mut.pdf"),prop.of.samples.with.mut,width=3,height=2)
cat("ED4a: written\n")

#-----------------------------------------------------------------------------------#
# Fig ED4b | Mutation burden in HSCs vs HPCs
#
# Burden is the sum of VAF across all retained mutations in a sample. Restricted
# to the individuals for whom both cell types were sampled.
#-----------------------------------------------------------------------------------#

sum_of_vaf_df<-Map(list=mito_data,Exp_ID=names(mito_data),function(list,Exp_ID){
  mat<-filtered_vaf_matrix(list)
  keep<-!grepl("DEL|INS",rownames(mat)) & !rownames(mat)%in%exclude_muts
  as.data.frame(colSums(mat[keep,,drop=FALSE],na.rm=TRUE))%>%
    tibble::rownames_to_column(var="Sample")%>%
    mutate(exp_ID=Exp_ID)%>%
    dplyr::rename(sum_of_vaf=2)%>%
    dplyr::filter(Sample!="global")
})%>%
  dplyr::bind_rows()%>%
  left_join(mito_cn%>%dplyr::select(Sample,Cell_type,Phenotype),by="Sample")%>%
  dplyr::mutate(Cell_type=convert_vec[Cell_type])

individuals_with_both=c("18pcw","CB002","SX001","AX001","KX004")

#Wilcoxon per individual. Computed directly rather than with
#ggpubr::stat_compare_means(), which fails against the installed ggplot2
#(its internal create_p_label() is not found).
wilcox_labels<-sum_of_vaf_df%>%
  filter(exp_ID%in%individuals_with_both,!is.na(Cell_type),!is.na(sum_of_vaf))%>%
  group_by(exp_ID)%>%
  dplyr::summarise(label=tryCatch(paste0("p = ",signif(wilcox.test(sum_of_vaf~Cell_type)$p.value,2)),
                                  error=function(e) ""),
                   y=max(sum_of_vaf,na.rm=TRUE),.groups="drop")

sum_of_vaf_by_celltype<-sum_of_vaf_df%>%
  filter(exp_ID%in%individuals_with_both,!is.na(Cell_type))%>%
  ggplot(aes(x=Cell_type,y=sum_of_vaf))+
  geom_violin(aes(fill=Cell_type),alpha=0.2)+
  scale_color_manual(values=c("#1a80bb","#ea801c"))+
  scale_fill_manual(values=c("#1a80bb","#ea801c"))+
  geom_jitter(aes(col=Cell_type),width=0.2,alpha=0.5)+
  geom_text(data=wilcox_labels,aes(x=1.5,y=y,label=label),size=2.6,fontface="italic",inherit.aes=FALSE)+
  facet_grid(cols=vars(factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])),scales="free",space="free")+
  theme_classic()+
  my_theme+
  labs(y="Mutation burden\n(sum of VAF)")+
  theme(axis.title.x=element_blank(),legend.position="none")

ggsave(filename=paste0(ed4_dir,"ExtDataFig4b.sum_of_vaf_by_celltype.pdf"),sum_of_vaf_by_celltype,width=7,height=2.2)
cat("ED4b: written\n")

#-----------------------------------------------------------------------------------#
# Fig ED4c and ED4d | Burden by driver status, and by clonal expansion
#
# Nuclear driver calls live in data/blood_adult/annot_files_filtered/ rather than
# inside mito_data.Rds. Each file holds $mat (mutations annotated with Gene, CDS,
# Protein) plus the NV/NR count matrices. treemut::assign_to_tree() re-assigns
# each nuclear mutation to a tree branch, so a driver can be traced to every
# sample descended from that branch.
#
# Restricted to the four elderly individuals - the only ones with enough cells
# carrying drivers, or belonging to clonal expansions, to compare.
#-----------------------------------------------------------------------------------#

IDs_to_include=paste0("KX00",c(7,8,4,3))
mut_set=c("DNMT3A","TET2","ASXL1")

#Attach the nuclear calls and re-assign them to tree branches
mito_data<-Map(list=mito_data,exp_ID=names(mito_data),function(list,exp_ID) {
  annot_muts_file<-paste0(root_dir,"/data/blood_adult/annot_files_filtered/annotated_muts_filt_",exp_ID,".Rds")
  if(!file.exists(annot_muts_file)) return(list)
  list$nDNA_mats<-readRDS(annot_muts_file)
  if(nrow(list$nDNA_mats$mat)==0) return(list)

  tree_samples<-list$tree.ultra$tip.label[list$tree.ultra$tip.label!="Ancestral"]
  if(any(!tree_samples%in%colnames(list$nDNA_mats$NV))) {
    drop_tips<-tree_samples[!tree_samples%in%colnames(list$nDNA_mats$NV)]
    list$tree.ultra<-drop.tip(list$tree.ultra,tip=drop_tips)
    list$tree<-drop.tip(list$tree,tip=drop_tips)
    tree_samples<-tree_samples[tree_samples%in%colnames(list$nDNA_mats$NV)]
  }
  #The appended zero-count/depth-10 column represents the ancestral outgroup
  mtr<-as.matrix(cbind(list$nDNA_mats$NV[,tree_samples,drop=F],
                       matrix(0,ncol=1,nrow=nrow(list$nDNA_mats$NV),dimnames=list(NULL,"Ancestral"))))
  dep<-as.matrix(cbind(list$nDNA_mats$NR[,tree_samples,drop=F],
                       matrix(10,ncol=1,nrow=nrow(list$nDNA_mats$NR),dimnames=list(NULL,"Ancestral"))))
  res<-treemut::assign_to_tree(tree=list$tree.ultra,mtr=mtr,dep=dep)
  list$nDNA_mats$mat$node<-res$tree$edge[res$summary$edge_ml,2]
  return(list)
})

#Cut the tree at 100 mutations of molecular time; each resulting clade (some are
#singletons) gets its mean mutation burden and any driver on its ancestral branch
have_nDNA<-names(mito_data)[sapply(mito_data,function(x) !is.null(x$nDNA_mats))]
sov_comparison_df<-dplyr::bind_rows(Map(list=mito_data[have_nDNA],exp_ID=have_nDNA,function(list,exp_ID) {
  ec_df<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off=100,min_clonal_fraction=0,min_samples=1)
  #vapply, not sapply: a donor with no clades would otherwise yield an empty list
  #and make mean_sov a list column, which then cannot be row-bound with the rest
  ec_df$mean_sov<-if(nrow(ec_df)) vapply(ec_df$nodes,function(node) {
    Samples<-getTips(list$tree.ultra,node=node)
    mean(sum_of_vaf_df%>%filter(Sample%in%Samples)%>%pull(sum_of_vaf),na.rm=TRUE)
  },numeric(1)) else numeric(0)
  if(!nrow(ec_df)) return(NULL)
  left_join(ec_df,list$nDNA_mats$mat%>%dplyr::select(node,mut_ref,Gene,CDS,Protein),by=c("nodes"="node"))%>%
    mutate(exp_ID=exp_ID,.before=1)
}))

#' Violin comparison with a Wilcoxon p-value per individual
#'
#' ggpubr::stat_compare_means() is not used: it fails against the installed
#' ggplot2 (its internal create_p_label() is not found).
burden_comparison_plot<-function(df,xlab_angle=0) {
  labs_df<-df%>%filter(!is.na(mean_sov))%>%group_by(exp_ID)%>%
    dplyr::summarise(label=tryCatch(paste0("p = ",signif(wilcox.test(mean_sov~type)$p.value,2)),
                                    error=function(e) ""),
                     y=max(mean_sov,na.rm=TRUE),.groups="drop")
  ggplot(df,aes(x=type,y=mean_sov,col=type,fill=type))+
    geom_violin(alpha=0.1)+
    geom_jitter(width=0.1,height=0,alpha=0.3)+
    scale_color_brewer(palette="Set1")+
    scale_fill_brewer(palette="Set1")+
    geom_text(data=labs_df,aes(x=1.5,y=y,label=label),size=3,fontface="italic",inherit.aes=FALSE)+
    facet_grid(~exp_ID)+
    theme_classic()+
    my_theme+
    theme(legend.position="none",axis.text.x=element_text(angle=xlab_angle))+
    labs(x="",y="Mean mtDNA mutation burden\n (sum of VAF)")
}

if(nrow(sov_comparison_df)) {
  #ED4c - clades whose ancestral branch carries a DTA driver vs those that do not
  sov_mut_vs_wt_df<-sov_comparison_df%>%
    mutate(exp_ID=factor(exp_ID,levels=ref_df%>%filter(Age>0)%>%arrange(Age)%>%pull(Sample)))%>%
    tidyr::replace_na(list(Gene=""))%>%
    mutate(type=ifelse(Gene%in%mut_set,"mut","wt"))%>%
    filter(exp_ID%in%IDs_to_include)
  ggsave(paste0(ed4_dir,"ExtDataFig4c.sov_mut_vs_wt.pdf"),burden_comparison_plot(sov_mut_vs_wt_df),width=3.3,height=2.5)
  cat("ED4c: written\n")

  #ED4d - clades with more than one sample (an expansion) vs singletons
  exp_vs_singleton_df<-sov_comparison_df%>%
    mutate(exp_ID=factor(exp_ID,levels=ref_df%>%filter(Age>0)%>%arrange(Age)%>%pull(Sample)))%>%
    filter(exp_ID%in%IDs_to_include)%>%
    mutate(type=ifelse(n_samples>1,"Clonal\nexpansion",ifelse(n_samples==1,"Singleton",NA)))%>%
    filter(!is.na(type))
  ggsave(paste0(ed4_dir,"ExtDataFig4d.sov_exp_vs_singleton.pdf"),burden_comparison_plot(exp_vs_singleton_df),width=3.3,height=2.5)
  cat("ED4d: written\n")
} else {cat("ED4c/d: no nuclear annotation files found - skipped\n")}

#-----------------------------------------------------------------------------------#
# Fig ED4e | Rate of mutation acquisition across datasets
#
# Mutation burden (sum of VAF) per individual against age, faceted by tissue,
# with dataset-specific slopes from a mixed model carrying Age:dataset
# interaction terms - i.e. each tissue is allowed its own accumulation rate, with
# individual as a random effect. Ribbons are 95% intervals on the predicted mean,
# obtained by parametric bootstrap (lme4::bootMer).
#
# Panel code follows mtDNA_mutations_comparator_tissues.Rmd, chunks
# mutation_rate_by_tissue_lmer and mutation_rate_by_tissue_and_age_lmer_plots2
# (the faceted mutburden_by_age_and_tissue_plot).
#-----------------------------------------------------------------------------------#

all_cohorts_ed4e<-c("KY","HL","SO","PR","LM","NW","lymph","blood")
name_conversion_vec=c("Normal blood","Lymphoid","MPN","Colon","IBD-affected\nColon",
                      "MUTYH-mutant\nColon","Endometrium","Bronchial\nepithelium")
names(name_conversion_vec)=c("blood","lymph","NW","HL","SO","PR","LM","KY")
tissue_cols<-c("#96e97c","#fd81c8","#145a6a","#65e6f9","#781486","#5f70cc","#aea2eb","#1d6d1f")
names(tissue_cols)<-name_conversion_vec

vaf_cut_off=0.03 #only count mutations above this VAF towards the burden
n_boot=1000      #bootstrap replicates for the confidence ribbons

ref_df_nonblood<-readxl::read_excel(paste0(root_dir,"/data/metadata/non_blood_metadata.xlsx"))

cohort_data<-lapply(all_cohorts_ed4e,function(d) {
  f<-paste0(root_dir,"/data/nonblood/mito_mutation_data_",d,".RDS")
  if(file.exists(f)) readRDS(f) else NULL
})
names(cohort_data)<-all_cohorts_ed4e
cohort_data<-cohort_data[!sapply(cohort_data,is.null)]

#Per-sample burden. Heteroplasmic oocyte mutations are excluded as well as the
#usual artefacts: a cell inheriting one starts with a head start in burden that
#has nothing to do with its age.
all_sum_of_vaf_df<-dplyr::bind_rows(Map(dataset_mito_data=cohort_data,dataset=names(cohort_data),
                                        function(dataset_mito_data,dataset) {
  CN_correlating_muts_ds<-dataset_mito_data$CN_correlating_muts
  dplyr::bind_rows(Map(list=dataset_mito_data,Exp_ID=names(dataset_mito_data),function(list,Exp_ID){
    if(is.null(list$matrices)) return(NULL)
    cn_muts<-if(!is.null(list$CN_correlating_muts)) list$CN_correlating_muts else CN_correlating_muts_ds
    CN_removal<-list$matrices$implied_mutCN>mutCN_cutoff |
      (matrix(!rownames(list$matrices$vaf)%in%cn_muts,ncol=1)%*%
         matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
    keep<-!grepl("DEL|INS",rownames(list$matrices$vaf)) &
      !rownames(list$matrices$vaf)%in%c(exclude_muts,list$het_oocyte_muts)
    mat<-(list$matrices$vaf*(list$matrices$vaf>vaf_cut_off)*list$matrices$SW*CN_removal)[keep,,drop=FALSE]
    as.data.frame(colSums(mat,na.rm=TRUE))%>%
      tibble::rownames_to_column(var="Sample")%>%
      mutate(exp_ID=Exp_ID)%>%
      dplyr::rename(sum_of_vaf=2)
  }))%>%mutate(dataset=dataset,.before=1)
}))

#Within a clonal expansion the samples are not independent observations of the
#mutation rate, so all but one member of each expanded clade is dropped.
set.seed(42)
all_samples_to_drop<-unlist(lapply(cohort_data,function(dataset_mito_data) {
  lapply(names(dataset_mito_data),function(this_exp_ID) {
    tr<-dataset_mito_data[[this_exp_ID]]$tree.ultra
    if(is.null(tr)) return(NULL)
    mut_burden<-get_mut_burden(tr)[1]
    ecn<-get_expanded_clade_nodes(tree=tr,min_samples=2,height_cut_off=mut_burden/2,min_clonal_fraction=0)
    unlist(lapply(ecn$nodes,function(node) {
      clade_samples<-getTips(tree=tr,node)
      sample(clade_samples,size=length(clade_samples)-1)
    }))
  })
}))

min_samples_per_individual=8
all_sum_of_vaf_summary_df<-dplyr::bind_rows(lapply(names(cohort_data),function(this_dataset) {
  all_sum_of_vaf_df%>%
    filter(dataset==this_dataset & !Sample%in%all_samples_to_drop)%>%
    group_by(exp_ID)%>%
    dplyr::summarise(n=n(),mean=mean(sum_of_vaf,na.rm=T),median=median(sum_of_vaf,na.rm=T),.groups="drop")%>%
    mutate(dataset=this_dataset,.before=1)%>%
    filter(n>=min_samples_per_individual)%>%
    left_join(ref_df_nonblood%>%filter(Cohort==this_dataset),by=c("exp_ID"="ID"))
}))%>%filter(!is.na(Age))

cat("ED4e: ",nrow(all_sum_of_vaf_summary_df)," individuals across ",
    length(unique(all_sum_of_vaf_summary_df$dataset))," datasets\n",sep="")

if(nrow(all_sum_of_vaf_summary_df)>10) {
  #Age:dataset interaction only - this parameterisation gives a slope for every
  #dataset including blood, rather than contrasts against a reference tissue
  lmer2<-lmerTest::lmer(mean~Age:dataset+(1|exp_ID),data=all_sum_of_vaf_summary_df)
  print(summary(lmer2))

  pred_df<-expand.grid(Age=seq(min(all_sum_of_vaf_summary_df$Age),max(all_sum_of_vaf_summary_df$Age),length.out=100),
                       dataset=unique(all_sum_of_vaf_summary_df$dataset))
  pred_df$fit<-predict(lmer2,newdata=pred_df,re.form=NA) #re.form=NA -> population-level prediction

  boot_ci<-lme4::bootMer(lmer2,FUN=function(x) predict(x,newdata=pred_df,re.form=NA),
                         nsim=n_boot,use.u=FALSE,type="parametric")
  pred_df$lower<-apply(boot_ci$t,2,quantile,0.025)
  pred_df$upper<-apply(boot_ci$t,2,quantile,0.975)

  #One facet per tissue rather than all slopes overlaid - the same model, easier
  #to read when the slopes are similar
  lab<-function(d) factor(name_conversion_vec[d],levels=name_conversion_vec)
  mutburden_by_age_and_tissue_plot<-ggplot()+
    geom_ribbon(data=pred_df%>%mutate(dataset=lab(dataset)),
                aes(x=Age,ymin=lower,ymax=upper,fill=dataset),alpha=0.1)+
    geom_line(data=pred_df%>%mutate(dataset=lab(dataset)),
              aes(x=Age,y=fit,colour=dataset))+
    geom_point(data=all_sum_of_vaf_summary_df%>%mutate(dataset=lab(dataset)),
               aes(x=Age,y=mean,colour=dataset),alpha=0.6)+
    scale_colour_manual(values=tissue_cols)+
    scale_fill_manual(values=tissue_cols)+
    facet_grid(~dataset)+
    theme_classic()+
    my_theme+
    labs(y="Mean mutation burden (sum of VAF)")+
    theme(legend.position="none")

  ggsave(filename=paste0(ed4_dir,"ExtDataFig4e.mutburden_by_age_and_tissue_plot.pdf"),
         mutburden_by_age_and_tissue_plot,width=7,height=2)
  cat("ED4e: written\n")
} else {cat("ED4e: too few individuals with age data - skipped\n")}

cat("\nExtended Data Fig. 4 panels written to",ed4_dir,"\n")


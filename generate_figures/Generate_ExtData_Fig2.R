#-----------------------------------------------------------------------------------#
# Generate_ExtData_Fig2.R
#
# Extended Data Fig. 2 | Mitochondrial DNA copy number.
#
#   a  Density ridge plot of mtDNA copy number across samples from each tissue
#      (log scale), coloured by tissue.
#   b  Tissue-specific fixed-effect terms from a linear mixed effects regression
#      of copy number on tissue, with individual as a random effect. Error bars
#      are 95% confidence intervals. The regression is fitted on log copy number;
#      coefficients are back-transformed to absolute copy number for display.
#   c  Median mtDNA copy number per individual against age, coloured by tissue.
#
# Panel code follows mtDNA_mutations_comparator_tissues.Rmd (the current
# cross-tissue analysis), chunks copy_number_analysis_ridge_plot,
# copy_number_by_age, mean_cn_regression_create_df / mean_cn_regression and
# lmer2_visualization.
#
# Note: that notebook derives the sample list from a fully preprocessed
# all_mito_datasets object. Only the per-donor sample identifiers are needed
# here, so this script reads the tip labels straight out of each cohort's data
# object instead - far quicker, and identical in result.
#-----------------------------------------------------------------------------------#

cran_packages=c("dplyr","ggplot2","ggridges","readxl","tibble","lme4","lmerTest")
for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

options(stringsAsFactors = F)

#Set these file paths before running the script
source(here::here("config.R")) #sets root_dir, genomeFile, project paths and my_theme
plots_dir=paste0(root_dir,"/plots/")
dir.create(paste0(plots_dir,"Extended_Data_Figure_02"),showWarnings=FALSE,recursive=TRUE)
ed2_dir=paste0(plots_dir,"Extended_Data_Figure_02/")

nonblood_ref_file=paste0(root_dir,"/data/metadata/non_blood_metadata.xlsx")
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R")) #supplies plot_tree

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))+
  theme(legend.key.height=unit(3,"mm"),legend.title = element_text(size=8))

#-----------------------------------------------------------------------------------#
# Cohorts, tissue labels and colours (as in the cross-tissue notebook)
#-----------------------------------------------------------------------------------#

all_cohorts<-c("KY","HL","SO","PR","LM","NW","lymph","blood")

#Maps the dataset key (generally the initials of the scientist who generated the
#data) to the tissue label used in the figures
name_conversion_vec=c("Normal blood","Lymphoid","MPN","Colon","IBD-affected\nColon",
                      "MUTYH-mutant\nColon","Endometrium","Bronchial\nepithelium")
names(name_conversion_vec)=c("blood","lymph","NW","HL","SO","PR","LM","KY")

tissue_cols<-c("#96e97c","#fd81c8","#145a6a","#65e6f9","#781486","#5f70cc","#aea2eb","#1d6d1f")
names(tissue_cols)<-name_conversion_vec

#Dataset key -> the tissue name used in the copy number table
translate=data.frame(dataset=c("KY","HL","SO","PR","LM","NW","lymph","EM"),
                     al_ref=c("lung_organoid","colon","colon_ibd","muty_mutant",
                              "endometrium","blood_MPN","immune","blood"))

#-----------------------------------------------------------------------------------#
# Data
#-----------------------------------------------------------------------------------#

#Per-sample mtDNA copy number from whole-genome coverage
mito_cn=read.csv(paste0(root_dir,"/data/whole_genome_coverage_pileup_and_bedtools_annotated.csv"),header=T)%>%
  dplyr::rename("pileup_median_mtDNA_coverage"=pilup_median_mtDNA_coverage) #This column has a typo

#Individual-level metadata (ages)
ref_df<-readxl::read_excel(nonblood_ref_file)

#Which samples belong to which donor and cohort. Taken from the phylogeny tip
#labels, which define the set of colonies included in the analysis.
cohort_file<-function(cohort) {
  if(cohort=="blood") paste0(root_dir,"/data/nonblood/mito_mutation_data_blood.RDS")
  else paste0(root_dir,"/data/nonblood/mito_mutation_data_",cohort,".RDS")
}
sample_summary_info<-lapply(all_cohorts,function(dataset) {
  f<-cohort_file(dataset)
  if(!file.exists(f)) {cat("  no data for cohort",dataset,"- skipped\n"); return(NULL)}
  d<-readRDS(f)
  dplyr::bind_rows(Map(ind=d,exp_ID=names(d),function(ind,exp_ID) {
    data.frame(dataset=dataset,exp_ID=exp_ID,SampleID=ind$tree$tip.label)
  }))
})%>%dplyr::bind_rows()
cat("Samples:",nrow(sample_summary_info),"across",length(unique(sample_summary_info$dataset)),"cohorts\n")

#Copy number joined to cohort/tissue, restricted to samples in the analysis
mito_cn_by_tissue<-mito_cn%>%
  filter(Tissue%in%translate$al_ref)%>%
  right_join(sample_summary_info,by=c("Sample"="SampleID"))%>%
  mutate(Tissue=factor(name_conversion_vec[dataset],levels=name_conversion_vec))%>%
  filter(!is.na(bedtools_mtDNA_genomes))%>%
  mutate(log_mitoCN=log(bedtools_mtDNA_genomes))

#-----------------------------------------------------------------------------------#
# Fig ED2a | Copy number distribution by tissue
#-----------------------------------------------------------------------------------#

mitoCN_ridges_plot<-mito_cn_by_tissue%>%
  ggplot(aes(y=Tissue,x=bedtools_mtDNA_genomes,fill=Tissue))+
  ggridges::geom_density_ridges(linewidth=0.3)+
  scale_fill_manual(values=tissue_cols)+
  scale_x_log10()+ #copy number is log-normally distributed
  theme_classic()+
  my_theme+
  theme(legend.position="none",axis.title.y=element_blank())+
  labs(x="mtDNA copy number")

ggsave(filename=paste0(ed2_dir,"ExtDataFig2a.mitoCN_ridges_plot.pdf"),mitoCN_ridges_plot,width=3,height=2)

#-----------------------------------------------------------------------------------#
# Fig ED2b | Tissue-specific fixed effects from the mixed model
#
# Copy number is log-normal, so the regression is fitted on log copy number and
# the coefficients are exponentiated for display. Individual is a random effect
# because samples are clustered within donors. Age was tested as a covariate in
# the notebook and dropped: it increased AIC, i.e. no convincing age effect.
#-----------------------------------------------------------------------------------#

mito_cn_lmer_df<-mito_cn_by_tissue%>%
  left_join(ref_df%>%dplyr::select(ID,Age,Smoking_years)%>%filter(!duplicated(.)),by=c("exp_ID"="ID"))

mito_cn.lmer<-lme4::lmer(log_mitoCN~Tissue+(1|exp_ID),data=mito_cn_lmer_df)
print(summary(mito_cn.lmer))
lmer.CIs<-confint(mito_cn.lmer)

#Coefficients are relative to the reference tissue, so add the intercept back to
#each to express them on the same absolute scale
intercept_value<-summary(mito_cn.lmer)$coefficients["(Intercept)","Estimate"]
tissue_lmer_exp_coefs_plot<-as.data.frame(summary(mito_cn.lmer)$coefficients)%>%
  tibble::rownames_to_column(var="Tissue")%>%
  left_join(as.data.frame(lmer.CIs)%>%tibble::rownames_to_column(var="Tissue"),by="Tissue")%>%
  dplyr::select(Tissue,Estimate,lowerCI=`2.5 %`,upperCI=`97.5 %`)%>%
  mutate(Estimate=ifelse(grepl("Tissue",Tissue),Estimate+intercept_value,Estimate),
         lowerCI=ifelse(grepl("Tissue",Tissue),lowerCI+intercept_value,lowerCI),
         upperCI=ifelse(grepl("Tissue",Tissue),upperCI+intercept_value,upperCI))%>%
  filter(Tissue!="Age")%>%
  mutate(Tissue=ifelse(Tissue=="(Intercept)","Normal blood",gsub("Tissue","",Tissue)))%>%
  mutate(Tissue=factor(Tissue,levels=name_conversion_vec))%>%
  filter(!is.na(Tissue))%>%
  mutate_at(c("Estimate","lowerCI","upperCI"),function(x) {exp(x)})%>% #back to absolute copy number
  ggplot(aes(y=Tissue,x=Estimate,xmin=lowerCI,xmax=upperCI,col=Tissue))+
  geom_point(size=0.4)+
  geom_errorbar(width=0.2)+
  scale_x_continuous(limits=c(0,2700))+
  scale_color_manual(values=tissue_cols)+
  theme_classic()+
  my_theme+
  labs(x="Tissue-specific mtDNA copy number\nfixed effect coefficients")+
  theme(legend.position="none")

ggsave(filename=paste0(ed2_dir,"ExtDataFig2b.tissue_lmer_exp_coefs_plot.pdf"),tissue_lmer_exp_coefs_plot,width=3,height=2)

#-----------------------------------------------------------------------------------#
# Fig ED2c | Median copy number per individual against age
#-----------------------------------------------------------------------------------#

median_cn_by_age_plot<-mito_cn_by_tissue%>%
  group_by(Tissue,exp_ID)%>%
  dplyr::summarise(Tissue=Tissue[1],median_cn=median(bedtools_mtDNA_genomes),.groups="drop")%>%
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

ggsave(filename=paste0(ed2_dir,"ExtDataFig2c.median_cn_by_age_plot.pdf"),median_cn_by_age_plot,width=3.3,height=2)

#-----------------------------------------------------------------------------------#
# Fig ED2d | mtDNA copy number on the 8 pcw foetal phylogeny
#
# Copy number is shown as a bar plot beneath the tree, on a log scale because it
# spans an order of magnitude. Cmean (Abouheif's phylogenetic autocorrelation)
# tests whether copy number tracks the phylogeny; the pre-computed result is read
# from data/mito_cn_phylo.Rds, which the blood notebook caches because
# phyloSignal is slow.
#-----------------------------------------------------------------------------------#

mito_cn_phylo_file<-paste0(root_dir,"/data/mito_cn_phylo.Rds")
if(file.exists(mito_cn_phylo_file)) {
  mito_data<-readRDS(paste0(root_dir,"/data/mito_data.Rds"))
  mito_cn_phylo<-readRDS(mito_cn_phylo_file)
  #The published panel shows the 8 post-conception week fetus. NB the two objects
  #key this donor differently - "8 pcw" in mito_cn_phylo, "8pcw" in mito_data -
  #so match on the de-spaced name rather than a literal.
  key<-function(nms) nms[gsub(" ","",nms)=="8pcw"][1]
  phylo_key<-key(names(mito_cn_phylo)); data_key<-key(names(mito_data))
  ed2d_donor<-data_key
  if(!is.na(phylo_key) && !is.na(data_key)) {
    phy<-mito_cn_phylo[[phylo_key]]
    cmean<-phy$phyloSignal_res$stat["dt","Cmean"]; pval<-phy$phyloSignal_res$pvalue["dt","Cmean"]
    pdf(paste0(ed2_dir,"ExtDataFig2d.",ed2d_donor,"_mitoCN_on_phylogeny.pdf"),width=7,height=3)
    tree<-plot_tree(mito_data[[ed2d_donor]]$tree.ultra,cex.label=0,bars=log10(phy$mito_cn),title="8 pcw")
    mtext(sprintf("Log10 MitoCN, Cmean = %.3f (p-value = %.3f)",cmean,pval),side=1,line=2,cex=0.6,adj=0)
    dev.off()
    cat("ED2d: written (Cmean =",round(cmean,3),", p =",round(pval,3),")\n")
  } else {cat("ED2d: donor",ed2d_donor,"not found - skipped\n")}
} else {cat("ED2d: data/mito_cn_phylo.Rds not found - skipped\n")}

#-----------------------------------------------------------------------------------#
# Fig ED2e | Copy number in bulk-sorted HSC and HPC populations
#
# In vitro clonal expansion perturbs mtDNA copy number, so bulk-sorted HSC
# (CD34+CD38-) and HPC (CD34+CD38+) populations were sequenced at low coverage to
# measure copy number in unmanipulated cells. Panel code from
# full_analysis_scripts/mitoCN_bulkbloodpops.R.
#
# NB that script reads "data/mitoCN_metadata.csv"; the file is actually
# MitoCN_metadata.csv, which works only on a case-insensitive filesystem.
#-----------------------------------------------------------------------------------#

mitoCN_metadata_file<-paste0(root_dir,"/data/MitoCN_metadata.csv")
if(file.exists(mitoCN_metadata_file)) {
  mitoCN_metadata<-readr::read_csv(mitoCN_metadata_file,col_select=1:5,show_col_types=FALSE)
  dataset_levels=c("Foetal","Cord blood","Adult")
  cell_type_replace=c("HSC pool","CD34+CD38+\nHPC pool")
  names(cell_type_replace)=c("HSC pool","CD34plusCD38plus haematopoietic progenitors")

  data_clean<-left_join(mitoCN_metadata,mito_cn,by=c("PD_number"="Sample"))%>%
    dplyr::filter(!is.na(Tissue))%>%
    dplyr::select(PDID=PD_number,ID,Dataset,`Cell type`,"mtDNA copy number"=bedtools_mtDNA_genomes)%>%
    dplyr::mutate(`Cell type`=cell_type_replace[`Cell type`],Dataset=factor(Dataset,levels=dataset_levels))

  #HPC vs HSC comparison within each age group. Computed directly rather than with
  #ggpubr::stat_compare_means(), which is broken against the installed ggplot2
  #(its internal create_p_label() is not found).
  wilcox_labels<-data_clean%>%
    dplyr::filter(!is.na(`mtDNA copy number`))%>%
    dplyr::group_by(Dataset)%>%
    dplyr::summarise(label=tryCatch(paste0("p = ",signif(wilcox.test(`mtDNA copy number`~`Cell type`)$p.value,2)),
                                    error=function(e) ""),.groups="drop")

  mitoCN_age_facet<-data_clean%>%
    ggplot(aes(x=`Cell type`,y=`mtDNA copy number`))+
    geom_boxplot(linewidth=0.25,outlier.shape=NA)+
    geom_jitter(aes(col=ID),width=0.1,size=0.75,alpha=0.75)+
    scale_y_continuous(limits=c(0,2600))+
    facet_grid(~Dataset)+
    theme_classic()+
    my_theme+
    geom_text(data=wilcox_labels,aes(x=1.5,y=2500,label=label),size=2,inherit.aes=FALSE)+
    theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())+
    labs(col="Sample ID")

  ggsave(filename=paste0(ed2_dir,"ExtDataFig2e.mitoCN_HSC_vs_HPC.pdf"),mitoCN_age_facet,width=4,height=2.5)
  cat("ED2e: written\n")
} else {cat("ED2e: data/MitoCN_metadata.csv not found - skipped\n")}

#-----------------------------------------------------------------------------------#
# Fig ED2f | Mutation VAF against 1/mtDNA copy number
#
# A mutation whose VAF rises as copy number falls is usually an artefact: a fixed
# number of contaminating or NUMT-derived reads makes up a larger fraction of a
# smaller mtDNA pool. Testing 1/copy_number against VAF therefore identifies
# these, and the 20 with the strongest evidence are shown.
#
# Panel code follows full_analysis_scripts/Compile_mitochondrial_data.R
# (~line 210). The per-mutation regressions take a few minutes, so the result is
# cached; delete the cache file to force a recomputation.
#-----------------------------------------------------------------------------------#

cn_corr_cache<-paste0(root_dir,"/data/CN_correlation_pval_df.Rds")
blood_ref_file<-paste0(root_dir,"/data/Samples_metadata_ref.csv")

#Tidy VAF table for the blood cohorts. Only shearwater-passing, non-indel calls
#with a non-zero VAF are tested.
df_tidy_full<-dplyr::bind_rows(Map(list=mito_data,exp_ID=names(mito_data),function(list,exp_ID) {
  if(is.null(list)) return(NULL)
  as.data.frame(list$matrices$vaf*list$matrices$SW)%>%
    tibble::rownames_to_column(var="mut_ref")%>%
    mutate(rho_vals=list$rho_vals)%>%
    dplyr::select(-any_of("global"))%>%
    tidyr::gather(key="Sample",value="vaf",-mut_ref,-rho_vals)%>%
    dplyr::filter(!grepl("DEL|INS",mut_ref))%>%
    dplyr::filter(vaf>0)%>%
    mutate(exp_ID=exp_ID)
}))

mito_cn_unique<-mito_cn%>%dplyr::filter(!duplicated(Sample))

test_muts<-df_tidy_full%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
  left_join(mito_cn_unique,by="Sample")%>%
  filter(!is.na(bedtools_mtDNA_genomes))%>%
  pull(mut_ref)%>%unique()

if(file.exists(cn_corr_cache)) {
  pval_df<-readRDS(cn_corr_cache)
  cat("ED2f: loaded cached regressions for",nrow(pval_df),"mutations\n")
} else {
  cat("ED2f: regressing",length(test_muts),"mutations against copy number (cached afterwards)\n")
  pval_df<-dplyr::bind_rows(lapply(test_muts,function(test_mut) {
    temp<-df_tidy_full%>%
      filter(mut_ref==test_mut)%>%
      mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
      left_join(mito_cn_unique,by="Sample")%>%
      mutate(implied_mut_CN=vaf*bedtools_mtDNA_genomes)%>%
      dplyr::select(Sample,exp_ID,mut_ref,vaf,bedtools_mtDNA_genomes,implied_mut_CN)
    usable<-function(d) nrow(d%>%filter(vaf!=0 & !is.na(bedtools_mtDNA_genomes)))>2
    if(!usable(temp)) return(data.frame(mut_ref=test_mut,coef=NA,pval=NA,max_vaf=NA))

    #A mutation seen almost entirely in one individual is assessed within that
    #individual; otherwise high implied copy number calls are dropped, as those
    #are the genuine mutations rather than the artefacts being looked for.
    max_prop<-max(table(temp$exp_ID))/length(temp$exp_ID)
    which_max_prop<-names(table(temp$exp_ID))[which.max(table(temp$exp_ID))]
    if(max_prop>0.95) {
      temp<-temp%>%filter(exp_ID==which_max_prop); assessment_type<-"single_sample"
    } else {
      temp<-temp%>%filter(implied_mut_CN<=25); assessment_type<-"multi_sample"
      if(!usable(temp)) return(data.frame(mut_ref=test_mut,coef=NA,pval=NA,max_vaf=NA))
    }

    lm.summary<-summary(lm(1/bedtools_mtDNA_genomes~vaf,
                           data=temp%>%filter(vaf!=0 & !is.na(bedtools_mtDNA_genomes))))
    if(!"vaf"%in%rownames(lm.summary$coefficients)) return(data.frame(mut_ref=test_mut,coef=NA,pval=NA,max_vaf=NA))
    vaf_coefs<-lm.summary$coefficients['vaf',]
    data.frame(mut_ref=test_mut,sum_of_vaf=sum(temp$vaf),assessment_type=assessment_type,
               sample_assessed=which_max_prop,coef=vaf_coefs[1],pval=vaf_coefs[4],
               r2=lm.summary$r.squared,max_vaf=temp%>%dplyr::slice_max(vaf,n=1)%>%pull(vaf))
  }))
  saveRDS(pval_df,file=cn_corr_cache)
}

#Positive coefficient only (VAF rising as copy number falls), then Benjamini-Hochberg
pval_df_filt<-pval_df%>%dplyr::filter(!is.na(pval)&!is.nan(pval)&coef>0)
pval_df_filt$qval<-p.adjust(pval_df_filt$pval,method="BH")
cat("ED2f: ",sum(pval_df_filt$qval<1e-2)," copy-number-correlating mutations at 1% FDR\n",sep="")

top20<-pval_df_filt%>%dplyr::filter(qval<1e-2)%>%dplyr::slice_min(order_by=qval,n=20)%>%pull(mut_ref)%>%sort()

if(length(top20)) {
  #Lymphocyte samples are excluded, as in the blood notebook, leaving the 12
  #individuals the Paired palette is sized for
  blood_ref_df<-read.csv(blood_ref_file)%>%filter(Dataset!="Lymphocyte")
  Individual_cols<-RColorBrewer::brewer.pal(12,"Paired")
  names(Individual_cols)<-blood_ref_df$Sample[order(blood_ref_df$Age)]

  CN_correlation_plot_20<-df_tidy_full%>%
    filter(mut_ref%in%top20)%>%
    mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
    left_join(mito_cn_unique,by="Sample")%>%
    mutate(implied_mut_CN=vaf*bedtools_mtDNA_genomes)%>%
    dplyr::select(exp_ID,Sample,mut_ref,vaf,bedtools_mtDNA_genomes,implied_mut_CN)%>%
    filter(implied_mut_CN<8 & vaf!=0)%>% #low implied copy number - the artefact regime
    mutate(mut_ref=gsub("MT_","",mut_ref))%>%
    ggplot(aes(x=vaf,y=1/bedtools_mtDNA_genomes,col=exp_ID))+
    geom_point(alpha=0.25,size=0.2)+
    geom_smooth(aes(x=vaf,y=1/bedtools_mtDNA_genomes),col="black",linewidth=0.5,method="lm",inherit.aes=FALSE)+
    facet_wrap(~mut_ref,ncol=10)+
    scale_color_manual(values=Individual_cols)+
    scale_x_log10()+
    scale_y_log10()+
    theme_bw()+
    my_theme+
    theme(strip.text.x=element_text(size=5,margin=unit(c(1,0,1,0),"mm")),
          axis.text.x=element_text(angle=90),
          legend.key.height=unit(3,"mm"))+
    labs(x="Mutation VAF",y="1/mtDNA genomes",col="Individual")+
    guides(colour=guide_legend(override.aes=list(size=2,alpha=1)))

  ggsave(filename=paste0(ed2_dir,"ExtDataFig2f.CN_correlation_plot_20.pdf"),CN_correlation_plot_20,height=3,width=7)
  cat("ED2f: written\n")
} else {cat("ED2f: no copy-number-correlating mutations found - skipped\n")}

cat("\nExtended Data Fig. 2 panels written to",ed2_dir,"\n")

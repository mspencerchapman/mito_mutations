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
genomeFile="~/Documents/Reference_files/hs37d5.fa"
root_dir="~/R_work/mito_mutations_blood"
source(paste0(root_dir,"/data/mito_mutations_blood_functions.R"))

#Set the key file paths using the root dir
tree_file_paths = list.files(paste0(root_dir,"/data/tree_files"),pattern=".tree",full.names = T)
ref_file=paste0(root_dir,"/data/Samples_metadata_ref.csv")
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
  left_join(sample_level_metadata,by="Sample")%>%
  left_join(ref_df%>%dplyr::rename("exp_ID"=Sample))
mito_cn$exp_ID<-factor(mito_cn$exp_ID,levels=ref_df$Sample[order(ref_df$Age)]) #Make the exp_ID a factor, with levels increasing by individual age
mito_cn<-left_join(mito_cn,phenotype_data_EM)

#Now import the mitochondrial mutation data
mito_data_file=paste0(root_dir,"/data/mito_data.Rds")
mito_data<-readRDS(mito_data_file)
CN_correlating_muts<-readRDS(paste0(root_dir,"/data/CN_correlation.RDS"))
exclude_muts=c("MT_302_A_C","MT_311_C_T","MT_567_A_C","MT_574_A_C","MT_16181_A_C","MT_16182_A_C","MT_16183_A_C","MT_16189_T_C") #These are the ones to exclude for the endometrial analysis

#----------------MITOCHONDRIAL COVERAGE STATISTICS----------------
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
ggsave(filename=paste0(figures_dir,"Supp_Figure_01/Coverage_plot.pdf"),coverage.plot,width=4,height=3)

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
ggsave(filename=paste0(figures_dir,"Supp_Figure_01/Coverage_median_plot.pdf"),median.coverage.plot,width=7,height=2)

median.coverage.ridges.plot<-mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  ggplot(aes(x=pilup_median_mtDNA_coverage,y=exp_ID))+
  ggridges::geom_density_ridges(linewidth=0.2,fill="lightblue")+
  scale_x_log10(labels=scales::label_comma())+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = margin(1,0,1,0, "mm")))+
  labs(x="Median mitochondrial DNA coverage",y="")
ggsave(filename=paste0(figures_dir,"Supp_Figure_01/median.coverage.ridges.plot.pdf"),median.coverage.ridges.plot,width=3.3,height=2.5)

uniformity_ridges.plot<-mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  ggplot(aes(x=perc_over0.8,y=exp_ID))+
  ggridges::geom_density_ridges(linewidth=0.2,fill="lightblue")+
  scale_x_log10()+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = margin(1,0,1,0, "mm")))+
  labs(x="Coverage uniformity\n(proportion of mtDNA genome with coverage >80% of mean)",y="")
ggsave(filename=paste0(figures_dir,"Supp_Figure_01/Uniformity_perc80_ridges_plot.pdf"),uniformity_ridges.plot,width=3.3,height=2.5)

#Print coverage summary statistics
mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  group_by(exp_ID)%>%
  summarise(samples=n(),mean_mtDNA_coverage=mean(bedtools_mtDNA_coverage),Under_1000X=sum(bedtools_mtDNA_coverage<1000))

mito_cn%>%
  dplyr::filter(Study%in%ref_df$Canapps.project & Sample%in% unlist(lapply(mito_data,function(list) list$tree$tip.label)))%>%
  summarise(mean_mtDNA_coverage=mean(bedtools_mtDNA_coverage))

#----------------MITOCHONDRIAL COPY NUMBER STATISTICS----------------
#Vector to convert the foetal cell type codes into basic 'HSC'/ 'Progenitor' types
convert_vec=c("HSC","HSC","HPC","HPC","HPC","HPC")
names(convert_vec)=c("H","HSC","C","M","HSPC","Progenitor")

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
ggsave(filename=paste0(figures_dir,"Additional_plots/mito_cn_by_pheno.pdf"),mito_cn_by_pheno,width=4,height=2)

mito_cn_by_pheno_all_HSCprog<-mito_cn%>%
  filter(!is.na(Phenotype))%>%
  mutate(Cell_type=ifelse(Cell_type=="Progenitor","HPC",Cell_type))%>%
  dplyr::select(exp_ID,Cell_type,Sample,bedtools_mtDNA_genomes,Phenotype)%>%
  ggplot(aes(x=factor(Cell_type),y=bedtools_mtDNA_genomes,col=factor(Cell_type)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.3,height=0,alpha=0.1,size=0.1)+
  ggpubr::stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=1.0,label.y=3.7,method="wilcox.test",size=1.8,fontface="italic")+
  scale_color_brewer(palette = "Set2")+
  scale_y_log10(limits=c(8,20000),breaks=c(10,100,1000,10000))+
  facet_grid(~factor(Phenotype,levels=pheno_levels))+
  theme_bw()+
  labs(x="",col="Founder cell\nphenotype",y="mtDNA copies per cell")+
  my_theme+
  theme(axis.text.x = element_text(angle=90),strip.text.x = element_text(size=6,margin = unit(c(1,0,1,0),"mm")),legend.position="none",legend.title = element_text(size=7))

#Generate the table of average CN by cell type - these values are overlaid onto the above plot
mito_cn%>%
  filter(!is.na(Phenotype))%>%
  mutate(Cell_type=ifelse(Cell_type=="Progenitor","HPC",Cell_type))%>%
  dplyr::select(exp_ID,Cell_type,Sample,bedtools_mtDNA_genomes,Phenotype)%>%
  group_by(Cell_type,Phenotype)%>%
  summarise(median_CN=median(bedtools_mtDNA_genomes))

ggsave(filename=paste0(figures_dir,"Figure_01/mito_cn_by_pheno_HSCprog.pdf"),mito_cn_by_pheno_all_HSCprog,width=4,height=2)

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
mtDNA.copy.number.logscale<-mito_cn%>%
  mutate(Cell_type=convert_vec[Cell_type])%>%
  filter(Sample%in%unlist(lapply(mito_data,function(list) list$tree$tip.label))&
           !is.na(bedtools_mtDNA_genomes)&
           !is.na(exp_ID))%>%
  ggplot(aes(x=factor(Sample,levels=mito_cn%>%dplyr::filter(Study%in%studies_to_include)%>%arrange(bedtools_mtDNA_genomes)%>%pull(Sample)%>%unique()),
             y=bedtools_mtDNA_genomes))+
  geom_point(size = 0.7, stroke = 0, shape = 16,col="black")+
  scale_y_log10(limits=c(8,20000),breaks=c(10,100,1000,10000))+
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

ggsave(filename = paste0(figures_dir,"Figure_01/mtDNA.copy.number.logscale.pdf"),plot=mtDNA.copy.number.logscale,height=2,width=7)

#Same plot, but now colour points by the colony phenotype
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

ggsave(filename = paste0(figures_dir,"Figure_01/mtDNA.copy.number.logscale.by.phenotype.pdf"),plot=mtDNA.copy.number.logscale.by.phenotype,height=2,width=7)

#Look in more details at the 18 pcw copy number by cell type & location
convert_vec2=c("HSC","MEP","CMP"); names(convert_vec2)=c("H","M","C")
convert_vec3=c("fBM","fLiver");names(convert_vec3)<-c("F","L")
foetal_18pcw_by_cell_type_loc<-mito_cn%>%
  filter(exp_ID=="18pcw")%>%
  mutate(Cell_type2=convert_vec[Cell_type])%>%
  mutate(Cell_type=convert_vec2[Cell_type],Tissue=convert_vec3[gsub("1|2","",Tissue.y)])%>%
  mutate(Cell_type_loc=paste(Tissue,Cell_type,sep="\n"))%>%
  ggplot(aes(x=factor(Cell_type_loc,levels=apply(expand.grid(convert_vec3,convert_vec2),1,paste,collapse="\n")),
             y=bedtools_mtDNA_genomes,
             col=Cell_type2))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.15,height=0,size=0.3,alpha=0.5)+
  scale_color_brewer(palette="Set2")+
  theme_bw()+
  my_theme+
  scale_y_log10(limits=c(8,20000),breaks=c(10,100,1000,10000))+
  theme(legend.title = element_text(size=7),legend.text = element_text(size=6))+
  labs(x="Cell type & location",y="mtDNA copies per cell",col="Founder cell\nphenotype")
ggsave(filename=paste0(figures_dir,"Figure_01/foetal_18pcw_by_cell_type_loc.pdf"),foetal_18pcw_by_cell_type_loc,width=3,height=2)

##--------------Plot mitochondrial CN by age--------------
median_mtDNA_CN_by_age<-mito_cn%>%
  filter(Sample%in%unlist(lapply(mito_data,function(list) list$tree$tip.label))&
           !is.na(bedtools_mtDNA_genomes)&
           !is.na(exp_ID)&
           Cell_type=="HSC")%>%
  group_by(exp_ID)%>% 
  summarise(mean=mean(bedtools_mtDNA_genomes,na.rm=T),median=median(bedtools_mtDNA_genomes,na.rm=T),sd=sd(bedtools_mtDNA_genomes,na.rm=T))%>%
  left_join(ref_df,by=c("exp_ID"="Sample"))%>%
  ggplot(aes(x=Age,y=median,col=factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])))+
  geom_smooth(aes(x=Age,y=median),method="lm",inherit.aes = F,col="black",size=0.5)+
  geom_point(alpha=0.8)+
  scale_color_manual(values = Individual_cols)+
  theme_bw()+
  scale_y_continuous(limits=c(0,1000))+
  labs(y="Median mitochondrial genomes per cell",col="Individual")+
  my_theme+
  theme(legend.key.height = unit(3,"mm"),legend.title=element_text(size=7),legend.text=element_text(size=6))
ggsave(filename = paste0(figures_dir,"Supp_Figure_01/median_mtDNA_CN_by_age.pdf"),plot = median_mtDNA_CN_by_age,width=3.5,height=2)

CN_by_age_summary<-mito_cn%>%
  filter(Sample%in%unlist(lapply(mito_data,function(list) list$tree$tip.label))&
           !is.na(bedtools_mtDNA_genomes)&
           !is.na(exp_ID)&
           Cell_type=="HSC")%>%
  group_by(exp_ID)%>% 
  summarise(mean=mean(bedtools_mtDNA_genomes,na.rm=T),median=median(bedtools_mtDNA_genomes,na.rm=T),sd=sd(bedtools_mtDNA_genomes,na.rm=T))%>%
  left_join(ref_df,by=c("exp_ID"="Sample"))
lm.CN_by_age<-lm(median~Age,data=CN_by_age_summary)
summary(lm.CN_by_age)
confint(lm.CN_by_age)

##---------------Test mitochondrial copy number for phylogenetic signal---------------

if(!file.exists(paste0(root_dir,"/data/mito_cn_phylo.Rds"))) {
  mito_cn_phylo<-lapply(mito_data,function(list) {
    tree.no.ancestral=drop.tip(list$tree.ultra,"Ancestral")
    mito_cn_tree<-sapply(tree.no.ancestral$tip.label,function(sampleID) mito_cn%>%filter(Sample==sampleID&Study%in%studies_to_include)%>%slice_head(n=1)%>%pull(bedtools_mtDNA_genomes))
    tree.4d<-phylobase::phylo4d(tree.no.ancestral,tip.data=mito_cn_tree[tree.no.ancestral$tip.label])
    barplot.phylo4d(tree.4d)
    res<-phyloSignal(tree.4d)
    mut.crlg <- phyloCorrelogram(tree.4d)
    return(list(mito_cn=mito_cn_tree,tree4d=tree.4d,phyloSignal_res=res,phyloCor_res=mut.crlg))
  })
  saveRDS(mito_cn_phylo,file=paste0(root_dir,"/data/mito_cn_phylo.Rds"))
} else {
  mito_cn_phylo<-readRDS(paste0(root_dir,"/data/mito_cn_phylo.Rds"))
}

#Plot the mitochondrial copy number onto the trees
par(mfrow=c(6,2))
Map(list=mito_data,mito_cn_phy=mito_cn_phylo,exp_ID=names(mito_data),function(list,mito_cn_phy,exp_ID){
  tree=plot_tree(list$tree.ultra,bars = log(mito_cn_phy$mito_cn),title=exp_ID)
})

## ----------------------ANALYSIS OF MITOCHONDRIAL MUTATION BURDENS----------------------

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

## ----------------------Print basic stats of total mutation numbers called----------------------

vaf_cut_off<-0
rho_cut_off<-0
sum(unlist(lapply(mito_data,function(list) {
  if(is.null(list)){stop(return(NULL))}
  vaf.filt<-(list$matrices$vaf*list$matrices$SW)[!grepl("DEL|INS",rownames(list$matrices$vaf))&
                                                   list$rho_vals>rho_cut_off&
                                                   !rownames(list$matrices$vaf)%in%exclude_muts,]
  total_pos=rowSums(vaf.filt>0,na.rm = T)
  return(sum(total_pos>0,na.rm=T))
})))

#Total number of SNVs identified by shearwater - once CN correlating muts excluded
sum(unlist(lapply(mito_data,function(list) {
  if(is.null(list)){stop(return(NULL))}
  CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
  vaf.filt<-(list$matrices$vaf*list$matrices$SW*CN_correlating_mut_removal_mat)[!grepl("DEL|INS",rownames(list$matrices$vaf))&
                                                                                  list$rho_vals>rho_cut_off&
                                                                                  !rownames(list$matrices$vaf)%in%exclude_muts,]
  total_pos=rowSums(vaf.filt>0,na.rm = T)
  return(sum(total_pos>0,na.rm=T))
})))

#Total number of SNVs identified by shearwater - once CN correlating muts excluded and only N1 muts included
sum(unlist(lapply(mito_data,function(list) {
  if(is.null(list)){stop(return(NULL))}
  CN_correlating_mut_removal_mat=list$matrices$implied_mutCN>mutCN_cutoff|(matrix((!rownames(list$matrices$vaf)%in%CN_correlating_muts),ncol=1)%*%matrix(rep(1,ncol(list$matrices$vaf)),nrow=1))
  vaf.filt<-(list$matrices$vaf*list$matrices$SW*(list$matrices$ML_Sig=="N1")*CN_correlating_mut_removal_mat)[!grepl("DEL|INS",rownames(list$matrices$vaf))&
                                                                                                               list$rho_vals>rho_cut_off&
                                                                                                               !rownames(list$matrices$vaf)%in%exclude_muts,]
  total_pos=rowSums(vaf.filt>0,na.rm = T)
  return(sum(total_pos>0,na.rm=T))
})))

Individual_burden_means_df<-df_tidy%>%
  dplyr::filter(Sample!="global"&
                  #ML_Sig=="N1"&
                  !mut_ref%in%exclude_muts)%>%
  #mutate(Sample=factor(Sample,levels=sample_order))%>%
  group_by(Sample)%>%
  dplyr::summarise(exp_ID=unique(exp_ID),n_mut=sum(vaf>vaf_cut_off))%>%
  tidyr::complete(Sample,fill=list(n_mut=0))%>%
  mutate(exp_ID=sapply(Sample,function(SampleID) {mito_cn$exp_ID[mito_cn$Sample==SampleID][1]}))%>%
  mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
  group_by(exp_ID)%>%
  summarise(mean=mean(n_mut))

samples_with_mut_df<-dplyr::bind_rows(lapply(seq(0.01,0.99,0.01),function(cut_off) {
  cut_off_df<-df_tidy%>%
    dplyr::filter(Sample!="global"&
                    #ML_Sig=="N1"&
                    !mut_ref%in%exclude_muts)%>%
    dplyr::mutate(Sample=factor(Sample,levels=sample_order))%>%
    group_by(Sample)%>%
    dplyr::summarise(n_mut=sum(vaf>cut_off))%>%
    tidyr::complete(Sample,fill=list(n_mut=0))%>%
    mutate(exp_ID=sapply(Sample,function(SampleID) {mito_cn$exp_ID[mito_cn$Sample==SampleID][1]}))%>%
    group_by(exp_ID)%>%
    dplyr::summarise(mean=mean(n_mut),n_with_mut=sum(n_mut>0),n_samp=n())%>%
    mutate(prop_with_mut=n_with_mut/n_samp,cut_off=cut_off)
  return(cut_off_df)
}))

prop.of.samples.with.mut<-samples_with_mut_df%>%
  dplyr::mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
  ggplot(aes(x=cut_off,y=prop_with_mut,col=factor(exp_ID,levels = ref_df$Sample[order(ref_df$Age)])))+
  #geom_point()+
  geom_line(alpha=0.6)+
  theme_bw()+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  scale_color_brewer(palette="Paired")+
  labs(x="Cut-off",
       y=str_wrap("Proportion of samples with at least 1 mutation with VAF > cut-off",width=40),
       col="Individual")+
  my_theme+
  theme(legend.key.height = unit(3,"mm"),legend.title = element_text(size=7))

ggsave(filename=paste0(figures_dir,"Figure_01/Sample_proportions_with_mut.pdf"),prop.of.samples.with.mut,width=3,height=2)

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

sum_of_vaf_df%>%
  left_join(mito_cn,by=c("Sample","exp_ID"))%>%
  dplyr::mutate(exp_ID=gsub("8 pcw","8pcw",exp_ID))%>%
  ggplot(aes(x=factor(Cell_type,levels = c("HSC","Progenitor")),col=Cell_type,y=sum_of_vaf))+
  ggpubr::stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=1.0,label.y=3.7,method="wilcox.test",size=1.8,fontface="italic")+
  geom_point(alpha=0.1)+
  facet_grid(~factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)]),drop=T)+
  geom_boxplot()+
  theme_bw()+
  my_theme

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
ggsave(filename=paste0(figures_dir,"Figure_01/mean_sum_of_vaf_by_age.pdf"),plot=mean_sum_of_vaf_by_age,height=2,width=3)

lm.mean_mutburden_by_age=lm(mean~Age,data=sum_of_vaf_summary%>%
                              left_join(ref_df,by=c("exp_ID"="Sample")))
summary(lm.mean_mutburden_by_age)$coefficients
summary(lm.mean_mutburden_by_age)$r.squared
confint(lm.mean_mutburden_by_age)
  
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
  geom_hline(aes(yintercept=median),linetype=2,data=sum_of_vaf_summary%>%mutate(median=ifelse(median==0,0.001,median)),col="red",size=0.5)+
  geom_text(aes(x=0,y=median_pos,label=paste0("tilde(x) == ",round(median,3))),size=2,nudge_x=+115,nudge_y=+0.5,data=sum_of_vaf_summary%>%mutate(median_pos=ifelse(median==0,0.001,median)),parse=T)+
  labs(x="Sample",y="Mutation burden\n(sum of VAF)")+
  theme(panel.spacing.x=unit(1, "mm"))
ggsave(filename = paste0(figures_dir,"Figure_01/sum_of_vaf_log.pdf"),plot=sum_of_vaf_plot_log,height=2,width=7)

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
ggsave(filename = paste0(figures_dir,"Figure_01/sum_of_vaf_log_by_celltype.pdf"),plot=sum_of_vaf_plot_log_by_celltype,height=2.2,width=7)

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
# ----------------------MUTATIONAL PROFILES----------------------
#-----------------------------------------------------------------------------------#

# read in the trinucleotide context
mtdna_trinuc_freq <- readRDS(paste0(root_dir,"/data/mtDNA_trinuc_freqs_coding_dloop_heavy_light.Rds"))

# set regions for coding and dloop regions
coding_region <- 577:16023
d_loop_region <- c(1:576,16024:16569)

trinucleotide_plot(mutations = df_tidy%>%
                     filter(Sample!="global"&
                              ML_Sig=="N1"&
                              !mut_ref%in%exclude_muts)%>%
                     separate(mut_ref,into=c("chr","pos","ref","mut"))%>%
                     dplyr::rename("donor"=exp_ID)%>%
                     mutate(pos=as.numeric(pos)),
                   analysis_region = "all_mtDNA",
                   analysis_type = "obs_exp",
                   file_name = paste0(plots_dir,"Mutational_signature.pdf"))

trinucleotide_plot(mutations = df_tidy%>%
                     filter(Sample!="global"&
                              ML_Sig=="N1"&
                              !mut_ref%in%exclude_muts&
                              vaf<0.01 & vaf>0.005)%>%
                     separate(mut_ref,into=c("chr","pos","ref","mut"))%>%
                     dplyr::rename("donor"=exp_ID)%>%
                     mutate(pos=as.numeric(pos)),
                   analysis_region = "all_mtDNA",
                   analysis_type = "obs_exp",
                   file_name = paste0(plots_dir,"Mutational_signature_0.005_0.01.pdf"))

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
#--------ANALYSIS OF MUTATION BURDEN ACCORDING TO CLONE EXPANSION/ SINGLETON--------
#-----------------------------------------------------------------------------------#

#Do clonal expansions have higher mutation burdens compared to singleton colonies?
#Need to take into account that different cells within a clade are not independent, therefore take an average
#This may slightly alter the distribution of mutation burdens as some of the variation is smoothed out (might effect the k-s test statistic)
#Therefore, could alternatively randomly pick a single sample

exp_vs_singleton_sov<-Map(list=mito_data[3:10],exp_ID=names(mito_data)[3:10],function(list,exp_ID) {
  
  #Using the 'get_expanded_clade_nodes' function to effectively cut across the tree at 100 mutations of molecular time
  #Some of the clades will be singletons
  ec_df<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = 100,min_clonal_fraction = 0,min_samples=1)
  
  ec_df$mean_sov<-sapply(ec_df$nodes,function(node) {
    Samples<-getTips(list$tree.ultra,node=node)
    mean_sov_in_clade<-sum_of_vaf_df%>%filter(Sample%in%Samples)%>%pull(sum_of_vaf)%>%mean(na.rm=T)
    return(mean_sov_in_clade)
  })
  
  return(ec_df%>%mutate(exp_ID=exp_ID,.before=1))
  
})%>%dplyr::bind_rows()

exp_vs_singleton_sov%>%
  mutate(exp_ID=factor(exp_ID,levels=ref_df%>%filter(Age>0)%>%arrange(Age)%>%pull(Sample)))%>%
  mutate(type=ifelse(n_samples==1,"Sing.","Exp."))%>%
  ggplot(aes(x=type,y=mean_sov,col=type,fill=type))+
  geom_violin(alpha=0.1)+
  geom_jitter(width=0.1,height=0,alpha=0.3)+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  facet_grid(~exp_ID)+
  theme_classic()+
  theme(legend.position="none")+
  stat_compare_means(method = "wilcox.test",label = "p.format",fontface="italic")+
  labs(x="",y="Mean mtDNA mutation burden\n (sum of VAF)")

all.ks.tests<-lapply(unique(exp_vs_singleton_sov$exp_ID),function(this_exp_ID) {
  ks<-ks.test(x=exp_vs_singleton_sov%>%filter(n_samples==1 & exp_ID==this_exp_ID)%>%pull(mean_sov),
          y=exp_vs_singleton_sov%>%filter(n_samples>1 & exp_ID==this_exp_ID)%>%pull(mean_sov))
  return(ks)
})

p.adjust(sapply(all.ks.tests,function(res) res$p.value),method = "BH")

#-----------------------------------------------------------------------------------#
# --------ANALYSIS OF MUTATION BURDEN ACCORDING TO DRIVER STATUS mut/ wt--------------------
#-----------------------------------------------------------------------------------#

#Could also do sub-analysis of selection here.. would probably do dN/dS per VAF bin looking at the means across the clone??

#1. Add basic nDNA driver mutation info to the object ----
mito_data<-Map(list=mito_data,exp_ID=names(mito_data),function(list,exp_ID) {
  annot_muts_file<-paste0(root_dir,"/data/blood_adult/annot_files_filtered/annotated_muts_filt_",exp_ID,".Rds")
  if(file.exists(annot_muts_file)) {
    cat(paste("Importing annotated mutations dataframe for",exp_ID),sep="\n")
    list$nDNA_mats<-readRDS(annot_muts_file)
    
    if(nrow(list$nDNA_mats$mat)==0) {
      stop(return(list))
    }
    
    tree_samples<-list$tree.ultra$tip.label[-which(list$tree.ultra$tip.label=="Ancestral")]
    
    if(any(!tree_samples%in%colnames(list$nDNA_mats$NV))) {
      cat("Dropping samples not found in the count matrix.")
      drop_tips<-tree_samples[!tree_samples%in%colnames(list$nDNA_mats$NV)]
      list$tree.ultra<-drop.tip(list$tree.ultra,tip = drop_tips)
      list$tree<-drop.tip(list$tree,tip = drop_tips)
      tree_samples<-tree_samples[tree_samples%in%colnames(list$nDNA_mats$NV)]
    }
    
    mtr<-as.matrix(cbind(list$nDNA_mats$NV[,tree_samples,drop=F],matrix(0,ncol = 1,nrow = nrow(list$nDNA_mats$NV),dimnames = list(rownames(list$nDNA_mats$NV),"Ancestral"))))
    dep<-as.matrix(cbind(list$nDNA_mats$NR[,tree_samples,drop=F],matrix(10,ncol = 1,nrow = nrow(list$nDNA_mats$NR),dimnames = list(rownames(list$nDNA_mats$NV),"Ancestral"))))
    res<-treemut::assign_to_tree(tree=list$tree.ultra,mtr=mtr,dep=dep)
    list$nDNA_mats$mat$node<-res$tree$edge[res$summary$edge_ml,2]
    return(list)
  } else {
    return(list)
  }
})

#2. Generate a similar 'sum of vaf' dataframe, but now incorporating driver information----
sov_comparison_df<-Map(list=mito_data[3:10],exp_ID=names(mito_data)[3:10],function(list,exp_ID) {
  cat(exp_ID,sep="\n")
  #Using the 'get_expanded_clade_nodes' function to effectively cut across the tree at 100 mutations of molecular time
  #Some of the clades will be singletons
  ec_df<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = 100,min_clonal_fraction = 0,min_samples=1)
  
  ec_df$mean_sov<-sapply(ec_df$nodes,function(node) {
    Samples<-getTips(list$tree.ultra,node=node)
    mean_sov_in_clade<-sum_of_vaf_df%>%filter(Sample%in%Samples)%>%pull(sum_of_vaf)%>%mean(na.rm=T)
    return(mean_sov_in_clade)
  })
  
  left_join(ec_df,list$nDNA_mats$mat%>%dplyr::select(node,mut_ref,Gene,CDS,Protein),by=c("nodes"="node"))%>%
    mutate(exp_ID=exp_ID,.before=1)
  
})%>%dplyr::bind_rows()

# Split samples by whether they harbour a DTA driver mutation or not
mut_set=c("DNMT3A","TET2","ASXL1")
IDs_to_include=paste0("KX00",c(7,8,4,3))

sov_mut_vs_wt_df<-sov_comparison_df%>%
  mutate(exp_ID=factor(exp_ID,levels=ref_df%>%filter(Age>0)%>%arrange(Age)%>%pull(Sample)))%>%
  replace_na(list(Gene=""))%>%
  mutate(type=ifelse(Gene%in%mut_set,"mut","wt"))%>%
  filter(exp_ID%in%IDs_to_include)

sov_mut_vs_wt<-sov_mut_vs_wt_df%>%
  ggplot(aes(x=type,y=mean_sov,col=type,fill=type))+
  geom_violin(alpha=0.1)+
  geom_jitter(width=0.1,height=0,alpha=0.3)+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  facet_grid(~exp_ID)+
  theme_classic()+
  my_theme+
  theme(legend.position="none")+
  ggpubr::stat_compare_means(method = "wilcox.test",label = "p.format",fontface="italic",vjust = +0.5,size=3)+
  labs(x="",y="Mean mtDNA mutation burden\n (sum of VAF)")

ggsave(filename = paste0(rebuttal_figs_dir,"sov_mut_vs_wt.pdf"),sov_mut_vs_wt,width=3.3,height=2.5)

mut_vs_wt.lmer<-lmerTest::lmer(mean_sov~type+(1|exp_ID),data=sov_mut_vs_wt_df)
summary(mut_vs_wt.lmer)
confint(mut_vs_wt.lmer)

mut_vs_wt_list<-Map(list=mito_data[5:8],exp_ID=names(mito_data)[5:8],function(list,exp_ID) {
  cat(exp_ID,sep="\n")
  #Using the 'get_expanded_clade_nodes' function to effectively cut across the tree at 100 mutations of molecular time
  #Some of the clades will be singletons
  ec_df<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = 100,min_clonal_fraction = 0,min_samples=1)%>%
    left_join(list$nDNA_mats$mat%>%dplyr::select(node,mut_ref,Gene,CDS,Protein),by=c("nodes"="node"))%>%
    mutate(exp_ID=exp_ID,.before=1)
  
  mut_nodes<-ec_df%>%filter(n_samples>1)%>%pull(nodes)
  wt_nodes<-ec_df%>%filter(n_samples==1)%>%pull(nodes)

  mut_samples<-unlist(lapply(mut_nodes,function(node) sample(x=getTips(list$tree.ultra,node=node),size=1)))
  wt_samples<-unlist(lapply(wt_nodes,function(node) sample(x=getTips(list$tree.ultra,node=node),size=1)))
  
  return(list(mut=mut_samples,wt=wt_samples))
})

all_mut<-unlist(lapply(mut_vs_wt_list,function(x) x$mut))
all_wt<-unlist(lapply(mut_vs_wt_list,function(x) x$wt))

IDs_to_include=paste0("KX00",c(7,8,4,3))

exp_vs_singleton_df<-sov_comparison_df%>%
  mutate(exp_ID=factor(exp_ID,levels=ref_df%>%filter(Age>0)%>%arrange(Age)%>%pull(Sample)))%>%
  filter(exp_ID%in%IDs_to_include)%>%
  mutate(type=ifelse(n_samples>1,"Clonal\nexpansion",ifelse(n_samples==1,"Singleton",NA)))%>%
  filter(!is.na(type))

sov_exp_vs_singleton<-exp_vs_singleton_df%>%
  ggplot(aes(x=type,y=mean_sov,col=type,fill=type))+
  geom_violin(alpha=0.1)+
  geom_jitter(width=0.1,height=0,alpha=0.3)+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  facet_grid(~exp_ID)+
  theme_classic()+
  my_theme+
  theme(legend.position="none",axis.text.x=element_text(angle=0))+
  ggpubr::stat_compare_means(method = "wilcox.test",label = "p.format",fontface="italic",vjust = +0.5,size=3)+
  labs(x="",y="Mean mtDNA mutation burden\n (sum of VAF)")

ggsave(filename = paste0(rebuttal_figs_dir,"sov_exp_vs_singleton.pdf"),sov_exp_vs_singleton,width=3.3,height=2.5)

exp_vs_singleton_df%>%group_by(type)%>%dplyr::summarise(n=n())
exp_vs_singleton.lmer<-lmerTest::lmer(mean_sov~type+(1|exp_ID),data=exp_vs_singleton_df)
summary(exp_vs_singleton.lmer)
confint(exp_vs_singleton.lmer)

#-----------------------------------------------------------------------------------#
# ------------------------------SELECTION ANALYSIS----------------------------------
#-----------------------------------------------------------------------------------#

library(dndscv)
mtref_rda_path=ifelse(Sys.info()['sysname']=="Darwin",paste0(root_dir,"/data/mtref.rda"),"/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/mtref.rda")
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
                          refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
  
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
  mutate(impact=ifelse(impact%in%c("Stop_loss","Nonsense"),"Truncating",impact))%>% #rename stop loss and nonsense mutations
  tidyr::unite(col="mut_ref",chr,pos,ref,mut,sep="_",remove=F)

my_comparisons <- list(c("Missense", "Synonymous"), c("Synonymous", "Truncating"))

### GENERATE FIG. 2B
mutation_category_by_vaf_violin_plot<-complete.annotated.mutation.table%>%
  dplyr::filter(impact!="Non-Coding")%>%
  dplyr::mutate(tissue=factor(tissue,levels=c("Foetal blood","Cord blood","Adult blood")))%>%
  ggplot(aes(impact, vaf, color = impact)) +
  geom_violin(aes(fill=impact),alpha=0.2) +
  geom_jitter(width=0.25,size = .1, alpha = 0.1) +
  theme_classic() + 
  ggsci::scale_color_lancet(palette = "lanonc") +
  ggsci::scale_fill_lancet(palette = "lanonc",) +
  scale_y_log10() + 
  facet_wrap(~ tissue, scales='free_x', nrow = 1) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons, size = 2) +
  labs(x="",y="Variant allele fraction",col="Impact",fill="Impact") + 
  my_theme+
  theme(legend.position="none")

ggsave(filename = paste0(rebuttal_figs_dir,"mutation_category_by_vaf_violin_plot.pdf"),mutation_category_by_vaf_violin_plot,width=3,height=2.5)


### Generate FIG. 2A ---------
dndscv_by_tissue<-lapply(c("Foetal blood","Cord blood","Adult blood"),function(tissue.id) {
  
  cat(tissue.id,sep="\n")
  tissue_info<-df_tidy%>%
    dplyr::filter(Tissue==tissue.id)%>%
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

globaldnds_res<-Map(dndsout=dndscv_by_tissue,tissue=c("Foetal blood","Cord blood","Adult blood"),function(dndsout,tissue) {
  dndsout$globaldnds%>%
    filter(!is.na(name) & complete.cases(.))%>%
    mutate(name=rename_vec[name])%>%
    mutate(tissue=tissue)
})%>%dplyr::bind_rows()

tissue_order=c("Foetal blood","Cord blood","Adult blood")
nb.cols <- length(tissue_order)
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)

max_dnds_value<-2
stats_to_include=c("Missense","Truncating")
global_dnds_blood<-globaldnds_res%>%
  filter(tissue%in%tissue_order & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         tissue=factor(tissue,levels=tissue_order))%>%
  ggplot(aes(x=tissue, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), size = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Tissue") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6))+
  scale_color_manual(values = mycolors[1:2], name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', linewidth = 0.5) +
  my_theme+
  #facet_wrap(~tissue,nrow=1)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.title=element_blank(),
        legend.position="none")

ggsave(filename = paste0(rebuttal_figs_dir,"global_dnds_blood.pdf"),global_dnds_blood,width=1.5,height=2.5)

### Generate Extended Data FIG. 2C ---------
gene_order=c(paste0("MT-ND",1:6),"MT-ND4L","MT-CYB",paste0("MT-CO",1:3),"MT-ATP6","MT-ATP8")
gene_order<-gene_order[gene_order%in%target_genes]

selcv_res_by_tissue<-Map(dndsout=dndscv_by_tissue,tissue=tissue_order,function(dndsout,tissue) {
  temp<-dndsout$sel_cv%>%
    dplyr::select(gene_name,wmis_cv,wnon_cv,qmis_cv,qtrunc_cv)
  
  dplyr::bind_rows(temp%>%dplyr::select(gene_name,"dNdS"=wmis_cv,"qval"=qmis_cv)%>%mutate(type="Missense"),
                   temp%>%dplyr::select(gene_name,"dNdS"=wnon_cv,"qval"=qtrunc_cv)%>%mutate(type="Nonsense"))%>%
    mutate(tissue=tissue)
})%>%dplyr::bind_rows()

dnds_heatmap_by_gene<-selcv_res_by_tissue%>%
  mutate(dnds_if_signif=ifelse(qval<0.1,paste0(round(dNdS,2),"*"),""),
         gene_name=factor(gene_name,levels=rev(gene_order)))%>%
  ggplot(aes(y=gene_name,x=tissue,fill=dNdS,label=dnds_if_signif))+
  geom_tile()+
  facet_grid(~type)+
  geom_text(aes(col=type),size=1.6)+
  scale_color_manual(values=c("black","white"),guide="none")+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),legend.position="none")+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(n=11,name = "RdYlGn"),limits=c(0,2))+
  labs(x="",y="")

ggsave(filename = paste0(rebuttal_figs_dir,"dnds_heatmap_by_gene.pdf"),dnds_heatmap_by_gene,width=1.7,height=2.5)

#-----------------------------------------------------------------------------------#
##----------------------Enhance the tidy dataframe with 'VAF bin' info---------------
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

#-----------------------------------------------------------------------------------#
##----------------------Run global dNdS by VAF group---------------
#-----------------------------------------------------------------------------------#
dndscv_by_vaf_adultblood<-lapply(new_VAF_groups,function(VAF_bin) {
  cat(VAF_bin,sep="\n")
  tissue_vaf_info<-df_tidy%>%
    dplyr::filter(Tissue=="Adult blood" & VAF_group==VAF_bin)%>%
    separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>% 
    dplyr::filter(!duplicated(.))%>% #only count each mutation once per individual in each VAF group
    arrange(sampleID,pos)
  
  mtdna.dndsout <- dndscv(tissue_vaf_info, gene_list=target_genes, 
                          refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
  return(mtdna.dndsout)
})

globaldnds_res_by_vaf<-Map(dndsout=dndscv_by_vaf_adultblood,VAF_bin=new_VAF_groups,function(dndsout,VAF_bin) {
  dndsout$globaldnds%>%
    filter(!is.na(name) & complete.cases(.))%>%
    mutate(name=rename_vec[name])%>%
    mutate(VAF_group=VAF_bin)
})%>%dplyr::bind_rows()

tissue_order=c("Foetal blood","Cord blood","Adult blood")
nb.cols <- length(tissue_order)
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)

max_dnds_value<-1.7

stats_to_include=c("Missense","Truncating")
global_dnds_by_vaf_plot<-globaldnds_res_by_vaf%>%
  filter(VAF_group%in%VAF_groups$labels[2:6] & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         VAF_group=factor(VAF_group,levels=VAF_groups$labels[2:6]))%>%
  ggplot(aes(x=VAF_group, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), size = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mycolors[1:2], name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', size = 0.5) +
  theme_classic()+
  my_theme+
  #facet_wrap(~tissue,nrow=1)+
  theme(legend.position="none",axis.text.x=element_text(angle=90))

ggsave(filename = paste0(rebuttal_figs_dir,"global_dnds_by_vaf_adult_blood.pdf"),global_dnds_by_vaf_plot,width=2,height=2.5)

#-----------------------------------------------------------------------------------#
##----------------------Get gene-level dNdS info by VAF group---------------
#-----------------------------------------------------------------------------------#
selcv_res_by_vaf<-Map(dndsout=dndscv_by_vaf_adultblood,VAF_bin=new_VAF_groups,function(dndsout,VAF_bin) {
  temp<-dndsout$sel_cv%>%
    dplyr::select(gene_name,wmis_cv,wnon_cv,qmis_cv,qtrunc_cv)
  
  dplyr::bind_rows(temp%>%dplyr::select(gene_name,"dNdS"=wmis_cv,"qval"=qmis_cv)%>%mutate(type="Missense"),
                   temp%>%dplyr::select(gene_name,"dNdS"=wnon_cv,"qval"=qtrunc_cv)%>%mutate(type="Nonsense"))%>%
        mutate(VAF_group=VAF_bin)
})%>%dplyr::bind_rows()%>%
  mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels))

dnds_heatmap_by_vaf_by_gene<-selcv_res_by_vaf%>%
  mutate(dnds_if_signif=ifelse(qval<0.1,paste0(round(dNdS,2),"*"),""),
         gene_name=factor(gene_name,levels=rev(gene_order)))%>%
  filter(VAF_group%in%VAF_groups$labels[2:6])%>%
  ggplot(aes(y=gene_name,x=VAF_group,fill=dNdS,label=dnds_if_signif))+
  geom_tile()+
  facet_grid(~type)+
  geom_text(aes(col=type),size=1.6)+
  scale_color_manual(values=c("black","white"),guide="none")+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(n=11,name = "RdYlGn"),limits=c(0,2))+
  labs(x="Variant allele fraction",y="")

ggsave(filename = paste0(rebuttal_figs_dir,"dnds_heatmap_by_vaf_by_gene.pdf"),dnds_heatmap_by_vaf_by_gene,width=4,height=2.5)
#-----------------------------------------------------------------------------------#
##--------Split analysis by whether young or old individual---------------
#-----------------------------------------------------------------------------------#

dndscv_by_young_or_old<-lapply(c("Young adult","Older adult"),function(tissue.id) {
  
  cat(tissue.id,sep="\n")
  if(tissue.id=="Young adult") {
    exp_IDs_to_include<-c("KX001","KX002","SX001")
  } else {
    exp_IDs_to_include<-c("KX003","KX004","KX007","KX008")
  }
  
  tissue_info<-df_tidy%>%
    dplyr::filter(exp_ID%in%exp_IDs_to_include)%>%
    separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>%
    dplyr::filter(!duplicated(.))%>%
    arrange(sampleID,pos)
  
  mtdna.dndsout <- dndscv(tissue_info, gene_list=target_genes, 
                          refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
  return(mtdna.dndsout)
})

global_dnds_by_young_or_old_status<-Map(dndsout=dndscv_by_young_or_old,status=c("Young adult","Older adult"),function(dndsout,status) {
  
    dndsout$globaldnds%>%
      filter(!is.na(name) & complete.cases(.))%>%
      mutate(name=rename_vec[name])%>%
    mutate(mut_status=status)
  
})%>%dplyr::bind_rows()

max_dnds_value<-2
stats_to_include=c("Missense","Truncating")
global_dnds_young_or_old_plot<-global_dnds_by_young_or_old_status%>%
  filter(name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh))%>%
  ggplot(aes(x=mut_status, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), size = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Tissue") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6))+
  scale_color_manual(values = mycolors[1:2], name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', linewidth = 0.5) +
  my_theme+
  #facet_wrap(~tissue,nrow=1)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.title=element_blank(),
        legend.position="none")

ggsave(filename = paste0(rebuttal_figs_dir,"global_dnds_blood.pdf"),global_dnds_blood,width=1.5,height=2.5)

dndscv_by_vaf_and_young_or_old_status<-lapply(c("Young adult","Older adult"),function(tissue.id) {
  
  if(tissue.id=="Young adult") {
    exp_IDs_to_include<-c("KX001","KX002","SX001")
  } else {
    exp_IDs_to_include<-c("KX003","KX004","KX007","KX008")
  }
  
  dndscv_by_vaf<-lapply(new_VAF_groups,function(VAF_bin) {
    cat(VAF_bin,sep="\n")
    tissue_vaf_info<-df_tidy%>%
      dplyr::filter(exp_ID%in%exp_IDs_to_include & VAF_group==VAF_bin)%>%
      separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
      mutate(pos=as.numeric(pos))%>%
      dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>% 
      dplyr::filter(!duplicated(.))%>% #only count each mutation once per individual in each VAF group
      arrange(sampleID,pos)
    
    mtdna.dndsout <- dndscv(tissue_vaf_info, gene_list=target_genes, 
                            refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
    return(mtdna.dndsout)
  })
  names(dndscv_by_vaf)<-new_VAF_groups
  return(dndscv_by_vaf)
})

global_dnds_by_vaf_and_young_or_old_status<-Map(list1=dndscv_by_vaf_and_young_or_old_status,status=c("Young adult","Older adult"),function(list1,status) {
  
  Map(dndsout=list1,VAF_bin=names(list1),function(dndsout,VAF_bin) {
    dndsout$globaldnds%>%
      filter(!is.na(name) & complete.cases(.))%>%
      mutate(name=rename_vec[name])%>%
      mutate(VAF_group=VAF_bin)
  })%>%dplyr::bind_rows()%>%
    mutate(mut_status=status)
  
})%>%dplyr::bind_rows()

global_dnds_by_vaf_and_young_or_old_status_plot<-global_dnds_by_vaf_and_young_or_old_status%>%
  filter(VAF_group%in%VAF_groups$labels[2:6] & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         VAF_group=factor(VAF_group,levels=VAF_groups$labels[2:6]),
         mut_status=factor(mut_status,levels=c("Young adult","Older adult")))%>%
  ggplot(aes(x=VAF_group, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), size = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Variant allele fraction") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mycolors[1:2], name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', size = 0.5) +
  theme_classic()+
  my_theme+
  facet_wrap(~mut_status,nrow=1)+
  theme(legend.position="top",axis.text.x=element_text(angle=90))

ggsave(filename = paste0(rebuttal_figs_dir,"global_dnds_by_vaf_and_young_or_old_status_plot.pdf"),global_dnds_by_vaf_and_young_or_old_status_plot,width=3.3,height=2.5)


#---------------------------------------------------------------------------------#
##--------Split analysis by whether DNMT3A/ TET2/ ASXL1 mutant or NOT---------------
#-----------------------------------------------------------------------------------#

#Generate a similar 'sum of vaf' dataframe, but now incorporating driver information
mut_vs_wt_list<-Map(list=mito_data[5:8],exp_ID=names(mito_data)[5:8],function(list,exp_ID) {
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
  
  mut_samples<-unlist(lapply(mut_nodes,function(node) getTips(list$tree.ultra,node=node)))
  wt_samples<-unlist(lapply(wt_nodes,function(node) getTips(list$tree.ultra,node=node)))
  
  return(list(mut=mut_samples,wt=wt_samples))
})

all_mut<-unlist(lapply(mut_vs_wt_list,function(x) x$mut))
all_wt<-unlist(lapply(mut_vs_wt_list,function(x) x$wt))

dndscv_by_vaf_and_mut_status<-lapply(list(mut=all_mut,wt=all_wt),function(sample_set) {
  dndscv_by_vaf<-lapply(new_VAF_groups,function(VAF_bin) {
    cat(VAF_bin,sep="\n")
    tissue_vaf_info<-df_tidy%>%
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
  names(dndscv_by_vaf)<-new_VAF_groups
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
rename_mut_status_vec=c("DNMT3A/ TET2/ ASXL1 mut","Wild-type")
names(rename_mut_status_vec)<-c("mut","wt")

global_dnds_by_mut_status_plot<-global_dnds_by_mut_status%>%
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
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', size = 0.5) +
  theme_classic()+
  my_theme+
  facet_wrap(~mut_status,nrow=1)+
  theme(legend.position="top",axis.text.x=element_text(angle=90))

ggsave(filename = paste0(rebuttal_figs_dir,"global_dnds_by_mut_status_plot.pdf"),global_dnds_by_mut_status_plot,width=3.5,height=2.5)


selcv_res_by_mut_status<-Map(list1=dndscv_by_vaf_and_mut_status,status=c("mut","wt"),function(list1,status) {
  
  Map(dndsout=list1,VAF_bin=new_VAF_groups,function(dndsout,VAF_bin) {
    temp<-dndsout$sel_cv%>%
      dplyr::select(gene_name,wmis_cv,wnon_cv,qmis_cv,qtrunc_cv)
    
    dplyr::bind_rows(temp%>%dplyr::select(gene_name,"dNdS"=wmis_cv,"qval"=qmis_cv)%>%mutate(type="Missense"),
                     temp%>%dplyr::select(gene_name,"dNdS"=wnon_cv,"qval"=qtrunc_cv)%>%mutate(type="Nonsense"))%>%
      mutate(VAF_group=VAF_bin)
  })%>%dplyr::bind_rows()%>%
    mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels),mut_status=status)
  
})%>%dplyr::bind_rows()

dnds_heatmap_by_mut_status_by_gene<-selcv_res_by_mut_status%>%
  mutate(dnds_if_signif=ifelse(qval<0.1,paste0(sprintf(dNdS, fmt = '%#.1f'),"*"),ifelse(dNdS>2,sprintf(dNdS, fmt = '%#.1f'),"")),
         dNdS=ifelse(dNdS>2,2,dNdS),
         gene_name=factor(gene_name,levels=rev(gene_order)),
         mut_status=rename_mut_status_vec[mut_status])%>%
  filter(VAF_group%in%VAF_groups$labels[2:6])%>%
  ggplot(aes(y=gene_name,x=VAF_group,fill=dNdS,label=dnds_if_signif))+
  geom_tile()+
  facet_grid(mut_status~type)+
  geom_text(col="white",size=1.5)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(n=11,name = "RdYlGn"))+
  labs(x="Variant allele fraction",y="")+
  theme(legend.position = "none")

ggsave(filename = paste0(rebuttal_figs_dir,"dnds_heatmap_by_mut_status_by_gene.pdf"),dnds_heatmap_by_mut_status_by_gene,width=2,height=3.5)

#-----------------------------------------------------------------------------------#
##--------Split analysis by whether in clonal expansion or NOT---------------
#-----------------------------------------------------------------------------------#

#Generate a similar 'sum of vaf' dataframe, but now incorporating driver information
mut_vs_wt_list<-Map(list=mito_data[5:8],exp_ID=names(mito_data)[5:8],function(list,exp_ID) {
  cat(exp_ID,sep="\n")
  #Using the 'get_expanded_clade_nodes' function to effectively cut across the tree at 100 mutations of molecular time
  #Some of the clades will be singletons
  ec_df<-get_expanded_clade_nodes(list$tree.ultra,height_cut_off = 100,min_clonal_fraction = 0,min_samples=1)%>%
    left_join(list$nDNA_mats$mat%>%dplyr::select(node,mut_ref,Gene,CDS,Protein),by=c("nodes"="node"))%>%
    mutate(exp_ID=exp_ID,.before=1)
  
  mut_nodes<-ec_df%>%filter(n_samples>1)%>%pull(nodes)
  wt_nodes<-ec_df%>%filter(n_samples==1)%>%pull(nodes)
  
  mut_samples<-unlist(lapply(mut_nodes,function(node) getTips(list$tree.ultra,node=node)))
  wt_samples<-unlist(lapply(wt_nodes,function(node) getTips(list$tree.ultra,node=node)))
  
  return(list(mut=mut_samples,wt=wt_samples))
})

all_mut<-unlist(lapply(mut_vs_wt_list,function(x) x$mut))
all_wt<-unlist(lapply(mut_vs_wt_list,function(x) x$wt))

dndscv_by_vaf_and_expansion_status<-lapply(list(mut=all_mut,wt=all_wt),function(sample_set) {
  dndscv_by_vaf<-lapply(new_VAF_groups,function(VAF_bin) {
    cat(VAF_bin,sep="\n")
    tissue_vaf_info<-df_tidy%>%
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
  names(dndscv_by_vaf)<-new_VAF_groups
  return(dndscv_by_vaf)
})

global_dnds_by_expansion_status<-Map(list1=dndscv_by_vaf_and_expansion_status,status=c("mut","wt"),function(list1,status) {
  
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
rename_mut_status_vec=c("Clonal expansion","Singleton")
names(rename_mut_status_vec)<-c("mut","wt")

global_dnds_by_expansion_status_plot<-global_dnds_by_expansion_status%>%
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
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', size = 0.5) +
  theme_classic()+
  my_theme+
  facet_wrap(~mut_status,nrow=1)+
  theme(legend.position="top",axis.text.x=element_text(angle=90))

ggsave(filename = paste0(rebuttal_figs_dir,"global_dnds_by_expansion_status_plot.pdf"),global_dnds_by_expansion_status_plot,width=3.5,height=2.5)

#Extract and combine the gene level analysis from the different dndscv output objects
selcv_res_by_expansion_status<-Map(list1=dndscv_by_vaf_and_expansion_status,status=c("mut","wt"),function(list1,status) {
  
  Map(dndsout=list1,VAF_bin=new_VAF_groups,function(dndsout,VAF_bin) {
    temp<-dndsout$sel_cv%>%
      dplyr::select(gene_name,wmis_cv,wnon_cv,qmis_cv,qtrunc_cv)
    
    dplyr::bind_rows(temp%>%dplyr::select(gene_name,"dNdS"=wmis_cv,"qval"=qmis_cv)%>%mutate(type="Missense"),
                     temp%>%dplyr::select(gene_name,"dNdS"=wnon_cv,"qval"=qtrunc_cv)%>%mutate(type="Nonsense"))%>%
      mutate(VAF_group=VAF_bin)
  })%>%dplyr::bind_rows()%>%
    mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels),mut_status=status)
  
})%>%dplyr::bind_rows()

dnds_heatmap_by_expansion_status_by_gene<-selcv_res_by_expansion_status%>%
  mutate(dnds_if_signif=ifelse(qval<0.1,paste0(sprintf(dNdS, fmt = '%#.1f'),"*"),ifelse(dNdS>2,sprintf(dNdS, fmt = '%#.1f'),"")),
         dNdS=ifelse(dNdS>2,2,dNdS),
         gene_name=factor(gene_name,levels=rev(gene_order)),
         mut_status=rename_mut_status_vec[mut_status])%>%
  filter(VAF_group%in%VAF_groups$labels[2:6])%>%
  ggplot(aes(y=gene_name,x=VAF_group,fill=dNdS,label=dnds_if_signif))+
  geom_tile()+
  facet_grid(mut_status~type)+
  geom_text(col="white",size=1.5)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90))+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(n=11,name = "RdYlGn"))+
  labs(x="Variant allele fraction",y="")

ggsave(filename = paste0(rebuttal_figs_dir,"dnds_heatmap_by_expansion_status_by_gene.pdf"),dnds_heatmap_by_expansion_status_by_gene,width=2.6,height=3.5)



#-----------------------------------------------------------------------------------#
#--------VAF DISTRIBUTION ANALYSIS--------
#-----------------------------------------------------------------------------------#

## 1. Re-define the VAF 'bins' that will be used for comparing VAF distributions ----
#Can set this on a log scale, or as decile bins
#boundaries<-seq(0,1,0.05)
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
muts_per_VAF_per_sample_mat<-Map(exp_ID=names(mito_data)[5:8],list=mito_data[5:8],function(exp_ID,list) {
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
  mut_vs_wt_list<-Map(list=mito_data[5:8],exp_ID=names(mito_data)[5:8],function(list,exp_ID) {
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

ggsave(filename = paste0(rebuttal_figs_dir,"VAF_dist_mut_vs_wt.pdf"),VAF_dist_mut_vs_wt,width=3.3,height=2.5)


#Repeat analysis but dividing by singleton vs member of clonal expansion
nboot=100
all_boot_res_clone_vs_singleton<-lapply(1:nboot,function(i) {
  cat(i,sep="\n")
  mut_vs_wt_list<-Map(list=mito_data[5:8],exp_ID=names(mito_data)[5:8],function(list,exp_ID) {
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

ggsave(filename = paste0(rebuttal_figs_dir,"VAF_dist_exp_vs_singleton.pdf"),VAF_dist_exp_vs_singleton,width=3.3,height=2.5)

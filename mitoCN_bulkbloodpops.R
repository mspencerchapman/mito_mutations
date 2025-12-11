library(dplyr)
library(tidyr)
library(ggplot2)

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


root_dir="~/R_work/mito_mutations_blood/"
plots_dir=paste0(root_dir,"rebuttal_plots/")
mitoCN_metadata<-readr::read_csv(paste0(root_dir,"data/mitoCN_metadata.csv"),col_select = 1:5)
table(mitoCN_metadata$`Cell type`)
#Import the mitochondrial copy number data
mito_cn=read.csv(paste0(root_dir,"data/whole_genome_coverage_pileup_and_bedtools_annotated.csv"),header=T)

sum(mitoCN_metadata$PD_number%in%mito_cn$Sample)

dataset_levels=c("Foetal","Cord blood","Adult")
cell_type_replace=c("HSC pool","CD34+CD38+\nHPC pool")
names(cell_type_replace)=c("HSC pool","CD34plusCD38plus haematopoietic progenitors")
data_clean<-left_join(mitoCN_metadata,mito_cn,by=c("PD_number"="Sample"))%>%
  dplyr::filter(!is.na(Tissue))%>%
  dplyr::select(PDID=PD_number,ID,Dataset,`Cell type`,"mtDNA copy number"=bedtools_mtDNA_genomes)%>%
  dplyr::mutate(`Cell type`=cell_type_replace[`Cell type`],Dataset=factor(Dataset,levels=dataset_levels))

#Plot with individual age group as facets
limits=c(0,2600)
mitoCN_age_facet<-data_clean%>%
  ggplot(aes(x=`Cell type`,y=`mtDNA copy number`))+
  geom_boxplot(linewidth=0.25,outlier.shape=NA)+
  geom_jitter(aes(col=ID),width=0.1,size=0.75,alpha=0.75)+
  scale_y_continuous(limits=limits)+
  facet_grid(~Dataset)+
  theme_classic()+
  my_theme+
  ggpubr::stat_compare_means(method="wilcox.test",label="p.format",size=2,vjust=+0)+
  theme(axis.text.x=element_text(angle=90),axis.title.x = element_blank())+
  labs(col="Sample ID")

ggsave(filename=paste0(plots_dir,"mitoCN_age_facet.pdf"),mitoCN_age_facet,width=4,height = 2.5)

#Plot with individual age group as facets
mitoCN_celltype_facet_noCB<-data_clean%>%
  filter(Dataset!="Cord blood")%>%
  ggplot(aes(x=`Dataset`,y=`mtDNA copy number`))+
  geom_boxplot(linewidth=0.25,outlier.shape = NA)+
  geom_jitter(aes(col=ID),width=0.1,size=0.75,alpha=0.75)+
  scale_y_continuous(limits=limits)+
  facet_grid(~`Cell type`)+
  theme_classic()+
  my_theme+
  ggpubr::stat_compare_means(method="wilcox.test",label="p.format",size=2,vjust=+0)+
  theme(axis.text.x=element_text(angle=90),axis.title.x = element_blank())+
  labs(col="Sample ID")

ggsave(filename=paste0(plots_dir,"mitoCN_celltype_facet_noCB.pdf"),mitoCN_celltype_facet_noCB,width=4,height = 2.5)

#Plot with cell type as facets
mitoCN_celltype_facet<-data_clean%>%
  ggplot(aes(x=`Dataset`,y=`mtDNA copy number`))+
  geom_boxplot(linewidth=0.25,outlier.shape = NA)+
  geom_jitter(aes(col=ID),width=0.1,size=0.75,alpha=0.75)+
  scale_y_continuous(limits=limits)+
  facet_grid(~`Cell type`)+
  theme_classic()+
  my_theme+
  ggpubr::stat_compare_means(label="p.format",size=2,vjust=+0.2)+
  theme(axis.text.x=element_text(angle=90),axis.title.x = element_blank())+
  labs(col="Sample ID")


ggsave(filename=paste0(plots_dir,"mitoCN_bulks_by_celltype.pdf"),mitoCN_celltype_facet,width=4,height = 2.5)

#Summarise by individual age group & cell type
data_clean%>%
  group_by(Dataset,`Cell type`)%>%
  dplyr::summarise(n=n(),mean_CN=mean(`mtDNA copy number`),median_CN=median(`mtDNA copy number`),sd_CN=sd(`mtDNA copy number`))%>%
  dplyr::mutate(sem=sd_CN/sqrt(n))

#Do between group comparisons
#1. Set up the groups
cats<-apply(expand_grid(unique(data_clean$Dataset),unique(data_clean$`Cell type`)),1,paste0,collapse="_")
comparisons<-expand.grid(cats,cats)%>%filter(Var1!=Var2)

#2. Do pairwise comparisons between groups
comparisons$p.value<-sapply(1:nrow(comparisons),function(i) {
  Dataset1=stringr::str_split(comparisons$Var1[i],pattern="_")[[1]][1]
  CellType1=stringr::str_split(comparisons$Var1[i],pattern="_")[[1]][2]
  Dataset2=stringr::str_split(comparisons$Var2[i],pattern="_")[[1]][1]
  CellType2=stringr::str_split(comparisons$Var2[i],pattern="_")[[1]][2]
  
  test.res<-wilcox.test(x=data_clean%>%filter(Dataset==Dataset1 & `Cell type`==CellType1)%>%pull(`mtDNA copy number`),
         y=data_clean%>%filter(Dataset==Dataset2 & `Cell type`==CellType2)%>%pull(`mtDNA copy number`))
  
  test.res$p.value
  
})

#Do ANOVA within the HPC compartment
summary(aov(data=data_clean%>%filter(`Cell type`=="CD34+CD38+\nHPC pool"),`mtDNA copy number`~Dataset))

#Do ANOVA within the HSC compartment
summary(aov(data=data_clean%>%filter(`Cell type`=="HSC pool"),`mtDNA copy number`~Dataset))

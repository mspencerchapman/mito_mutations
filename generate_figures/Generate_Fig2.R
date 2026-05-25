#-----------------------------------------------------------------------------------#
# --------Load packages (and install if they are not installed yet)-------------------
#-----------------------------------------------------------------------------------#
cran_packages=c("devtools","ape","stringr","dplyr","tidyr","ggplot2","gridExtra","phylosignal","ggpmisc")
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
# (therefore the mutational signature normalization will not work)
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

#-----------------------------------------------------------------------------------#
### Generate FIG. 2A ---------
#-----------------------------------------------------------------------------------#

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
    mutate(tissue=tissue,.before=1)
})%>%dplyr::bind_rows()

readr::write_csv(globaldnds_res,file = paste0(root_dir,"/tables/global_dnds_estimates.csv"))

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

ggsave(filename = paste0(plots_dir,"Figure_02/Fig2a.global_dnds_blood.pdf"),global_dnds_blood,width=1.5,height=2.5)

#-----------------------------------------------------------------------------------#
## Generate FIG. 2B----
#-----------------------------------------------------------------------------------#

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

ggsave(filename = paste0(plots_dir,"Figure_01/Fig2b.mutation_category_by_vaf_violin_plot.pdf"),mutation_category_by_vaf_violin_plot,width=3,height=2.5)


#-----------------------------------------------------------------------------------#
### Generate FIG. 2C ---------
#-----------------------------------------------------------------------------------#

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

ggsave(filename = paste0(plots_dir,"Figure_02/Fig2c.dnds_heatmap_by_gene.pdf"),dnds_heatmap_by_gene,width=1.7,height=2.5)

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
### Generate FIG. 2D ---------
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

ggsave(filename = paste0(plots_dir,"Figure_02/Fig2d.global_dnds_by_vaf_adult_blood.pdf"),global_dnds_by_vaf_plot,width=2,height=2.5)

#-----------------------------------------------------------------------------------#
### Generate FIG. 2E ---------
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

ggsave(filename = paste0(plots_dir,"Figure_02/Fig2e.dnds_heatmap_by_vaf_by_gene.pdf"),dnds_heatmap_by_vaf_by_gene,width=4,height=2.5)

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

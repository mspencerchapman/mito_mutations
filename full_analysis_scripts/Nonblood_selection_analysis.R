################################################################################################################################################
##                                                                                                                      
##  ANNOTATE THE MTDNA VARIANTS USING DNDSCV AND VISUALISE VAF PLOTS ACROSS TISSUES
##                                                                                                                      
##  Date: 4 DECEMBER 2024                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("BiocManager", "reshape2", "ggrepel", "readr", "stringr", "tidyverse", "dndscv", "ggsci", "ggpubr", "ape", "ggplot2", "gridExtra", "phylosignal")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))


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

#####################################################################################
# READ IN & PROCESS THE DATA
#####################################################################################

my_working_directory<-"~/Mounts/Lustre2/Mitochondria_study/nonblood"
R_function_files = list.files("~/Mounts/Lustre2/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="~/Mounts/Lustre2/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)
root_dir="~/R_work/mito_mutations_blood/"
plots_dir=paste0(root_dir,"rebuttal_plots/")

#Import the mitochondrial copy number data
mito_cn=read.csv(paste0(root_dir,"data/whole_genome_coverage_pileup_and_bedtools_annotated.csv"),header=T)
mtref_rda_path=paste0(root_dir,"data/mtref.rda")

#Import sample metadata
ref_df<-readxl::read_excel(paste0(root_dir,"data/non_blood_metadata.xlsx"))

translate=data.frame(dataset=c("KY","HL","SO","PR","LM","NW","lymph"),
                     al_ref=c("lung_organoid","colon","colon_ibd","muty_mutant","endometrium","blood_MPN","immune"))

muty_samples=c("PD44887","PD44888","PD44889","PD44890","PD44891")

##CHOOSE DATASET AND START ANALYSIS
data.list = list()
datasets = c("KY","HL","SO","PR","LM","NW","lymph","blood")
for (dataset in datasets){
  
  if(dataset=="blood") {
    dataset_mito_data<-readRDS("../mito_data.Rds")
    CN_correlating_muts<-readRDS("../CN_correlation.RDS")
  } else {
    CN_correlating_muts=readRDS(paste0("CN_correlation_",dataset,".RDS"))
    dataset_mito_data=readRDS(paste0("mito_mutation_data_",dataset,".RDS"))
  }
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
    
    if(dataset=="lymph") {list$tree$tip.label<-unique(list$sample_shearwater_calls$Sample)}
    
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
      dplyr::filter(!grepl("DEL|INS",mut_ref)&
                      #!mut_ref%in%list$het_oocyte_muts&
                      !mut_ref%in%exclude_muts)%>%
      dplyr::filter(vaf>=vaf_cut_off)%>%
      mutate(exp_ID=exp_ID)
    return(comb_tidy)
  }))
  
  #Filter the CN-correlating mutations - due to mis-mapping of nuclear reads.
  #There may be some genuine mutations at these sites, in which case the implied mutation copy number will be much higher than ~2 (here an arbitrary threshold of 8 is applied)
  df_tidy<-df_tidy%>%
    left_join(mito_cn)%>%
    mutate(implied_mut_CN=vaf*bedtools_mtDNA_genomes)%>%
    dplyr::filter(!(mut_ref%in%CN_correlating_muts & implied_mut_CN<=mutCN_cutoff & Sample!="global"))
  
  data.list[[dataset]] <- df_tidy
}

all.data.list = bind_rows(data.list)
all.data.list<-all.data.list[!is.na(all.data.list$Tissue),]
dim(all.data.list)

# rename data
all.data.list[all.data.list$Tissue == "lung_organoid", "Tissue"] <- "Lung"
all.data.list[all.data.list$Tissue == "endometrium_infertility", "Tissue"] <- "Endometrium"
all.data.list[all.data.list$Tissue == "endometrium", "Tissue"] <- "Endometrium"
all.data.list[all.data.list$Tissue == "colon_ibd", "Tissue"] <- "Colon"
all.data.list[all.data.list$Tissue == "colon", "Tissue"] <- "Colon"
all.data.list[all.data.list$Tissue == "muty_mutant", "Tissue"] <- "Colon"
all.data.list[all.data.list$Tissue == "immune", "Tissue"] <- "Immune"
all.data.list[all.data.list$Tissue == "blood_MPN", "Tissue"] <- "MPN"
all.data.list[all.data.list$Tissue %in% c("blood_emily","blood_foetal","blood_emily_hiseq"), "Tissue"] <- "blood"

#Split the blood individuals into 'foetal', 'cord blood' and 'adult'
foetal_IDs<-c("8 pcw","18 pcw")
cordblood_IDs<-c("CB001","CB002")
adultblood_IDs<-c(paste0("KX00",c(1:4,7,8)),"SX001","AX001")
all.data.list[all.data.list$exp_ID %in% foetal_IDs, "Tissue"] <- "Foetal blood"
all.data.list[all.data.list$exp_ID %in% cordblood_IDs & all.data.list$Tissue=="blood", "Tissue"] <- "Cord blood"
all.data.list[all.data.list$exp_ID %in% adultblood_IDs & all.data.list$Tissue=="blood", "Tissue"] <- "Adult blood"

#####################################################################################
# READ IN & PROCESS THE DATA
#####################################################################################

input.dir <- paste0(root_dir,"dnds_tables/")

#exclude ND6 from analysis since its on the other strand
target_genes <- c("MT-CYB", "MT-ND5", "MT-ND2", "MT-ND4", "MT-ND1", "MT-CO3", "MT-ATP6","MT-ND3", "MT-ATP8", "MT-ND4L", "MT-CO2", "MT-CO1")
#target_genes <- c("MT-CYB", "MT-ND5", "MT-ND2", "MT-ND4", "MT-ND1", "MT-CO3", "MT-ATP6","MT-ND3", "MT-ATP8", "MT-ND4L", "MT-CO2", "MT-CO1","MT-ND5")

valid.tissues <- unique(all.data.list$Tissue)
all.mutation.table <- list()
all.annotated.mutation.table <- list()

# iterate over each variant file
for (i in 1:length(valid.tissues)){
  
  # get the study id
  tissue.id <- valid.tissues[i]
  print(tissue.id)
  
  # read in the mtdna variant file and remove patient id
  mtdna.variant.data <- all.data.list[all.data.list$Tissue == tissue.id, ]
  mtdna.variant.data <- mtdna.variant.data %>% separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")
  mtdna.variant.data$pos<-as.numeric(mtdna.variant.data$pos)
  
  # next, reorder the columns to make them compatible with dnds input format
  mtdna.variant.data <- mtdna.variant.data[,c("Sample", "chr", "pos", "ref", "mut", "vaf")]
  
  # remove duplicated variants
  mtdna.variant.data <- mtdna.variant.data[!duplicated(mtdna.variant.data),]
  
  # order the mutations within each sample
  mtdna.variant.data <- mtdna.variant.data[order(mtdna.variant.data$Sample, mtdna.variant.data$pos),] 
  all.mutation.table[[i]] <- mtdna.variant.data
  
  # run dnds with the mtDNA variants to annotate them (the selection analysis isn't actually used here)
  mtdna.dndsout <- dndscv(mtdna.variant.data, gene_list=target_genes, 
                          refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
  
  # get the results with the annotated variants
  annotated.mtdna.variants <- mtdna.dndsout$annotmuts
  annotated.vaf.mutation.table <- merge(annotated.mtdna.variants, mtdna.variant.data, by.x = c("sampleID", "chr", "pos", "ref", "mut"), by.y = c("Sample", "chr", "pos", "ref", "mut"), all.y = T)
  annotated.vaf.mutation.table[is.na(annotated.vaf.mutation.table$impact), "impact"] <- "Non-Coding"
  annotated.vaf.mutation.table$tissue <- tissue.id
  all.annotated.mutation.table[[tissue.id]] <- annotated.vaf.mutation.table
}

all.mutation.table <- bind_rows(all.mutation.table)
complete.annotated.mutation.table <- bind_rows(all.annotated.mutation.table)
complete.annotated.mutation.table <- complete.annotated.mutation.table[complete.cases(complete.annotated.mutation.table$vaf),]

complete.annotated.mutation.table<-complete.annotated.mutation.table%>%
  left_join(all.data.list%>%dplyr::select(Sample,exp_ID),by=c("sampleID"="Sample"),relationship="many-to-many")%>%
  dplyr::rename("patientID"=exp_ID)%>%
  filter(!duplicated(.))%>%
  arrange(patientID,pos)%>% # order the dataframe
  mutate(impact=ifelse(impact%in%c("Stop_loss","Nonsense"),"Truncating",impact))%>% #rename stop loss and nonsense mutations
  tidyr::unite(col="mut_ref",chr,pos,ref,mut,sep="_",remove=F)

#####################################################################################
# REVIEW THE CODING CONSEQUENCES OF MUTATIONS USED FOR THE DRIFT ANALYSES
#####################################################################################

MPN_drift_muts<-c("MT_9804_G_A","MT_14111_T_C","MT_5136_T_C","MT_2270_A_G","MT_12565_T_C")
Blood_drift_muts<-c("MT_13552_G_A","MT_9151_A_G","MT_4232_T_C","MT_4487_A_G","MT_4965_A_G","MT_6379_T_C","MT_11790_T_C","MT_6180_G_A","MT_12396_T_C","MT_14035_T_C")

putative_MPN_muts=c("MT_11653_A_G","MT_13813_G_A","MT_8348_A_G","MT_6384_G_A","MT_14168_T_C","MT_11718_G_A","MT_5992_G_A","MT_14698_G_A")

no_coding_change_MPN_muts=data.frame(exp_ID=c("PD6646","PD9478","PD6629","PD5179","PD6646"),
                                     mut=c("MT_8348_A_G","MT_11653_A_G","MT_14168_T_C","MT_14698_G_A","MT_2270_A_G"))
readr::write_csv(x=no_coding_change_MPN_muts,file="~/Mounts/lustre2/Mitochondria_study/nonblood/Drift_ABC_MPN_nocoding/abc_muts.csv")

complete.annotated.mutation.table%>%
  dplyr::select(-sampleID,-vaf,-tissue,-patientID)%>%
  filter(mutation%in%c(Blood_drift_muts))%>%
  filter(!duplicated(.))%>%
  dplyr::select(chr,pos,ref,mut,gene,aachange,impact)%>%
  readr::write_csv(file=paste0(root_dir,"tables/blood_mutations_used_for_drift.csv"))

complete.annotated.mutation.table%>%
  dplyr::select(-sampleID,-vaf,-tissue,-patientID)%>%
  filter(mutation%in%c(MPN_drift_muts))%>%
  filter(!duplicated(.))%>%
  dplyr::select(chr,pos,ref,mut,gene,aachange,impact)%>%
  readr::write_csv(file=paste0(root_dir,"tables/MPN_mutations_used_for_drift.csv"))

complete.annotated.mutation.table%>%
  dplyr::select(-sampleID,-vaf,-tissue,-patientID)%>%
  filter(mutation%in%c(no_coding_change_MPN_muts$mut))%>%
  filter(!duplicated(.))%>%
  dplyr::select(chr,pos,ref,mut,gene,aachange,impact)%>%
  readr::write_csv(file=paste0(root_dir,"tables/MPN_nocoding_mutations_used_for_drift.csv"))

complete.annotated.mutation.table%>%
  #dplyr::select(-sampleID,-vaf,-tissue,-patientID)%>%
  filter(tissue=="MPN" & mutation%in%c(no_coding_change_MPN_muts$mut))%>%
  dplyr::distinct(chr,pos,ref,mut,.keep_all=T)

complete.annotated.mutation.table%>%
  dplyr::select(-sampleID,-tissue,-patientID)%>%
  filter(mutation%in%c(Blood_drift_muts))%>%group_by(mutation)%>%dplyr::summarise(max_vaf=max(vaf))

MPN_drift_muts%in%all.data.list$mut_ref

#####################################################################################
# VISUALIZE ALL ANNOTATED MTDNA MUTATIONS FIRST - PER GROUP
#####################################################################################

category.annotated.mutation.table <- complete.annotated.mutation.table
category.annotated.mutation.table$mutation <- paste(category.annotated.mutation.table$chr, category.annotated.mutation.table$pos, category.annotated.mutation.table$ref, 
                                                    category.annotated.mutation.table$mut, sep = "_")
category.annotated.mutation.table <- category.annotated.mutation.table[category.annotated.mutation.table$impact != "Non-Coding",]

# my_comparisons <- list(c("Missense", "Synonymous"), c("Synonymous", "Truncating"), c("Synonymous", "Non-Coding"), c("Missense", "Non-Coding"), c("Truncating", "Non-Coding"))
my_comparisons <- list(c("Missense", "Synonymous"), c("Synonymous", "Truncating"))

ggplot(category.annotated.mutation.table, aes(impact, vaf, color = impact)) +
  geom_violin(aes(fill=impact),alpha=0.2) +
  geom_jitter(width=0.3,size = .1, alpha = 0.1) +
  theme_classic() + 
  scale_color_lancet(palette = "lanonc") +
  scale_fill_lancet(palette = "lanonc") +
  scale_y_log10() + 
  facet_wrap(~ tissue, scales='free_x', nrow = 1) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  stat_compare_means(comparisons = my_comparisons, size = 2) +
  labs(x="",
       y="VAF") + 
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.75),
        strip.text = element_text(face="plain", size=6, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(color = "black", size = 0, angle = 0, hjust = .5, vjust = 0.5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 6, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 8, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 8, angle = 90, hjust = .5, vjust = .5, face = "plain", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 6, hjust = .5, face = "plain"), 
        legend.position = "bottom")



#####################################################################################
# RUN DNDSCV ON ALL THE MTDNA VARIANTS FROM EACH STUDY-----
#####################################################################################

sel_cv.complete <- list()
global.dnds.complete <- data.frame('name' = NA, "mle" = NA, "cilow" = NA, "cihigh" = NA, "tissue" = NA, stringsAsFactors = F)

# iterate over each variant file
for (i in 1:length(valid.tissues)){
  
  # get the study id
  tissue.id <- valid.tissues[i]
  print(tissue.id)
  
  #Get mtdna variant data in format for dnds, & only include each variant from each individual once
  mtdna.variant.data <- all.data.list%>%
    filter(Tissue==tissue.id)%>%
    separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")%>%
    mutate(pos=as.numeric(pos))%>%
    dplyr::select("sampleID"=exp_ID,chr,pos,ref,mut)%>% #Make compatible with dnds input format
    dplyr::filter(!duplicated(.))%>% # remove duplicated variants
    arrange(sampleID,pos) # order the mutations within each sample
  
  # run dnds with the mtDNA variants to annotate them
  mtdna.dndsout <- dndscv(mtdna.variant.data, gene_list=target_genes, 
                          refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
  
  # get the results with the annotated variants
  annotated.mtdna.variants <- mtdna.dndsout$annotmuts
  write.table(annotated.mtdna.variants, paste0(input.dir, "/", tissue.id, "_all_mtDNA_variants_dNdS_annotated.txt"),  sep = "\t", quote = F, row.names = F, col.names = T)
  
  # add results to list
  sel_cv <- mtdna.dndsout$sel_cv
  sel_cv$tissue <- tissue.id
  sel_cv.complete[[i]] <- sel_cv 
  
  # get the global dnds values
  mtdna.global.dnds <- mtdna.dndsout$globaldnds
  mtdna.global.dnds$tissue <- tissue.id
  
  # store them as well
  global.dnds.complete <- rbind(global.dnds.complete, mtdna.global.dnds)
  
}

sel_cv.complete <- bind_rows(sel_cv.complete)

# processing the output dataframe
rename_vec=c("Missense","Nonsense","Truncating","Overall")
names(rename_vec)=c("wmis","wnon","wtru","wall")
global.dnds.complete <- global.dnds.complete%>%
  filter(!is.na(name) & complete.cases(.))%>%
  mutate(name=rename_vec[name])

#####################################################################################
# VISUALIZE THE RESULTS EQUIVALENT TO MARTINCORENA ET AL. 2017 FIGURE 1C -----------
#####################################################################################

nb.cols <- length(unique(valid.tissues))
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)

max_dnds_value<-2
tissue_order=c("Adult blood","MPN","Immune","Cord blood","Foetal blood","Colon","Lung","Endometrium")
stats_to_include=c("Missense","Truncating")
tissue_dnds_plot<-global.dnds.complete%>%
  filter(tissue%in%tissue_order & name%in%stats_to_include)%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh),
         tissue=factor(tissue,levels=tissue_order))%>%
  ggplot(aes(x=name, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), size = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Tissue") +
  theme_classic() + 
  ylim(c(0,max_dnds_value))+
  scale_color_manual(values = mycolors[3:4], name="Mutation type") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', size = 0.5) +
  dnds_theme+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank())+
  facet_wrap(~tissue,nrow=1)

ggsave(plot=tissue_dnds_plot,filename=paste0(plots_dir,"tissue_dnds_plot.pdf"), width = 7, height = 2.5)


#####################################################################################
# RUN DNDSCV ON A PER VAF CATEGORY -------------------------------------------------
#####################################################################################

all.mutation.table = all.data.list
all.mutation.table = all.mutation.table %>% separate("mut_ref", c("chr", "pos", "ref", "mut"), "_")
all.mutation.table$pos<-as.numeric(all.mutation.table$pos)
all.mutation.table$PID <-all.mutation.table$exp_ID
all.mutation.table = all.mutation.table[,c("Sample", "chr", "pos", "ref", "mut", "vaf", "Tissue","PID")]

# remove duplicated variants
all.mutation.table <- all.mutation.table[!duplicated(all.mutation.table),]
all.mutation.table <- all.mutation.table[order(all.mutation.table$Sample, all.mutation.table$pos),] 
all.sample.mutation.table <- all.mutation.table[complete.cases(all.mutation.table),]

# divide mutations according to VAF
all.sample.mutation.table$vaf_cat <- "NA"
all.sample.mutation.table[all.sample.mutation.table$vaf < 0.05, "vaf_cat"] <- "<5%"
all.sample.mutation.table[all.sample.mutation.table$vaf < 0.1 & all.sample.mutation.table$vaf >= 0.05, "vaf_cat"] <- "5-10%"
all.sample.mutation.table[all.sample.mutation.table$vaf < 0.2 & all.sample.mutation.table$vaf >= 0.1, "vaf_cat"] <- "10-20%"
all.sample.mutation.table[all.sample.mutation.table$vaf < 0.5 & all.sample.mutation.table$vaf >= 0.2, "vaf_cat"] <- "20-50%"
all.sample.mutation.table[all.sample.mutation.table$vaf >= 0.5, "vaf_cat"] <- ">=50%"

table(all.sample.mutation.table$Tissue, all.sample.mutation.table$vaf_cat)

# <1% >=50%  1-5% 10-20% 20-50% 5-10%
# Blood Transplant 10234   588  7328    742    761  1016
# Colon            19463   198  1111    263    219   325
# Endometrium       1186   134   445     31     47    49
# Immune             401   199   205     72     99    74
# Liver             8246    23  1590     23     13    86

# 
# <5% >=50% 10-20% 20-50% 5-10%
# Colon       352   254    288    271   363
# Endometrium  81   155     68     70    79
# Immune       48   179     72     93    59
# Lung        127   279    162    180   183
# MPN         393   286    330    347   329

#####################################################################################
# ITERATE OVER EACH GROUP AND VAF CATEGORY AND STORE VALUES FOR IT
#####################################################################################

tissues <- valid.tissues
vaf.cat <- c("<5%", "5-10%", "10-20%", "20-50%", ">=50%")

vaf.cat.sel_cv.complete <- list()
vaf.cat.global.dnds.complete <- data.frame('name' = NA, "mle" = NA, "cilow" = NA, "cihigh" = NA, 
                                           "tissues" = NA, "vaf_cat" = NA, stringsAsFactors = F)


for(tissue in tissues){
  
  # subset data
  tissue.data <- all.sample.mutation.table[all.sample.mutation.table$Tissue == tissue,]
  print(tissue)
  
  sel_cv.complete <- list()
  global.dnds.complete <- data.frame('name' = NA, "mle" = NA, "cilow" = NA, "cihigh" = NA, 
                                     "tissues" = NA, "vaf_cat" = NA, stringsAsFactors = F)
  
  # cat <- vaf.cat[1]
  for (cat in vaf.cat){
    
    # subset data according to vaf category
    vaf.group.data <- tissue.data%>%
      filter(vaf_cat==cat)%>%
      #dplyr::select(Sample)%>%
      distinct(chr,pos,ref,mut,PID,.keep_all = T)
    
    print(paste0(cat, " - ", nrow(vaf.group.data)))
    
    if(nrow(vaf.group.data) > 35){
      
      # run dnds with the mtDNA variants to annotate them
      mtdna.dndsout <- dndscv(vaf.group.data, gene_list=target_genes, 
                              refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)
      
      # add results to list
      sel_cv <- mtdna.dndsout$sel_cv
      sel_cv$vaf_cat <- cat
      sel_cv.complete[[cat]] <- sel_cv   
      
      # get the global dnds values
      mtdna.global.dnds <- mtdna.dndsout$globaldnds
      mtdna.global.dnds$tissues <- tissue
      mtdna.global.dnds$vaf_cat <- cat
      mtdna.global.dnds <- mtdna.global.dnds[complete.cases(mtdna.global.dnds),]
      mtdna.global.dnds[mtdna.global.dnds$name == "wmis","name"] <- "Missense"
      mtdna.global.dnds[mtdna.global.dnds$name == "wnon","name"] <- "Nonsense"
      mtdna.global.dnds[mtdna.global.dnds$name == "wtru","name"] <- "Truncating"
      mtdna.global.dnds[mtdna.global.dnds$name == "wall","name"] <- "Overall"

      # store them as well
      global.dnds.complete <- rbind(global.dnds.complete, mtdna.global.dnds)
      
    } else {
      
      next
    }
    
  }
  
  sel_cv.complete <- bind_rows(sel_cv.complete) 
  sel_cv.complete$tissues <- tissue
  
  global.dnds.complete <- global.dnds.complete[-1,]
  vaf.cat.sel_cv.complete[[tissue]] <- sel_cv.complete
  vaf.cat.global.dnds.complete <- rbind(vaf.cat.global.dnds.complete, global.dnds.complete)
}

# combine
vaf.cat.sel_cv.complete <- bind_rows(vaf.cat.sel_cv.complete)

nb.cols <- length(unique(vaf.cat.global.dnds.complete$vaf_cat))
mycolors <- colorRampPalette(ggsci::pal_lancet(palette = "lanonc")(9))(nb.cols)

max_dnds_value<-2.5
blood_tissues<-c("Cord blood","Adult blood","MPN")
blood_dnds_by_vaf<-vaf.cat.global.dnds.complete%>%
  filter(!name %in% c("Overall","Nonsense") & !is.na(name) & tissues%in%blood_tissues)%>%
  mutate(vaf_cat=factor(vaf_cat,levels=vaf.cat),
         tissues=factor(tissues,levels=blood_tissues))%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh))%>%
  ggplot(aes(x=vaf_cat, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), size = 0.3, color="darkgrey",linetype=1) +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Variant allele frequency") +
  theme_classic() + 
  facet_wrap(~ tissues, nrow = 1) +
  ylim(0, max_dnds_value) +
  scale_color_manual(values = mycolors[3:4], name = "") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "") + 
  geom_hline(yintercept=1, linetype=2, col = 'darkgrey', size = 0.5) +
  coord_flip()+
  dnds_theme

ggsave(plot = blood_dnds_by_vaf, filename = paste0(plots_dir,"blood_dnds_by_vaf.pdf"),width = 7, height = 2.5)

max_dnds_value<-3.5
nonblood_tissues<-c("Lung","Colon","Endometrium","Immune")
nonblood_dnds_by_vaf<-vaf.cat.global.dnds.complete%>%
  filter(!name %in% c("Overall","Nonsense") & !is.na(name) & tissues%in%nonblood_tissues)%>%
  mutate(vaf_cat=factor(vaf_cat,levels=vaf.cat),
         tissues=factor(tissues,levels=nonblood_tissues))%>%
  mutate(cihigh=ifelse(cihigh>max_dnds_value,max_dnds_value,cihigh))%>%
  ggplot(aes(x=vaf_cat, y = mle, ymin = cilow, ymax = cihigh, color = name, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), size = 0.3, color="darkgrey",linetype=1) +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Variant allele frequency") +
  theme_classic() + 
  facet_wrap(~ tissues, nrow = 1) +
  ylim(0, max_dnds_value) +
  scale_color_manual(values = mycolors[3:4], name = "") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "") + 
  geom_hline(yintercept=1, linetype=2, col = 'darkgrey', size = 0.5) +
  coord_flip()+
  dnds_theme

ggsave(plot = nonblood_dnds_by_vaf, filename = paste0(plots_dir,"nonblood_dnds_by_vaf.pdf"),width = 7, height = 2.5)

# only plot the adult blood and the missense and nonsense mutations
subset.vaf.cat.global.dnds.complete <- vaf.cat.global.dnds.complete[vaf.cat.global.dnds.complete$name != "Truncating",]

ggplot(subset.vaf.cat.global.dnds.complete, aes(x=vaf_cat, y = mle, ymin = cilow, ymax = cihigh, color = vaf_cat, shape = name)) +
  geom_linerange(position= position_dodge2(width=0.75), size = 0.5, color="darkgrey") +
  geom_point(position=position_dodge2(width=0.75), size = 1.5) +
  labs(y = "Genome-wide dN/dS", x = "Variant allele frequency") +
  theme_classic() + 
  facet_wrap(~ tissues, nrow = 1) +
  scale_color_manual(values = mycolors) +
  ylim(0, 3) +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Mutation type") + 
  geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey', size = 0.5) +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.75),
        strip.text = element_text(face="plain", size=6, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1),
        axis.text.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = 0.5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 6, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 8, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 8, angle = 90, hjust = .5, vjust = .5, face = "plain", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(color = "black", size = 10, hjust = .5, face = "plain"), 
        legend.position = "top") + guides(color = "none")
ggsave("/Users/mp34/team154_campbell/plots/mtDNA/vaf_cat_dnds_mtDNA_multitissue_wo_<1_updated.pdf", width = 10, height = 2.5, dpi = 600)

#####################################################################################
# VISUALIZE THE DNDS COEFFICIENTS ACROSS GENES
#####################################################################################

vaf.cat.sel_cv.complete <- vaf.cat.sel_cv.complete[complete.cases(vaf.cat.sel_cv.complete),]

# replace nas with 1s
vaf.cat.sel_cv.complete[is.na(sel_cv.complete$qmis_cv), "qmis_cv"] <- 1
vaf.cat.sel_cv.complete[is.na(vaf.cat.sel_cv.complete$qtrunc_cv), "qtrunc_cv"] <- 1

# make numeric
vaf.cat.sel_cv.complete$wmis_cv <- as.numeric(as.character(vaf.cat.sel_cv.complete$wmis_cv))
vaf.cat.sel_cv.complete$wnon_cv <- as.numeric(as.character(vaf.cat.sel_cv.complete$wnon_cv))

# add a column for the colouring
vaf.cat.sel_cv.complete$mis_sig <- FALSE
vaf.cat.sel_cv.complete[vaf.cat.sel_cv.complete$qmis_cv <= 0.05, "mis_sig"] <- TRUE

# add a column for the colouring
vaf.cat.sel_cv.complete$non_sig <- FALSE
vaf.cat.sel_cv.complete[vaf.cat.sel_cv.complete$qtrunc_cv <= 0.05, "non_sig"] <- TRUE

vaf.cat.sel_cv.complete$vaf_cat <- factor(vaf.cat.sel_cv.complete$vaf_cat, levels = c("<5%", "5-10%", "10-20%", "20-50%", ">=50%"))


for (tissue in tissues){
  
  sel.group.data <- vaf.cat.sel_cv.complete[vaf.cat.sel_cv.complete$tissues == tissue,]
  
  # create matrices for missense mutations
  sel_cv.matrix <- as.data.frame((reshape2::acast(sel.group.data[,c("gene_name", "vaf_cat", "wmis_cv", "qmis_cv")], gene_name ~ vaf_cat, value.var = "wmis_cv", fill = "NA")))
  qmis.matrix <- as.data.frame((reshape2::acast(sel.group.data[,c("gene_name", "vaf_cat", "wmis_cv", "qmis_cv")], gene_name ~ vaf_cat, value.var = "qmis_cv", fill = "NA")))
  
  for (i in 1:ncol(sel_cv.matrix)){
    
    sel_cv.matrix[,i] <- as.numeric(as.character(sel_cv.matrix[,i]))
    
  }
  
  sel_cv.matrix <- as.matrix(t(sel_cv.matrix))
  qmis.matrix <- as.matrix(t(qmis.matrix))
  col_fun = colorRamp2(c(0, 1, 2), c("darkgreen", "white", "#B63679FF"))
  ht1 <- Heatmap(sel_cv.matrix, name = "Missense dN/dS", col = col_fun,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(qmis.matrix[i, j] < 0.05)
                     grid.text(sprintf("%.2f", sel_cv.matrix[i, j]), x, y, gp = gpar(fontsize = 10, fontface = "plain", col = "black"))
                 },
                 cluster_rows = F, 
                 cluster_columns = FALSE,
                 show_row_names = TRUE, 
                 row_names_side = "left",
                 show_row_dend = FALSE,
                 show_column_names = TRUE,
                 column_names_rot = 45,
                 row_names_gp = gpar(fontsize = 9, fontface = "plain"), 
                 column_names_gp = gpar(fontsize = 9, fontface = "plain"), 
                 rect_gp = gpar(col = "white", lwd = 1), border_gp = gpar(col = "black", lwd = 2),
                 column_title_gp = gpar(fontsize = 15, fontface = "plain"),
                 heatmap_legend_param = list(color_bar = "continous",
                                             at = c(0, 0.5, 1, 1.5, 2),
                                             title = "Missense dN/ds"), 
                 border = T)
  
  
  # pdf(paste0("/Users/mp34/team154_campbell/plots/mtDNA/", tissue, "_heatmap_wo_<1_missenseMLE_dnds.pdf"), width=10, height = 2.5, pointsize=0.1)
  # print(ht1)
  # dev.off()
  
  # create matrices for missense mutations
  sel_cv.non.matrix <- as.data.frame((reshape2::acast(sel.group.data[,c("gene_name", "vaf_cat", "wnon_cv", "qtrunc_cv")], gene_name ~ vaf_cat, value.var = "wnon_cv", fill = "NA")))
  qnon.matrix <- as.data.frame((reshape2::acast(sel.group.data[,c("gene_name", "vaf_cat", "wnon_cv", "qtrunc_cv")], gene_name ~ vaf_cat, value.var = "qtrunc_cv", fill = "NA")))
  
  for (i in 1:ncol(sel_cv.non.matrix)){
    
    sel_cv.non.matrix[,i] <- as.numeric(as.character(sel_cv.non.matrix[,i]))
    
  }
  
  sel_cv.non.matrix <- as.matrix(t(sel_cv.non.matrix))
  qnon.matrix <- as.matrix(t(qnon.matrix))
  col_fun = colorRamp2(c(0, 1, 2), c("darkgreen", "white", "#B63679FF"))
  ht2 <- Heatmap(sel_cv.non.matrix, name = "Nonsense dN/dS", col = col_fun,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(qnon.matrix[i, j] < 0.05)
                     grid.text(sprintf("%.2f", sel_cv.non.matrix[i, j]), x, y, gp = gpar(fontsize = 10, fontface = "plain", col = "black"))
                 },
                 cluster_rows = F, 
                 cluster_columns = FALSE,
                 show_row_names = TRUE, 
                 row_names_side = "left",
                 show_row_dend = FALSE,
                 show_column_names = TRUE,
                 column_names_rot = 45,
                 row_names_gp = gpar(fontsize = 9, fontface = "plain"), 
                 column_names_gp = gpar(fontsize = 9, fontface = "plain"), 
                 rect_gp = gpar(col = "white", lwd = 1), border_gp = gpar(col = "black", lwd = 2),
                 column_title_gp = gpar(fontsize = 15, fontface = "plain"),
                 heatmap_legend_param = list(color_bar = "continous",
                                             at = c(0, 0.5, 1, 1.5, 2),
                                             title = "Nonsense dN/ds"), 
                 border = T)
  
  
  # pdf(paste0("/Users/mp34/team154_campbell/plots/mtDNA/", tissue, "_heatmap_wo_<1_nonsenseMLE_dnds.pdf"), width=10, height = 2.5, pointsize=0.1)
  # print(ht2)
  # dev.off()
  
  # make a list out of both heatmaps
  ht_list <- ht2 %v%  ht1
  
  # plot and save the list as a complete figure
  pdf(paste0("/Users/mp34/team154_campbell/plots/mtDNA/", tissue, "_heatmap_wo_<1_dnds_combined_updated.pdf"), width=10, height = 5, pointsize=0.1)
  heatmap <- draw(ht_list, merge_legend = T,
                  heatmap_legend_side = "right",
                  gap = unit(0.25, "cm"), padding = unit(c(5,10,5,10), "mm"))
  dev.off()
  
}








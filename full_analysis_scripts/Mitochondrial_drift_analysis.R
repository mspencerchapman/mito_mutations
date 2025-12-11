library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

root_dir="~/R_work/mito_mutations_blood/"
figures_dir=paste0(root_dir,"figures/")
source(paste0(root_dir,"data/mito_mutations_blood_functions.R"))
mito_data_file=paste0(root_dir,"data/mito_data.Rds")
ref_file=paste0(root_dir,"data/Samples_metadata_ref.csv")

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

###Read in the necessary data for this analysis
mito_data<-readRDS(mito_data_file)
ref_df=read.csv(ref_file)%>%filter(Dataset!="Lymphocyte")
Individual_cols=RColorBrewer::brewer.pal(12,"Paired")
names(Individual_cols)<-ref_df$Sample[order(ref_df$Age)]

exposures<-read.csv(paste0(root_dir,"data/mutational_signatures/exposures.csv"))
components<-read.csv(paste0(root_dir,"data/mutational_signatures/components.csv"))
mutation_profiles_mat=read.table(paste0(root_dir,"data/mutational_signatures/trinuc_mut_mat.txt"))
sig_ref<-readRDS(file=paste0(root_dir,"data/mutational_signatures/sig_ref_file.Rds"))

#Clean up the row and column names
sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G"); ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
colnames(components)<-c(paste(full_vec,"H",sep="_"),paste(full_vec,"L",sep="_"))
rownames(components)<-paste0("N",0:(nrow(components)-1))
rownames(exposures)<-paste0("N",0:(nrow(exposures)-1))
colnames(exposures)<-gsub("^X","",colnames(exposures))
colnames(exposures)<-gsub("\\.pcw"," pcw",colnames(exposures))


VAF_groups=unique(sapply(stringr::str_split(colnames(exposures),pattern="_"),function(vec) paste(vec[2:3],collapse="_")))
new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
names(new_VAF_groups)<-VAF_groups

mut_numbers=data.frame(SampleID=rownames(mutation_profiles_mat),nmuts=rowSums(mutation_profiles_mat))

#Fill in the rare bins with too low mutation numbers - assume all are N1, as these are high VAF mutations
# (in samples with adequate mutations for extraction, high VAF mutations are all N1)
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
  geom_bar(stat="identity",fill="#FDBF6F",col="black",linewidth=0.25)+
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
  geom_bar(stat="identity",fill="#FDBF6F",col="black",linewidth=0.25)+
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

#-----------------------------------------------------------------------------------#
##----------------------VISUALIZE WRIGHT-FISHER MODELS OF DRIFT----------------------
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

#How many 'high VAF' (>50%) mutaitons are there depending on the generation?
gens_record%>%
  filter(VAF_group==">50%")%>%
  #mutate(total_gens=factor(total_gens,levels=seq(10,1500,10)))%>%
  dplyr::group_by(total_gens)%>%
  summarise(n=n())%>%
  ggplot(aes(x=total_gens,y=n))+
  scale_x_continuous(limits=c(0,1500))+
  geom_point()+
  geom_line()+
  theme_classic()+
  my_theme

#One notable feature is that these simulations suggest that high VAF mutations were acquired early in life
#whereas low VAF mutations were invariably acquired very recently.
acquisition_time_by_VAF<-gens_record%>%
  mutate(VAF_group=factor(VAF_group,levels=VAF_groups$labels))%>%
  filter(total_gens==1250)%>%
  ggplot(aes(x=as.numeric(gen),y=VAF_group))+
  ggridges::geom_density_ridges(panel_scaling = F)+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90),strip.text.x = element_text(size=7,margin = unit(c(1,0,1,0),"mm")))+
  labs(x="Time of mutation acquisition\n(WF Generations)",y="Density")
ggsave(filename = paste0(figures_dir,"Supp_Figure_04/Acquisition_time_by_VAF.pdf"),acquisition_time_by_VAF,width=7,height=2.5)

#Plot same data but as a ridges plot
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

ggsave(filename = paste0(figures_dir,"Supp_Figure_04/acquisition_time_by_VAF_ridges.pdf"),acquisition_time_by_VAF_ridges,width=3.3,height=2.5)

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

#-----------------------------------------------------------------------------------#
##----------------------PERFORM THE ABC----------------------
#-----------------------------------------------------------------------------------#

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

#-----------------------------------------------------------------------------------#
##----Now a version of the ABC with only synonymous/ non-coding mutation------------
#-----------------------------------------------------------------------------------#

# This to ensure that non-synonymous mutations are not introducing a bias in drift estimates due to their positive selection
CN_correlating_muts<-readRDS(paste0(root_dir,"/data/CN_correlation.RDS"))
length(CN_correlating_muts)
exclude_muts=c("MT_16183_A_C") #This is a specific artefact primarily in the 8pcw foetus due to a germline SNP and consequent frequent PCR slippage

#Define the VAF bins to use for signature extraction
bins=2^c(-10:0)

#Generate tidy data frame of samples, mutations and vafs
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
    dplyr::filter(!grepl("DEL|INS",mut_ref) & !mut_ref%in%exclude_muts)%>%
    dplyr::filter(vaf>0)%>%
    mutate(exp_ID=exp_ID)
  return(comb_tidy)
}))

###Annotate these using dndscv
library(dndscv)
mtref_rda_path=paste0(root_dir,"data/mtref.rda")
#All mtDNA genes
target_genes <- c("MT-CYB", "MT-ND5", "MT-ND2", "MT-ND4", "MT-ND1", "MT-CO3", "MT-ATP6","MT-ND3", "MT-ATP8", "MT-ND4L", "MT-CO2", "MT-CO1","MT-ND5")

dndscv_input<-df_tidy%>%tidyr::separate(col=mut_ref,into = c("chr","pos","ref","alt"),sep = "_")%>%
  dplyr::select("sampleID"=Sample,chr,pos,ref,alt)

dndscv_output <- dndscv(dndscv_input, gene_list=target_genes, 
                        refdb = mtref_rda_path, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf)

mtDNA.annotated<-dndscv_output$annotmuts%>%
  tidyr::unite(col = "mut_ref",chr,pos,ref,mut,sep = "_")%>%
  dplyr::rename("Sample"=sampleID)

df_tidy_annotated<-left_join(df_tidy,mtDNA.annotated,by=c("Sample","mut_ref"))%>%
  tidyr::replace_na(replace=list(impact="Non-coding"))
table(df_tidy_annotated$impact)
df_tidy_syn_noncoding_only<-df_tidy_annotated%>%filter(impact=="Synonymous"|impact=="Non-coding")

# get the results with the annotated variants
#Define the VAF 'bins' that will be used for comparing VAF distributions
new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")
VAF_groups=data.frame(
  labels=new_VAF_groups,
  lower_limit=c(0,2^(-10:-1)),
  upper_limit=2^(-10:0)
)

syn_muts_distribution_df<-lapply(names(mito_data),function(this_exp_ID) {
  print(this_exp_ID)
  exp_bins<-lapply(1:length(bins),function(i) {
    lower_limit=ifelse(i==1,0,bins[i-1])
    upper_limit=bins[i]
    print(paste(lower_limit,"to",upper_limit))
    bin_exp_muts<-df_tidy_syn_noncoding_only%>%
      dplyr::filter(exp_ID==this_exp_ID & vaf>lower_limit & vaf<= upper_limit)%>%
      dplyr::distinct(exp_ID,mut_ref,.keep_all = T)%>%
      nrow()
    
    data.frame(exp_ID=this_exp_ID,VAF_range=VAF_groups$labels[i],Signature="N1",abs_muts=bin_exp_muts)
  })
  return(exp_bins)
})%>%dplyr::bind_rows()

#Bind together mutation sets & annotate the trinucleotide context
real_syn_mut_distribution_plot_normalized<-left_join(syn_muts_distribution_df,nsamp_df)%>%
  mutate(muts_per_samp=abs_muts/n_samp,VAF_range=factor(VAF_range,levels=new_VAF_groups))%>%
  ggplot(aes(x=VAF_range,y=muts_per_samp))+
  geom_bar(stat="identity",fill="#FDBF6F",col="black",linewidth=0.25)+
  geom_line(aes(group=exp_ID))+
  theme_bw()+
  my_theme+
  theme(axis.text.x = element_text(size=5,angle=90),strip.text.x = element_text(size=6))+
  #facet_grid(rows=vars(factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)])),scales="free")+
  facet_wrap(~factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)]),scales="fixed",dir="v",ncol=6)+
  labs(x="VAF range",y="Number of mutations not affecting coding sequence detected per sample")
ggsave(filename=paste0(figures_dir,"Figure_04/synonymous_noncoding_distribution_plot.pdf"),real_mut_distribution_plot_normalized,width=7,height=2.5)

sumstats.syn.data<-left_join(syn_muts_distribution_df,nsamp_df)%>%
  mutate(muts_per_samp=abs_muts/n_samp)%>%
  dplyr::select(exp_ID,VAF_range,muts_per_samp)%>%
  pivot_wider(names_from="VAF_range",values_from="muts_per_samp")%>%
  replace(is.na(.), 0)

#-----------------------------------------------------------------------------------#
##----------------------REPEAT ABC WITH NO MISSENSE MUTATIONS----------------------
#-----------------------------------------------------------------------------------#

#This section uses summary statistics collected from many simulations and compares the results to those from the
#data stastics for each individual (as collected above).
#The actual simulations are run in the separate script "Mito_VAF_distribution_simulations.R"

library(abc)
all_sumstats=readRDS(file=paste0(root_dir,"data/Drift_ABC/VAF_distribution_ABC_simulation_sumstats_combined.Rds"))

new_VAF_groups=c("<0.1%","0.1-0.2%","0.2-0.4%","0.4-0.8%","0.8-1.6%","1.6-3.1%","3.1-6.2%","6.2-12.5%","12.5-25%","25-50%",">50%")

#Exclude the lowest VAF categories from abc that may be unreliably captured by sequencing
VAF_categories_to_include_in_ABC=new_VAF_groups[3:11]

mtDNA_CN=600
abc_syn_res_df<-lapply(1:nrow(sumstats.syn.data),function(i) {
  Exp_ID=sumstats.syn.data[i,1]
  res<-abc(target = as.numeric(sumstats.syn.data[i,VAF_categories_to_include_in_ABC]),
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

abc_synres_plot_mut_rate<-abc_syn_res_df%>%
  filter(Age>=0)%>%
  ggplot(aes(x=forcats::fct_reorder(exp_ID,Age),y=`muts_per_mitochondria_per_generation_50%`,
             ymin=`muts_per_mitochondria_per_generation_2.5%`,
             ymax=`muts_per_mitochondria_per_generation_97.5%`))+
  geom_smooth(method="lm",col="black",linewidth=0.5)+
  geom_point(alpha=0.75,size=0.5)+
  geom_errorbar(width=0.3,alpha=0.5)+
  theme_classic()+
  scale_y_continuous(limits=c(1e-4,8e-4),breaks=seq(2e-4,1e-3,2e-4))+
  labs(x="Individual",y="Mutations per mitochondria\nper generation")+
  my_theme+
  coord_flip()

ggsave(filename=paste0(figures_dir,"Figure_04/abc_synres_plot_mut_rate.pdf"),abc_synres_plot_mut_rate,width=1.75,height=2)

##Now do linear mixed effects modelling using all values of the posterior distribution and using individual as a random effect
abc_syn_res_df<-lapply(1:nrow(sumstats.syn.data),function(i) {
  Exp_ID=sumstats.syn.data[i,1]
  res<-abc(target = as.numeric(sumstats.syn.data[i,VAF_categories_to_include_in_ABC]),
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
lme.gens_by_age.syn<-lme4::lmer(total_generations~Age+(1|exp_ID),data=abc_syn_res_df)
summary(lme.gens_by_age.syn)
mtDNA_CN*365/lme.gens_by_age.syn@beta[2]
mtDNA_CN*365/confint(lme.gens_by_age.syn)["Age",]

#Plot these results
abc_synres_plot_total_generations<-abc_syn_res_df%>%
  dplyr::filter(!exp_ID%in%c("8pcw","18pcw"))%>%
  mutate(exp_ID=factor(exp_ID,levels=ref_df$Sample[order(ref_df$Age)]))%>%
  ggplot(aes(x=Age,y=total_generations))+
  geom_point(aes(col=exp_ID),alpha=0.05,size=0.25)+
  scale_y_continuous(limits=c(0,1700))+
  geom_abline(slope=lme.gens_by_age.syn@beta[2],intercept = lme.gens_by_age.syn@beta[1],linetype=1)+
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

ggsave(filename=paste0(figures_dir,"Figure_04/abc_synres_plot_total_generations.pdf"),abc_synres_plot_total_generations,width=2,height=2)


#------------------------------------------------------------------------------------------------#
##-------------- Driver mutation-specific analysis for the HSPCs cells------------------------------
#------------------------------------------------------------------------------------------------#

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


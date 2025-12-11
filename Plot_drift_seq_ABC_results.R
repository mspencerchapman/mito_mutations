library(dplyr)
library(readr)
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


exp_df<-data.frame(dirs=c("/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/Mitochondrial_ABC/Drift_ABC_normal_individual/",
                  "/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/Mitochondrial_ABC/Drift_ABC_normal_SEQ/",
                  "/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_MPN_individual/",
                  "/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_MPN/",
                  "/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_MPN_nocoding_individual/",
                  "/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_MPN_nocoding/",
                  "/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_CML_individual/",
                  "/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_CML/"),
           type=rep(c("normal","MPN","MPN_nocoding","CML"),each=2),
           abc_type=rep(c("individual","sequential"),times=4))


for(j in 1:nrow(exp_df)) {
  type=exp_df$type[j]
  abc_type=exp_df$abc_type[j]
  abc_dir=exp_df$dirs[j]
  
  plots_dir=paste0(abc_dir,"plots/")
  muts_file<-list.files(path=abc_dir,pattern="abc_muts")
  
  setwd(abc_dir)
  abc_muts<-readr::read_csv(muts_file)
  posterior_files<-list.files(path="output",pattern="posterior_table")
  
  prior=data.frame(generation_time=10^runif(1e4,min=-1,max=2.7),starting_vaf=runif(1e4))%>%
    mutate(mut="Prior",order=0)
  all_posts<-lapply(abc_muts$mut,function(mut) {readRDS(paste0("output/posterior_table_",mut,".Rds"))})
  
  comb_posts<-Map(post=all_posts,mut=abc_muts$mut,order=1:length(all_posts),function(post,mut,order) {
    as.data.frame(post)%>%dplyr::select(starting_vaf,generation_time)%>%mutate(mut=mut,order=order,.before=1)
  })%>%dplyr::bind_rows()%>%
    dplyr::bind_rows(prior)%>%
    mutate(mut=factor(mut,levels=c("Prior",abc_muts$mut)))
  
  #Print the posterior parameter summary stats (i.e. median and 95% CI)
  mtDNA_CN_used_in_simulations=600
  comb_posts%>%
    group_by(factor(order))%>%
    dplyr::summarise(lowerCI=quantile(generation_time,0.025),median=median(generation_time),upperCI=quantile(generation_time,0.975))%>%
    mutate(across(where(is.numeric), list(drift_param=function(x) x*mtDNA_CN_used_in_simulations)))
  
  seq_posts<-comb_posts%>%
    ggplot(aes(x=generation_time))+
    geom_histogram(bins=100,fill="pink",alpha=0.5,linewidth=0.1,col="gray")+
    scale_x_log10()+
    facet_grid(rows=vars(mut))+
    theme_classic()+
    my_theme+
    theme(strip.text.y=element_text(angle=0))+
    labs(x="Generation time (days)",y="Count")
  
  ggsave(filename = paste0(plots_dir,"seq_posts_",type,"_",abc_type,".pdf"),seq_posts,width=3.5,height=5)
  
  if(abc_type=="sequential") {
    mut_cols=colorRampPalette(RColorBrewer::brewer.pal(n=8,name="YlOrRd"))(length(all_posts)+1)
  } else if (abc_type=="individual") {
    mut_cols=colorRampPalette(RColorBrewer::brewer.pal(n=12,name="Paired"))(length(all_posts)+1)
  }
  
  names(mut_cols)<-levels(comb_posts$mut)
  seq_posts_ridges<-comb_posts%>%
    ggplot(aes(x=generation_time,y=mut,fill=mut))+
    ggridges::geom_density_ridges(alpha=0.9,linewidth=0.3,col="black",scale=3)+
    scale_fill_manual(values=mut_cols)+
    scale_x_log10(breaks=c(0.1,1,10,100),labels=c(0.1,1,10,100))+
    theme_classic()+
    my_theme+
    theme(strip.text.y=element_text(angle=0),legend.position = "none")+
    labs(x="Generation time (days)",y="Count")
  
  ggsave(filename = paste0(plots_dir,"seq_posts_ridges_",type,"_",abc_type,".pdf"),seq_posts_ridges,width=3.3,height=2)
}


#Plot the final posteriors for normal, MPN and CML drift rates against each other

plots_dir="/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/plots/"

drift_rate_comparison<-lapply(c("normal","MPN","CML"), function(setting) {
  
  abc_dir<-exp_df%>%filter(type==setting & abc_type=="sequential")%>%pull(dirs)
  muts_file<-list.files(path=abc_dir,pattern="abc_muts")
  
  setwd(abc_dir)
  abc_muts<-readr::read_csv(muts_file)
  posterior_files<-list.files(path="output",pattern="posterior_table")
  
  
  all_posts<-lapply(abc_muts$mut,function(mut) {readRDS(paste0("output/posterior_table_",mut,".Rds"))})
  prior=data.frame(generation_time=10^runif(1e4,min=-1,max=2.7),starting_vaf=runif(1e4))%>%
    mutate(mut="Prior",order=factor(0,levels=0:length(all_posts)))
  
  final_post_order<-nrow(abc_muts)
  
  comb_posts<-Map(post=all_posts,mut=abc_muts$mut,order=1:length(all_posts),function(post,mut,order) {
    as.data.frame(post)%>%dplyr::select(starting_vaf,generation_time)%>%mutate(mut=mut,order=factor(order,levels=0:length(all_posts)),.before=1)
  })%>%dplyr::bind_rows()%>%
    dplyr::bind_rows(prior)%>%
    mutate(mut=factor(mut,levels=c("Prior",abc_muts$mut)))
  
  final_post<-comb_posts%>%filter(order==final_post_order)%>%mutate(setting=setting)
  
  return(final_post)
  
})%>%dplyr::bind_rows()

prior=data.frame(generation_time=10^runif(1e4,min=-1,max=2.7),starting_vaf=runif(1e4))%>%
  mutate(mut="Prior",
         order=factor(0),
         setting="prior")

normal_disease_comparison_ridges_plot<-drift_rate_comparison%>%
  dplyr::bind_rows(prior)%>%
  mutate(setting=factor(setting,levels=c("prior","normal","MPN","CML")))%>%
  ggplot(aes(x=generation_time,y=setting,fill=setting))+
  ggridges::geom_density_ridges(linewidth=0.3)+
  scale_fill_brewer(palette = "RdPu")+
  scale_x_log10(breaks=c(0.1,1,10,100),labels=c(0.1,1,10,100))+
  theme_classic()+
  #ggridges::theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  my_theme+
  theme(strip.text.y=element_text(angle=0),legend.position = "none")+
  labs(x="Generation time (days)",y="Density")

ggsave(filename = paste0(plots_dir,"normal_disease_comparison_ridges_plot.pdf"),normal_disease_comparison_ridges_plot,width=2.5,height=2.2)
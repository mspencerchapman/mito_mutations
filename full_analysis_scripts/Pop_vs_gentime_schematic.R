##Plot schematic of population size against generation time given the drift information
library(ggplot2)
library(dplyr)
root_dir="~/R_work/mito_mutations_blood/"
figures_dir=paste0(root_dir,"figures/")

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


drift_parameter_ml=15960
drift_parameter_lowerCI=13660
drift_parameter_upperCI=19190

population_size_estimate_1=760
population_size_estimate_2=150
xrange=seq(45,5000,10)
text_size=1.5

pop_vs_gentime_schematic<-data.frame(x=xrange,y=drift_parameter_ml/xrange,ymin=drift_parameter_lowerCI/xrange,ymax=drift_parameter_upperCI/xrange)%>%
  ggplot(aes(x=x,y=y,ymin=ymin,ymax=ymax))+
  geom_line()+
  annotate(geom="text",label="Combinations consistent with drift rate",x=800,y=30,col="black",angle=313,size=text_size+0.5)+
  geom_ribbon(alpha=0.2)+
  annotate(geom="rect",xmin = 0,xmax=population_size_estimate_2,ymin = drift_parameter_ml/population_size_estimate_1,ymax=drift_parameter_ml/population_size_estimate_2,fill="#CCCCCC50")+
  annotate(geom="rect",xmin = population_size_estimate_2,xmax=population_size_estimate_1,ymin = 0,ymax=drift_parameter_ml/population_size_estimate_1,fill="#CCCCCC50")+
  geom_path(data=data.frame(x=c(0,population_size_estimate_1,population_size_estimate_1),y=c(drift_parameter_ml/population_size_estimate_1,drift_parameter_ml/population_size_estimate_1,0)),aes(x=x,y=y),col="red",linetype=2,inherit.aes = F)+
  annotate(geom="text",label="Pop. size = mtDNA CN",x=population_size_estimate_1+120,y=7,col="red",angle=270,size=text_size)+
  annotate(geom="text",label=paste0(round(drift_parameter_ml/population_size_estimate_1)," days"),x=60,y=(drift_parameter_ml/population_size_estimate_1)+5,col="red",size=text_size)+
  geom_path(data=data.frame(x=c(0,population_size_estimate_2,population_size_estimate_2),y=c(drift_parameter_ml/population_size_estimate_2,drift_parameter_ml/population_size_estimate_2,0)),aes(x=x,y=y),col="blue",linetype=2,inherit.aes = F)+
  annotate(geom="text",label="Pop. size = # of mitochondria",x=population_size_estimate_2+25,y=10,col="blue",angle=270,size=text_size)+
  annotate(geom="text",label=paste0(round(drift_parameter_ml/population_size_estimate_2)," days"),x=60,y=(drift_parameter_ml/population_size_estimate_2)+30,col="blue",size=text_size)+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  my_theme+
  labs(x="Effective mitochondrial\npopulation size",y="Generation time (days)")

ggsave(filename = paste0(figures_dir,"Figure_04/pop_vs_gentime_schematic.pdf"),pop_vs_gentime_schematic,width=1.75,height=2)

##This script is for submitting jobs on the farm in a way that "triages" the jobs to different queues based on the cell population size parameter
##To do this, the cell population parameter needs to be chosen outside of the simulation script and supplied as a parameter to the script

for(id in 1101:10000) {
  cat(id,sep="\n")
  log_total_cell_pop=runif(1,min=3,max=4.9)#2e3
  
  #Set the queue to submit to based on the cell population, which is the primary predictor of job time
  queue=ifelse(log_total_cell_pop>4.4,"basement",ifelse(log_total_cell_pop>3.9,"long","normal"))
  
  #combine the command to submit the job
  command=paste0("bsub -o $PWD/log_files/log.%J -e $PWD/err_files/err.%J -q ",queue," -G team78-grp -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J 'mito_dNdS' Rscript mito_dNdS_simulations.R -n ",id," -s 1 -p ",log_total_cell_pop," -o output")
  
  #review this command to check it makes sense
  cat(command,sep="\n")
  
  #submit it outside of R
  system(command)  
}
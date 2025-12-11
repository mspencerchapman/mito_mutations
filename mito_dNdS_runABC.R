#Import the parameters
param_files=list.files(path=".",pattern="parameters")
param_data_list=lapply(param_files,readr::read_delim,show_col_types = FALSE)
param_data_df<-Reduce(rbind,param_data_list)
hist(param_data_df$log_total_cell_pop)

#Import the summary stats
ss_files=list.files(path=".",pattern="summary")
length(ss_files)
ss_data_list=lapply(ss_files,readr::read_delim,show_col_types = FALSE)
ss_data_df<-Reduce(rbind,ss_data_list)

#match up the two data sets
ss_data_df<-ss_data_df[match(ss_data_df$run_id,param_data_df$run_id),]
param_data_df<-param_data_df[match(param_data_df$run_id,ss_data_df$run_id),]


library(abc)
full_param_cols=c("cell_fitness_param_Th","cell_fitness_param_h", "cell_fitness_param_n","fitness_threshold" ,"gamma_shape",
             "gamma_rate","mito_CN","total_mitochondrial_generations","n_cell_gen","n_mito_gen","log_total_cell_pop","log_mito_mutation_rate","prop_of_ns_under_selection","non_synonymous_prob")
full_ss_cols=c("mut_burden_quant_0.","mut_burden_quant_10.","mut_burden_quant_20.","mut_burden_quant_30.","mut_burden_quant_40.",
          "mut_burden_quant_50.","mut_burden_quant_60.","mut_burden_quant_70.","mut_burden_quant_80.","mut_burden_quant_90.",
          "mut_burden_quant_100.","total_per_cell_<1%","total_per_cell_1-5%","total_per_cell_5-10%","total_per_cell_10-20%",
          "total_per_cell_20-50%","total_per_cell_>50%","dNdS_<1%","dNdS_1-5%","dNdS_5-10%","dNdS_10-20%","dNdS_20-50%","dNdS_>50%")
param_cols=full_param_cols
ss_cols=full_ss_cols[19:23]
target=c(1.3,1.4,1.2,1,0.8)
abc_res=abc(target = target,param=as.matrix(param_data_df[,param_cols]),sumstat = as.matrix(ss_data_df[,ss_cols]),method = "rejection",tol = 0.05)



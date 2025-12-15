#Thie script goes through how to run the sequential ABC to infer mtDNA drift based on the data of varying heteroplasmy levels across a clonal expansion with a known phylogeny.

#Structure
#(1) Set the order of mutations for the sequential ABC, based on the order within the abc_muts dataframe in the scripts
#(2) Run the 20,000 simulations of mtDNA drift through the phylogeny for the 1st mutation. Here the prior is taken from the starting prior distribution.
# This is done with the script: Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_simulations.R
ABC_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_MPN_nocoding
SCRIPTS_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_MPN
MUTS_PATH=abc_muts_MPN_nocoding.csv

cd $ABC_DIR
mkdir -p log_files
mkdir -p output
QUEUE=normal
MEM=16000
N=1
bsub -o log_files/log.%J -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n10 -J mitoSIM Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_simulations.R -m $MUTS_PATH -j $N
#(3) Extract the summary statistics of each simulation. This is done in 200 separate batches of 100 simulations (as it is a slow step; particularly the 'phylosignal' function to infer phylogenetic signal)
# This is done with the script: 
cd $ABC_DIR
QUEUE=normal
MEM=16000
N=2
#This runs via a wrapper script where the LSB_JOBINDEX is fed into the R script as the batch index
bsub -o log_files/log.%J -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J ss_batch[1-200] bash ${SCRIPTS_DIR}/extract_sumstats_wrapper.sh $MUTS_PATH $N 

#bsub -o log_files/log.%J -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J mitoSIM Rscript Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_extract_sumstats.R -m $MUTS_PATH -j $N 

#(4) Combine the summary statistics into a single set and run the ABC to generate an updated posterior
cd $ABC_DIR
QUEUE=normal
MEM=16000
N=2
bsub -o log_files/log.%J -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J mito_comb Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_combine_sumstats.R -m $MUTS_PATH -j $N -b 200

#(5) Sequentially run through the remaining mutations updating the posterior each time (using the posterior from the previous mutation as the new prior)
cd $ABC_DIR
mkdir -p log_files
QUEUE=normal
MEM=16000
N=5
JID1=`uuidgen`
JID2=`uuidgen`
JID3=`uuidgen`
bsub -o log_files/log.${JID1} -e log_files/log.${JID1} -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n10 -J $JID1 Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_simulations.R -m $MUTS_PATH -j $N
bsub -o log_files/log.${JID2}.%I -w "done($JID1)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J "$JID2[1-200]" bash ${SCRIPTS_DIR}/extract_sumstats_wrapper.sh $MUTS_PATH $N 
bsub -o log_files/log.%J -w "done($JID2)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J mito_comb Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_combine_sumstats.R -m $MUTS_PATH -j $N -b 200


########################################################################
#### HERE ARE THE SCRIPTs FOR THE SEQUENTIAL NORMAL ABC
########################################################################

#Structure
#(1) Set the order of mutations for the sequential ABC, based on the order within the abc_muts dataframe in the scripts
#(2) Run the 20,000 simulations of mtDNA drift through the phylogeny for the 1st mutation. Here the prior is taken from the starting prior distribution.
# This is done with the script: Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_simulations.R
ABC_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/Mitochondrial_ABC/Drift_ABC_normal_SEQ
SCRIPTS_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/Mitochondrial_ABC/Drift_ABC_normal_SEQ
MUTS_PATH=abc_muts_normal.csv

cd $ABC_DIR
mkdir -p log_files
mkdir -p output
QUEUE=normal
MEM=16000
N=7
bsub -o log_files/log.%J -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n10 -J mitoSIM Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_normal_SEQ_simulations.R -m $MUTS_PATH -j $N
#(3) Extract the summary statistics of each simulation. This is done in 200 separate batches of 100 simulations (as it is a slow step; particularly the 'phylosignal' function to infer phylogenetic signal)
# This is done with the script: 
cd $ABC_DIR
QUEUE=normal
MEM=16000
N=6
#This runs via a wrapper script where the LSB_JOBINDEX is fed into the R script as the batch index
bsub -o log_files/log.%J -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J ss_batch[1-200] bash ${SCRIPTS_DIR}/extract_sumstats_wrapper.sh $MUTS_PATH $N 

#bsub -o log_files/log.%J -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J mitoSIM Rscript Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_extract_sumstats.R -m $MUTS_PATH -j $N 

#(4) Combine the summary statistics into a single set and run the ABC to generate an updated posterior
cd $ABC_DIR
QUEUE=normal
MEM=16000
N=6
bsub -o log_files/log.%J -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J mito_comb Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_normal_SEQ_combine_sumstats.R -m $MUTS_PATH -j $N -b 200

#(5) Sequentially run through the remaining mutations updating the posterior each time (using the posterior from the previous mutation as the new prior)
cd $ABC_DIR
mkdir -p log_files
QUEUE=normal
MEM=16000
N=10
JID1=`uuidgen`
JID2=`uuidgen`
JID3=`uuidgen`
bsub -o log_files/log.${JID1} -e log_files/log.${JID1} -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n10 -J $JID1 Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_normal_SEQ_simulations.R -m $MUTS_PATH -j $N
bsub -o log_files/log.${JID2}.%I -w "done($JID1)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J "$JID2[1-200]" bash ${SCRIPTS_DIR}/extract_sumstats_wrapper.sh $MUTS_PATH $N 
bsub -o log_files/log.%J -w "done($JID2)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J mito_comb Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_normal_SEQ_combine_sumstats.R -m $MUTS_PATH -j $N -b 200



####INDIVIDUAL ABC for MPN
ABC_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_MPN_individual
SCRIPTS_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_MPN_individual
MUTS_PATH=abc_muts_MPN.csv

cd $ABC_DIR
mkdir -p log_files
mkdir -p output
QUEUE=normal
MEM=16000

for N in {1..5};do
    JID1=`uuidgen`
    JID2=`uuidgen`
    JID3=`uuidgen`
    bsub -o log_files/log.${JID1} -e log_files/log.${JID1} -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n10 -J $JID1 Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_simulations.R -m $MUTS_PATH -j $N
    bsub -o log_files/log.${JID2}.%I -w "done($JID1)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J "$JID2[1-200]" bash ${SCRIPTS_DIR}/extract_sumstats_wrapper.sh $MUTS_PATH $N 
    bsub -o log_files/log.%J -w "done($JID2)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J mito_comb Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_combine_sumstats.R -m $MUTS_PATH -j $N -b 200
done

####INDIVIDUAL ABC for MPN nocoding
ABC_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_MPN_nocoding_individual
SCRIPTS_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_MPN_individual
MUTS_PATH=abc_muts_MPN_nocoding.csv

cd $ABC_DIR
mkdir -p log_files
mkdir -p output
QUEUE=normal
MEM=16000

for N in {1..5};do
    JID1=`uuidgen`
    JID2=`uuidgen`
    JID3=`uuidgen`
    bsub -o log_files/log.${JID1} -e log_files/log.${JID1} -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n10 -J $JID1 Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_simulations.R -m $MUTS_PATH -j $N
    bsub -o log_files/log.${JID2}.%I -w "done($JID1)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J "$JID2[1-200]" bash ${SCRIPTS_DIR}/extract_sumstats_wrapper.sh $MUTS_PATH $N 
    bsub -o log_files/log.%J -w "done($JID2)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J mito_comb Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_MPN_SEQ_combine_sumstats.R -m $MUTS_PATH -j $N -b 200
done


####INDIVIDUAL ABC for normal
ABC_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/Mitochondrial_ABC/Drift_ABC_normal_individual
SCRIPTS_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/Mitochondrial_ABC/Drift_ABC_normal_individual
MUTS_PATH=abc_muts_normal.csv

cd $ABC_DIR
mkdir -p log_files
mkdir -p output
QUEUE=normal
MEM=16000

for N in {1..10};do
    JID1=`uuidgen`
    JID2=`uuidgen`
    JID3=`uuidgen`
    bsub -o log_files/log.${JID1} -e log_files/log.${JID1} -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n10 -J $JID1 Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_normal_SEQ_simulations.R -m $MUTS_PATH -j $N
    bsub -o log_files/log.${JID2}.%I -w "done($JID1)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J "$JID2[1-200]" bash ${SCRIPTS_DIR}/extract_sumstats_wrapper.sh $MUTS_PATH $N 
    bsub -o log_files/log.%J -w "done($JID2)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J mito_comb Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_normal_SEQ_combine_sumstats.R -m $MUTS_PATH -j $N -b 200
done



####INDIVIDUAL ABC for CML
ABC_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_CML
SCRIPTS_DIR=/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/Drift_ABC_CML
MUTS_PATH=abc_muts_CML.csv

cd $ABC_DIR
mkdir -p log_files
mkdir -p output
QUEUE=normal
MEM=16000


N=3
JID1=`uuidgen`
JID2=`uuidgen`
JID3=`uuidgen`
bsub -o log_files/log.${JID1} -e log_files/log.${JID1} -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n10 -J $JID1 Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_CML_SEQ_simulations.R -m $MUTS_PATH -j $N
bsub -o log_files/log.${JID2}.%I -w "done($JID1)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J "$JID2[1-200]" bash ${SCRIPTS_DIR}/extract_sumstats_wrapper.sh $MUTS_PATH $N 
bsub -o log_files/log.%J -w "done($JID2)" -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J mito_comb Rscript ${SCRIPTS_DIR}/Mitochondrial_drift_through_phylo_ABC_CML_SEQ_combine_sumstats.R -m $MUTS_PATH -j $N -b 200


    
    



bsub -o log_files/log.${JID2}.%I -e log_files/log.%J -q $QUEUE -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" -M${MEM} -n1 -J "$JID2[1-200]" bash ${SCRIPTS_DIR}/extract_sumstats_wrapper.sh $MUTS_PATH $N 
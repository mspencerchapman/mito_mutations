# Mitochondrial mutation calling pipeline
## Devised by Andrew Lawson

## Step 1: Set the variables for the tissue


```bash
ALL_PROJECT_NUMBERS=3044,3045,3224
TISSUE_NAME='blood_CML'
ROOT_DIR='/lustre/scratch125/casm/team268im/al28/mtDNA/'
SCRIPT_FOLDER='/lustre/scratch126/casm/team154pc/ms56/my_programs/mtDNA_code/'

module load ISG/sanger-samtools-refpath
module load R/4.4.0

```

## Step 2: Generate a sample list

This must be a tab-separated file ('.tsv') with project/ sample info for all samples in the 'tissue' (or study). \
The first column is the project, and the second the sample name. \
There should not be any column names. This can be generated manually, otherwise the following code block will generate a suitable list of all samples in the listed projects. \
(Note that if the list should only contain a subset of samples in any of the given projects it will need to be edited manually)

```bash

IFS=',' read -r -a PROJECT_ARRAY <<< "$ALL_PROJECT_NUMBERS"

for PROJECT_NUMBER in "${PROJECT_ARRAY[@]}"; do
    PROJECT_SAMPLES=$(ls /nfs/cancer_ref01/nst_links/live/${PROJECT_NUMBER})
    for SAMPLE in $PROJECT_SAMPLES; do
      echo -e "${PROJECT_NUMBER}\t$SAMPLE"
    done
done>${ROOT_DIR}/sample_lists/${TISSUE_NAME}.tsv

```

## Step 3: Extract the mitochondrial bam files and generate index files

Generate the mtDNA bams

```bash

#Make and go into the Data directory for the tissue
mkdir -p ${ROOT_DIR}/Data/${TISSUE_NAME}
cd ${ROOT_DIR}/Data/${TISSUE_NAME}

#Set off the wrapper script for generating the mtDNA bam files
#This will check if the file is already present in the folder and only run it for that sample if not there
Rscript ${SCRIPT_FOLDER}/mtDNA_bam_generation_GRCh38.R ${ROOT_DIR}/sample_lists/${TISSUE_NAME}.tsv

```

Then index the bam files (run this once all the bams have been generated)
This runs the bedtools index command across all the generated bams.

```bash

cd ${ROOT_DIR}/Data/${TISSUE_NAME}
Rscript ${SCRIPT_FOLDER}/mtDNA_bai_generation.R

```

## Step 4a: Create the pileup files

This wrapper script will call the script 'mtDNA_genome_wide_pileup_GRCh38.R' in the SCRIPTS_FOLDER. \
This script uses the 'bam2R' function from deepSNV to generate a pileup across all sites in the mtDNA genome, outputting it as a csv file. \
It must be run in the Data folder for the tissue. \
This can take a few hours for each sample, so it is run in the long queue.

```bash

cd ${ROOT_DIR}/Data/${TISSUE_NAME}

#This will check if the file is already present in the folder and only generate if not there
Rscript ${SCRIPT_FOLDER}/mtDNA_genome_wide_pileup_wrapper_GRCh38.R ${ROOT_DIR}/sample_lists/${TISSUE_NAME}.tsv

```

## Step 4b-1: Run bedtools to generate coverage histograms for all samples

This files are later used to calculate the mtDNA copy number.
Script must be run from the sample list folder. Note that this is a wrapper script that will set off LSF jobs for each file independently using the 'bedtools genomecov' command.
It will check to see if the output files already exist before setting off the job to generate new files.
Note that this script currently contains hard file paths to Andrew's file system.

```bash
#Make the folder for the output if it does already exist
mkdir -p ${ROOT_DIR}/Analysis/bedtools_coverage/${TISSUE_NAME}

cd ${ROOT_DIR}/sample_lists
Rscript ${SCRIPT_FOLDER}/mtDNA_bedtools_coverage_histogram.R ${TISSUE_NAME}.tsv

```

Once complete, run this code to strip off the LSF job summary at the end of the output.
This generates a new file with '.final' appended to the file names, and removes the original file.
Any failed job iles will also be removed so that the job may be re-attempted.

```bash

cd ${ROOT_DIR}/Analysis/bedtools_coverage/${TISSUE_NAME}

GOOD=$(grep -lr "Successfully completed" .)
for i in $GOOD; do
    sed '/^$/Q' "$i" > "$i.final"
    rm $i
done &
BAD=$(grep -lr "Exit" .)
for i in $BAD; do
    rm $i
done &

```

## Step 4b-2: Get bedtools summaries

This files are later used to calculate the mtDNA copy number.
Note that this script currently has a fixed path to the chromosome lengths file: /lustre/scratch125/casm/team268im/al28/bed_ref/GRCh38_chrom_sizes_no_gaps.csv (can easily be replaced)
This step runs across all samples without submitting any jobs to LSF.

```bash

cd ${ROOT_DIR}/Analysis/bedtools_coverage/${TISSUE_NAME}

Rscript ${SCRIPT_FOLDER}mtDNA_bedtools_coverage_summary_grCh38.R

```

## Step 5a: Run verifyBamID

This uses the downloaded verifyBamID program in Andrew's lustre space: /lustre/scratch125/casm/team268im/fa8/119/VerifyBamId-farm5/VerifyBamID

As detailed below, it also uses multiple fixed path reference files in Andrew's lustre space for various paths used in the verifyBamID command: \
--UDPath /lustre/scratch125/casm/team268im/al28/bed_ref/striped_ref_files/1000g.phase3.100k.b38.vcf.gz.dat.UD
--BedPath /lustre/scratch125/casm/team268im/al28/bed_ref/striped_ref_files/1000g.phase3.100k.b38.vcf.gz.dat.bed
--MeanPath /lustre/scratch125/casm/team268im/al28/bed_ref/striped_ref_files/1000g.phase3.100k.b38.vcf.gz.dat.mu
--Reference /lustre/scratch125/casm/team268im/al28/bed_ref/GRCh38_full_analysis_set_plus_decoy_hla_genome.fa'

The wrapper will submit an LSF job for each sample in the 'normal' queue.

```bash
mkdir -p ${ROOT_DIR}/Analysis/verify_bam_id/${TISSUE_NAME}
cd ${ROOT_DIR}/Analysis/verify_bam_id/${TISSUE_NAME}

Rscript ${SCRIPT_FOLDER}mtDNA_verifyBamID_wrapper_GRCh38.R ${ROOT_DIR}/sample_lists/${TISSUE_NAME}.tsv

```

## Step 5b: Run haplocheck

First generate vcf files to be used by haplocheck using mutserve. Needs to be run from the 'sample_lists' folder.
This uses 'mutserve' which runs in java. \
The mutserve java file is here: /lustre/scratch125/casm/team268im/al28/software/mutserve.jar \
The relevant reference fasta path is here: /lustre/scratch125/casm/team268im/al28/software/rCRS.fasta


```bash
module load openjdk-17.0.8.1_1

mkdir -p ${ROOT_DIR}/Analysis/haplocheck/${TISSUE_NAME}/logs
cd ${ROOT_DIR}/sample_lists/

Rscript ${SCRIPT_FOLDER}mtDNA_haplocheck_mutserve_wrapper.R ${TISSUE_NAME}.tsv

```

Next, run the haplocheck program downloaded to Andrew's lustre space here: /lustre/scratch125/casm/team268im/al28/software/haplocheck
The wrapper will submit an LSF job for each sample in the 'normal' queue.

```bash
cd ${ROOT_DIR}/Analysis/haplocheck/${TISSUE_NAME}
Rscript ${SCRIPT_FOLDER}mtDNA_haplocheck_contamination_wrapper.R

```

## Step 6: Update the coverage table

As currently structured, there is one large master coverage/ sample info table located here: $ROOT_DIR/whole_genome_coverage_pileup_and_bedtools_annotated.csv

If you add a new study, you need to update this file with the new sample info.

This currently includes a manual steps of: \
(1) downloading the project info spreadsheets from CanApps
(2) convert these to .csv files with the blank lines at the top removed
(3) move these to a folder on the HPC farm for parsing into R and updating the spreadsheet

```bash
CANAPPS_DIR='/lustre/scratch126/casm/team154pc/ms56/Mitochondria_study/nonblood/canapps_info/'

Rscript ${SCRIPT_FOLDER}mtDNA_genome_cov_table_add_new_projects.R $TISSUE_NAME $CANAPPS_DIR/$TISSUE_NAME
```

Initially, a lot of this data is 'NA'. Therefore you need to fill it up from the available data using a separate script.
This iterates through each sample (each row in the dataframe) and adds: \
1. mtDNA coverage from the 'pileup' step, this is then converted into an implied mtDNA copy number by dividing by the SeqX value \
2. the summary output from the bedtools genomecov summaries \ 
3. the levels of sample contamination as assessed by the 'haplocheck' algorithm \ 
4. stats from the ascat summary (if run) \
5. an estimate of the mtDNA copy number from bedtools, using the average coverage across the autosomes

```bash
Rscript ${SCRIPT_FOLDER}mtDNA_genome_coverage_table_generation.R
```

## Step 7: Generate the global VAF heatmaps

This sets of a single farm job (using script mtDNA_global_vaf_heatmap_farm.R) that runs across donors.

```bash
mkdir -p ${ROOT_DIR}/Analysis/global_vaf_heatmaps/${TISSUE_NAME}
Rscript ${SCRIPT_FOLDER}mtDNA_global_vaf_heatmap_farm_wrapper.R $TISSUE_NAME
```

## Step 8: Assess for 'high VAF somatic' or germline mtDNA mutations

Extract the high VAF mutations
- classify as germline or 'high VAF heteroplasmic'

This is later used to assign a haplotype to each donor
(There is a separate script that runs this across all tissues sequentially)

Deopending on the number of donors/ samples, this script may take a few hours.

```bash
Rscript ${SCRIPT_FOLDER}mtDNA_germline_high_vaf_somatic_mutation_extraction_SINGLE_TISSUE.R $TISSUE_NAME
```

## Step 9: Assign likely haplotypes to each individual based on the germline mtDNA mutations

Assign each donor in the study to a specific haplotype.
This uses the reference info stored as a fixed path here: mtDNA_haplotype_markers_mitomap_2021-01-15.txt

```bash
Rscript ${SCRIPT_FOLDER}mtDNA_germline_haplotype_assignment_SINGLE_TISSUE.R $TISSUE_NAME
```

## Step 10: Run shearwater

This step runs the shearwater over all individuals from a specific individual.
The error profile for 'normal' has to be pre-set from a defined set of individuals (the 'normal' panel). This is ~102 individuals.

```bash
mkdir -p ${ROOT_DIR}/Analysis/shearwater_stringent_exclusion/${TISSUE_NAME}/logs
Rscript ${SCRIPT_FOLDER}mtDNA_shearwater_farm_wrapper_SINGLE_TISSUE.R $TISSUE_NAME
```


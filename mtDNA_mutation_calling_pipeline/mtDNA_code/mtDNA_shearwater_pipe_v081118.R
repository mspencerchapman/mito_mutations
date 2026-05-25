# Inigo Martincorena - Feb 2017
#
# Modified by Federico Abascal - April 2017, to:
#   -make the preparation of bam files faster
#   -realign BAMs around INDELs
#   -combine the preparation of bams and the running of shearwater into one single script
#   -concatenate all shearwater results into a single file of potential mutations
# Post-processing and filtering of the potential mutations is to be done by a different
#   script
#
# Modified by Andrew Lawson - July 2017, to:
#   -update environment links to local copies of files
#   -delete duplicated row naming 
#   -fixed bugs in wait2finish (job_table before assignment, pending jobs not considered, issue with escaping)
#   -changed version of RScript used for running shearwater to 3.4.0 (No GenomicRanges library installed for 3.3.0)
#   -TODO: Proper handling of wait2finish (check status of job is DONE or EXIT instead)
#   -changed shearwater error and output files to go to job_id file rather than bsub file
#   -check .out rather than .err file for success of Shearwater jobs
#   -changed output of BAM realignment to dedicated folders
#   -renamed shearwater output folder to be consistent with mutation annotation script
#   -samples_with_no_bam subsetting changed from - to !
#   -check whether realign is set to true and whether any bam files are marked for realignment before submitting jobs

# Example to run the pipeline including indel realignment: 
#    nohup Rscript shearwater_pipe.R my_analysis_name true fg_samples_or_tumor.file bg_samples_or_normal.file bait_file.bed output

# Input arguments:
# 1. Analysis name, eg PD12345_oeso_kk
# 2. Either true or false, depending on whether (GATK's) indel realignment has to be done before executing Shearwater
# 3. Name of file with the list of samples to call mutations in ("tumor"). Format of the file:
#    tab-delimited list of projectID and sampleID, no header.
# 4. Name of file with the list of samples to be used as reference ("normals"). Format of the file:
#    tab-delimited list of projectID and sampleID, no header.
# 5. Bait capture bed file [Default: "/lustre/scratch112/sanger/casm/cgp/im3/Projects/Oesophagus/Bait_set_ALL_EXONS_and_TP53.filt.nochrprefix.bed"]
# 6. Output directory [Default: current directory]


## 1. Environment

args <- commandArgs(TRUE)
analysis_name <- args[1]                #e.g. "PD12345_oeso_kk"
realign       <- as.character(args[2])  #e.g. true
samples1_file <- args[3]                #"tumor" list of samples
samples2_file <- args[4]                #"normal" list of samples
if (length(args)>4) { baits_bed <- args[5] } else { baits_bed <- "/lustre/scratch125/casm/team268im/al28/bladder_old/Bait_set_ALL_EXONS_and_TP53.filt.nochrprefix.bed" }
if (length(args)>5) { outdir    <- args[6] } else { outdir    <- getwd() }
numsegments_per_job <- 1
shearwater_fun <- "/lustre/scratch125/casm/team268im/al28/bladder_manuscript/al28_code/wrapper_mtDNA_shearwaterML_multiplebams_v081118.R"
realign_script <- "/lustre/scratch125/casm/team268im/al28/bladder_targeted/al28_code/realign_around_indels.pl"

cat("Analysis_name = ", analysis_name,"\n")
cat("realign       = ", realign,"\n")
cat("samples1_file = ", samples1_file,"\n")
cat("samples2_file = ", samples2_file,"\n")

# Libraries
library("GenomicRanges")
library("rtracklayer")

# Loading input files
samples1 <- read.table(samples1_file, header=0, sep="\t", stringsAsFactors=F)[,1:2]; samples1$ref = F
if (is.na(as.numeric(samples1[1,1]))) { samples1 = samples1[-1,]; samples1[,1] = as.numeric(samples1[,1]) }
samples2 <- read.table(samples2_file, header=0, sep="\t", stringsAsFactors=F)[,1:2]; samples2$ref = T
if (is.na(as.numeric(samples2[1,1]))) { samples2 = samples2[-1,]; samples2[,1] = as.numeric(samples2[,1]) }
rownames(samples1) <- samples1[,2];		
rownames(samples2) <- samples2[,2];
b = read.table(baits_bed, header=0, sep="\t", stringsAsFactors=F)
baits = GRanges(b[,1], IRanges(b[,2],b[,3]))


## 2. Creating a directory with all required symbolic links and launching GATK's indel realignment

setwd(outdir);

if (!file.exists("RDat")) {
    system("mkdir RDat")
}

if (file.exists("BAM_files")){
    cat("Warning: ", outdir,"/BAM_files directory already exists\n",sep="")
} else {
	system("mkdir BAM_files")
}
setwd("BAM_files")
sample_table = rbind(samples1, samples2)
sample_table <- unique(sample_table);
colnames(sample_table) = c("pid","sampleID","ref")
sample_table$path = sample_table$bai = sample_table$bam = NA
remove_these <- vector(mode="numeric",length=0); # <-- In case some sample / bam is not found
counter <- 1;

# For all the projects, get the list of "*sample.dupmarked.bam" files:
bamfilesF <- vector();
baifilesF <- vector();
pathsF    <- vector();
for(pid in unique(sample_table$pid)) {
	fullfiles        <- list.files(path = paste("/nfs/cancer_ref01/nst_links/live/",pid,"/",sep=""), pattern = "\\.bam$", recursive = TRUE, full.names = TRUE)
	fullfiles        <- grep("pindel",fullfiles,invert=T,value=T)
	fullfiles        <- grep("brm",   fullfiles,invert=T,value=T)
	files            <- sapply(strsplit(fullfiles, split="\\/"), tail, 1) 
	samples          <- sapply(strsplit(files,     split="\\."), head, 1) 
	names(files)     <- samples
	names(fullfiles) <- samples
	bamfilesF        <- c(bamfilesF, files    [which(samples %in% sample_table$sampleID)])
	pathsF           <- c(pathsF,    fullfiles[which(samples %in% sample_table$sampleID)])
	fullfiles        <- list.files(path = paste("/nfs/cancer_ref01/nst_links/live/",pid,"/",sep=""), pattern = "\\.bai$", recursive = TRUE, full.names = TRUE)
	fullfiles        <- grep("pindel",fullfiles,invert=T,value=T)
	fullfiles        <- grep("brm",   fullfiles,invert=T,value=T)
	files            <- sapply(strsplit(fullfiles, split="\\/"), tail, 1) 
	samples          <- sapply(strsplit(files, split="\\."), head, 1) 
	names(fullfiles) <- samples
	baifilesF        <- c(baifilesF, fullfiles[which(samples %in% sample_table$sampleID)])
}

for(i in 1:length(bamfilesF)) {
	sample_table[which(sample_table$sampleID == names(bamfilesF)[i]),"bam"] <- bamfilesF[i]
	sample_table[which(sample_table$sampleID == names(bamfilesF)[i]),"bai"] <- baifilesF[i]
	sample_table[which(sample_table$sampleID == names(bamfilesF)[i]),"path"] <- pathsF[i]
}

samples_with_no_bam <- rownames(sample_table[is.na(sample_table$bam),])
if(length(samples_with_no_bam) > 0) {
	cat("Removing ", length(samples_with_no_bam), " samples for which a bam file was not found - check project id!\n")
	for(j in c(1:length(samples_with_no_bam))) {	cat("    Missing bam for: ", samples_with_no_bam[j], "\n");		}
	sample_table <- sample_table[!(rownames(sample_table) %in% samples_with_no_bam),]
}

# If there is another shearwater_pipe.R process working with bams on this directory, wait til it finishes
if(file.exists("SEMAPHORE_INDEL_REALIGNMENT_ONGOING")) {
	cat("Semaphore directory found... Waiting until removed by a parallel process\n");
	cat("If there is not other process working on it, just remove that directory and rerun the script\n");
	cat("Do: rm -r SEMAPHORE_INDEL_REALIGNMENT_ONGOING\n");
	while(file.exists("SEMAPHORE_INDEL_REALIGNMENT_ONGOING")) {
		Sys.sleep(200)
	}
}

# Link BAMs to the current directory:
cat("Linking BAM and BAI files to the BAM_files directory...\n")
sample_table$realignable = NA
for(j in c(1:nrow(sample_table))) {
	if(file.exists(sample_table[j,"bam"])) {
		cat("Warning: The BAM file ", sample_table[j,"bam"], " already exists in the BAM_files directory. Not linking nor realigning it!\n")
		sample_table[j,"realignable"] = F
	} else {
		system(sprintf("ln -s %s .",sample_table[j,"path"]))
		system(sprintf("ln -s %s .",sample_table[j,"bai" ]))
		sample_table[j,"realignable"] = T
	}
}

# Function to wait for jobs to finish (im3)
wait2finish = function(bjob_vec) {
	#counter_jobs <- length(bjob_vec)
	Sys.sleep(30)
	if (!is.null(bjob_vec)) {
		cont = 1
		while (cont) {
			bjobs_table = system("bjobs -w | awk '{print $7}'", intern=T)
			if(identical(bjobs_table, character(0))) {
				cat("Number of jobs running: 0\n")
				break;
			}
			job_table <- read.table(text=bjobs_table, header=0, quote="\"",na.strings="-",fill=TRUE)
			running_jobs = unique(as.vector(job_table[-1,]))
			cont = any(bjob_vec %in% running_jobs)
			if(length(intersect(running_jobs,bjob_vec)) > 0) {
				#counter_jobs <- length(intersect(running_jobs,bjob_vec)) 
				cat("Number of jobs running: ", length(intersect(running_jobs,bjob_vec)) , "; Jobs finished: ",length(bjob_vec)-length(intersect(running_jobs,bjob_vec)) ,"\n")
			}
			Sys.sleep(30)
		}
	}
}


## 3. If realign is "true", then run indel realignment for all these bam files
#     Invoke the perl script and wait until all jobs are finished
#     After indel realignment is finished, replace the original files with the realigned ones
if(realign == "true" & any(sample_table$realignable)) {
	if (file.exists(paste(outdir,"/BAM_realigned",sep=""))){
    		cat("Warning: ", outdir,"/BAM_realigned directory already exists\n",sep="")
	} else {
        	system(paste("mkdir ",outdir,"/BAM_realigned",sep=""))
		system(paste("mkdir ",outdir,"/BAM_realigned/log",sep=""))
	}
	dir.create("SEMAPHORE_INDEL_REALIGNMENT_ONGOING") # Turn red the semaphore
	sample_table$job <- NA
	job_prefix = paste(sample(letters,5,TRUE),collapse="")
	job_count  = 1
	job_vec    = vector(length=nrow(sample_table[which(sample_table$realignable==TRUE),]))
	cat("Preparing to run ", length(job_vec), " indel realignment jobs...\n")
	for(sample in rownames(sample_table)) {
		if(sample_table[sample,"realignable"] == TRUE) { # So we only realign it if it's the first time we get the bam from nfs
			bam_file <- sample_table[sample,"bam"]
			job_id   <- paste(job_prefix,job_count,sep="")
			log_path <- paste(outdir,"/BAM_realigned/log/",job_id,sep="")
			cat("   Submitting job: ", sprintf("bsub -M 40000 -R \"select[mem>40000] rusage[mem=40000]\" -J %s -e %s.err -o %s.err %s %s \n", 
					job_id, log_path, log_path, realign_script, bam_file))
			system(sprintf("bsub -M 40000 -R \"select[mem>40000] rusage[mem=40000]\" -J %s -e %s.err -o %s.err %s %s", 
					job_id, log_path, log_path, realign_script, bam_file))
			sample_table$job <- job_id
			job_vec[job_count] <- job_id
			job_count <- job_count + 1
		}
	}
	
	
	save.image(paste(outdir,"/RDat/",analysis_name,"_realignment.Rdat",sep=""))
	# Wait for the jobs to finish:
	cat("Be patient, indel realignment can take long (a few hours)\n");
	wait2finish(job_vec)
	
	# Check all works finished properly:
	failed_jobs <- vector()
	for(job in job_vec) {
		out_lines <- system(paste("grep Success ",outdir,"/BAM_realigned/log/",job,".err",sep=""),intern=T)
		if(length(out_lines)==0) {
			cat("Job ", job, " did not finish properly\n");
			out_lines <- system(paste("grep -i exit ",outdir,"/BAM_realigned/log/",job,".err",sep=""),intern=T)
			cat("   Exit message: ", out_lines, "\n")
			failed_jobs[job] <- job
		}
	}
	if(length(failed_jobs) > 0) {	
		cat("There were ", length(failed_jobs), " failed jobs\n")	
	} else 						{	
		cat("All indel realignment jobs finished properly\n")		
	}
	
	# Move all the realigned bams
	old_bam_files <- sample_table$bam
	old_bai_files <- paste(sample_table$bam,".bai",sep="")
	new_bam_files <- paste(sample_table$bam,".out.bam",sep="")
	new_bai_files <- paste(sample_table$bam,".out.bai",sep="")
	file.rename(new_bam_files,paste(outdir,"/BAM_realigned/",old_bam_files,sep=""))
	file.rename(new_bai_files,paste(outdir,"/BAM_realigned/",old_bai_files,sep=""))
	system("rm -r SEMAPHORE_INDEL_REALIGNMENT_ONGOING") # Turn green the semaphore
}

# Prepare lists of bam files to run shearwater
setwd(outdir);
if (file.exists(analysis_name)){
    cat("Warning: ", outdir,"/",analysis_name," directory already exists. Content will be overwritten\n",sep="")
} else {
	system(sprintf("mkdir %s",analysis_name))
}
setwd(analysis_name);

if (realign == "true") {
    write.table(paste("../BAM_realigned/",sample_table[which(sample_table$ref==FALSE),"bam"],"\n",collapse="",sep=""),"tumor.list.tsv", col.names=F,row.names=F,sep="\t",quote=F)
    write.table(paste("../BAM_realigned/",sample_table[which(sample_table$ref==TRUE),"bam"],"\n",collapse="",sep=""),"normal.list.tsv", col.names=F,row.names=F,sep="\t",quote=F)
} else {
    write.table(paste("../BAM_files/",sample_table[which(sample_table$ref==FALSE),"bam"],"\n",collapse="",sep=""),"tumor.list.tsv", col.names=F,row.names=F,sep="\t",quote=F)
    write.table(paste("../BAM_files/",sample_table[which(sample_table$ref==TRUE), "bam"],"\n",collapse="",sep=""),"normal.list.tsv",col.names=F,row.names=F,sep="\t",quote=F)
}

## 4. Running Shearwater

entry_start = seq(from=1, to=length(baits), by=numsegments_per_job)
entry_end = pmin(entry_start+numsegments_per_job-1, length(baits))
cmds = NULL
job_prefix = paste(sample(letters,5,TRUE),collapse="")
job_count  = 1
job_vec    = vector(length=length(entry_start))
dir.create("log")
cat("Preparing to run ", length(job_vec), " shearwater jobs...\n")
for (h in 1:length(entry_start)) {
	job_id   <- paste(job_prefix,job_count,sep="")
	log_path <- paste("log/",job_id,sep="")
	cmds[length(cmds)+1] = sprintf("bsub -J %s -q normal -R 'select[mem>=5000] rusage[mem=5000]' -M5000 -e %s.err -o %s.out /software/R-3.5.3/bin/Rscript %s %s %0.0f %0.0f %s %s ./shearwater_temp", job_id, log_path, log_path, shearwater_fun, baits_bed, entry_start[h], entry_end[h], "tumor.list.tsv", "normal.list.tsv")
	cat("   Submitting job: ", cmds[length(cmds)], "\n");

	system(cmds[length(cmds)])
	job_vec[job_count] <- job_id
	job_count <- job_count + 1
}
save.image(paste(outdir,"/RDat/",analysis_name,"_shearwater.Rdat",sep=""))
# Wait for the jobs to finish:
cat("Be patient, shearwater can take long\n");
wait2finish(job_vec)

# Check all works finished properly:
failed_jobs <- vector()
for(job in job_vec) {
	out_lines <- system(paste("grep Success log/",job,".out",sep=""),intern=T)
	if(length(out_lines)==0) {
		cat("Job ", job, " did not finish properly\n");
		out_lines <- system(paste("grep -i exit log/",job,".out",sep=""),intern=T)
		cat("   Exit message: ", out_lines, "\n")
		failed_jobs[job] <- job
	}
}
if(length(failed_jobs) > 0) {	
	cat("There were ", length(failed_jobs), " failed jobs\n")	
} else 						{	
	cat("All shearwater jobs finished properly\n")		
}

save.image(paste(outdir,"/RDat/",analysis_name,"_pipeline_end.Rdat",sep=""))

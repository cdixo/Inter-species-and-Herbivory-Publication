#!/bin/bash

#SBATCH --job-name=RNAseq_gunzip_to_fastq
#SBATCH --time=1:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=2
#SBATCH --mem=1gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444


cd $SLURM_SUBMIT_DIR #this allows you to say "run on files in this folder from which I am sending the job"

#cd /fs/project/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir
#you can also use the above to change the directory to the directory in which you want the following command to be done in if you would rather that

find . -mindepth 1 -type f -print -exec mv {} . \; #this pulls out all of the gunzipped folders from inside their standard folders named by sample.  It then pulls them out to the parent folder location (the location you are currently in)
# without the above step pulling out all gunzipped folders to the forefront, the below command will not work

gunzip *fq.gz
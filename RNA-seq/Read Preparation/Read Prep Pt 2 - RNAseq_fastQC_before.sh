#!/bin/bash

#SBATCH --job-name=RNAseq_clean_fastqc_before
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444


cd $SLURM_SUBMIT_DIR
#cd /fs/project/PAS1444/labrusca_RNAseq_rawdata/211013_Fiorella_GSL-FCC-2418

# -q = quiet (surpress printing updates to the console screen)
# -t = sprecifies the number of files it can simultaneously work on; allocates 250MB to each so account for that in the RAM request 
# -f = specifies what format you are inputing (in this case we say fastq)
/users/PAS1444/li10917/miniconda2/envs/rnaseq/bin/fastqc -t 4 -f fastq *.fq

mkdir fastQC_output_before

mv *fastqc.html fastQC_output_before
mv *fastqc.zip fastQC_output_before


#this took 2.5hr for 112 files on 1/27/22

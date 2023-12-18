#!/bin/bash

#SBATCH --job-name=ATR_box1.1_try4rev.sh_9-19-22
#SBATCH --time=40:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=60gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444

#move to where we need to work in
cd /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/GREM4transcriptevidence/211013_Fiorella_GSL-FCC-2418

#concatenate (combine) all of the rev reads together
zcat *_R2_*.fastq.gz | gzip > GREM4RevTranscripts.fq.gz

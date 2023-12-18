#!/bin/bash

#SBATCH --job-name=FA_funannotate_annotate_try6.1_submission1_12-2-22
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=100gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444

singularity exec /fs/project/PAS1444/software/funannotate/funannotate_mask.sif \
bash /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/functionalannotation/FA_funannotate_annotate_try6.1.sh

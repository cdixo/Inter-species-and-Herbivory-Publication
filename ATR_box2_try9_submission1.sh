#!/bin/bash

#SBATCH --job-name=ATR_box2_try9.sh_9-25-22
#SBATCH --time=164:00:00
#SBATCH --nodes=1
#SBATCH --mem=117gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444

singularity exec /fs/project/PAS1444/software/funannotate/funannotate_mask.sif \
bash /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/ATR_box2_try9.sh
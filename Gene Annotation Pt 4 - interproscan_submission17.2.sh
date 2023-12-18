#!/bin/bash

#SBATCH --job-name=interproscan_submission17.2_10-6-22
#SBATCH --time=105:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=100gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444


module load java/12.0.2

#Running InterProScan, program installed locally, upon the funannotate update AA.fasta output file to provide it with functional annotation (and hopefully AED scores for each specific annotation so I can run an AED test on it)
/fs/project/PAS1444/GeneAnnoWorkDir/interproscan2_17_22/interproscan-5.57-90.0/interproscan.sh \
-cpu 23 \
-T /fs/scratch/PAS1444/tempinterproscanrunfiles \
-appl CDD,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,SFLD,SMART,SUPERFAMILY,TIGRFAM \
-iprlookup \
-goterms \
-t p \
-f TSV, XML, JSON, GFF3 \
-i /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/functionalannotation/practice/9_25_22_ATR_box2/update_results/Vitis_labrusca.proteins.fa \
-d /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/functionalannotation/FA_interproscan_submission17.2_output
#-cpu = how many cores to run with
#-t = input is nucleotides (n); input is amino acids (p)			
#-T = listing the temporary directory for the run
#-o = output file name [may only provide a '-o' or a '-d', not both]
#-d = output directory name [may only provide a '-o' or a '-d', not both]

#!/bin/bash

#SBATCH --job-name=CoCo_counts_try10_PN.sh_3-1-23
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444


# Load the module you intend run (need for samtools) (with this specific version of conda which is very necessary) 
# samtools is used within the CoCo program
module load python/3.6-conda5.2

# Activate the RNAseq environment which contains samtools
source activate /users/PAS1444/cullendixon/.conda/envs/RNAseq



# Set the placeholder variables
coco_correct_gtf_file=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo_PN/Vitis_vinifera.PN40024.v4.56.cococorrectannoout.gtf   # genome annotation file which is the output from the 'coco correct_annotation' command
bam_input_files=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/STAR_PN/STAR_align/bam/ 																						# input file directory (the real deal .bam files) in .bam format
output_location=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo_PN/																														# output location (just location, not the name of the output file, that comes in the for loop)


# Run the 'coco correct_count' command to get the count output file for use in DESeq2 DEG identication and downstream analysis
# Inserted in a 'for' loop
for file in "$bam_input_files"*.bam
do
  /users/PAS1444/cullendixon/software/CoCo/coco/bin/coco correct_count \
	-c both \
	-p \
	-t 22 \
	"$coco_correct_gtf_file" \
	${file} \
	${file}_cococountout.count
done

#This script actually places the output files at location '/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/STAR/STAR_align/bam'
#You must therefore mv all of the output files to location '/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo/' by use of the command line.


# Options:
# -c = Decide whether to consider only uniquely mapped reads (uniqueOnly) or both uniquely and multimapped reads
#                       (both) to produce read count values for genes. Default: both
# -s = Strandedness. Default: both
# -p = Use this option if your dataset is a paired-end read dataset. Used by featureCounts.
# -t = Thread.  How many threads the program can run on. Default: 1
# -r = Use this option to dictate that only the raw count output file should be produced and not also the CPM & TPM files


#---------------------------------------------------------------------------------------------
# Location of Actual .bam Files
#/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/STAR/STAR_align/bam

# Location of Copy of Three of the Actual .bam Files to Use for Practice
#/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo/prac
#!/bin/bash

#SBATCH --job-name=STAR_align_all_try4_PN.sh_2-28-23
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=43
#SBATCH --mem=160gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444


# Load the module you intend to run's preresquisite module
module load gnu/10.3.0

# Load the module you actually intend to run
module load star/2.7.9a

# Manually create these directories
# 'STAR_align/bam' at location '/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/STAR'


# Set the placeholder variables
sample_list=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/sampleswFILLER.txt                 # the txt file which lists which samples (transcriptome samples) which the above command should look for to then feed into STAR; this file contains the term 'FILLER' at the end - see elec. lab NB on 2/28/23 for explanation
sample=($(cat "$sample_list"))                                                                   # dictates the samples (transcriptome samples) which STAR will align, but needs an input txt file to do so
indexed_ref_genome_dir=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/STAR_PN/STAR_index      # input directory of the indexed reference genome from the last step ('STAR_index_try2.sh')
ref_gtf=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/Vitis_vinifera.PN40024.v4.56.gtf       # genome annotation file
output_dir=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/STAR_PN/STAR_align/bam              # output directory which will contain a finalized, sorted .bam .bam file


# Run STAR to align the transcripts to the index reference genome, but on all files this time, therefore, use this loop
for sample in "${sample[@]}"; do
  R1=$(ls *"$sample"*_1.clean.fq)
  R2=$(ls *"$sample"*_2.clean.fq)
  slurm_log=slurm-align-"$sample"-%j.out
  echo "Submitting star_align.sh for $sample"
  STAR --genomeDir "$indexed_ref_genome_dir" \
    --readFilesIn "$R1" "$R2" \
    --sjdbGTFfeatureExon "$ref_gtf" \
    --runThreadN "$SLURM_CPUS_ON_NODE" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$output_dir/$sample"_   
done


# R1=$(ls *"$sample"*_1.clean.fq)     # dictating the fwd read demarcater for the FASTQ files
# R2=$(ls *"$sample"*_2.clean.fq)     # dictating the rev read demarcater for the FASTQ files

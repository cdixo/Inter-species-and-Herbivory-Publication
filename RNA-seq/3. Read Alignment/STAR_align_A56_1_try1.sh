#!/bin/bash

#SBATCH --job-name=STAR_align_A56_1_try1.sh_11-2-22
#SBATCH --time=0:10:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=76gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444

# Load the module you intend to run's preresquisite module
module load gnu/10.3.0

# Load the module you actually intend to run
module load star/2.7.9a

# Manually create these directories
# 'STAR_align_one/bam' at location '/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/STAR'


# Set the placeholder variables
sample=A56_1
indexed_ref_genome_dir=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/STAR/STAR_index  # input directory of the indexed reference genome from the last step ('STAR_index_try2.sh')
ref_gtf=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/Vitis_labrusca_annoV2.gtf       # genome annotation file
output_dir=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/STAR/STAR_align_A56_1/bam      # output directory which will contain a finalized, sorted .bam .bam file
R1=$(ls *"$sample"*_1.clean.fq)                                                           # dictating the fwd read demarcater for the FASTQ files
R2=$(ls *"$sample"*_2.clean.fq)                                                           # dictating the rev read demarcater for the FASTQ files


# Run STAR to align the transcripts to the index reference genome
STAR --genomeDir "$indexed_ref_genome_dir" \
   --readFilesIn "$R1" "$R2" \
   --sjdbGTFfeatureExon "$ref_gtf" \
   --runThreadN "$SLURM_CPUS_ON_NODE" \
   --outSAMtype BAM SortedByCoordinate \
   --outFileNamePrefix "$output_dir/$sample"_ \

# --genomeDir = dictate the location of your indexed genome for input
# --readFilesIn = dictate the location (names) of your fwd and rev transcriptome reads
# --runThreadN = dictates the amount of computing resources it can use to run with
# --outSAMtype = dictates the output file format; we are using BAM here; also asking it to sort it by coordinate so we can use it for downstream analysis
# --outFileNamePrefix = allows us to dictate the output files' prefixes on their output names
#
# --sjdbGTFfeatureExon - dictates that you will list the .gtf file herein, and then it will only use the lines that read 'exon' in the 3rd column
# --sjdbGTFtagExonParentTranscript = dictates that you will list the .gtf file herein, and then it will only use the lines that contain 'transcript_id' (which is all of them)
#     Rachel is using this as her option to feed in the gene annotation (.gtf file),
#     The manual says says "string: GTF attribute name for parent transcript ID (default "transcript id"
#     works for GTF files)" - althought I am not sure exactly what the means.  Another passage says this 
#     "In addition to the aforementioned options, for GFF3 formatted annotations you need to use
#     --sjdbGTFtagExonParentTranscript Parent. In general, for --sjdbGTFfile files STAR only processes lines 
#     which have --sjdbGTFfeatureExon (=exon by default) in the 3rd field (column). The exons are assigned to
#     the transcripts using parent-child relationship defined by the --sjdbGTFtagExonParentTranscript 
#     (=transcript id by default) GTF/GFF attribute."
#     While I was originally going to use this one now I am not.
# --sjdbGTFfile = dictate your genome annotation file and then follow it up with another option to dictate how you want it to pull the rows for use (see manual pg 27)
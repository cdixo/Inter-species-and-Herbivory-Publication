#!/bin/bash

#SBATCH --job-name=STAR_index_try2_PN1_2-28-23
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=76gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444



# Load the module
module load gnu/10.3.0
module load star/2.7.9a

# Set the placeholder variables
input_fasta=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/Vitis_vinifera.PN40024.v4.dna.toplevel.fa      # reference genome FASTA file (input)
ref_gtf=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/Vitis_vinifera.PN40024.v4.56.gtf	 				 # genome annotation file
output_dir=/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/STAR_PN/STAR_index      						 # Output dir for the indexed genome


# Run STAR to index the reference genome
STAR --runMode genomeGenerate \
	--genomeFastaFiles "$input_fasta" \
	--sjdbGTFtagExonParentTranscript "$ref_gtf" \
	--runThreadN "$SLURM_CPUS_PER_TASK" \
	--sjdbOverhang maxReadLength-1 \
	--genomeSAindexNbases 13 \
	--genomeDir "$output_dir" \
	

# --runMode = Tells STAR to index the genome (tells it what task we are trying to accomplish?)
# --runThreadN = How many threads it can use to accomplish its job
# --genomeDir = dictates the output directory (where you want it to place the finalized indexed genome)
# --genomeFastaFiles = dictates the input file (the genome (a fasta file))
# --sjdbGTFtagExonParentTranscript
# --genomeSAindexNbases = Used to apply the correct scale to the analysis, see manual.
# --sjdbOverhang = "Specifies the length of the genomic sequence around the annotated junction 
#					to be used in constructing the splice junctions database. Ideally, this length 
#					should be equal to the ReadLength-1, where ReadLength is the length of the reads. 
#					For instance, for Illumina 2x100b paired-end reads, the ideal value is 100-1=99. In 
#					case of reads of varying length, the ideal value is max(ReadLength)-1. In most cases, 
#					the default value of 100 will work as well as the ideal value." - Manual
#
#					In looking at 5 qc reports, the reads all ran from 36-150bps in size.  Rachel had previously
#					been using 13 only, but I don't understand why she would be doing that based
#					on what I am reading.  I am going to mention this to her and see what she was
#					thinking with her 13 she used. - Cullen 10/31/22.

#### Initial Notes - 
# '~' = /users/PAS1444/cullendixon
# All output files will be placed at the location from which the job was submitted of the submission file

#### Change Names to Match Formatting - 
# First, you have to change 'labrusca.fasta' to your genome name; labrusca is prefix of database
#system "~/miniconda3/envs/RNAseq/bin/hisat2-build -p 2 /fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/labrusca.V1.fasta labrusca"; 
# 	I turned this off because I don't think I really need to change this name over any longer as of 11/1/22 - CWD


#### .gff3 File to .gtf File Conversion - 
# When quantifying gene expression, stringtie needs a .gtf file - if you do not have .gft but rather a more conventional .gff3 file, you need to convert the .gff file to a .gft file first

#Turn on the python module
module load python/3.7-2019.10

#Turn on the RNAseq environment within which 'gffread' resides (activate the environment)
source activate RNAseq

#Run the command which creates a .gtf file from the .gff file
gffread -T /fs/ess/PAS1444/GeneAnnoWorkDir/funannotate/funannotate_output/9_25_22_ATR_box2/update_results/Vitis_labrusca.gff3 -o /fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/Vitis_labrusca_annoV2.gtf


#Remember to now back on out of the RNAseq environment and miniconda and python
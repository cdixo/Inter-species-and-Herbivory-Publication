######    Run by hand in the command line of OSC    ######

## This time we are now trying to use an altogether different program called AGAT
## This program should give us what we want - transcript, CDS, exon, AND gene feature lines
## To run AGAT, since the conda installer didn't work, I had to download a .sif file via
## apptainer.  This apptainer .sif file can be found here - 
##					/users/PAS1444/cullendixon/software/AGAT_via_apptainer/agat_0.8.0--pl5262hdfd78af_0.sif


#### .gff3 File to .gtf File Conversion - 
# When quantifying gene expression, stringtie needs a .gtf file - if you do not have .gtf but rather a more conventional .gff3 file, you need to convert the .gff file to a .gtf file first

#Run the singularity (apptainer) .sif file that contains the AGAT program
singularity run /users/PAS1444/cullendixon/software/AGAT_via_apptainer/agat_0.8.0--pl5262hdfd78af_0.sif


#Run the command which creates a .gtf file from the .gff file using the .sif file of AGAT
agat_convert_sp_gff2gtf.pl \
	--in /fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/Vitis_labrusca.gff3 \
	--out /fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo/Vitis_labrusca_annoV2.AGAT_try3.10.gtf

# Options:
# --in = input .gff file
# --out = output .gtf file


#Exit the singularity (apptainer) session
exit


#See what the first 50 lines look like of the file to see if they contain the word 'gene' now in the 3rd column
head -n 50 \
/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo/Vitis_labrusca_annoV2.AGAT_try3.10.gtf \
> /fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo/Vitis_labrusca_annoV2.AGAT_try3.10.gtf.50lines





###  Did AGAT's conversion script capture all of the lines that contained 'gene'?  ###

# Query the .gff file for how many times the word 'gene' occurs in its 3rd column - the feature column
awk -F '\t' '$3=="gene"' /fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/Vitis_labrusca.gff3 | wc -l
# Result = 37,444

# Query the .gtf file for how many times the word 'gene' occurs in its 3rd column - the feature column
awk -F '\t' '$3=="gene"' /fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo/Vitis_labrusca_annoV2.AGAT_try3.10.gtf | wc -l
# Result = 37,444





###  Rename the new .gtf file to something more reasonable  ###

cp \
/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo/Vitis_labrusca_annoV2.AGAT_try3.10.gtf \
> /fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo/Vitis_labrusca_annoV2.plusgenes.gtf


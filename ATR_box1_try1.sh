#Zach's original code \/
#zcat <SET1_READ_1>.fq.gz <SETn_READ_1>.fq.gz | gzip > <OME>_1.fq.gz

#My code \/ --------------------------------------------------------
#(I am running all of this in the command line)

#move to where we need to work in
cd /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/GREM4transcriptevidence/211013_Fiorella_GSL-FCC-2418

#pull all of the contents of the various folders at this location to the location we are at now (the above location)
find . -mindepth 1 -type f -print -exec mv {} . \;

#see the next commands to run (which need to be run with a job script as they require much memory) as 'ATR_box1.1_try3.sh'.
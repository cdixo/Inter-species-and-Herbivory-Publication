source /fs/project/PAS1444/software/funannotate/source.sh

funannotate update -f /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/labrusca.V1.masked.fasta \
-g /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/labrusca.V1.gff \
--species "Vitis labrusca" --cpus 28 \
--max_intronlen 5000 \
-l /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/GREM4transcriptevidence/211013_Fiorella_GSL-FCC-2418/GREM4FwdTranscripts.fq.gz -r /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/GREM4transcriptevidence/211013_Fiorella_GSL-FCC-2418/GREM4RevTranscripts.fq.gz \
-o /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/funannotate_output/9_25_22_ATR_box2
#the '--species' line may need to match the output when you type in the command 'funannotate -species' for the V. labrusca work that was previously done, but perhaps not, idk, just a note

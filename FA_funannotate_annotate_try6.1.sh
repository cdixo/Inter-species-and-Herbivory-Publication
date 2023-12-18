
source /fs/project/PAS1444/software/funannotate/source.sh

funannotate annotate \
--cpus 26 \
--iprscan /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/functionalannotation/FA_interproscan_submission17.2_output/Vitis_labrusca.proteins.fa.xml \
--input /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/functionalannotation/practice/9_25_22_ATR_box2 \
--fasta /fs/project/PAS1444/GeneAnnoWorkDir/funannotate/labrusca.V1.masked.fasta \
--species "Vitis labrusca" \
--busco_db /fs/project/PAS1444/databases/funannotate/embryophyta \
--out /fs/ess/PAS1444/GeneAnnoWorkDir/funannotate/functionalannotation/12_2_22_attempts

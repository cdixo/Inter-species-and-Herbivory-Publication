for file in Uniq_genes_GREM_GeneList.txt
do
   ../DEGs/seqtk/seqtk subseq Vitis_labrusca.cds-transcripts_noT.fa ${file} > ${file}.cdsseqs.fa
done

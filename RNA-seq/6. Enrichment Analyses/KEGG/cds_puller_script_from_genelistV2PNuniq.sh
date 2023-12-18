for file in Uniq_genes_PN_GeneList.txt
do
   ../DEGs/seqtk/seqtk subseq PN40024.v4.1.REF.cds_noT.fa ${file} > ${file}.cdsseqs.fa
done

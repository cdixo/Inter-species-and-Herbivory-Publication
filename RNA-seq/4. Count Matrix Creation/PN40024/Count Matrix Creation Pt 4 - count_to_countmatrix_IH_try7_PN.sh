###   Creating a Count Matrix from All of the Individual .count Files output by 'coco correct_count'   ###

# All of the following was run at location '/fs/ess/PAS1444/IH_Timecourse_V1_RNA-seq/work_dir/CoCo_PN' on 3/2/23 - CWD 3/2/23
# Sort the file by the first column (gene name) which is important because for some odd reason all of the count files don't have their
#       genes sorted in the first column in the same order.
for file in *.count
do
  (tail -n +3 ${file} | sort -k 1 ) > ${file}.sorted
done
#          ^ skip the first two lines (its a reverse pyschology head (tail?) command)

# Add a header to the column over the count column for use later when we cut them out
for file in *.count.sorted
do
  awk 'NR==1 {print "\t" "\t" FILENAME "\t" "\t"}; 1' ${file} > temp && mv temp ${file}.headed # the "; 1" portion ensures that it just doesn't rewrite the entire file as just the filename.  Without this, this it deletes all text within and just prints the tab + file name.
done

# Combine the various sorted count files into one large file
paste *.sorted.headed | # this says take all of the .headed files just made above and combine them into one file, moving from left to right, with each file being placed in their own successive columns

# Cut out the columns we want in the count matrix file ('gene name/id' & 'count')
#       IMPORTANT: THIS \/ IS THE ORDER OF THE INFORMATION IN THE .count FILES:
#             gene_id gene_name count cpm tpm
cut -f1,3,8,13,18,23,28,33,38,43,48,53,58,63,68,73,78,83,88,93,98,103,108,113,118,123,128,133,138,143,148,153,158,163,168,173,178,183,188,193,198,203,208,213,218,223,228,233,238,243,248,253,258,263,268,273,278 \
    > countmatrix_IH_PN.txt

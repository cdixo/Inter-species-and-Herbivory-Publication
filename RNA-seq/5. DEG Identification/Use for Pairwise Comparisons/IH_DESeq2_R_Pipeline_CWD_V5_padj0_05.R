###### Environmental Setup Steps #####

library(DESeq2)
library(gplots)
library(amap)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)
library(ggrepel)
library(data.table)
library(writexl)
library(readr)
library(dplyr)
library(ggplot2)

#   On/Off Switches   #
#       - Analysis 1
#       - Analysis 2


##Set you working directory:
setwd("C:/Users/dixon.778/OneDrive - The Ohio State University/School/College - Graduate School/19-20/Gschwend Lab/AAA   Science/Coding/Aim 3 IH Timecourse V1/Pairwise Comparisons")


##Feed in your data

#GREM4
#data <- read.table("countmatrix_IH.txt", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")

#PN
#data <- read.table("countmatrix_IH_PN.txt", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")


#Inter-species Comparisons
#data <- read.table("countmatrix_IH_collinearity_GREM4_PN_sorted.txt", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")
data <- read.table("countmatrix_IH_orthology_GREM4_PN_sorted.txt", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")



## Print to screen the number of lines in the dataset and the number of columns
dim(data)


##Set Your Parameters of the Run
lfcvalue <- 1  #this creates a variable which will set the log2foldchange value where inserting a '0'=lfc of 1, '1'=lfc of 2 (what we want)
padjvalue <- 0.05
pvaluevalue <- 0.05







###### Analysis 1: All Samples vs. All Conditions vs. All Species Analysis ######

#### Preparatory Work ####

if(FALSE) {

## read in group information in for population; must be done before any of the following analysis
sample <- read.table("Condition.Matrix.IH.nodashorrep.txt", header=T, row.names=1, com='',
                     quote='', check.names=F, sep="\t", colClasses="factor")
sample <- sample[match(colnames(data), rownames(sample)),, drop=F]
sample_rowname <- rownames(sample)


## Round the data using the round function to get it into integers and no decimal places (since the data may or many not be in integer format if you use CoCo for example)
##        the 'DESeqDataSetFromMatrix' command needs integer data to run
data <- round(data)


## Dictate the format of the incoming data into a new variable for use below
# 'design = ~ ' should be EXACTLY the name of the condition matrix condition column; you can set this to whatever you want, like treatment or state, but must be written exactly the same
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data, colData = sample,  design = ~ condition)




#### Run DESeq2 ####


## Run DESeq2 from the above-formatted dataset to produce DESeq2 output (dds)
dds <- DESeq(ddsFullCountTable)


## Create an Output File of the DESeq2 Results

# Normalize the DESeq2 output data (dds) and convert to a file-printable format
normalized_counts <- counts(dds, normalized=TRUE)

# DO NOT normalize the DESeq2 output data (dds) and convert to a file-printable format
unnormalized_counts <- counts(dds, normalized=FALSE)

# Sort the DESeq2 Output by Median Absolute Deviation on both objects
#       Normalized
normalized_counts_mad <- apply(  #apply takes a data frame or matrix and provides output in the form of a vector, list, or array
  normalized_counts,  #object to work on
  1, #a value between 1 and 2 to define where to apply the function where 1=manipulate rows and 2=manipulate columns and c(1,2)=manipulate rows and columns
  mad #a function to apply; mad=median absolute deviation; mean; median; sum; min; max
)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

#       Un-Normalized
unnormalized_counts_mad <- apply(  #apply takes a data frame or matrix and provides output in the form of a vector, list, or array
  unnormalized_counts,  #object to work on
  1, #a value between 1 and 2 to define where to apply the function where 1=manipulate rows and 2=manipulate columns and c(1,2)=manipulate rows and columns
  mad #a function to apply; mad=median absolute deviation; mean; median; sum; min; max
)
unnormalized_counts <- unnormalized_counts[order(unnormalized_counts_mad, decreasing=T), ]

# Create a output excel file for both Normalized and Un-Normalized (raw) DESeq2 Output
write_xlsx(as.data.frame(normalized_counts), paste("All_sample_comparison.coco.count.matrix.DESeq2.Normalized.xlsx", sep=""))
write_xlsx(as.data.frame(unnormalized_counts), paste("All_sample_comparison.coco.count.matrix.DESeq2.Unnormalized.xlsx", sep=""))




#### Create Graphics ####


## Conduct an rlog Transformation of the DESeq2 Output Data (dds)
##      This is the data which I will use for plotting in graphics
rlogdds <- rlog(dds,     #this functions rlog transforms the DESeq2 output data (dds)
                blind=FALSE   #this should be left as 'FALSE' if planning on using for downstream analysis
)
#The above step will take quite a while to run FYI


## Set the Coloration Pattern of Future Graphics
hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(100)



### Create Pearson Correlation Graphic ###

#Reformat the 'rlogdds' created just above into a matrix which is readable to the 'pearson_cor' command
rlogddsMatrix <- assay(rlogdds)
rlogddsMatrix <- rlogddsMatrix[order(normalized_counts_mad, decreasing=T), ]

#Create the Pearson Correlation
pearson_cor <- as.matrix(cor(rlogddsMatrix, method="pearson"))

#Create a File to Output the Pearson Correlation Data
write_xlsx(as.data.frame(pearson_cor), "All_sample_comparison.coco.count.matrix.DESeq2.rlog.pearsonsmatrix.xlsx")

#Run 'hcluster' to Obtain Hierarchical Clustering Data for the Graphic
hc <- hcluster(t(rlogddsMatrix), method="pearson")

#Create the Pearson Correlation Graphic
#     The coloration pattern below is based on the exact order of the 'pearson_cor' object
#     and was created by hand to match the order by using the condition matrix to determine
#     which samples were to get which color (condition)
pdf("HierarchialClusteringCorrelationMap.pdf", pointsize=10, wi = 14, he = 14)
par(cex.main=2.0) #this sets the font size for the below title
heatmap.2(pearson_cor, 
          Rowv = as.dendrogram(hc), 
          key = T, #determines if you want the corner key graph to be displayed
          keysize = 1, #sets the size of the key by a ratio
          #lhei = 14, #sets the height by physical number; never got to work
          #lwid = 14, #sets the height by physical number; never got to work
          key.par = list(cex=1.1), #sets the font size of the key
          symm = T, 
          trace = "none",
          col = hmcol, 
          margins = c(15,15), #11,11
          main = "Pearson Correlation of All Samples", #title here
          cexRow = 1.2,  #sets font size for rows
          cexCol = 1.2,  #sets font size for columns
          ColSideColors = c(    # grouping column-variables into different
            rep("lightslateblue", 1),   # GREM_C_30 A10_2_cococountout.count      (color, total number of each variable)
            rep("maroon", 1),    # GREM_H_1 A11_1_cococountout.count
            rep("sienna3", 1),    # GREM_C_1  A12_1_cococountout.count
            rep("cadetblue1", 1), # GREM_H_0  A1_2_cococountout.count
            rep("chartreuse1", 1), # GREM_H_4 A13_3_cococountout.count
            rep("coral1", 1), # GREM_C_4  A14_1_cococountout.count
            rep("cadetblue1", 1), # GREM_H_0  A15_1_cococountout.count
            rep("khaki1", 1), # GREM_H_30 A16_1_cococountout.count
            rep("lightslateblue", 1),   # GREM_C_30 A17_2_cococountout.count
            rep("maroon", 1),    # GREM_H_1 A18_1_cococountout.count
            rep("sienna3", 1),    # GREM_C_1  A19_1_cococountout.count	
            rep("chartreuse1", 1), # GREM_H_4 A20_1_cococountout.count	
            rep("coral1", 1), # GREM_C_4  A21_1_cococountout.count	
            rep("khaki1", 1), # GREM_H_30 A2_1_cococountout.count	
            rep("cadetblue1", 1), # GREM_H_0  A22_1_cococountout.count	
            rep("khaki1", 1), # GREM_H_30 A23_1_cococountout.count	
            rep("lightslateblue", 1),   # GREM_C_30 A24_2_cococountout.count	
            rep("maroon", 1),    # GREM_H_1 A25_1_cococountout.count	
            rep("sienna3", 1),    # GREM_C_1  A26_1_cococountout.count	
            rep("chartreuse1", 1), # GREM_H_4 A27_1_cococountout.count	
            rep("coral1", 1), # GREM_C_4  A28_1_cococountout.count	
            
            rep("deeppink", 1), # PN_H_0	A29_1_cococountout.count	
            rep("blueviolet", 1), # PN_H_30  A30_1_cococountout.count	
            rep("brown1", 1), # PN_C_30	A31_1_cococountout.count
            rep("lightslateblue", 1),   # GREM_C_30 A3_1_cococountout.count	
            rep("darkolivegreen1", 1), # PN_H_1	A32_1_cococountout.count	
            rep("darkorchid2", 1), # PN_C_1	A33_1_cococountout.count	
            rep("mistyrose", 1), # PN_H_4	A34_1_cococountout.count	
            rep("lightsalmon", 1), # PN_C_4	A35_1_cococountout.count	
            rep("deeppink", 1), # PN_H_0	A36_1_cococountout.count		
            rep("blueviolet", 1), # PN_H_30  A37_1_cococountout.count
            rep("brown1", 1), # PN_C_30	A38_1_cococountout.count
            rep("darkolivegreen1", 1), # PN_H_1	A39_1_cococountout.count		
            rep("darkorchid2", 1), # PN_C_1	A40_1_cococountout.count	
            rep("mistyrose", 1), # PN_H_4	A41_1_cococountout.count	
            rep("lightsalmon", 1), # PN_C_4	A42_1_cococountout.count	
            rep("maroon", 1),    # GREM_H_1 A4_2_cococountout.count		
            rep("deeppink", 1), # PN_H_0	A43_1_cococountout.count	
            rep("blueviolet", 1), # PN_H_30  A44_1_cococountout.count	
            rep("brown1", 1), # PN_C_30	A45_2_cococountout.count
            rep("darkolivegreen1", 1), # PN_H_1	A46_1_cococountout.count	
            rep("darkorchid2", 1), # PN_C_1	A47_1_cococountout.count	
            rep("mistyrose", 1), # PN_H_4	A48_1_cococountout.count	
            rep("lightsalmon", 1), # PN_C_4	A49_1_cococountout.count	
            rep("deeppink", 1), # PN_H_0	A50_1_cococountout.count		
            rep("blueviolet", 1), # PN_H_30  A51_1_cococountout.count	
            rep("sienna3", 1),    # GREM_C_1  A5_1_cococountout.count	
            rep("brown1", 1), # PN_C_30	A52_1_cococountout.count	
            rep("darkolivegreen1", 1), # PN_H_1	A53_1_cococountout.count	
            rep("darkorchid2", 1), # PN_C_1	A54_1_cococountout.count		
            rep("mistyrose", 1), # PN_H_4	A55_1_cococountout.count		
            rep("lightsalmon", 1), # PN_C_4	A56_1_cococountout.count	
            
            rep("chartreuse1", 1), # GREM_H_4 A6_2_cococountout.count		
            rep("coral1", 1), # GREM_C_4  A7_2_cococountout.count		
            rep("cadetblue1", 1), # GREM_H_0  A8_1_cococountout.count	
            rep("khaki1", 1) # GREM_H_30 A9_1_cococountout.count
                ),
          RowSideColors = c(    # grouping row-variables into different
            rep("lightslateblue", 1),   # GREM_C_30 A10_2_cococountout.count      (color, total number of each variable)
            rep("maroon", 1),    # GREM_H_1 A11_1_cococountout.count
            rep("sienna3", 1),    # GREM_C_1  A12_1_cococountout.count
            rep("cadetblue1", 1), # GREM_H_0  A1_2_cococountout.count
            rep("chartreuse1", 1), # GREM_H_4 A13_3_cococountout.count
            rep("coral1", 1), # GREM_C_4  A14_1_cococountout.count
            rep("cadetblue1", 1), # GREM_H_0  A15_1_cococountout.count
            rep("khaki1", 1), # GREM_H_30 A16_1_cococountout.count
            rep("lightslateblue", 1),   # GREM_C_30 A17_2_cococountout.count
            rep("maroon", 1),    # GREM_H_1 A18_1_cococountout.count
            rep("sienna3", 1),    # GREM_C_1  A19_1_cococountout.count	
            rep("chartreuse1", 1), # GREM_H_4 A20_1_cococountout.count	
            rep("coral1", 1), # GREM_C_4  A21_1_cococountout.count	
            rep("khaki1", 1), # GREM_H_30 A2_1_cococountout.count	
            rep("cadetblue1", 1), # GREM_H_0  A22_1_cococountout.count	
            rep("khaki1", 1), # GREM_H_30 A23_1_cococountout.count	
            rep("lightslateblue", 1),   # GREM_C_30 A24_2_cococountout.count	
            rep("maroon", 1),    # GREM_H_1 A25_1_cococountout.count	
            rep("sienna3", 1),    # GREM_C_1  A26_1_cococountout.count	
            rep("chartreuse1", 1), # GREM_H_4 A27_1_cococountout.count	
            rep("coral1", 1), # GREM_C_4  A28_1_cococountout.count	
            
            rep("deeppink", 1), # PN_H_0	A29_1_cococountout.count	
            rep("blueviolet", 1), # PN_H_30  A30_1_cococountout.count	
            rep("brown1", 1), # PN_C_30	A31_1_cococountout.count
            rep("lightslateblue", 1),   # GREM_C_30 A3_1_cococountout.count	
            rep("darkolivegreen1", 1), # PN_H_1	A32_1_cococountout.count	
            rep("darkorchid2", 1), # PN_C_1	A33_1_cococountout.count	
            rep("mistyrose", 1), # PN_H_4	A34_1_cococountout.count	
            rep("lightsalmon", 1), # PN_C_4	A35_1_cococountout.count	
            rep("deeppink", 1), # PN_H_0	A36_1_cococountout.count		
            rep("blueviolet", 1), # PN_H_30  A37_1_cococountout.count
            rep("brown1", 1), # PN_C_30	A38_1_cococountout.count
            rep("darkolivegreen1", 1), # PN_H_1	A39_1_cococountout.count		
            rep("darkorchid2", 1), # PN_C_1	A40_1_cococountout.count	
            rep("mistyrose", 1), # PN_H_4	A41_1_cococountout.count	
            rep("lightsalmon", 1), # PN_C_4	A42_1_cococountout.count	
            rep("maroon", 1),    # GREM_H_1 A4_2_cococountout.count		
            rep("deeppink", 1), # PN_H_0	A43_1_cococountout.count	
            rep("blueviolet", 1), # PN_H_30  A44_1_cococountout.count	
            rep("brown1", 1), # PN_C_30	A45_2_cococountout.count
            rep("darkolivegreen1", 1), # PN_H_1	A46_1_cococountout.count	
            rep("darkorchid2", 1), # PN_C_1	A47_1_cococountout.count	
            rep("mistyrose", 1), # PN_H_4	A48_1_cococountout.count	
            rep("lightsalmon", 1), # PN_C_4	A49_1_cococountout.count	
            rep("deeppink", 1), # PN_H_0	A50_1_cococountout.count		
            rep("blueviolet", 1), # PN_H_30  A51_1_cococountout.count	
            rep("sienna3", 1),    # GREM_C_1  A5_1_cococountout.count	
            rep("brown1", 1), # PN_C_30	A52_1_cococountout.count	
            rep("darkolivegreen1", 1), # PN_H_1	A53_1_cococountout.count	
            rep("darkorchid2", 1), # PN_C_1	A54_1_cococountout.count		
            rep("mistyrose", 1), # PN_H_4	A55_1_cococountout.count		
            rep("lightsalmon", 1), # PN_C_4	A56_1_cococountout.count	
            
            rep("chartreuse1", 1), # GREM_H_4 A6_2_cococountout.count		
            rep("coral1", 1), # GREM_C_4  A7_2_cococountout.count		
            rep("cadetblue1", 1), # GREM_H_0  A8_1_cococountout.count	
            rep("khaki1", 1) # GREM_H_30 A9_1_cococountout.count
          ) #this closes the rowside color function
)   #this closes the heatmap.2 function
legend("topright",
       #inset = -0.05,
        title = "Conditions",
       text.font = 14,
        legend = c("GREM_C_30",
                        "GREM_H_1",
                        "GREM_C_1",
                        "GREM_H_0",
                        "GREM_H_4",
                        "GREM_C_4",
                        "PN_H_0",
                        "PN_H_30",
                        "PN_C_30",
                        "PN_H_1",
                        "PN_C_1",
                        "PN_H_4",
                        "PN_C_4"
             ), #this closes the legend details
             fill = c("lightslateblue",
                                      "maroon",
                                      "sienna3",
                                      "cadetblue1",
                                      "chartreuse1",
                                      "coral1",
                                      "deeppink",
                                      "blueviolet",
                                      "brown1",
                                      "darkolivegreen1",
                                      "darkorchid2",
                                      "mistyrose",
                                      "lightsalmon"
             ) #this closes the fill details
             ) #this closes legend
             #) #this closes the heatmap.2 function
dev.off()

} #closes the on/off switch



###### Analysis 2: Pair-wise Analysis with DEseq2 ######

#### Preparatory Work ####

if(TRUE) {

## First you have to open up the count matrix and condition matrix files and make them yourself by copying and pasting the columns you need
   ## You cannot simply write a code to do this because we have different numbers of reps for each tissue type so simply pulling out, for example, every 6 columns won't work

# Insert the Names of the Pairs of Data You Wish to Analyze Here \/ make a nice list for a for loop
#allname <- list("leaf_root", "preveraisonberry_veraisonberry", "veraisonberry_matureberry")
#     ^ example from the transcriptome by tissue analysis
#                 \/ must be written exactly as done in the '...Cond.Input.File.txt' files

#Do just one
#allname <- list("GREM_H_30_v_GREM_H_1")

#Do all of them intra GREM4
#allname <- list("GREM_C_1_v_GREM_H_1","GREM_C_1_v_GREM_H_0","GREM_H_0_v_GREM_H_4","GREM_H_0_v_GREM_H_30","GREM_H_1_v_GREM_C_1","GREM_H_1_v_GREM_H_0","GREM_H_1_v_GREM_H_30","GREM_C_4_v_GREM_H_0","GREM_C_30_v_GREM_H_0","GREM_C_30_v_GREM_H_30","GREM_H_0_v_GREM_H_1","GREM_H_1_v_GREM_H_4","GREM_H_4_v_GREM_H_0","GREM_H_4_v_GREM_H_1","GREM_H_4_v_GREM_H_30","GREM_H_30_v_GREM_H_0","GREM_H_30_v_GREM_H_1","GREM_H_30_v_GREM_H_4","GREM_H_30_v_GREM_C_30","GREM_C_4_v_GREM_H_4","GREM_H_4_v_GREM_C_4","GREM_H_0_v_GREM_C_30","GREM_H_0_v_GREM_C_1","GREM_H_0_v_GREM_C_4")

#Do all non-problematic GREM4
#allname <- list("GREM_C_1_v_GREM_H_0","GREM_H_0_v_GREM_H_4","GREM_H_0_v_GREM_H_30","GREM_H_1_v_GREM_H_0","GREM_C_4_v_GREM_H_0","GREM_C_30_v_GREM_H_0","GREM_C_30_v_GREM_H_30","GREM_H_0_v_GREM_H_1","GREM_H_1_v_GREM_H_4","GREM_H_4_v_GREM_H_0","GREM_H_4_v_GREM_H_1","GREM_H_4_v_GREM_H_30","GREM_H_30_v_GREM_H_0","GREM_H_30_v_GREM_H_4","GREM_H_30_v_GREM_C_30","GREM_C_4_v_GREM_H_4","GREM_H_4_v_GREM_C_4","GREM_H_0_v_GREM_C_30","GREM_H_0_v_GREM_C_1","GREM_H_0_v_GREM_C_4")

#Problematic but completed GREM4 (don't have any up or down regulated genes)

    #No Up-regulated Genes
        #"GREM_C_1_v_GREM_H_1"

    #No Down-regulated Genes
        #"GREM_H_1_v_GREM_C_1"

    #No Up or Down-regulated Genes
        #"GREM_H_1_v_GREM_H_30"
        #"GREM_H_30_v_GREM_H_1"

  
  

#Do all of them intra PN
#allname <- list("PN_H_30_v_PN_H_0","PN_H_30_v_PN_H_1","PN_C_30_v_PN_H_0","PN_H_1_v_PN_H_0","PN_H_1_v_PN_H_30","PN_C_1_v_PN_H_0","PN_H_4_v_PN_H_0","PN_H_30_v_PN_H_4","PN_C_4_v_PN_H_0","PN_H_0_v_PN_H_30","PN_H_4_v_PN_H_30","PN_H_0_v_PN_C_30","PN_H_0_v_PN_H_1","PN_H_1_v_PN_H_4","PN_H_0_v_PN_C_1","PN_H_0_v_PN_H_4","PN_H_4_v_PN_H_1","PN_H_0_v_PN_C_4","PN_H_30_v_PN_C_30","PN_H_1_v_PN_C_1","PN_H_4_v_PN_C_4","PN_C_30_v_PN_H_30","PN_C_1_v_PN_H_1","PN_C_4_v_PN_H_4")

#Do the select list of intra PN
#allname <- list("PN_H_30_v_PN_H_0","PN_H_1_v_PN_H_0","PN_H_4_v_PN_H_0","PN_H_30_v_PN_C_30","PN_H_1_v_PN_C_1","PN_H_4_v_PN_C_4","PN_H_1_v_PN_H_30","PN_H_4_v_PN_H_30","PN_H_4_v_PN_H_1","PN_C_30_v_PN_H_0","PN_C_1_v_PN_H_0","PN_C_4_v_PN_H_0")

#Do all non-problematic PN
#allname <- list("PN_H_30_v_PN_H_0","PN_H_1_v_PN_H_0","PN_H_4_v_PN_H_0","PN_H_4_v_PN_H_30","PN_H_4_v_PN_H_1","PN_C_30_v_PN_H_0","PN_C_1_v_PN_H_0","PN_C_4_v_PN_H_0")

#Problematic but completed PN (don't have any up or down regulated genes)

    #No Up-regulated Genes
        #None

    #No Down-regulated Genes
        #"PN_H_30_v_PN_C_30"
        #"PN_H_1_v_PN_C_1"
        #"PN_H_4_v_PN_C_4"
  
        #"PN_H_1_v_PN_H_30"
  
  
  
  
#Do all of them inter
allname <- list("GREM_H_0_v_PN_H_0","GREM_H_30_v_PN_H_30","GREM_H_1_v_PN_H_1","GREM_H_4_v_PN_H_4")

#Do literally all of them
#allname <- list("GREM_C_1_v_GREM_H_0","GREM_C_1_v_GREM_H_1","GREM_C_4_v_GREM_H_0","GREM_C_4_v_GREM_H_4","GREM_C_30_v_GREM_H_0","GREM_C_30_v_GREM_H_30","GREM_H_0_v_GREM_H_1","GREM_H_0_v_GREM_H_4","GREM_H_0_v_GREM_H_30","GREM_H_0_v_PN_H_0","GREM_H_1_v_GREM_C_1","GREM_H_1_v_GREM_H_0","GREM_H_1_v_GREM_H_4","GREM_H_1_v_GREM_H_30","GREM_H_1_v_PN_H_1","GREM_H_4_v_GREM_C_4","GREM_H_4_v_GREM_H_0","GREM_H_4_v_GREM_H_1","GREM_H_4_v_GREM_H_30","GREM_H_4_v_PN_H_4","GREM_H_30_v_GREM_C_30","GREM_H_30_v_GREM_H_0","GREM_H_30_v_GREM_H_1","GREM_H_30_v_GREM_H_4","GREM_H_30_v_PN_H_30","PN_C_1_v_PN_H_0","PN_C_4_v_PN_H_0","PN_C_30_v_PN_H_0","PN_H_0_v_GREM_H_0","PN_H_0_v_PN_H_1","PN_H_0_v_PN_H_4","PN_H_0_v_PN_H_30","PN_H_1_v_GREM_H_1","PN_H_1_v_PN_C_1","PN_H_1_v_PN_H_0","PN_H_1_v_PN_H_4","PN_H_1_v_PN_H_30","PN_H_4_v_GREM_H_4","PN_H_4_v_PN_C_4","PN_H_4_v_PN_H_0","PN_H_4_v_PN_H_1","PN_H_4_v_PN_H_30","PN_H_30_v_GREM_H_30","PN_H_30_v_PN_C_30","PN_H_30_v_PN_H_0","PN_H_30_v_PN_H_1","PN_H_30_v_PN_H_4","GREM_H_0_v_GREM_C_30","GREM_H_0_v_GREM_C_1","GREM_H_0_v_GREM_C_4")
#  Experimental Sample /\       /\ Control Sample      !!!!!!!!!!!!

  

#Do a subset
#allname <- list("PN_H_4_v_PN_H_30","PN_H_4_v_PN_H_1","PN_C_30_v_PN_H_0","PN_C_1_v_PN_H_0","PN_C_4_v_PN_H_0")




for (p in allname){
   print (p)

## Read in the Count Matrix Input File (the data which you will do the comparisons on):
data <- read.table(paste (p,".countmatrix_IH.txt",sep=""), header=T, row.names=1, com='', quote='', check.names=F, sep="\t")


## Round the data using the round function to get it into integers and no decimal places (since the data may or many not be in integer format if you use CoCo for example)
data <- round(data)


## Read in the Condition Input File (the comparison you want to make):
sample <- read.table(paste (p,".Cond.Input.File.txt",sep=""), header=T, row.names=1, com='',
                       quote='', check.names=F, sep="\t", colClasses="factor")
sample <- sample[match(colnames(data), rownames(sample)),, drop=F]
sample_rowname <- rownames(sample)
  



#### Variability Check with PCA Plot ####


## Read Files In to DESeq2-usable Variables with Specific Format Parameters Dictated
# (DESeqDataSetFromMatrix is a DESeq2 parameter we are currently setting)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data, 
                                            colData = sample,  
                                            design = ~ condition)
head(assay(ddsFullCountTable), 3)


## Pre-filter the dataset (data) so that low read genes are culled from the analysis as they are likely
#       uninformative.  See 'http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering'
#       I am not pre-filtering as it only speeds up the analysis
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]


## Normalize the 'ddsFullCountTable' input count data so that is can be more visually amiable to plotting
#       Using the 'rlog' transformation
rlogddsFullCountTable <- rlog(ddsFullCountTable, blind = FALSE)
head(assay(rlogddsFullCountTable), 3)


## PCA plot with rlog normalized data - the program will ONLY allow this to be run on normalized data

#Create the PCA Plot (dataset that is)
PCAplotoutputdata <- plotPCA(rlogddsFullCountTable,
                        ntop = 30000,  #how many of the top most expressed genes should be included in the analysis? (default appears to be 500 of them)
                        intgroup = "condition",  #how should you sort the samples for the PCA plot? (by condition in my case)
                        returnData = TRUE  #return the output data that would originally be graphed as data and not a graphical output (FALSE produces a graphic)
                        )

#Plot the PCA Data with Sample Labels
pdf(paste(p,".rlogPCAPlotwlabels.pdf",sep=""),wi = 8, he = 6)
print(ggplot(PCAplotoutputdata,
       aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 5) +
  geom_label_repel(aes(label = name)) +
  ggtitle(paste(p," PCA Plot of Sample Variance After rlog Transformation",sep="")) +
  theme(text = element_text(size = 10)
))
dev.off()

#Plot the PCA Data without Sample Labels
pdf(paste(p,".rlogPCAPlot.pdf",sep=""),wi = 8, he = 6)
print(ggplot(PCAplotoutputdata,
             aes(x = PC1, y = PC2, color = condition)) +
        geom_point(size = 5) +
        #geom_label_repel(aes(label = name)) +
        ggtitle(paste(p," PCA Plot of Sample Variance After rlog Transformation",sep="")) +
        theme(text = element_text(size = 10)
        ))
dev.off()




#### Run DESeq2 ####


## Run DESeq2
dds <- DESeq(ddsFullCountTable)


## Create an Output File of the DESeq2 Results

# Normalize the DESeq2 output data (dds) and convert to a file-printable format
normalized_counts <- counts(dds, normalized=TRUE)

# DO NOT normalize the DESeq2 output data (dds) and convert to a file-printable format
unnormalized_counts <- counts(dds, normalized=FALSE)

# Sort the DESeq2 Output by Median Absolute Deviation on both objects
#       Normalized
normalized_counts_mad <- apply(  #apply takes a data frame or matrix and provides output in the form of a vector, list, or array
                               normalized_counts,  #object to work on
                               1, #a value between 1 and 2 to define where to apply the function where 1=manipulate rows and 2=manipulate columns and c(1,2)=manipulate rows and columns
                               mad #a function to apply; mad=median absolute deviation; mean; median; sum; min; max
                               )
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

#       Un-Normalized
unnormalized_counts_mad <- apply(  #apply takes a data frame or matrix and provides output in the form of a vector, list, or array
                                unnormalized_counts,  #object to work on
                                1, #a value between 1 and 2 to define where to apply the function where 1=manipulate rows and 2=manipulate columns and c(1,2)=manipulate rows and columns
                                mad #a function to apply; mad=median absolute deviation; mean; median; sum; min; max
                                )
unnormalized_counts <- unnormalized_counts[order(unnormalized_counts_mad, decreasing=T), ]

# Create a output excel file for both Normalized and Un-Normalized (raw) DESeq2 Output
write_xlsx(as.data.frame(normalized_counts), paste(p,".coco.count.matrix.DESeq2.Normalized.xlsx", sep=""))
write_xlsx(as.data.frame(unnormalized_counts), paste(p,".coco.count.matrix.DESeq2.Unnormalized.xlsx", sep=""))




#### Prepare DESeq2 Output for Downstream Analysis ####


## Create the naming conventions based on the file name needed for use throughout the rest of this script
#                               This \/ is where you dictate how to read in your input data (so that is to say, what is your divider in your file name)
tmpname <- matrix(unlist(strsplit(p,"_v_",2)),ncol = 2)
sampleA = tmpname[1,1]
sampleB = tmpname[1,2]  #these two extract sample names and give them to the next function which creates the variable contrastv


## Produce Summary of Differential Expression Results as an Object
contrastV <- c("condition", sampleA, sampleB) #in quotes put whatever your column header title is for the condition matrix
res <- results(dds,  contrast=contrastV, lfcThreshold = lfcvalue, alpha = padjvalue) #res is a necessary input list file that DEseq2 must use for analysis
head(res, 10) #display the first 10 lines of the object
summary(res)


## Add Additional Informative Information to 'res' which DESeq2 Does Not Naturally Include

# Calculate the Average Expression for Each Gene in Sample A
baseA <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleA]
if (is.vector(baseA)){
  baseMeanA <- as.data.frame(baseA)
} else {
  baseMeanA <- as.data.frame(rowMeans(baseA))
}
colnames(baseMeanA) <- sampleA

# Calculate the Average Expression for Each Gene in Sample B
baseB <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleB]
if (is.vector(baseB)){
  baseMeanB <- as.data.frame(baseB)
} else {
  baseMeanB <- as.data.frame(rowMeans(baseB))
}
colnames(baseMeanB) <- sampleB

# Add the Average Expression For Each Gene for Samples A and B Calculated Just Above to 'res'
res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))

# Add the Gene IDs in Again but with a Column Header of ID this time
res <- cbind(ID=rownames(res), as.data.frame(res))

# Not sure what this does, it does not appear to affect the contents of 'res' in any way
res$baseMean <- rowMeans(cbind(baseA, baseB))


## change padj as NA to 1
res$padj[is.na(res$padj)] <- 1


## Order the listing by p-value and padj value
res <- res[order(res$pvalue),]
res <- res[order(res$padj),]
#head(res)


## Write Files with the Newly Created 'res' object (creates a spreadsheet that outputs expression data for every single gene in human-readable format, whether a DEG or not)
file_base <- paste(paste(p,".coco.count.matrix.DESeq2.results",sep=""))
write_xlsx(as.data.frame(res), paste(file_base,".xlsx", sep=""))
write_csv(as.data.frame(res), paste(file_base,".csv", sep=""))
#^you need to write it as a .csv as well for the GSEA analysis downstream

## The below totally works, just saving for posterity
#file_base <- paste(paste(p,".coco.count.matrix.DESeq2",sep=""), sep=".")
#file_base1 <- paste(file_base, "results", sep=".")
#write_xlsx(as.data.frame(res), paste(file_base1, ".xlsx"))




#### Identify DEGs ####


## Identify Differentially Expressed Gene Identification Parameters: padj<0.1
#padj = p-value adjusted; used to adjust the p-value you wish to use as the cut off (the stat. sig. ie 0.05)
res_de <- subset(res, 
                 res$padj<padjvalue, # THIS IS WHERE YOU SET YOUR padj VALUE which you use for your padj cutoff for DEGs (which you previously calculated using 'results' command from DESeq2 but here is where you actually only pull out those DEGs whereas in 'results' it still outputs all of the genes, DEGs or not, but it makes the determination based on your parameters fed in such as log2foldchange and padj)
                 select=c('ID',
                           sampleA,
                           sampleB,
                           'log2FoldChange',
                           'padj')
                 )


## Write output files containing DEGs
# First, write a file listing genes that are up regulated DEGs
# 	    Determine foldchange (FC) threshold/cutoff: foldchange > 1 for up-regulated DEGs
#	      ('>=1' is actually equal to saying 'anything greater than or equal to 2')
res_de_up <- subset(res_de, res_de$log2FoldChange>=lfcvalue)

file_base <- paste(sampleA,"_higherThan_",sampleB, sep="") 
write_xlsx(as.data.frame(res_de_up), paste(file_base, ".xlsx", sep=""))

# Second, write a file listing genes that are down regulated DEGs
#	      Determine foldchange (FC) threshold/cutoff: foldchange > 1 for down-regulated DEGs
#	      This number below should match the number given above to make the analysis symmetrical
# 	    ('(-1)*1' is actually equal to saying 'anything less than or equal to -2')
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*lfcvalue)

file_base <- paste(sampleA,"_lowerThan_",sampleB, sep="") 
write_xlsx(as.data.frame(res_de_dw), paste(file_base, ".xlsx", sep=""))

# Combine the above up and down regulated DEG files files into one file that explicitly states in which condition each DEG is up or down regulated in in respect to the other condition (gene1 is higher in conditionA than conditionB)
res_de_up_id = data.frame(ID=res_de_up$ID, 
                          type=paste(sampleA,"_higherThan_", sampleB, sep="."))
res_de_dw_id = data.frame(ID=res_de_dw$ID, 
                          type=paste(sampleA,"_lowerThan_", sampleB, sep="."))
de_id = rbind(res_de_up_id, res_de_dw_id)

file_base <- paste(p,".coco.count.matrix.all.DEGs",sep="")
write_xlsx(as.data.frame(de_id), paste(file_base, ".xlsx", sep=""), col_names = FALSE)




#### Create Graphics ####


## Conduct an rlog Transformation of the DESeq2 Output Data (dds)
##      While these objects are made it appears that none of these objects are ever used again though...
rlogdds <- rlog(dds,     #this functions rlog transforms the DESeq2 output data (dds)
                blind=FALSE   #this should be left as 'FALSE' if planning on using for downstream analysis
)

#Not sure what this is doing...
#rlogMat <- assay(rlogdds)

#Not particularly sure what this is doing either, well, it's sorting the above, but what is the above for anyway?
#rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]



### Fig 1: Volcano Plot ###    
#It would appear that volcano plots do not take normally transformed data, nor did Rachel do that either
pdf(paste(p,".VolcanoPlot.pdf",sep=""), wi = 8, he = 8)
print(EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'pvalue',  #was previously 'pvalue', want to leave it this way so the column is scaled by p-value (because if it is scaled by padj then it is just a flat pancake of dots on the x-axis)
                      #xlim = c(-8, 8), #this sets the limits of the x axis of the graph, not the dimensions of the entire graph!
                      title = p, #this tells the program to just use the p variable from in the for loop for the title (which is just the first half of the input file name)
                      pCutoff = padjvalue, #define your cutoff value for your chosen statistical significance method (likely going to be p-adjusted (padj), likely going to be 0.05)
                      pCutoffCol = 'padj',  #Column name of statistical significance values to be used as the cut-off. A typical usage situation would be to pass nominal [un-adjusted] p-values as 'y', but adjusted p-values as pCutoffCol. In this way, a plot is generated via - log10(unadjusted p-value), but cut-offs based on adjusted p-values.
                      FCcutoff = lfcvalue, #change your fold change cutoff here
                      selectLab = ""  #dictate here if you want your DEGs to be labeled with their gene name or not (to not print them, leave this on; to print them, hash it off)
))
dev.off()
ggsave(paste(p,".VolcanoPlot.png",sep=""), width = 34, height = 22, unit = "cm", dpi = 300)



### Fig 2: Heatmap Creation (of top 20 DE-ed Genes between conditions) ###
#       This graphic does need normally transformed data


##Create a list of the top 20 up and down regulated DEGs (should 20 exist, if not, the list will compensate accordingly)
res_de_up_top20_id <- as.vector(head(res_de_up$ID,20))
res_de_dw_top20_id <- as.vector(head(res_de_dw$ID,20))


##Combine the top 20 up and down regulated DEGs into one list of 40 total genes (should 40 exist, if not, the list will compensate accordingly)
red_de_top20_updown_id <- c(res_de_up_top20_id, res_de_dw_top20_id)


##Create the Top 20 DEG Heatmap Graphic
pdf(paste(p,".top20updown.DEG.pdf",sep=""),width = 10, height = 18)
print(pheatmap(assay(rlogdds)[red_de_top20_updown_id,], #feeding in the rlog transformed data and then dictating to only use the top 20 up and down regulated DEGs (based on gene ID) and not all of them in the dataset
         main = paste(p, " Top 20 Up and Down Regulated Genes", sep=""),   #name main title of graphic
         cluster_row = T, 
         scale = "row",   #sets how the graphic should be scaled (row, column, none)
         cellwidth = 20,   #set cell width
         cellheight = 25,   #set cell height
         annotation_col = sample, #how to subset the samples for the condition legend bar
         display_numbers = FALSE, #don't display numbers in the colored boxes
         legend = TRUE,   #turns on or off the legend [TRUE or FALSE]
         legend_breaks = c(-1, 0, 1),  #dictate the numbers you want to display on the legend
         legend_labels = c("Low", "Medium", "High"),
         show_colnames = T,   #print the column names (T or F)
         fontsize_row = 12,   #font size for rows
         fontsize_column = 12, #font size for columns
         #width = 30,   #width of overall graphic [appears to not work after internet search]
         #height = 50   #height of overall graphic [appears to not work after internet search]
         )
)
dev.off()

} #Closes the name loop

} #closes the on/off switch


##### END #####




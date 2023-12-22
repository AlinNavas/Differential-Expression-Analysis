path = "C:/Users/alin2/Downloads"

folder_name = "brca_tcga_pan_can_atlas_2018.tar.gz"

folder = paste(path, folder_name, sep = "/")

untar(folder)

# Go to new path

new_dir = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )

setwd(new_dir)

# Read patient data, RNA seq data, count number
data_patient  = read.delim("data_clinical_patient.txt")
data_Rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")
data_cna =read.delim("data_cna.txt")

# Select row with erbb2
erbb2_row_id <- which(data_cna[,1] == "ERBB2")
cat(erbb2_row_id)

#Remove gene names and create matrix
assay_rna = as.matrix(data_Rnaseq[,-c(1,2)])
assay_cna = as.matrix(data_cna[,-c(1,2)])

# Build an empty metadata
metadata = matrix(0, dim(assay_rna)[2],1)

#Match the RNASeq patient ids with the CNA ids and the Patient Data ids.
#create a list of patient_ids
pat_id_rna = colnames(assay_rna)
pat_id_cna = colnames(assay_cna)

#create a index for overlapping patients 
## masterfile is rna therefore on left, like left join

filter = is.element(pat_id_rna,pat_id_cna)

#create a subset of rna patients
rna_subset =  assay_rna[,filter]
rna_subset_pat_id = colnames(rna_subset)

#create metadata of erbb2 as a filter 

for (i in 1:length(rna_subset_pat_id)){
  idx = which(pat_id_cna == rna_subset_pat_id[i])
  metadata[i,1] = 1*(data_cna[erbb2_row_id,idx] > 0 )
}
metadata[is.na(metadata)] =0

# checking number of positive values
sum(metadata)

#Label metadata column as a filter
colnames(metadata)= "ERBB2"


library(DESeq2)
# Build DESeq Object

assay_rna[is.na(assay_rna)] = 0  # Impute with zeros the NA
assay_rna[assay_rna<0] = 0

#The count data for the DESeqDataSet is provided by rounding the values in the assay matrix. The round function is used to ensure that the count data consists of integers.
# The column data (metadata) is provided by the metadata matrix,
# The design formula specifies the model that will be fitted to the data.

dds <- DESeqDataSetFromMatrix(countData = round(assay_rna),
                              colData = metadata,
                              design = ~ ERBB2)
# Filter
#counts(dds) >= 10: This creates a logical matrix where each element is TRUE if the corresponding count is greater than or equal to 10, and FALSE otherwise.
#This calculates the sum of TRUE values (counts greater than or equal to 10) for each row. The result is a numeric vector representing the count of values meeting the condition for each row.
#>= smallestGroupSize: This part checks if the row sum is greater than or equal to a threshold value smallestGroupSize. It's a way to ensure that there are a sufficient number of counts meeting the condition in each row
#This will create a new DESeqDataSet (dds_filtered) with only the rows that satisfy the count criterion. Adjust smallestGroupSize based on your analysis requirements.

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Normalize

dds <- DESeq(dds)

# Get Results

res <- results(dds)

# Summary

summary(res)
rownames(res) = data_Rnaseq[keep,1]

# Significantly Differentially Expressed are filtered

# res$padj refers to the adjusted p-values (corrected for multiple testing) obtained from some statistical analysis
# So, signif will be a vector containing the indices of the elements in res$padj that have an adjusted p-value less than 0.05. 

signif = which(res$padj<0.0005)
deg = res[signif,]



# Separate them into upregulation and downregulation
# dup will contain only those rows from the original data frame deg where the values in the second column are greater than 0.

dup = deg[deg[,2]>0.,]
dup_sorted <- dup[order(-dup$log2FoldChange), ]
dup_sorted[1:10,]

ddown = deg[deg[,2]<0.,]
ddown_sorted <- ddown[order(ddown$log2FoldChange), ]
ddown_sorted[1:10,]

# For Pathway Enrichment we need Entrez IDs
entrez_all = data_Rnaseq[keep[signif],2]
entrez_up = data_Rnaseq[keep[signif[deg[,2]>0.]],2]
entrez_down = data_Rnaseq[keep[signif[deg[,2]<0.]],2]

# Pathway Enrichment

library(clusterProfiler)

up_paths = enrichKEGG(gene = entrez_up, organism = 'hsa', pvalueCutoff = 0.05)
head(up_paths)


down_paths = enrichKEGG(gene = entrez_down, organism = 'hsa', pvalueCutoff = 0.05)
head(down_paths)


# Transform the data to visualize
# vst(): This function performs a variance stabilizing transformation on the count data. The variance stabilizing transformation is useful for making the variance of the data more uniform across the range of expression levels.
# blind=FALSE: The blind argument is set to FALSE, indicating that the design information present in the dds object will be used during the transformation. This is typically done when you have a design matrix in your DESeqDataSet and want to account for the experimental design during the transformation.

rld <- vst(dds, blind=FALSE)

# Do Principal Components Analysis
pc = prcomp(assay(rld))

# Plot 
# pc$rotation[,1]: This extracts the scores of the first principal component from the PCA result stored in the pc object.
# In PCA, the rotation component of the PCA result contains the loadings, which represent the weights assigned to each variable (feature) in the dataset for each principal component.
# This indexing syntax is used to extract the values in the first column of the pc$rotation matrix. In other words, it extracts the loadings of the first principal component.
# pch = 19: The pch argument specifies the type of point used in the plot. pch = 19 corresponds to filled circles.
# col = 1 + (metadata): The col argument specifies the color of each point in the plot. It looks like you're assigning colors based on the values in the metadata vector. The 1 + (metadata) is used to shift the color indices, possibly to avoid having points with color 0 (if metadata has values starting from 0).

plot(pc$rotation[,1], pc$rotation[,2], col = 1+(metadata), pch = 19)

# The pca plot might not show much but it is good to show it.

hist(assay[1000,], breaks = 50)

# this shows regularised assay

hist(assay(rld)[1000,], breaks = 50)

browseKEGG(down_paths, 'hsa05165')

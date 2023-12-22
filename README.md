Gene expression analysis and interpretation for HER2 +ve breast cancer

Patient data, RNA-seq data, and copy number alteration (CNA) data are loaded into R.
The row index of the ERBB2 gene in the CNA data is identified.
Patient IDs from RNA-seq and CNA data are matched, and a subset of RNA-seq data is created based on the overlapping patients.
Metadata is generated for ERBB2 status based on CNA data.
DESeq2 objects are created, and differential expression analysis is performed.
Rows with low counts are filtered out, and normalization using DESeq is applied.
Results and differentially expressed genes are obtained.
KEGG pathway enrichment analysis is performed for upregulated and downregulated genes.

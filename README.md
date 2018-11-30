# stat_genomics

In this work, I analyze high throughput sequencing data generated from four cell-lines in mice, involved in different stages of differentiation of blood cells. RNA-Seq data was generated using two protocols along with ATAC-Seq data for each of the cell lines. Two biological replicates were used for each protocol. I analyzed four pairs: HSC vs EBT, EBT vs CMP, CMP vs CFU-E and CFU-E vs EBT. The RNA-Seq and ATAC-Seq data were analyzed for each of these pairs, to characterize differentially expressed genes and chromatin patterns, using DE-Seq2 (and Limma) and Homer respectively. The results are discussed in the context of known biochemical function of the genes.


The analysis of the differential gene expression is based on the RNA-seq data and ATAC-seq data provided by Dr. Hardison’s lab, at Pennsylvania State University. All four of the cell lines were isolated from the bone marrow cells of 20-25 C57BL6J mouse femur and tibia as described in their protocol: https://www.encodeproject.org/biosamples/ENCBS209VOV/


I analyzed four pairs of cell lines for this project: HSC vs Erythroblast, CMP vs CFUE, CMP vs CFU-E, and CFUE vs Erythroblast. For RNA-seq data we worked with the tsv files and for ATAC-seq data we used BAM files. I used two technical replicates, namely ScriptSeq and TotalScript for every cell line in the RNA-Seq. Each technical replicate contains two biological replicates. 


I took the expected counts of every sample and filtered out the genes with less than an average of 10 counts per sample. For the exploratory data analysis, we first performed a regularized log transform using DESeq2. We then visualized the data with density plots (Appendix A) to check for unimodality, heatmaps (Appendix I) showing the agreement between replicates, and principal components plots to check data quality by checking that biological replicates were clustering (Appendices B,H). We normalized our data using Trimmed Mean of M-values (TMM) to eliminate the library size effect. We then assessed the normalization quality with MA-plots within replicates (Appendix C).


For each pair of cell lines that we analyzed, I used DESeq2 and Limma for differential expression analysis. With DESeq2, I first estimated the dispersion parameters of the normalized data. Then I fitted the data in a negative binomial GLM and tested for differential expression across the two treatment groups for each gene. I controlled the false discovery rate (FDR) at 0.01 by adjusting p-values using the Benjamini-Hochberg procedure. With Limma, I first transformed our data to log2-counts per million using voom. I then fit a one way ANOVA model to the transformed data using empirical Bayes methods to compute moderated test statistics. I used the same FDR cut-off, so as to make results comparable.


I filtered the bam files for reads with mapq greater than 10, unique mapping and both ends paired. I used these reads to find broad peak regions in Homer (style: Histone) for each of the replicates. I then removed all peaks which overlapped with mm10 blacklist regions (see references) for mm10 as identified by ENCODE and modENCODE consortia. 


I performed a tss centric analysis for peaks identified from each sample, counted number of tags that map near tss and obtained a raw counts matrix. This matrix was input to DESeq2 to differential expression analysis. To be conservative, I set the log2-fold change ratio cut-off at 2 to be considered differentially expressed and controlled FDR at 0.01.


I compared log2-fold change ratios between genes identified as differentially expressed in RNA-Seq analysis, and differentially chromatin accessible, in the ATAC-Seq analysis. I plot these log2-fold ratios in a scatter plot and check if differentially expressed genes also have differential chromatin accessibility and vice-versa. I also pay attention to whether this difference is positive or negative (upregulated or downregulated) in both these analyses.

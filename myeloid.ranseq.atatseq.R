
library(downloader)
ids = scan("encodeIds.txt",sep="\t",what="character")
cell.lines = c()
cell.ids = c()
for(i in 1:length(ids)){
  if(substr(ids[i],nchar(ids[i]),nchar(ids[i])) == ")"){
    cell.lines = c(cell.lines,unlist(strsplit(ids[i]," "))[1])
    next}
  cell.ids= c(cell.ids,unlist(strsplit(ids[i]," ",fixed=T))[2])
}

#Adding replicates to the names
cell.lines = paste(rep(paste(cell.lines, rep(1:3, times=4),sep="_"),each=2),1:2,sep="_")


#Downloading data; Replace file-tag if need to download alternative format
url = "https://www.encodeproject.org/files/<id>/@@download/<id>.tsv"
for(i in 1:length(cell.ids)){
  url_i = gsub("<id>",cell.ids[i],x=url,fixed=T)
  download.file(url_i,paste(cell.lines[i],"tsv",sep="."),method="wget")
}


#EDA of HSC replicates
setwd("/data/bzk18/Project/Misc/Spring17/STAT555/Project/tsv")
hsc1_1 = read.table("HSC_1_1.tsv",sep="\t",header=T)
hsc1_2 = read.table("HSC_1_2.tsv",sep="\t",header=T)
hsc2_1 = read.table("HSC_2_1.tsv",sep="\t",header=T)
hsc2_2 = read.table("HSC_2_2.tsv",sep="\t",header=T)

colnames(hsc1_1)
hsc = data.frame(hsc1_1$expected_count,hsc1_2$expected_count,hsc2_1$expected_count,hsc2_2$expected_count )
rownames(hsc) = hsc1_1$gene_id
colnames(hsc) = c("HSC1_1","HSC1_2","HSC2_1","HSC2_2")
dim(hsc)
cat("Total Number of genes = ", nrow(hsc),"\n")
cat("Percent genes expressed in each HSC sample = ",sapply(1:ncol(hsc), function(i) round(mean(hsc[,i]!=0)*100,1) ),"\n")


setwd("/data/bzk18/Project/Misc/Spring17/STAT555/Project/tsv")
ebt1_1 = read.table("Erythroblast_1_1.tsv",sep="\t",header=T)
ebt1_2 = read.table("Erythroblast_1_2.tsv",sep="\t",header=T)
ebt2_1 = read.table("Erythroblast_2_1.tsv",sep="\t",header=T)
ebt2_2 = read.table("Erythroblast_2_2.tsv",sep="\t",header=T)

ebt = data.frame(ebt1_1$expected_count,ebt1_2$expected_count,ebt2_1$expected_count,ebt2_2$expected_count )
rownames(ebt) = ebt1_1$gene_id
colnames(ebt) = c("EBT1_1","EBT1_2","EBT2_1","EBT2_2")
dim(ebt)
cat("Percent genes expressed in each EBT sample = ",sapply(1:ncol(ebt), function(i) round(mean(ebt[,i]!=0)*100,1) ),"\n")

sample = data.frame(hsc,ebt)
head(sample)
cat("Percent genes not expressed in any sample = ",round(sum(rowSums(sample)==0)/nrow(sample)*100,1),"\n")
#Filtering out genes not expressed in any sample
sample = sample[rowSums(sample)!=0,]
cat("Number of genes in our sample now is ",nrow(sample),"\n")
head(sample)


#Creating DESEqDataSet and estimating size factors for each column
library(DESeq2)
library(geneplotter)
sample=as.matrix(data.frame(hsc,ebt))
sample = round(sample)
#Filtering genes which don't have atleast <10 reads> on each sample
head(sample)
sample = sample[rowSums(sample)>=80,]
head(sample)
spikeins = grep(pattern="gSpikein",rownames(sample))
condition <- factor(rep(c("HSC","EBT"),each=4))
dds <- DESeqDataSetFromMatrix(sample, DataFrame(condition), ~ condition)

#RLogTransformation for clustering samples in heatmap
rlg = rlogTransformation(dds,blind=T)
#Plotting density plots 
df = reshape2::melt(assay(rlg),value.name ="values" )
df$Var1=NULL
df$condn = colData(dds)[df$Var2,1]
print(ggplot(df, aes(x = values, colour = Var2, fill = Var2)) + ylim(c(0,0.5)) +
  geom_density(alpha = 0.2, size = 1.25) + facet_wrap( ~condn) +
  xlab(expression(log[2](count + k))))

boxplot(assay(rlg),names=colnames(assay(rlg)),notch=T,ylab="log2(counts+k)")

#Plotting Heatmap
rlg_dist = dist(t(assay(rlg)))
rlg_mat = as.matrix(rlg_dist)
heatmap(rlg_mat)

#Plotting First 2 PCs 
plotPCA(rlg,intgroup="condition")


#Normalization using Relative Log Expressions
dds = estimateSizeFactors(dds,type="iterate",controlGenes=spikeins)
sizeFactors(dds)
#Trying TMM Normalization (I prefer this)
normfactors = calcNormFactors(counts(dds),method="TMM")
sizeFactors(dds) <- normfactors

#MA plot between all samples
library(ggplot2)
require(gridExtra)
MA.idx = t(combn(1:8, 2))
plots = list(rep(0,nrow(MA.idx)))
count=1
for( i in 1:length(MA.idx)){
  A = 0.5*(counts(dds,normalized=T)[,MA.idx[i,1]] + counts(dds,normalized=T)[,MA.idx[i,2]])
  M = counts(dds,normalized=T)[,MA.idx[i,1]] - counts(dds,normalized=T)[,MA.idx[i,2]]
  df = data.frame(A, M)
  plots[[count]] = ggplot(df, aes(x = A, y = M)) + geom_point(size = 1.5, alpha = 1/5) +
    geom_hline(yintercept = 0,color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")+ggtitle(paste(colnames(sample[i])," vs ",colnames(sample)[j]))
  count = count+1
}

ggsave("MA_plots.pdf", arrangeGrob(plots,nrow=nrow(MA.idx),ncol=ncol(MA.idx)))


#Estimate Dispersion Parameters for each gene
dds = estimateDispersions(dds)
plotDispEsts(dds)

#Differential Expression Analysis 
dds = nbinomWaldTest(dds)
#Setting Significance level
alpha = 0.01
#Adjusting P-values for multiple testing
#Make sure to set independent filtering to FALSE
ddres = results(dds,p.adjust="bonf",alpha=alpha,independentFiltering = F)
table(ddres$padj<alpha)/length(ddres$padj)
table(ddres$padj<alpha)
#Too many significiant P values?
#Checking if P-values follow correct distribution
hist(ddres$padj,xlab="Adjusted Pvalues",main="Hist of P.adj")
# DESeq2::plotMA(ddres,alpha=alpha)
#Checking for outliers
cat("Number of Outliers =",sum(is.na(ddres$pvalue)),"\n")
#Estimating Variance of the Null Model
library(fdrtool)
ddres2 = ddres
fdr.res = fdrtool(ddres2$stat,statistic="normal",plot=F)
print(fdr.res$param[1,"sd"])
#Null Model sd seems significantly larger than 1
#Adjusting New Pvalues
new.p.adj = p.adjust(fdr.res$pval,"fdr")
new.p.adj = qvalue(fdr.res$pval,0.01)$qvalues
table(new.p.adj<alpha)/length(new.p.adj)
table(new.p.adj<alpha)
#297 for bonferroni or 777 for BH seems reasonable
hist(new.p.adj,xlab="Adjusted Pvalues",main="Hist of P.adj")
# ddres$padj = new.p.adj
#Now the hist is U shaped!!!


#DE With Limma
source("http://www.bioconductor.org/biocLite.R")
library(limma)
library(statmod)
library(edgeR)
#Significance level
alpha= 0.01
es = data.frame(hsc,ebt)
#Filter genes with < <10> reads per sample
keep =(rowSums(es)>80)
es = es[keep,]
head(es)
nrow(es)
#Create a DGEList
dge= DGEList(counts=es,samples=colnames(es),group=rep(c("HSC","EBT"),each=4))
#Normalize with TMM
dge = calcNormFactors(dge,method="TMM")
design = model.matrix(~ dge$samples$group)
#Doing voom stuff
v = voom(dge,design=design,plot=T)
#Fit into linear models for contrasting between samples
fit = lmFit(v,design)
#Emperical estimation of t-statistics (This gives Pvalues!)
fit = eBayes(fit)
#Nice dataframe extracting function, can adjust Pvalues here
lmres = topTable(fit,coef=ncol(design),lfc=0, p.value=1,number=nrow(es),adjust.method="fdr")
table(lmres$adj.P.Val<alpha)
dge$p.adj = lmres$adj.P.Val
dge$lfc = lmres$logFC
#Plot to see if Pvalues follow correct distn
hist(lmres$adj.P.Val,xlab="Adjusted Pvalues",main="Hist of P.adj")


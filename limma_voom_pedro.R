####################################################################
## Análisis de los resultados de STAR y featureCounts usando      ## 
## limma + voom                                                   ##
##                                                                ##
##               Aligner: STAR                                    ##
##               Quantification: featureCounts                    ##
##               Differential expression: limma-voom              ##
##               Some graphical representations: ggplot2          ##
##                                                                ##
##                                                                ##
## Autores: Pedro de los Reyes Rodriguez pedro.reyes@ibvf.csic.es ##
##          Francisco J. Romero Campero                           ##
####################################################################

## Required packages
library(edgeR) #also load limma as a dependency
library(ggplot2)
library(ggrepel)
library(sva) #batch correction
library(dplyr)
library(pheatmap)
library(readr)
library(viridis) #Default Color Maps
library(MetBrewer) #cool palette
library(DESeq2) #just to calculate fpkms
library(reshape2) #reorganize data frames
library(tidyr) #tidy messy data


## Import the read count matrix data into R.
counts <- read.table(file="counts.txt", sep = "\t", header = T)
rownames(counts) <- counts$Geneid
head(counts)
counts["AT5G15840",]
counts["AT1G65480",]
counts["AT5G24470",]
counts["AT1G22770",]


## Get only the counts
samples <- c(paste("35SCO", 1:3, sep = "_"), paste("Col-0", 1:3, sep = "_"), 
             paste("SUC2", 1:3, sep = "_"), paste("co10", 1:3, sep = "_"))
sample.counts <- counts[,7:18]
rownames(sample.counts) <- counts$Geneid
head(sample.counts,2)
colnames(sample.counts) <- samples
# write.table(sample.counts, file = "sample_counts.txt", sep = "\t", quote = FALSE)

## Library sizes
barplot(colSums(sample.counts)*1e-6, names = colnames(sample.counts), 
        ylab="Library size (millions)",las = 2, cex.names = 0.3)


## Sample info. Same sample order in coldata than in count matrix
condition <- c(rep("x35SCO", 3),rep("col0",3), rep("suc2", 3), rep("co10",3))
type <- rep("paired-end",12)
coldata <- data.frame(condition, type)
rownames(coldata) <- colnames(sample.counts)
coldata

## Previsualizamos la similitud entre las réplicas con con el log2(counts+1) 
plot(log2(sample.counts[,1]+1),log2(sample.counts[,2]+1),pch=19,cex=0.7,xlab="35SCO_1",ylab="35SCO_2",cex.lab=1.25)
plot(log2(sample.counts[,2]+1),log2(sample.counts[,3]+1),pch=19,cex=0.7,xlab="35SCO_2",ylab="35SCO_3",cex.lab=1.25)
plot(log2(sample.counts[,1]+1),log2(sample.counts[,3]+1),pch=19,cex=0.7,xlab="35SCO_1",ylab="35SCO_3",cex.lab=1.25)

plot(log2(sample.counts[,4]+1),log2(sample.counts[,5]+1),pch=19,cex=0.7,xlab="wt_1",ylab="wt_2",cex.lab=1.25)
plot(log2(sample.counts[,5]+1),log2(sample.counts[,6]+1),pch=19,cex=0.7,xlab="wt_2",ylab="wt_3",cex.lab=1.25)
plot(log2(sample.counts[,4]+1),log2(sample.counts[,6]+1),pch=19,cex=0.7,xlab="wt_1",ylab="wt_3",cex.lab=1.25)

plot(log2(sample.counts[,7]+1),log2(sample.counts[,8]+1),pch=19,cex=0.7,xlab="SUC2_1",ylab="SUC2_2",cex.lab=1.25)
plot(log2(sample.counts[,8]+1),log2(sample.counts[,9]+1),pch=19,cex=0.7,xlab="SUC2_2",ylab="SUC2_3",cex.lab=1.25)
plot(log2(sample.counts[,7]+1),log2(sample.counts[,9]+1),pch=19,cex=0.7,xlab="SUC2_1",ylab="SUC2_3",cex.lab=1.25)

plot(log2(sample.counts[,10]+1),log2(sample.counts[,11]+1),pch=19,cex=0.7,xlab="co10_1",ylab="co10_2",cex.lab=1.25)
plot(log2(sample.counts[,11]+1),log2(sample.counts[,12]+1),pch=19,cex=0.7,xlab="co10_2",ylab="co10_3",cex.lab=1.25)
plot(log2(sample.counts[,10]+1),log2(sample.counts[,12]+1),pch=19,cex=0.7,xlab="co10_1",ylab="co10_3",cex.lab=1.25)


##########################################################
######### Batch correction using Combat-seq #############
########################################################

# devtools::install_github("zhangyuqing/sva-devel")
# https://rnabio.org/module-03-expression/0003/05/01/Batch-Correction/library(sva)
# tissue <- rep(c("root", "noroot", "noroot"),3)

# You need to input at least two parameters - a raw count matrix from RNA-Seq studies, 
# without any normalization or transformation, and a vector for batch separation
combat1 <- ComBat_seq(counts = as.matrix(sample.counts), batch = tissue, group = NULL)
# write.table(combat1, file = "combat1.txt", sep = "\t", quote = FALSE)

# may specify biological covariates, whose signals will be preserved
# in the adjusted data. If the user would like to specify one 
# biological variable, they may use the group parameter:
combat2 <- ComBat_seq(counts = as.matrix(sample.counts), batch = tissue, group = coldata$condition)
# write.table(combat2, file = "combat2.txt", sep = "\t", quote = FALSE)

## Run one of the below lines if you want to follow with batch corrected counts
# sample.counts <- as.data.frame(combat1)
# sample.counts <- as.data.frame(combat2)
## Important: check again the scatterplots between replicates

## Construimos un boxplot para comprobar que las distribuciones globales de las
## muestras son similares y comparables.
boxplot(log2(sample.counts+1),col=met.brewer(n=9, name="Cassatt2"),ylab="log2(counts + 1)",cex.lab=1.5)

##################################################################
########### Differential expression analysis ####################
################################################################

## Create DGEList object
d0 <- DGEList(sample.counts)
dim(d0)
## Calculate normalization factors. It doesn't normalize! 
##Just calculates normalization factors for use downstream
d0.norm <- calcNormFactors(d0) # tmm as default
# d0.norm <- calcNormFactors(d0, method = "TMM") 
# d0.norm <- calcNormFactors(d0, method = "upperquartile")

## Filter low expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0.norm), 1, max) < cutoff)
d <- d0.norm[-drop,] 
dim(d) #number of genes left

colnames(sample.counts)
genotype <- c(rep("x35SCO", 3),rep("col0",3), rep("suc2", 3), rep("co10",3))
group <- as.factor(genotype) #Create a new variable “group” as factor

plotMDS(d, col = as.numeric(group), labels = genotype) #Multidimensional scaling (MDS) plot


###### Voom transformation and calculation of variance weights
## Specify the model to be fitted. We do this before using voom since voom 
# uses variances of the model residuals (observed - fitted)
mm <- model.matrix(~0 + group)

# When operating on a DGEList-object, voom converts raw counts 
# to log-CPM values by automatically extracting library sizes 
# and normalisation factors from x itself. 
# Additional normalisation to log-CPM values can be specified 
# within voom using the normalize.method argument.

## ¿What is voom doing?
# Counts are transformed to log2 counts per million reads (CPM),
# where “per million reads” is defined based on the normalization 
# factors we calculated earlier
# A linear model is fitted to the log2 CPM for each gene, and the 
# residuals are calculated. A smoothed curve is fitted to the 
# sqrt(residual standard deviation) by average expression 
# (see red line in the plot). The smoothed curve is used to 
# obtain weights for each gene and sample that are passed 
# into limma along with the log2 CPMs.


# The read counts are processed by the voom function in limma to 
# convert them into log2 counts per million (logCPM) with associated 
# precision weights. If the data are very noisy, The logCPM values 
# can be normalized between samples by the voom function or can be 
# pre-normalized by adding normalization factors within edgeR (calcNormFactors).
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
nofilter.voom <- voom(d0, mm, plot = T)


# y <- voom(d, mm, plot = T, normalize="quantile")
# nofilter.voom <- voom(d0, mm, plot = T, normalize="quantile")
# Esto lo haria si veo que los boxplot de los counts son muy diferentes,
# por ejemplo si son muestras de diferentes experimentos. 

# OJO: Aquí nosotros aplicamos TMM normalization en vez de quantile, (uno u otro),
# eso lo hacemos con el paquete EdgeR cuando hacemos el calcnormFactors, 
# así se recomienda en el manual de limma: "It is usual to apply scale normalization 
# to RNA-seq read counts, and the TMM normalization method [33] in particular 
# has been found to perform well in comparative studies". 

# Boxplot after voom WITHOUT normalization
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
y.no.norm <- voom(d, mm, plot = T)
boxplot(y.no.norm$E, col=met.brewer(n=9, name="Hiroshige"), ylab="log2CPM voom",cex.lab=1.5, main= "No normalisation")

# Boxplot after voom WITH normalization (if any). Notice the difference
boxplot(y$E, col=met.brewer(n=9, name="Hiroshige"),ylab="log2CPM voom",cex.lab=1.5, main= "Normalisation")

# Plot density
reshaped.y <- melt(y$E) #just reshaping data frame to input to ggplot
p <- ggplot(aes(x=value, colour=Var2), data=reshaped.y)
p + geom_density() + xlab("Log2-cpm") + ylab("Density")

# Ajuste lineal
fit <- lmFit(y, mm)
head(coef(fit))

contrast.matrix <- makeContrasts(groupx35SCO-groupcol0, groupco10-groupcol0,
                                 groupsuc2-groupcol0, levels=colnames(coef(fit)))
head(contrast.matrix)

contrast.linear.fit <- contrasts.fit(fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

###########################################
#######Contraste 35SCO vs wt (Col0)########
###########################################
x35SCO.wt <- topTable(contrast.results, number=nrow(contrast.results),coef=1,sort.by="logFC")
# x35SCO.wt <- topTable(contrast.results, number=33602,coef=1,sort.by="logFC", p.value=0.05)#, lfc = 1)
head(x35SCO.wt)
x35SCO.wt["AT1G65480",]
x35SCO.wt["AT5G24470",]
nrow(x35SCO.wt)

fold.change.35SCO.wt <- x35SCO.wt$logFC
genes.ids.35SCO.wt <- rownames(x35SCO.wt)
p.value <- x35SCO.wt$adj.P.Val

activated.genes.35SCO <- genes.ids.35SCO.wt[fold.change.35SCO.wt > log2(2) & p.value < 0.05] 
repressed.genes.35SCO <- genes.ids.35SCO.wt[fold.change.35SCO.wt < -log2(2) & p.value < 0.05]

length(activated.genes.35SCO)
length(repressed.genes.35SCO)

# write.table(activated.genes.35SCO, file="tablas/fc_1.5/activated_genes_35SCO_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(repressed.genes.35SCO, file="tablas/fc_1.5/repressed_genes_35SCO_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(activated.genes.35SCO, file="tablas/activated_genes_35SCO_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(repressed.genes.35SCO, file="tablas/repressed_genes_35SCO_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)


##### Scatterplot 35SCO vs wt#####
voom.counts <- y$E
x35SCO.counts<- (voom.counts[,1] + voom.counts[,2] + voom.counts[,3])/3
col0.counts<- (voom.counts[,4] + voom.counts[,5] + voom.counts[,6])/3
suc2.counts<- (voom.counts[,7] + voom.counts[,8] + voom.counts[,9])/3
co10.counts<- (voom.counts[,10] + voom.counts[,11] + voom.counts[,12])/3

names(x35SCO.counts) <- rownames(voom.counts)
names(col0.counts) <- rownames(voom.counts)
names(suc2.counts) <- rownames(voom.counts)
names(co10.counts) <- rownames(voom.counts)
mean.counts <- matrix(c(x35SCO.counts, col0.counts, suc2.counts, co10.counts), ncol=4)
colnames(mean.counts) <- c("x35SCO", "Col0", "suc2", "co10")
rownames(mean.counts) <- rownames(voom.counts)
head(mean.counts)
mean.counts["AT5G15840",]

mean.counts.1.2 <- as.data.frame(mean.counts[,1:2])
mean.counts.2.3 <- as.data.frame(mean.counts[,2:3])
mean.counts.2.4 <- as.data.frame(mean.counts[,c(2,4)])


ggplot(data = as.data.frame(mean.counts.1.2), aes(x=x35SCO, y=Col0)) + geom_point() + theme_minimal()

# add a column of NAs
mean.counts.1.2$diffexpressed <- "NO"
x35SCO.wt$diffexpressed <- "NO"
# if log2Foldchange > 1  set as "UP" 
for (i in 1:nrow(x35SCO.wt))
{
  if (x35SCO.wt[i,"logFC"]>1)
  {
    up.gene <- rownames(x35SCO.wt[i,])
    mean.counts.1.2[up.gene,"diffexpressed"] <- "UP"
    x35SCO.wt[up.gene,"diffexpressed"] <- "UP"
  }
}
# if log2Foldchange < -1  set as "DOWN"
for(i in 1:nrow(x35SCO.wt))
{
  if (x35SCO.wt[i,"logFC"] < -1)
  {
    down.gene <- rownames(x35SCO.wt[i,])
    mean.counts.1.2[down.gene,"diffexpressed"] <- "DOWN"
    x35SCO.wt[down.gene,"diffexpressed"] <- "DOWN"
  }
}
# Set colors
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

p <- ggplot(data = as.data.frame(mean.counts.1.2), aes(x=x35SCO, y=Col0, col=diffexpressed)) + geom_point() + theme_minimal() +scale_color_manual(values=mycolors)
p + ggtitle("Scatterplot") + xlab("log2(CPM) 35S:CO") + ylab("log2(CPM) Col-0")


## Quick Volcano plot
gene.names <- rownames(contrast.results$coefficients)
volcanoplot(contrast.results, coef=1, names = gene.names, highlight = 5)

## MA plot
# It is used to check whether data are comparable among groups
# (i.e. normalization worked properly). This can be created 
# using a function available in limma
limma::plotMA(contrast.results, coef=1, main="35SCO vs WT")
plotWithHighlights(x35SCO.wt$AveExpr, x35SCO.wt$logFC,
                   status=x35SCO.wt$diffexpressed,
                   hl.col = c("red", "blue"), hl.cex = 0.5,
                   bg.col = "grey", xlab="Average log2-expression",
                   ylab= "log2 FC", main = "35SCO vs WT")
abline(h=0, col="darkred", lwd=2)

#################################################
###############Contraste SUC2 vs wt#############
################################################
suc2.wt <- topTable(contrast.results, number=nrow(contrast.results),coef=3,sort.by="logFC")
# suc2.wt <- topTable(contrast.results, number=33602,coef=2,sort.by="logFC", p.value = 0.05)
head(suc2.wt)
suc2.wt["AT5G15840",]
suc2.wt["AT5G24470",]

fold.change.suc2.wt <- suc2.wt$logFC
p.value <- suc2.wt$adj.P.Val
genes.ids.suc2.wt <- rownames(suc2.wt)

activated.genes.suc2 <- genes.ids.suc2.wt[fold.change.suc2.wt > log2(1.5) & p.value < 0.05 ]
repressed.genes.suc2 <- genes.ids.suc2.wt[fold.change.suc2.wt < -log2(1.5) & p.value < 0.05]

length(activated.genes.suc2)
length(repressed.genes.suc2)

# write.table(activated.genes.suc2, file="tablas/activated_genes_suc2_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(repressed.genes.suc2, file="tablas/repressed_genes_suc2_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

##### Scatterplot suc2 vs WT #####
ggplot(data = as.data.frame(mean.counts.2.3), aes(x=Col0, y=suc2)) + geom_point() + theme_minimal()

# add a column of NAs
mean.counts.2.3$diffexpressed <- "NO"
suc2.wt$diffexpressed <- "NO"
# if log2Foldchange > 1  set as "UP" 
for (i in 1:nrow(suc2.wt))
{
  if (suc2.wt[i,"logFC"]>1)
  {
    up.gene <- rownames(suc2.wt[i,])
    mean.counts.2.3[up.gene,"diffexpressed"] <- "UP"
    suc2.wt[up.gene,"diffexpressed"] <- "UP"
  }
}
# if log2Foldchange < -1  set as "DOWN"
for(i in 1:nrow(suc2.wt))
{
  if (suc2.wt[i,"logFC"] < -1)
  {
    down.gene <- rownames(suc2.wt[i,])
    mean.counts.2.3[down.gene,"diffexpressed"] <- "DOWN"
    suc2.wt[down.gene,"diffexpressed"] <- "DOWN"
  }
}
# Set colors
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

p <- ggplot(data = as.data.frame(mean.counts.2.3), aes(x=Col0, y=suc2, col=diffexpressed)) + geom_point() + theme_minimal() +scale_color_manual(values=mycolors)
p + ggtitle("Scatterplot") + xlab("log2(CPM) Col-0") + ylab("log2(CPM) SUC2")



## Quick Volcano plot 
gene.names <- rownames(contrast.results$coefficients)
volcanoplot(contrast.results, coef=2, names = gene.names, highlight = 5)

##MA plot
# It is used to check whether data are comparable among groups
# (i.e. normalization worked properly). This can be created 
# using a function available in limma
limma::plotMA(contrast.results, coef=3, main="SUC2 vs WT")
plotWithHighlights(suc2.wt$AveExpr, suc2.wt$logFC,
                   status=suc2.wt$diffexpressed,
                   hl.col = c("red", "blue"), hl.cex = 0.5,
                   bg.col = "grey", xlab="Average log2-expression",
                   ylab= "log2 FC", main = "SUC2 vs WT")
abline(h=0, col="darkred", lwd=2)

#################################################
###############Contraste co10 vs wt#############
################################################
co10.wt <- topTable(contrast.results, number=nrow(contrast.results),coef=2,sort.by="logFC")
# co10.wt <- topTable(contrast.results, number=33602,coef=2,sort.by="logFC", p.value = 0.05)
head(co10.wt)
co10.wt["AT5G15840",]
co10.wt["AT5G24470",]

fold.change.co10.wt <- co10.wt$logFC
p.value <- co10.wt$adj.P.Val
genes.ids.co10.wt <- rownames(co10.wt)

activated.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt > log2(2)]
repressed.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt < -log2(2)]
# activated.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt > log2(1.5) & p.value < 0.05 ]
# repressed.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt < -log2(1.5) & p.value < 0.05]

length(activated.genes.co10)
length(repressed.genes.co10)

# write.table(activated.genes.co10, file="tablas/activated_genes_co10_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(repressed.genes.co10, file="tablas/repressed_genes_co10_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

##### Scatterplot co10 vs WT #####
ggplot(data = as.data.frame(mean.counts.2.4), aes(x=Col0, y=co10)) + geom_point() + theme_minimal()

# add a column of NAs
mean.counts.2.4$diffexpressed <- "NO"
co10.wt$diffexpressed <- "NO"
# if log2Foldchange > 1  set as "UP" 
for (i in 1:nrow(co10.wt))
{
  if (co10.wt[i,"logFC"]>1)
  {
    up.gene <- rownames(co10.wt[i,])
    mean.counts.2.4[up.gene,"diffexpressed"] <- "UP"
    co10.wt[up.gene,"diffexpressed"] <- "UP"
  }
}
# if log2Foldchange < -1  set as "DOWN"
for(i in 1:nrow(co10.wt))
{
  if (co10.wt[i,"logFC"] < -1)
  {
    down.gene <- rownames(co10.wt[i,])
    mean.counts.2.4[down.gene,"diffexpressed"] <- "DOWN"
    co10.wt[down.gene,"diffexpressed"] <- "DOWN"
  }
}
# Set colors
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

p <- ggplot(data = as.data.frame(mean.counts.2.4), aes(x=Col0, y=co10, col=diffexpressed)) + geom_point() + theme_minimal() +scale_color_manual(values=mycolors)
p + ggtitle("Scatterplot") + xlab("log2(CPM) Col-0") + ylab("log2(CPM) co-10")



## Quick Volcano plot 
gene.names <- rownames(contrast.results$coefficients)
volcanoplot(contrast.results, coef=2, names = gene.names, highlight = 5)

##MA plot
# It is used to check whether data are comparable among groups
# (i.e. normalization worked properly). This can be created 
# using a function available in limma
limma::plotMA(contrast.results, coef=2, main="co10 vs WT")
plotWithHighlights(co10.wt$AveExpr, co10.wt$logFC,
                   status=co10.wt$diffexpressed,
                   hl.col = c("red", "blue"), hl.cex = 0.5,
                   bg.col = "grey", xlab="Average log2-expression",
                   ylab= "log2 FC", main = "co-10 vs WT")
abline(h=0, col="darkred", lwd=2)


##############################################
##########      Volcano plots       ###########
##############################################
## 35SCO vs WT
volcano <- ggplot(data=x35SCO.wt, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()

volcano + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") 

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
x35SCO.wt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
x35SCO.wt$diffexpressed[x35SCO.wt$logFC > log2(2) & x35SCO.wt$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
x35SCO.wt$diffexpressed[x35SCO.wt$logFC < -log2(2) & x35SCO.wt$P.Value< 0.05] <- "DOWN"

#Set colors
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

#Add a color in aes section
volcano <- ggplot(data=x35SCO.wt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal()

volcano + geom_vline(xintercept=c(-log2(2), log2(2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  scale_color_manual(values=mycolors) 

x35SCO.wt$agi <- NA
x35SCO.wt$agi[x35SCO.wt$diffexpressed != "NO"] <- rownames(x35SCO.wt)[x35SCO.wt$diffexpressed != "NO"]

ggplot(data=x35SCO.wt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=agi)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 50) +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-log2(2), log2(2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


## suc2 vs WT
volcano <- ggplot(data=suc2.wt, aes(x=logFC, y=-log10(P.Value))) + geom_point() + theme_minimal()

volcano + geom_vline(xintercept=c(-log2(2), log2(2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
suc2.wt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
suc2.wt$diffexpressed[suc2.wt$logFC > log2(2) & suc2.wt$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
suc2.wt$diffexpressed[suc2.wt$logFC < -log2(2) & suc2.wt$P.Value < 0.05] <- "DOWN"

#Set colors
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

#Add a color in aes section
volcano <- ggplot(data=suc2.wt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal()

volcano + geom_vline(xintercept=c(-log2(2), log2(2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  scale_color_manual(values=mycolors)

suc2.wt$agi <- NA
suc2.wt$agi[suc2.wt$diffexpressed != "NO"] <- rownames(suc2.wt)[suc2.wt$diffexpressed != "NO"]

ggplot(data=suc2.wt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=agi)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-log2(1.2), log2(1.2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


## co10 vs WT
volcano <- ggplot(data=co10.wt, aes(x=logFC, y=-log10(P.Value))) + geom_point() + theme_minimal()

volcano + geom_vline(xintercept=c(-log2(2), log2(2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
co10.wt$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
co10.wt$diffexpressed[co10.wt$logFC > log2(2) & co10.wt$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
co10.wt$diffexpressed[co10.wt$logFC < -log2(2) & co10.wt$P.Value < 0.05] <- "DOWN"

#Set colors
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

#Add a color in aes section
volcano <- ggplot(data=co10.wt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal()

volcano + geom_vline(xintercept=c(-log2(2), log2(2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  scale_color_manual(values=mycolors)

co10.wt$agi <- NA
co10.wt$agi[co10.wt$diffexpressed != "NO"] <- rownames(co10.wt)[co10.wt$diffexpressed != "NO"]

ggplot(data=co10.wt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=agi)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-log2(1.2), log2(1.2)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")



#########################################
#######        HEATMAP      #############
#########################################
# Draw heatmap for DEGS
head(voom.counts) #counts transformed by voom 

degs <- c(activated.genes.35SCO, repressed.genes.35SCO,
          activated.genes.suc2, repressed.genes.suc2,
          activated.genes.co10, repressed.genes.co10)
length(degs)
degs <- unique(degs)
length(degs)

degs.table <- voom.counts[degs,]
colnames(degs.table) <- colnames(voom.counts)
coldata

tiff("images/heatmap_degs.tiff", height = 4, width = 8,
     units = 'in', res=400, compression="lzw")
pheatmap(as.matrix(degs.table), cluster_rows = T, cluster_cols = T,
         scale = "row", clustering_distance_rows = "correlation", 
         clustering_method = "complete", #annotation_col = coldata, 
         main="DEGs",fontsize_col=14, fontsize_row = 4, color = magma(28),
         show_rownames = F)
dev.off()


## Heatmap for selected genes

# selected.genes <- read.table(file="selected_genes/list.txt", sep = "\t", header = T)
# head(selected.genes)
# 
# up.selected <- selected.genes$up
# length(up.selected)
# up.selected <- up.selected[1:22]
# up.table <- filtered.data[as.vector(up.selected),]
# tiff("selected_genes/genes_up_cluster.tiff", height = 4, width = 6, 
#      units = 'in', res=300, compression="lzw")
# pheatmap(as.matrix(up.table), cluster_rows = F, cluster_cols = T,
#          scale = "row", clustering_distance_rows = "correlation", 
#          clustering_method = "complete", #annotation_col = pheno.data, 
#          main="DEGs",fontsize_col=14, fontsize_row = 6, color = greenred(28))
# dev.off()
# 
# down.selected <- selected.genes$down
# length(down.selected)
# down.table <- filtered.data[as.vector(down.selected),]
# tiff("selected_genes/genes_down_cluster.tiff", height = 4, width = 6, 
#      units = 'in', res=300, compression="lzw")
# pheatmap(as.matrix(down.table), cluster_rows = F, cluster_cols = T,
#          scale = "row", clustering_distance_rows = "correlation", 
#          clustering_method = "complete", #annotation_col = pheno.data, 
#          main="DEGs",fontsize_col=14, fontsize_row = 6, color = greenred(28))
# dev.off()



##########################################
#######         PCA plot         #########
##########################################

## Draw PCA plot
# transpose the data and compute principal components
# pca.data <- prcomp(t(degs.table)) #Only for DEGs

# Log transform the data before performing PCA
# the log-transform is used to reduce the influence 
# of extreme values or outliers. 
pca.data <- prcomp(t(voom.counts))
nrow(pca.data$rotation)

# Calculate PCA component percentages
pca.data.perc <- round(100*pca.data$sdev^2/sum(pca.data$sdev^2),1)

# Extract 1 and 2 principle components and create a data frame with sample names, first and second principal components and group information
df.pca.data <- data.frame(PC1 = pca.data$x[,1], PC2 = pca.data$x[,2], sample = colnames(sample.counts), 
                          condition = c(rep("35SCO",3), rep("Col0",3), rep("suc2",3), rep("co10",3)))
head(df.pca.data)
# tissue <- c("root", "air", "air")
# tissue <- rep(tissue, 3)
# df.pca.data$tissue <- tissue

# color by sample
ggplot(df.pca.data, aes(PC1,PC2, color = sample))+
  geom_point(size=8)+ 
  labs(x=paste0("PC1 (",pca.data.perc[1],")"), y=paste0("PC2 (",pca.data.perc[2],")")) 


# color by condition/group
# tiff("images/PCA.tiff", height = 4, width = 6,
# units = 'in', res=600, compression="lzw")
ggplot(df.pca.data, aes(PC1,PC2, color = condition))+
  geom_point(size=8)+
  labs(x=paste0("PC1 (",pca.data.perc[1],")"), y=paste0("PC2 (",pca.data.perc[2],")"))+
  geom_text_repel(aes(label=sample),point.padding = 0.75) 
# stat_ellipse(geom = "polygon", aes(fill=condition), alpha=0.2, show.legend = FALSE, level=0.95)
# dev.off()



#########################################################
###### Comparing gene expression within samples     #####
#########################################################

# We need to convert counts to tpm, because this unit is normalized
# by fragment length.

# TPM is very similar to RPKM and FPKM. The only difference is the order 
# of operations. Here’s how you calculate TPM:
#   
# * Divide the read counts by the length of each gene in kilobases. 
#   This gives you reads per kilobase (RPK).
#   
# * Count up all the RPK values in a sample and divide this number by 1,000,000. 
#   This is your “per million” scaling factor. 
# 
# * Divide the RPK values by the “per million” scaling factor. This gives you TPM.
# 
# So you see, when calculating TPM, the only difference is that you normalize for gene 
# length first, and then normalize for sequencing depth second. However, the effects 
# of this difference are quite profound. When you use TPM, the sum of all TPMs in 
# each sample are the same. This makes it easier to compare the proportion of reads 
# that mapped to a gene in each sample. In contrast, with RPKM and FPKM, the sum of 
# the normalized reads in each sample may be different, and this makes it harder to 
# compare samples directly. Here’s an example. If the TPM for gene A in Sample 1 
# is 3.33 and the TPM in sample B is 3.33, then I know that the exact same proportion 
# of total reads mapped to gene A in both samples. This is because the sum of the TPMs 
# in both samples always add up to the same number (so the denominator required to 
# calculate the proportions is the same, regardless of what sample you are looking at.)
# With RPKM or FPKM, the sum of normalized reads in each sample can be different. Thus, 
# if the RPKM for gene A in Sample 1 is 3.33 and the RPKM in Sample 2 is 3.33, I would 
# not know if the same proportion of reads in Sample 1 mapped to gene A as in Sample 2. 
# This is because the denominator required to calculate the proportion could be different
# for the two samples.

# TPMs are convenient when you want to add samples over time and change groups
# and samples being visually compared. The results of that will not be as 
# robust to outliers as normalized counts (produced with TMM or another method, like voom), 
# but they're usually GOOD ENOUGH FOR VISUALIZATIONS.

## Computing TPM
kb.length <- counts$Length/1000
rpk <- sample.counts/kb.length
head(sample.counts)
head(rpk)

scaling.factors <- colSums(rpk)/1e6

tpm <- matrix(nrow = nrow(rpk), ncol = ncol(rpk))
for (i in 1:ncol(rpk))
{
  tpm[,i] <- rpk[,i]/scaling.factors[i]
}

colnames(tpm) <- colnames(sample.counts)
rownames(tpm) <- rownames(sample.counts)
head(tpm)

sample.counts["AT5G15840",]
sample.counts["AT1G65480",]
tpm["AT5G15840",]
tpm["AT1G65480",]
# write.table(tpm, file = "expression/log_TPM.txt", sep = "\t", quote = F, col.names = NA)

# IMPORTANT NOTE: The "mean fragment length" has to be taken in account for 
# PAIRED-END data. But It is not possible to estimate fragment length from single-end 
# sequencing data. https://gist.github.com/slowkow/c6ab0348747f86e2748b
# For paired-end sequencing data, the mean fragment length can be obtained 
# from Picard using InsertSizeMetrics.

# The mean fragment length is used to compute EFFECTIVE LENGTHS of feature in each 
# library, that is different to annotated feature length.
# https://www.biostars.org/p/253789/


## Plotting the expression of some genes
# Duda, qué usamos para representar la expresión de los genes?
# el resultado de la transformación de voom? que no es exactamente log2(cpm)
# O log2(tpm + 1)
# O el log2CPM +1

# nothing wrong with a negative value if it was generated by 
# log-transforming a positive value. The CPM value must have been 
# less than one.

## Compute CPM with cpm function (library size and norm factors)
cpm <- cpm(d) ## Here the norm factors are used to do TMM normalisation
cpm["AT5G15840",]
# write.table(log2(cpm+1),file = "expression/log_CPM.txt", sep = "\t", quote = FALSE, col.names = NA)

# E-values from voom are calculated in the following way ->
# t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
#Los valores de voom tb son log(cpm) pero calculado
# un poco distinto de como lo hace la función cpm (edgeR). https://support.bioconductor.org/p/59846/#59917
voom.counts["AT5G15840",]
# write.table(voom.counts,file = "expression/log_CPM_voom.txt", sep = "\t", quote = FALSE, col.names = NA)


# Here you can choose the gene to plot and the unit to plot
gene.to.plot <- "AT1G22770"
gene.alias <- "GI"
# unit.to.plot <- "log2 (TPM + 1)" #for TPM
# unit.to.plot <- "log2 (CPM + 1)" #for CPM
unit.to.plot <- "log2 normalized counts" #for voom counts


gene.expression <- data.frame(expression=c(voom.counts[gene.to.plot,]),
                              # gene.expression <- data.frame(expression=log2(c(tpm[gene.to.plot,]+1)),
                              # gene.expression <- data.frame(expression=log2(cpm[gene.to.plot,]+1),
                              genotype= genotype)

sd.35sco <- sd(subset(gene.expression, genotype == "x35SCO")$expression)
sd.col0 <- sd(subset(gene.expression, genotype == "col0")$expression)
sd.co10 <- sd(subset(gene.expression, genotype == "co10")$expression)
sd.suc2 <- sd(subset(gene.expression, genotype == "suc2")$expression)


mean.35sco <- mean(subset(gene.expression, genotype == "x35SCO")$expression)
mean.col0 <- mean(subset(gene.expression, genotype == "col0")$expression)
mean.co10 <- mean(subset(gene.expression, genotype == "co10")$expression)
mean.suc2 <- mean(subset(gene.expression, genotype == "suc2")$expression)


gene.mean <- data.frame(mean=c(mean.35sco, mean.suc2, mean.col0, mean.co10),
                        sd=c(sd.35sco, sd.suc2, sd.col0, sd.co10 ),
                        genotype=c("35SCO", "SUC2", "Col-0", "co10"))
head(gene.mean)

# tiff(paste0("barplots/", gene.alias,".tiff"), height = 8, width = 10,
#      units = 'in', res=300, compression="lzw")
ggplot(data = gene.mean, aes(x=genotype, y=mean, fill=genotype)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) + 
  # scale_fill_brewer(palette="Dark2")
  scale_fill_manual(values=met.brewer(n=4, name="Degas")) +
  ggtitle(gene.alias) +
  xlab("") + 
  ylab(unit.to.plot)
# dev.off()


#########################################################
############   Stripchart   #############################
#########################################################

# We can quickly look at grouped expression using stripchart. 
# We can use the normalised log expression values in the voom 
# object (y$E).

stripchart(y$E["AT1G22770",]~group, vertical = TRUE)

gene.to.plot <- "AT1G22770"
gene.alias <- "GI"
nice.col <- met.brewer(6,name="Degas")
stripchart(y$E[gene.to.plot,]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main=gene.alias)


###########################################################
############        Compute FPKMs            ##############
###########################################################

dds <- DESeqDataSetFromMatrix(countData = sample.counts,
                              colData = coldata,
                              design = ~ condition)
dds
mcols(dds)$basepairs <- counts$Length #Add fragment length
dds <- estimateSizeFactors(dds) #estimate the size factors using the "median ratio method"
gene.expression <- as.matrix(fpkm(dds)) #compute fpkm

head(gene.expression,5)

colnames(gene.expression) <- samples 

## Perform some checks
# wrong.genes <- subset(as.data.frame(gene.expression), co10_1 == 0 & co10_2 == 0 & co10_3==0)
# nrow(wrong.genes)
zero.genes <- subset(as.data.frame(gene.expression), co10_1 == 0 & co10_2 == 0 & co10_3==0)
nrow(zero.genes)
zero.genes[3000:3050,]
zero.genes[2000:2050,]
zero.genes[5000:5100,]

zero.genes <- subset(as.data.frame(gene.expression), "35SCO_1" == 0 & "35SCO_2" == 0 & "35SCO_3"==0)
nrow(zero.genes)
zero.genes[3000:3050,]
zero.genes[2000:2050,]
zero.genes[5000:5100,]

zero.genes <- subset(as.data.frame(gene.expression), "Col-0_1" == 0 & "Col-0_2" == 0 & "Col-0_3" ==0)
nrow(zero.genes)
zero.genes[3000:3050,]
zero.genes[2000:2050,]
zero.genes[5000:5100,]


## Some interesting genes
gene.expression["AT5G15840",] #co
gene.expression["AT1G22770",] #gi
gene.expression["AT5G24470",] #prr5
gene.expression["AT1G65480",] #ft
gene.expression["AT1G01060",] #lhy
gene.expression["AT1G06530",]

head(gene.expression)

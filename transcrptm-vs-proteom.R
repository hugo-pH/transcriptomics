##transcriptm-vs-proteom.R
##Author: Hugo Pineda
##Create the plot compairing proteome and transcriptome log ratios over two groups.

library(limma)
library(reshape2)
library(dplyr)
library(ggplot2)



#Read the transcriptomics data
design <- read.delim("./design.txt")
design <- filter(design, Dye=="Cy3")
data.rma.full <- read.delim("./Expression.RMA.txt")
annot <- data.rma.full[,c(1:13)]
data.rma <- as.matrix(data.rma.full[,-c(1:13)])
groups <- design$SampleID

#Select groups to perform the tests and calculate logFC. Change these groups if you want to compare other samples.
#Select sample 1 and 5 as aerobes
group1<-as.list(1, 5)
#Select sample 4 which correponds to the proteomics sample. The results do not apparently change if you select all the anaerobic samples (2, 4, 6)
group2<-as.list(4)

#Get the columns of the selected samples from data.rma
idx1 <- which(groups %in% group1)
idx2 <- which(groups %in% group2)
x <- data.rma[,c(idx1,idx2)]

# The following code will perform a t-test, with Empirical Bayes correction and multiple testing correction:
factorlevels <- c(rep(1,length(idx1)),rep(2,length(idx2)))
design.mat <- model.matrix(~ 0+factor(factorlevels))
colnames(design.mat) <- c("group1", "group2")
fit <- lmFit(x, design.mat)
contrast.matrix <- makeContrasts(group2-group1, levels=design.mat)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tblPval <- topTable(fit2, adjust="BH",number=nrow(data.rma))
#Now merge the results with the annot df. It will be useful later to have this information in the same df. 
#First create a row in both dataframes containing the row name.
annot$idx<-as.vector(rownames(annot), mode="numeric")
tblPval$idx<-as.vector(rownames(tblPval), mode="numeric")
#Merge both df by the new column
Pval.annot<-merge(annot, tblPval, by="idx")
#Sort the df in decreasing order by the corrected pvalue 
Pval.annot<-arrange(Pval.annot, adj.P.Val)



#Read proteome file. You have to provide the path and name to your proteome data.
#CAUTION HERE. Be sure that your data is separated by tabs or choose the correct value for the "sep" argument
#Do not modify the column names of the original file, otherwise it won't work
proteom_raw<-read.delim("PATH TO FILE HERE", sep = "\t")
#Get the GeneName from the ID column of proteome data which is in this form: "sw|P00331|ADH2_YEAST"
#Remove the first part before the GeneName "sw|P00331|"
proteom_raw$ID<-gsub(".*\\|", "", proteom_raw$ID)
#Remove the string after the GeneNanem "_YEAST"
proteom_raw$ID<-gsub("_YEAST", "", proteom_raw$ID)
#Create a new column containing the log of l.h column
proteom<-mutate(proteom_raw, logproteome = log2(l.h))
#Select columns from Pval.annot (GeneName, logFC)
transcrip<-select(Pval.annot, GeneName, logFC)
#Merge both df by GeneName
df.plot<-merge(transcrip, proteom, by.x="GeneName", by.y="ID")

#Plot the results
p<-ggplot(df.plot, aes(x=logproteome, y=logFC)) +
  #Add axis
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  #Add the labels for GeneNames
  geom_text(aes(label=GeneName), size=4)
p

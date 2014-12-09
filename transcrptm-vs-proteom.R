##transcriptm-vs-proteom.R
##Author: Hugo Pineda
##Create the plot compairing proteome and transcriptome log ratios over two groups.

library(limma)
library(reshape2)
library(dplyr)
library(ggplot2)
library(pathview)


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
#Select columns of interest
transcrip<-dplyr::select(Pval.annot, Name, GeneName, logFC)



#Read proteome file. You have to provide the path and name to your proteome data.
#CAUTION HERE. Be sure that your data is separated by tabs or choose the correct value for the "sep" argument
#The .csv should contain TWO COLUMNS, the first is the ID column from the excel file.
#The second column is "l.h" column from the excel file
#Column names should be "ID" and "l.h"
proteom_raw<-read.delim("PATH TO FILE", sep = "\t")
#Remove NAs from raw data
proteom<-filter(proteom_raw, l.h != "NA")
#Get the Uniprot ID from the ID column of proteome data which is in this form: "sw|P00331|ADH2_YEAST"
#This regular expression returns "P00331" from "sw|P00331|ADH2_YEAST"
proteom$uniprot<-gsub("(^.+\\|)(.+)(\\|.+$)", "\\2", proteom$ID)
#Map ID from Uniprot to Entrez. "id2eg" is a function from pathview package
proteom$entrez <- id2eg(ids = proteom$uniprot, category = gene.idtype.list[7], pkg.name = "org.Sc.sgd.db")[,2]
#Get the logarithm of l.h, remove columns returning "-Inf" and select columns of interest
proteom<-mutate(proteom, logproteome = log2(l.h))
proteom<-filter(proteom, logproteome != "-Inf" )
proteom<-dplyr::select(proteom, entrez, logproteome)



#Merge both df by GeneName
df.plot<-merge(transcrip, proteom, by.x="Name", by.y="entrez")

# nor.trans<-ad.test(df.plot$logFC)
# nor.prot<-ad.test(df.plot$logproteome)
# cor<-cor.test(df.plot$logFC, df.plot$logproteome)
# p.value<-cor["p.value"]

#Plot the results
  p<-ggplot(df.plot, aes(x=logproteome, y=logFC)) +
    scale_x_continuous(breaks=-6:6)+
    scale_y_continuous(breaks=-5:5)+
    
    #Add axis
    geom_hline(aes(yintercept=0)) +
    geom_vline(aes(xintercept=0)) +
    theme(axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18)) +
    geom_smooth(method="lm") +
    #Add the labels for GeneNames
    geom_text(aes(label=GeneName), size=8) +
#     annotate("text" , x=-5, y=c(2, 1.5), label=c("r = 0.66", "p.value = 2.2e-16"), size=8, , color="#013ADF") + 
    #add axis names
    xlab(expression(log2(frac(-O[2] , +O[2]))*" Proteins")) +
    ylab(expression(log2(frac(-O[2] , +O[2]))*" Transcripts"))
  
  p

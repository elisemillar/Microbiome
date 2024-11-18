#
#title: "Invert Taxa and Water Microbiome"
#author: "Elise Millar"
#date: "Aug 27th, 2020"
#output: html_document
#

install.packages("rlang")
install.packages("phytools")

#set working directory
setwd("C:/Users/Home/Desktop/Fall_Samples/2020 Aug 26 Millar Masters Data + Code")
path ="C:/Users/Home/Desktop/Fall_Samples/2020 Aug 26 Millar Masters Data + Code"
getwd()

# load libraries 
library("cluster")
library("phyloseq")
packageDescription("phyloseq")$Version
library("ggplot2")
library("plyr")
library("grid")
library("ape")
library("phangorn")
library("phytools")
library("vegan")


#set default to white
theme_set(theme_bw()) 

#create otu data frame
otufile = "2018-19_InvertMicrobiomeWWTP_seqv34_v1.csv"
otu_df = read.csv(otufile, row.names = 1)   
seqs = rownames(otu_df)
rownames(otu_df) = NULL

#create taxa data frame
taxfile = "2018-19_InvertMicrobiomeWWTP_taxav34_v1.csv"
tax_df = read.csv(taxfile, row.names = 1)

#rownames(tax_df)=tax_df$X.TAXONOMY #assign contents of first column to rownames
#tax_df=tax_df[,-1] #remove first column
all(seqs == rownames(tax_df))
rownames(tax_df) = NULL

#create metadata data frame
metfile = "2018-19_InvertMicrobiomeWWTP_meta_v1.txt"
met_df = read.delim(metfile) # read tab-delimited instead of comma
#met_df = read.csv(metfile)
met_df$Sample_ID = gsub('-', '.', as.character(met_df$Sample_ID)) #in case there are dashes
met_df$Sample_ID = gsub(' ', '.', met_df$Sample_ID) #in case there are spaces
met_df$Sample_ID = as.factor(met_df$Sample_ID) #in case
rownames(met_df) = met_df$Sample_ID   #make sure rownames in metadata data frame are same as row (or column, whichever is samples) in otu data frame

#create phyloseq object
data = phyloseq(otu_table(otu_df, taxa_are_rows = TRUE), # or FALSE if false
                  tax_table(as.matrix(tax_df)),
                  sample_data(met_df))
data #52193 ASVs, 389 samples

#Subset to only the samples (no wash, no negatives, no VR rainbow darter)
data <- subset_samples(data, Sample_Label=="Perlidae" | Sample_Label=="Baetidae" | Sample_Label=="Ephemerellidae" | 
                Sample_Label=="Heptageniidae" | Sample_Label=="Hydropsychidae" | Sample_Label=="Spiders" | 
                  Sample_Label=="Mussels" | Sample_Label=="Water" | Sample_Label=="K-WWTP" | Sample_Label=="W-WWTP")
data <- prune_taxa(taxa_sums(data) >0, data)
data #38033 ASVs, 389 samples

#Subset to only the samples (only inverts, no mussel, no water)
inverts <- subset_samples(data, Sample_Label=="Perlidae" | Sample_Label=="Baetidae" | Sample_Label=="Ephemerellidae" | 
                         Sample_Label=="Heptageniidae" | Sample_Label=="Hydropsychidae" | Sample_Label=="Spiders")
inverts <- prune_taxa(taxa_sums(inverts) >0, inverts)
inverts #24008 ASVs, 304 samples

#################################################### Beta Diversity ##################################################################

newlabelorder_Location = c("Upstream", "Waterloo_WWTP", "Downstream_Waterloo", "Kitchener_WWTP", "Downstream_Kitchener")
newlabelname_Location = c("Upstream", "W-WWTP", "Downstream Waterloo", "K-WWTP", "Downstream Kitchener")

newlabelorder_Sample = c("Spiders", "Mussels", "Baetidae", "Ephemerellidae", "Heptageniidae", "Perlidae", "Hydropsychidae", 
                         "Water", "W-WWTP", "K-WWTP")
newlabelname_Sample = c("Spiders", "Mussels", "Baetidae", "Ephemerellidae", "Heptageniidae", "Perlidae", "Hydropsychidae", 
                         "Water", "W-WWTP", "K-WWTP")

#colours
library(scales)
show_col(hue_pal()(10))

SampleColours = c(Baetidae="#F8766D", Heptageniidae="#A3A500", `K-WWTP`="#D89000", Ephemerellidae="#00B0F6", 
                  Mussels="#E76BF3", `W-WWTP`="#39B600", Hydropsychidae="#00BF7D", Spiders="#9590FF", Perlidae="#FF62BC", 
                  Water="#00BFC4") 

# Transform data to proportions for Bray-Curtis distances
trans.dat <- transform_sample_counts(data, function(otu) otu/sum(otu))

# Data ordination
bray.ord = ordinate(trans.dat, method = "PCoA", distance = "bray")
bray.ord

# Plot Bray-Curtis + ellipses
bray.plot = plot_ordination(trans.dat, bray.ord, color="Sample_Label", title="Bray-Curtis All Samples", axes = c(1,2))
bray.plot$data$Sample_Label <- as.character(bray.plot$data$Sample_Label)
bray.plot$data$Sample_Label <- factor(bray.plot$data$Sample_Label, levels=newlabelorder_Sample, labels=newlabelname_Sample)
bray.plot + scale_color_manual(values = SampleColours, name="Sample")  # + stat_ellipse(geom = "polygon", type="norm", alpha=0.2, aes(fill=Sample)) + scale_fill_manual(values=colourblind2)


############################################## Beta Diversity Matrices #############################################################

#Adonis (PERMANOVA) for Bray Curtis matrices
#Bray-Curtis distance matrix
bray.test = phyloseq::distance(trans.dat, "bray")

#Transformed dataframe
sample.df = data.frame(sample_data(trans.dat))

adonis(formula = bray.test~Sample_Label, data =  sample.df, permutations = 9999, method = "bray")

#Pairwise adonis for significant pairs
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pairwise.adonis(bray.test, sample_data(data)$Sample_Label, perm = 99999)

####################################################### DESeq2 for families #####################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")
packageVersion("DESeq2")

##REF and DSW (Deseq2 compares 2 things so compared by location; reference and downstream Waterloo)

REF_DSW = subset_samples(inverts, Location_vs_WWTP=="Downstream_Waterloo" | Location_vs_WWTP=="Upstream")
REF_DSW <- prune_taxa(taxa_sums(REF_DSW) >0, REF_DSW)
REF_DSW #17759 ASVs, 183 samples

deseq2.REF_DSW <- phyloseq_to_deseq2(REF_DSW, ~Location_vs_WWTP)
gm_mean <- function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans <- apply(counts(deseq2.REF_DSW), 1, gm_mean)
deseq2.REF_DSW <- estimateSizeFactors(deseq2.REF_DSW, geoMeans = geoMeans)
deseq2.REF_DSW <- DESeq(deseq2.REF_DSW, test="Wald", fitType="parametric")
deseq2.REF_DSW

res = results(deseq2.REF_DSW, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dat_ref_dsw)[rownames(sigtab), ], "matrix"))
head(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
  theme(legend.text=element_text(size=15, color="black")) + 
  theme(legend.title = element_text(size=15, color="black")) + facet_wrap(~Sample_Label, scales="free")



##REF and DSK (Deseq2 compares 2 things so compared by location; reference and downstream Kitchener)

dat_ref_dsk = subset_samples(dat_lessKK, Location != "Downstream_Waterloo")
dat_ref_dsk

dat_deseq2_ref_dsk <- phyloseq_to_deseq2(dat_ref_dsk, ~Location)
dat_deseq2_ref_dsk <- DESeq(dat_deseq2_ref_dsk, test="Wald", fitType="parametric")
dat_deseq2_ref_dsk

res = results(dat_deseq2_ref_dsk, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dat_ref_dsk)[rownames(sigtab), ], "matrix"))
head(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
  theme(legend.text=element_text(size=15, color="black")) + 
  theme(legend.title = element_text(size=15, color="black"))



##DSW and DSK (Deseq2 compares 2 things so compared by location; downstream Waterloo and downstream Kitchener)
dat_dsw_dsk = subset_samples(dat_lessKK, Location != "Upstream")
dat_dsw_dsk

dat_deseq2_dsw_dsk <- phyloseq_to_deseq2(dat_dsw_dsk, ~Location)
dat_deseq2_dsw_dsk <- DESeq(dat_deseq2_dsw_dsk, test="Wald", fitType="parametric")
dat_deseq2_dsw_dsk

res = results(dat_deseq2_dsw_dsk, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dat_dsw_dsk)[rownames(sigtab), ], "matrix"))
head(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))










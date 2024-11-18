#
#title: "Effluent Bacteria"
#author: "Elise Millar"
#date: "Jan 10th, 2022"
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

#Subset (no wash, no negatives, no VR rainbow darter, no effluents)
data <- subset_samples(data, Sample_Label=="Perlidae" | Sample_Label=="Baetidae" | 
                         Sample_Label=="Heptageniidae" | Sample_Label=="Hydropsychidae" | Sample_Label=="Spiders" | 
                         Sample_Label=="Mussels" | Sample_Label=="Water")
data <- prune_taxa(taxa_sums(data) >0, data)
data #35375 ASVs, 373 samples

#Downstream sites
ups <- subset_samples(rivdat, Location_vs_WWTP=="Upstream")
ups <- prune_taxa(taxa_sums(ups) >0, ups)
ups



effluent = c('A7P-90m', 'Acetobacterium', 'Acetonema', 'Alysiella', 'Azoarcus', 'CAG-352', 
'Candidatus_Caldatribacterium', 'Candidatus_Fritschea', 'Candidatus_Latescibacter', 'Candidatus_Microthrix', 
'Candidatus_Xiphinematobacter', 'Castellaniella', 'Citrobacter', 'Enhydrobacter', 'Galbitalea', 'JL-ETNP-Z34', 
'KCM-B-112', 'Kouleothrix', 'Lachnoclostridium_12', 'Laribacter', 'Methanoculleus', 'Methanolinea', 
'Methanosaeta', 'Methyloversatilis', 'Moraxella', 'Nubsella', 'OLB17', 'Pelotomaculum', 'Uruburuella', 'W5', 
'Zooshikella', 'Acidaminococcus', 'Actinomyces', 'BD1-7_clade', 'Bifidobacterium', 'Candidatus_Paenicardinium', 
'Chiayiivirga', 'Cloacibacterium', 'Flavitalea', 'Fusibacter', 'Lelliottia', 'Ottowia', 'Prevotellaceae_UCG-004', 
'Steroidobacter', 'SWB02', 'Acetoanaerobium', 'Aquimonas', 'AUTHM297', 'Bact-08', 'C1-B045', 
'Candidatus_Cloacimonas', 'Candidatus_Protochlamydia', 'Desulfobacter', 'Leptotrichia', 'Mesotoga', 
'Neochlamydia', 'Planctopirus', 'Proteiniclasticum', 'SC103', 'Succinivibrio', 'Thermovirga', 'Turneriella', 
'U29-B03', 'XBB1006')
effgen = subset_taxa(data, (Genus %in% effluent))
effgen <- prune_samples(sample_sums(effgen) >0, effgen)
efflist = psmelt(effgen)
efflist[1:109, c(2,3,5,13,33)]

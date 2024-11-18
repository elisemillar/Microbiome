#
#title: "Mussel Microbiome"
#author: "Elise Millar"
#date: "Aug 27th, 2020"
#output: html_document
#

install.packages("rlang")
packageVersion("rlang")
install.packages("phytools")

#set working directory
setwd("C:/Users/Home/Desktop/2020 Aug 26 Millar Masters Data + Code")
path ="C:/Users/Home/Desktop/2020 Aug 26 Millar Masters Data + Code"
getwd()

# load libraries 
library("cluster")
library("phyloseq")
packageDescription("phyloseq")$Version
library("ggplot2")
library("plyr")
library("dplyr")
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
data

########################################################
#Phylum - sequence read counts
negdat <- subset_samples(data, Sample_Label=="Mussels" | Sample_Label=="Neg-1" | 
                           Sample_Label=="Neg-2")
negdat <- prune_taxa(taxa_sums(negdat) >0, negdat)
negdat #10668 ASVs, 45 samples

dat.aglo = tax_glom(negdat, taxrank = "Phylum")
dat.dataframe = psmelt(dat.aglo)
dat.agr = aggregate(Abundance~Sample_ID+Phylum, data=dat.dataframe, FUN=mean)

ggplot(dat.agr, aes(x=reorder(Sample_ID, -Abundance), y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity")+ scale_fill_manual(values=phylaColours1) +
  labs(y= "Sequence Read Count", x= "Sample ID") + 
  theme(axis.text.x = element_text(angle = 90))

#################

#Subset to only ~Mussel~ samples
musdat <- subset_samples(data, Sample_Label=="Mussels")
musdat <- prune_taxa(taxa_sums(musdat) >0, musdat)
musdat #10664 ASVs, 43 samples

##################################################### Different Subsets #############################################################

#Samples containing WWTP Effluent Genera
effluent = c('Ottowia', 'Fusibacter', 'Bifidobacterium','Actinomyces', 'Prevotellaceae_UCG-004', 'Cloacibacterium', 
             'Chiayiivirga', 'Lelliottia', 'BD1-7_clade', 'Steroidobacter', 'Acidaminococcus', 
             'Candidatus_Paenicardinium')
effgen = subset_taxa(musdat, (Genus %in% effluent))
effgen <- prune_samples(sample_sums(effgen) >0, effgen)
efflist = psmelt(effgen)
efflist[1:30, c(1:4,33)]

#By Phylum (Ex. Cyanobacteria)
cyano = subset_taxa(musdat, Phylum=="Cyanobacteria")
cyano

#By Site (Ex. REF1)
REF1 <- subset_samples(musdat, Site_ID=="REF1")
REF1 <- prune_taxa(taxa_sums(REF1) >0, REF1)
REF1

################################################ How many taxonomic levels ###############################################################
#Excluding NA taxa

get_taxa_unique(musdat, "Phylum") #41
get_taxa_unique(musdat, "Class") #86
get_taxa_unique(musdat, "Order") #201
get_taxa_unique(musdat, "Family") #303
get_taxa_unique(musdat, "Genus") #724

#################################################### Abundance Numbers ############################################################################

##Total number of sequence reads for dataset (including unassigned)
sum(sample_sums(musdat)) #1852118
#Read sequence range for samples
sort(sample_sums(musdat)) #8654-63510

#Number of sequence reads per rank (Ex. Phylum)
top.taxa = sort(tapply(taxa_sums(musdat), tax_table(musdat)[, "Phylum"], sum), TRUE)

top.taxa #highest amount of reads to lowest
top.taxa[1:10] #top 10 taxa by read count
sum(top.taxa) #total reads assigned to taxa = 1829596

################
#Didn't need this code but could be useful 
#Reads per Species
taxa_sums(musdat)

##Number of clustered ASVs 
ntaxa(musdat)
tax_table(musdat)
top_taxa = table(tax_table(musdat)[, "Phylum"], exclude = NULL) #Ex. Phylum
top_taxa
sort(top_taxa, decreasing = TRUE)
head(sort(top_taxa, decreasing = TRUE))

#ASV counts per Species
asv_df <- t(otu_table(musdat))
colSums(asv_df != 0)
###############

###################################################### Alpha Diversity #########################################################################

#Sharokh's rarefy method
rare.depth <- min(sample_sums(musdat))
rare.depth #8654
dat_r <- rarefy_even_depth(musdat, sample.size = rare.depth, rngseed=1414) #rngseed = random value
dat_r

#Plot by Site_ID
alpha_plot = plot_richness(dat_r, x="Site_ID", measures=c("Shannon", "Simpson"))
newlabelorder_Site_ID = c("REF1", "REF2", "DSW1", "DSK1", "DSK2")
alpha_plot$data$Site_ID <- as.character(alpha_plot$data$Site_ID)
alpha_plot$data$Site_ID <- factor(alpha_plot$data$Site_ID, levels=newlabelorder_Site_ID)
alpha_plot

################################################# Alpha Diversity Metrics ###################################################################

alpha_diversity <- estimate_richness(dat_r, split=TRUE, measures=c("Shannon", "Simpson"))
alpha_div <- data.frame(alpha_diversity)
alpha_div

#Add Sample_ID column 
alpha_div["Sample_ID"] <- NA 
alpha_div$Sample_ID <- c(319:361) #exact range for Mussel samples in metadata

#Double-check
colnames(alpha_div)
alpha_div

#Merge metadata and alpha diversity numbers
alpha_met_df <- mutate(met_df, Sample_ID = c(219:280, 282:283, 285:361, 365:414, 415:570, 
                                             964:1007)) #ranges for ALL samples in metadata
alpha_met_df

ncol(alpha_div)
ncol(alpha_met_df)
as.character(alpha_met_df$Sample_ID)
as.character(alpha_div$Sample_ID)

alpha_div_merge <- merge(alpha_met_df, alpha_div, by = "Sample_ID") 
summary(alpha_div_merge)

#Confidence intervals
yo = subset(alpha_div_merge, select = c('Site_ID','Shannon'))
View(yo)

library(dplyr)

means = aggregate(.~ Site_ID, data = yo, mean)
View(means)
sdevs = aggregate(.~ Site_ID, data = yo, sd)
View(sdevs)

n <- 10
xbar <- 4.001276
s <- 1.3176467

margin <- qt(0.975,df=n-1)*s/sqrt(n)

lowerinterval <- xbar - margin
lowerinterval 
upperinterval <- xbar + margin
upperinterval 

################################################# Simpson and Shannon ANOVAs ####################################################

alpha_anova_shan <- aov(Shannon~Site_ID, data = alpha_div_merge)
summary(alpha_anova_shan)
alpha_anova_simp <- aov(Simpson~Site_ID, data = alpha_div_merge)
summary(alpha_anova_simp)

#Tukey Post Hoc
tukey=TukeyHSD(alpha_anova_shan, ordered=FALSE, conf.level=0.95)
TukeyHSD(alpha_anova_simp, ordered=FALSE, conf.level=0.95)

#letters of significance
install.packages("multcompView")
library("multcompView")
cld = multcompLetters4(alpha_anova_shan, tukey)
print(cld)

################################################### Linear Model ##############################################################

#Linear model (will use to account for uneven sample #s)
alpha_shan_lm <- lm(Shannon~Site_ID-1, data=alpha_div_merge) 
summary(alpha_shan_lm)
alpha_anova_simp_lm <- lm(Simpson~Site_ID-1, data=alpha_div_merge) 
summary(alpha_anova_simp_lm)

#################################################### Brittany's Code ########################################################

# I think this should produce 4 plots, but if not we can get them individually

plot(alpha_shan_lm)

#This code pulls residuals from the model and then plots the fitted vs residuals. 
#We are looking for no patterns in the plot (no U shaped or fanning pattern).

res <- resid(alpha_shan_lm)
plot(fitted(alpha_shan_lm), res)
qqnorm(res)
plot(density(res))

# lsmeans package for pairwise

install.packages("lsmeans")
library("lsmeans")

# use for pairwise comparison, Site_ID = specs?
a = lsmeans(alpha_shan_lm, "Site_ID")
pairs(a)

anova(alpha_shan_lm)

################################################# Shannon + Simpson Boxplots ##################################################################

#Site colours
SiteColours = c(REF1 = "#F8766D", REF2 = "#A3A500", DSW1 = "#00BF7D", DSK1 = "#00B0F6", DSK2 = "#E76BF3")

#Orders for boxplots
alpha_div_merge$Site_ID <- factor(alpha_div_merge$Site_ID, 
                                  levels=c("REF1", "REF2", "DSW1", "DSK1", "DSK2"))

alpha_div_merge$Location_vs_WWTP <- factor(alpha_div_merge$Location_vs_WWTP, 
                                           levels=c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener"), 
                                           labels=c("Upstream", "Downstream Waterloo", "Downstream Kitchener"))

#Shannon Diversity boxplot
ggplot(data=alpha_div_merge, aes(x=Site_ID, y=Shannon)) + 
  geom_boxplot(stat="boxplot", fill="gray95", color=SiteColours) + 
  facet_grid(~Location_vs_WWTP, scale="free") + xlab("Site") + ylab("Shannon Diversity") + 
  theme(text = element_text(size = 15))
# + geom_text(label = c("a", "a", "a", "b", "c"), aes(y = c(5.5, 6.1, 5.2, 5, 5.9), x = "Site"), size = 6)
# + theme(strip.text = element_text(size = 11, color="black")) + theme(axis.text.x=element_text(angle=45, hjust=1, size=12, color="black"), axis.text=element_text(size=12, color="black"), axis.title=element_text(size=18))
  

#Simpson Diversity boxplot
ggplot(data=alpha_div_merge, aes(x=Site_ID, y=Simpson)) + 
  geom_boxplot(stat="boxplot", fill="gray95", color=SiteColours) + 
  facet_grid(~Location_vs_WWTP, scale="free") + xlab("Site") + ylab("Simpson Diversity")
# + theme(strip.text = element_text(size = 11, color="black")) + theme(axis.text.x=element_text(angle=45, hjust=1, size=12, color="black"), axis.text=element_text(size=12, color="black"), axis.title=element_text(size=18))

################################################### Beta Diversity ##################################################################

newlabelorder_Location = c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener")
newlabelname_Location = c("Upstream", "Downstream Waterloo", "Downstream Kitchener")
newlabelorder_Site = c("REF1", "REF2", "DSW1", "DSK1", "DSK2")

#SiteColours = c(REF1 = "#DE8C00", REF2 = "#B79F00", DSW1 = "#00C08B", DSK1 = "#F564E3", DSK2 = "#FF64B0")
SiteColours = c(REF1 = "#F8766D", REF2 = "#A3A500", DSW1 = "#00BF7D", DSK1 = "#00B0F6", DSK2 = "#E76BF3")

# Transform data to proportions for Bray-Curtis distances
trans.dat <- transform_sample_counts(musdat, function(otu) otu/sum(otu))

# Data ordination
bray.ord = ordinate(trans.dat, method = "PCoA", distance = "bray")
bray.ord

# Plot Bray-Curtis + ellipses
bray.plot = plot_ordination(trans.dat, bray.ord, color="Site_ID", shape="Site_ID", 
            axes = c(1,2)) + scale_color_manual(values = SiteColours, name="Site", labels = c("REF1", "REF2", 
            "DSW1", "DSK1", "DSK2")) + scale_shape_manual(name = "Site",labels = c("REF1", "REF2", "DSW1", "DSK1", 
            "DSK2"), values = c(15,3,16,8,17)) + stat_ellipse(geom = "polygon", type="norm", alpha=0.2, 
            aes(fill=Site_ID)) + scale_fill_manual(values=SiteColours, name="Site") + geom_point(size= 2) + 
            theme(text = element_text(size = 15))
bray.plot$data$Site_ID <- as.character(bray.plot$data$Site_ID)
bray.plot$data$Site_ID <- factor(bray.plot$data$Site_ID, levels=newlabelorder_Site)
bray.plot 

############################################## Beta Diversity Matrices #############################################################

#Adonis (PERMANOVA) for Bray Curtis matrices
bray.test = phyloseq::distance(musdat, "bray")
sample.df = data.frame(sample_data(trans.dat))
adonis(formula = bray.test~Site_ID, data =  sample.df, permutations = 9999, method = "bray")

#Pairwise adonis to determine significant pairs
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pairwise.adonis(bray.test, sample_data(musdat)$Location_vs_WWTP, perm = 99999) #by Location
pairwise.adonis(bray.test, sample_data(musdat)$Site_ID, perm = 99999) #by Site

#Confidence intervals
bray.ci = bray.ord$vectors[,1:2]

means = aggregate(bray.ci~Site_ID, data = sample.df, mean)
View(means)
sdevs = aggregate(bray.ci~Site_ID, data = sample.df, sd)
View(sdevs)

n <- 10
xbar <- -0.01190460
s <- 0.11174954

margin <- qt(0.975,df=n-1)*s/sqrt(n)

lowerinterval <- xbar - margin
lowerinterval 
upperinterval <- xbar + margin
upperinterval 

############################################## SIMPER Analysis #############################################################

#Which species contribute most to Bray Curtis dissimilarity
prune.dat = prune_taxa(taxa_sums(dat.aglo) > 90, dat.aglo)

## Simper analysis, extract abundance matrix from a phyloseq object
musASV = as(otu_table(prune.dat), "matrix")

# transpose so we have the OTUs as columns
musASV = t(musASV)

# Coerce the object to a data.frame
musASV_scaled = as.data.frame(musASV)

# running the simper analysis on the dataframe and the variable of interest "Location"
mus_simper <- simper(musASV_scaled, sample_data(musdat)$Location_vs_WWTP)

# printing the top OTUs
print(mus_simper)
summ = summary(mus_simper)

# summary for each
upwat = summ$Upstream_Downstream_Waterloo
upkit = summ$Upstream_Downstream_Kitchener
watkit = summ$Downstream_Waterloo_Downstream_Kitchener

#filter to significant p values
upwat = subset(upwat, p < 0.05)
upkit = subset(upkit, p < 0.05)
watkit = subset(watkit, p < 0.05)

# make tax table
tax = tax_table(musdat)
write.csv(tax, "TaxTable.csv")
met = read.csv("TaxTable.csv", row.names = 1)

#match Genus to Species
upwat$Genus = met$Genus[match(rownames(upwat),rownames(met))]
upkit$Genus = met$Genus[match(rownames(upkit),rownames(met))]
watkit$Genus = met$Genus[match(rownames(watkit),rownames(met))]

# percentages
dat.trans = transform_sample_counts(musdat, function(x) x/sum(x))
prune.dat = prune_taxa(taxa_sums(dat.trans) > 0.02, dat.trans)
dat.dataframe = psmelt(prune.dat)
dat.agr = aggregate(Abundance~Location_vs_WWTP+OTU, data=dat.dataframe, FUN=mean)

dat.agr$Genus = met$Genus[match(dat.agr$OTU,rownames(met))]
View(dat.agr)

write.csv(dat.agr, "Abund.csv")

################################################## Relative Abundance Plots ########################################################################

#colours
library(scales)
show_col(hue_pal()(29)) #ggplot colours, can adjust to amount of taxa

#Bar colours
phylaColours1 = c(Acidobacteria = "#c2c134", Actinobacteria = "#e84bd3", Bacteroidetes = "#e5e736", 
                  Chloroflexi = "#b16af2", Cyanobacteria = "#43e07b", Dependentiae = "#7488ec", Firmicutes = "#eab42c", 
                  Fusobacteria = "#d67ad6", Nitrospirae = "#EF7F48", Patescibacteria = "#3ceabf", 
                  Planctomycetes = "#ee603a", Proteobacteria = "#8fd9e5", TA06 = "#00C088", 
                  Tenericutes = "#ec5d6c", Verrucomicrobia = "#c29937")

phylaColours2 = c(Actinobacteria = "#e84bd3", Bacteroidetes = "#e5e736", Cyanobacteria = "#43e07b", 
                  Firmicutes = "#eab42c", Proteobacteria = "#8fd9e5", Tenericutes = "#ec5d6c", 
                  Verrucomicrobia = "#c29937", `< 2% abund.` = "#A9A9A9")

########################################################
#Phylum - remove taxa in <2% abundance (normalize)
dat.aglo = tax_glom(musdat, taxrank = "Phylum")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
#prune.dat = prune_taxa(taxa_sums(dat.trans) > 0.02, dat.trans)
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~Site_ID+Location_vs_WWTP+Phylum, data=dat.dataframe, FUN=mean)

#Edit labels for plot
dat.agr$Site_ID= factor(dat.agr$Site_ID, 
                        levels=c("REF1", "REF2", "DSW1", "DSK1", "DSK2"))
dat.agr$Location_vs_WWTP= factor(dat.agr$Location_vs_WWTP, 
                                 levels=c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener"), 
                                 labels=c("Upstream", "Downstream Waterloo", "Downstream Kitchener"))

ggplot(dat.agr, aes(x=Site_ID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity") + 
  facet_grid(~Location_vs_WWTP, scale = "free") + 
  scale_fill_manual(values=phylaColours1) +
  labs(y= "Abundance >2%", x= "Site")

#Relative abundance proportions (Ex. Cyanobacteria)
with(dat.agr[dat.agr$Phylum == "Actinobacteria", ], table(Abundance, Site_ID))

#For the normalized data (>2%), what is in the lower 2%? (Abundance column) - added to plot caption
dat.aglo = tax_glom(musdat, taxrank = "Phylum")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.frame = psmelt(dat.trans)

dat.agr = aggregate(Abundance~Site_ID+Phylum, data=dat.frame, FUN=mean)

dat.frame
dat.frame[750:850, 1:4]
tail(dat.frame, 10)

#########################################################
#Phylum - display legend category for <2% abundant taxa
dat.aglo = tax_glom(musdat, taxrank = "Phylum")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~Sample_Label+Site_ID+Location_vs_WWTP+Phylum, data=dat.dataframe, FUN=mean)
dat.agr$Phylum <- as.character(dat.agr$Phylum)
#simple way to rename phyla with < 2% abundance
dat.agr$Phylum[dat.agr$Abundance < 0.02] <- "< 2% abund."

#Edit labels for plot
unique(dat.agr$Phylum) #what will show up in legend
dat.agr$Phylum <- factor(dat.agr$Phylum, levels = c("Actinobacteria", "Bacteroidetes", "Cyanobacteria", "Firmicutes", 
                                                    "Proteobacteria", "Tenericutes", "Verrucomicrobia", "< 2% abund."))
dat.agr$Site_ID= factor(dat.agr$Site_ID, 
                        levels=c("REF1", "REF2", "DSW1", "DSK1", "DSK2"))
dat.agr$Location_vs_WWTP= factor(dat.agr$Location_vs_WWTP, 
                                 levels=c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener"), 
                                 labels=c("Upstream", "Downstream Waterloo", "Downstream Kitchener"))

ggplot(dat.agr, aes(x=Site_ID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity") + 
  facet_grid(~Location_vs_WWTP, scale = "free") + 
  scale_fill_manual(values=phylaColours2) +
  labs(y= "Relative Abundance", x= "Site")

#Relative abundance proportions (Ex. Cyanobacteria)
with(dat.agr[dat.agr$Phylum == "Tenericutes", ], table(Abundance, Site_ID))

##################################################### Cyanobacteria ############################################################

#Subset to Cyanobacteria
cyano = musdat %>% subset_taxa(Phylum == "Cyanobacteria")
cyano

top.taxa = sort(tapply(taxa_sums(cyano), tax_table(cyano)[, "Genus"], sum), TRUE)
top.taxa[1:20]

genusColours = c("#e96f29","#6d83f0","#d7da37","#b06eec","#58d658","#e659cb","#99de4b","#f14b8a","#75e187","#cd86d5",
                 "#75ad31","#3295e9","#dca62f","#6a86d6","#d3e06c","#8e9aec","#9d9e3b","#4f8ed0","#ea655d","#65e9b3",
                 "#e1789c","#5ea655","#53a7e8","#dd8853","#4ae6dc","#cda95c","#40b186","#c3e08e")

#Genus
dat.aglo = tax_glom(cyano, taxrank = "Genus")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~Site_ID+Location_vs_WWTP+Genus, data=dat.dataframe, FUN=mean)

dat.agr$Site_ID= factor(dat.agr$Site_ID, 
                        levels=c("REF1", "REF2", "DSW1", "DSK1", "DSK2"))
dat.agr$Location_vs_WWTP= factor(dat.agr$Location_vs_WWTP, 
                                 levels=c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener"), 
                                 labels=c("Upstream", "Downstream Waterloo", "Downstream Kitchener"))

ggplot(dat.agr, aes(x=Site_ID, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity") + 
  facet_grid(~Location_vs_WWTP, scale = "free") + 
  labs(y="Abundance", x="Site") + 
  scale_fill_manual(values=genusColours) + theme(legend.position='none') #+ theme(legend.position="bottom")

# Relative abundance matrices (Ex. Cyanobium_PCC-6307)
with(dat.agr[dat.agr$Genus == "Cyanobium_PCC-6307", ], table(Abundance, Site))

####################################################### DESeq2 for genera #####################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")
packageVersion("DESeq2")

library(scales)
show_col(hue_pal()(13)) #ggplot colours, can adjust to amount of taxa

DColours = c(Fusobacteria = "#F8766D", Actinobacteria = "#E18A00", Bacteroidetes = "#BE9C00", 
             Cyanobacteria = "#8CAB00", `Deinococcus-Thermus` = "#24B700", Firmicutes = "#00BE70", 
             Nitrospirae = "#00C1AB", Patescibacteria = "#00BBDA", Planctomycetes = "#00ACFC", 
             Proteobacteria = "#8B93FF", Synergistetes = "#D575FE", Tenericutes = "#F962DD", Verrucomicrobia = "#FF65AC")

#####################################
#####Upstream and Downstream Waterloo 

REF_DSW = subset_samples(dat_r, Location_vs_WWTP=="Downstream_Waterloo" | Location_vs_WWTP=="Upstream")
REF_DSW <- prune_taxa(taxa_sums(REF_DSW) >0, REF_DSW)
REF_DSW #8476 ASVs, 30 samples

deseq2.REF_DSW <- phyloseq_to_deseq2(REF_DSW, ~Location_vs_WWTP)
deseq2.REF_DSW <- DESeq(deseq2.REF_DSW, test="Wald", fitType = "parametric")
deseq2.REF_DSW
plotMA(deseq2.REF_DSW)
plotDispEsts(deseq2.REF_DSW)

res = results(deseq2.REF_DSW, cooksCutoff = FALSE)
res
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(REF_DSW)[rownames(sigtab), ], "matrix"))
head(sigtab)

#plot
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

ggplot(sigtab, aes(x=log2FoldChange, y=Genus, color=Phylum)) + geom_point(size=4) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13)) +
  theme(legend.text=element_text(size=13, color="black")) + 
  theme(legend.title = element_text(size=13, color="black")) + ggtitle("Upstream vs Downstream Waterloo") + 
  scale_colour_manual(values=DColours) + xlim(-33, 33) +
  theme(legend.position="none")

#####Upstream and Downstream Kitchener

REF_DSK = subset_samples(dat_r, Location_vs_WWTP=="Downstream_Kitchener" | Location_vs_WWTP=="Upstream")
REF_DSK <- prune_taxa(taxa_sums(REF_DSK) >0, REF_DSK)
REF_DSK #3273 ASVs, 68 samples

deseq2.REF_DSK <- phyloseq_to_deseq2(REF_DSK, ~Location_vs_WWTP)
deseq2.REF_DSK <- DESeq(deseq2.REF_DSK, test="Wald", fitType="parametric")
deseq2.REF_DSK

res = results(deseq2.REF_DSK, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(REF_DSK)[rownames(sigtab), ], "matrix"))
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

ggplot(sigtab, aes(x=log2FoldChange, y=Genus, color=Phylum)) + geom_point(size=4) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13)) +
  theme(legend.text=element_text(size=13, color="black")) + 
  theme(legend.title = element_text(size=13, color="black")) + ggtitle("Upstream vs Downstream Kitchener")+ 
  scale_colour_manual(values=DColours) + xlim(-33, 33)+
  theme(legend.position="none")

#####Downstream Waterloo and Kitchener

DSW_DSK = subset_samples(dat_r, Location_vs_WWTP=="Downstream_Waterloo" | Location_vs_WWTP=="Downstream_Kitchener")
DSW_DSK <- prune_taxa(taxa_sums(DSW_DSK) >0, DSW_DSK)
DSW_DSK #2453 ASVs, 68 samples

deseq2.DSW_DSK <- phyloseq_to_deseq2(DSW_DSK, ~Location_vs_WWTP)
deseq2.DSW_DSK <- DESeq(deseq2.DSW_DSK, test="Wald", fitType="parametric")
deseq2.DSW_DSK

res = results(deseq2.DSW_DSK, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(DSW_DSK)[rownames(sigtab), ], "matrix"))
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

ggplot(sigtab, aes(x=log2FoldChange, y=Genus, color=Phylum)) + geom_point(size=4) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13)) +
  theme(legend.text=element_text(size=13, color="black")) + 
  theme(legend.title = element_text(size=13, color="black")) + ggtitle("Downstream Waterloo vs Downstream Kitchener")+ 
  scale_colour_manual(values=DColours) + xlim(-33, 33)+
  theme(legend.position="none")

######################################## lfcShrink for zero-inflation #########################################

dds <- phyloseq_to_deseq2(DSW_DSK, ~ Location_vs_WWTP)
dds <- DESeq(dds, test = "Wald", fitType = "local", sfType = "poscounts")
plotMA(dds)
plotDispEsts(dds)
res <- results(dds, cooksCutoff = FALSE)
res

df_res <- cbind(as(res, "data.frame"), as(tax_table(DSW_DSK)[rownames(res), ], "matrix"))
df_res <- df_res %>%
  rownames_to_column(var = "OTU") %>%
  arrange(padj)

(fdr_otu <- df_res %>%
    dplyr::filter(padj < 0.1))     #none with FDR p-val < 0.1

resultsNames(dds)

res_LFC <- lfcShrink(dds, coef=2, type="apeglm")   #can also consider alternative shrinkage estimator and plot most DE species
plotMA(res_LFC)

df_res_lfc <- cbind(as(res_LFC, "data.frame"), as(tax_table(DSW_DSK)[rownames(res_LFC), ], "matrix"))
df_res_lfc <- df_res_lfc %>%
  rownames_to_column(var = "OTU") %>%
  arrange(desc(log2FoldChange))

lfc_otu <- df_res_lfc %>%
  filter(abs(log2FoldChange) > 1) %>%
  droplevels() %>%
  arrange(desc(log2FoldChange))

ggplot(lfc_otu, aes(x = Genus, y = log2FoldChange, color = Genus)) +
  geom_point(size = 4) +
  labs(y = "\nLog2 Fold-Change for REF vs. DSW", x = "") +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dotted")



#
#title: "Heptageniidae Microbiome"
#author: "Elise Millar"
#date: "Aug 27th, 2020"
#output: html_document
#

install.packages("rlang")
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

#################

#Subset to only ~Heptageniidae~ samples
hepdat <- subset_samples(data, Sample_Label=="Heptageniidae")
hepdat <- prune_taxa(taxa_sums(hepdat) >0, hepdat)
hepdat #7533 ASVs, 50 samples

##################################################### Different Subsets #############################################################

#Samples containing WWTP Effluent Genera
effluent = c('Ottowia', 'Fusibacter', 'Bifidobacterium','Actinomyces', 'Prevotellaceae_UCG-004', 'Cloacibacterium', 
             'Chiayiivirga', 'Lelliottia', 'BD1-7_clade', 'Steroidobacter', 'Acidaminococcus', 
             'Candidatus_Paenicardinium')
effgen = subset_taxa(hepdat, (Genus %in% effluent))
effgen <- prune_samples(sample_sums(effgen) >0, effgen)
efflist = psmelt(effgen)
efflist[1:30, c('Sample','Genus','Abundance','Site_ID')]

#By Phylum (Ex. Cyanobacteria)
cyano = subset_taxa(hepdat, Phylum=="Cyanobacteria")
cyano

#By Site (Ex. REF1)
REF1 <- subset_samples(hepdat, Site_ID=="REF1")
REF1 <- prune_taxa(taxa_sums(REF1) >0, REF1)
REF1

################################################ How many taxonomic levels ###############################################################
#Excluding NA taxa

get_taxa_unique(hepdat, "Phylum") #31
get_taxa_unique(hepdat, "Class") #64
get_taxa_unique(hepdat, "Order") #169
get_taxa_unique(hepdat, "Family") #246
get_taxa_unique(hepdat, "Genus") #447

#################################################### Abundance Numbers ############################################################################

##Total number of sequence reads for dataset (including unassigned)
sum(sample_sums(hepdat)) #1578826
#Read sequence range for samples
sort(sample_sums(hepdat)) #18149-46092

#Number of sequence reads per rank (Ex. Phylum)
top.taxa = sort(tapply(taxa_sums(hepdat), tax_table(hepdat)[, "Phylum"], sum), TRUE)

top.taxa #highest amount of reads to lowest
top.taxa[1:10] #top 10 taxa by read count
sum(top.taxa) #total reads assigned to taxa = 1577369

################
#Didn't need this code but could be useful 
#Reads per Species
taxa_sums(hepdat)

##Number of clustered ASVs 
ntaxa(hepdat)
tax_table(hepdat)
top_taxa = table(tax_table(hepdat)[, "Phylum"], exclude = NULL) #Ex. Phylum
top_taxa
sort(top_taxa, decreasing = TRUE)
head(sort(top_taxa, decreasing = TRUE))

#ASV counts per Species
asv_df <- t(otu_table(hepdat))
colSums(asv_df != 0)
###############

###################################################### Alpha Diversity #########################################################################

#Sharokh's rarefy method
rare.depth <- min(sample_sums(hepdat))
rare.depth #18149
dat_r <- rarefy_even_depth(hepdat, sample.size = rare.depth, rngseed=1414) #rngseed = random value
dat_r

#Plot by Site_ID
alpha_plot = plot_richness(dat_r, x="Site_ID", measures=c("Shannon", "Simpson"))
newlabelorder_Site_ID = c("REF1", "REF2", "REF3", "DSW1", "DSW2", "DSW3", "DSK1", "DSK2", "DSK3", "DSK4")
alpha_plot$data$Site_ID <- as.character(alpha_plot$data$Site_ID)
alpha_plot$data$Site_ID <- factor(alpha_plot$data$Site_ID, levels=newlabelorder_Site_ID)
alpha_plot + xlab("Site")

################################################# Alpha Diversity Metrics ###################################################################

alpha_diversity <- estimate_richness(dat_r, split=TRUE, measures=c("Shannon", "Simpson"))
alpha_div <- data.frame(alpha_diversity)
alpha_div

#Add Sample_ID column 
alpha_div["Sample_ID"] <- NA 
alpha_div$Sample_ID <- c(365:414) #exact range for Heptageniidae samples in metadata

#Double-check
colnames(alpha_div)
alpha_div

#Merge metadata and alpha diversity numbers
alpha_met_df <- mutate(met_df, Sample_ID = c(219:280, 282:283, 285:361, 365:414, 415:570, 
                                             964:1005)) #ranges for ALL samples in metadata
alpha_met_df

ncol(alpha_div)
ncol(alpha_met_df)
as.character(alpha_met_df$Sample_ID)
as.character(alpha_div$Sample_ID)

alpha_div_merge <- merge(alpha_met_df, alpha_div, by = "Sample_ID") 
summary(alpha_div_merge)

################################################# Simpson and Shannon ANOVAs ####################################################

alpha_anova_simp <- aov(Simpson~Site_ID, data = alpha_div_merge)
summary(alpha_anova_simp)
alpha_anova_shan <- aov(Shannon~Site_ID, data = alpha_div_merge)
summary(alpha_anova_shan)

#Tukey Post Hoc
TukeyHSD(alpha_anova_simp, ordered=FALSE, conf.level=0.95)
TukeyHSD(alpha_anova_shan, ordered=FALSE, conf.level=0.95)

#Linear model (didn't use)
alpha_anova_simp_lm <- lm(Simpson~Site_ID-1, data=alpha_div_merge) 
summary(alpha_anova_simp_lm)
alpha_anova_shan_lm <- lm(Shannon~Site_ID-1, data=alpha_div_merge) 
summary(alpha_anova_shan_lm)

################################################# Shannon + Simpson Boxplots ##################################################################

#Site colours
SiteColours = c(REF1 = "#F8766D", REF2 = "#DE8C00", REF3 = "#B79F00", DSW1 = "#7CAE00", DSW2 = "#00BA38", DSW3 = "#00BFC4", 
                DSK1 = "#00B4F0", DSK2 = "#619CFF", DSK3 = "#C77CFF", DSK4 = "#FF64B0")

#Orders for boxplots
alpha_div_merge$Site_ID <- factor(alpha_div_merge$Site_ID, 
                                  levels=c("REF1", "REF2", "REF3", "DSW1", "DSW2", "DSW3", "DSK1", "DSK2", "DSK3", "DSK4"))

alpha_div_merge$Location_vs_WWTP <- factor(alpha_div_merge$Location_vs_WWTP, 
                                           levels=c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener"), 
                                           labels=c("Upstream", "Downstream Waterloo", "Downstream Kitchener"))

#Shannon Diversity boxplot
ggplot(data=alpha_div_merge, aes(x=Site_ID, y=Shannon)) + 
  geom_boxplot(stat="boxplot", fill="gray95", color=SiteColours) + 
  facet_grid(~Location_vs_WWTP, scales = "free_x", space = "free_x") + 
  expand_limits(y=c(0, 6)) + xlab("Site") + ylab("Shannon Diversity") + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
# + theme(strip.text = element_text(size = 11, color="black")) + theme(axis.text.x=element_text(angle=45, hjust=1, size=12, color="black"), axis.text=element_text(size=12, color="black"), axis.title=element_text(size=18))

#Simpson Diversity boxplot
ggplot(data=alpha_div_merge, aes(x=Site_ID, y=Simpson)) + 
  geom_boxplot(stat="boxplot", fill="gray95", color=SiteColours) + 
  facet_grid(~Location_vs_WWTP, scale="free") + xlab("Site") + ylab("Simpson Diversity")
# + theme(strip.text = element_text(size = 11, color="black")) + theme(axis.text.x=element_text(angle=45, hjust=1, size=12, color="black"), axis.text=element_text(size=12, color="black"), axis.title=element_text(size=18))


################################################### Beta Diversity ##################################################################

newlabelorder_Location = c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener")
newlabelname_Location = c("Upstream", "Downstream Waterloo", "Downstream Kitchener")
newlabelorder_Site = c("REF1", "REF2", "REF3", "DSW1", "DSW2", "DSW3", "DSK1", "DSK2", "DSK3", "DSK4")

SiteColours = c(REF1 = "#F8766D", REF2 = "#DE8C00", REF3 = "#B79F00", DSW1 = "#7CAE00", DSW2 = "#00BA38", DSW3 = "#00BFC4", 
                DSK1 = "#00B4F0", DSK2 = "#619CFF", DSK3 = "#C77CFF", DSK4 = "#FF64B0")

# Transform data to proportions for Bray-Curtis distances
trans.dat <- transform_sample_counts(hepdat, function(otu) otu/sum(otu))

# Data ordination
bray.ord = ordinate(trans.dat, method = "PCoA", distance = "bray")
bray.ord

# Plot Bray-Curtis + ellipses
bray.plot = plot_ordination(trans.dat, bray.ord, color="Site_ID", title="Bray-Curtis Heptageniidae", axes = c(1,2)) + 
  scale_color_manual(values = SiteColours, name="Site") 
bray.plot$data$Site_ID <- as.character(bray.plot$data$Site_ID)
bray.plot$data$Site_ID <- factor(bray.plot$data$Site_ID, levels=newlabelorder_Site)
bray.plot + stat_ellipse(geom = "polygon", type="norm", alpha=0.2, aes(fill=Site_ID)) + 
  scale_fill_manual(values=SiteColours, name="Site")

############################################## Beta Diversity Matrices #############################################################

#Adonis (PERMANOVA) for Bray Curtis matrices
bray.test = phyloseq::distance(trans.dat, "bray")
sample.df = data.frame(sample_data(trans.dat))
adonis(formula = bray.test~Site_ID, data =  sample.df, permutations = 9999, method = "bray")


#Pairwise adonis to determine significant pairs
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pairwise.adonis(bray.test, sample_data(hepdat)$Location_vs_WWTP, perm = 99999) #by Location
pairwise.adonis(bray.test, sample_data(hepdat)$Site_ID, perm = 99999) #by Site

################################################## Relative Abundance Plots ########################################################################

#colours
library(scales)
show_col(hue_pal()(29)) #ggplot colours, can adjust to amount of taxa. Didn't use but handy to have.

#Bar colours
phylaColours1 = c(Acidobacteria = "#49e3d9", Actinobacteria = "#e45ad1", Bacteroidetes = "#a8e14b", 
                 Chloroflexi = "#9e73f2", Cyanobacteria = "#decf3f", Firmicutes = "#e66aa1", Gemmatimonadetes = "#57ca98", 
                 Nitrospirae = "#eb6071", Patescibacteria = "#51a6ea", Planctomycetes = "#e77543", 
                 Proteobacteria = "#528cd6", RsaHf231 = "#cf9c45", Synergistetes = "#90a342", 
                 Tenericutes = "#b289df", Verrucomicrobia = "#cde086")

phylaColours2 = c(Acidobacteria = "#c2c134", Actinobacteria = "#e84bd3", Bacteroidetes = "#e5e736", 
                  Chloroflexi = "#b16af2", Cyanobacteria = "#43e07b", Dependentiae = "#7488ec", Firmicutes = "#eab42c", 
                  Fusobacteria = "#d67ad6", Gemmatimonadetes = "#31abb2", 
                  Nitrospirae = "#ea5da0", Patescibacteria = "#3ceabf", Planctomycetes = "#ee603a", 
                  Proteobacteria = "#8fd9e5", RsaHf231 = "#a19d2e", Synergistetes = "#e5d869", TA06 = "#df8538", 
                  Tenericutes = "#ec5d6c", Verrucomicrobia = "#c29937", `< 2% abund.` = "#A9A9A9")

########################################################
#Phylum - remove taxa in <2% abundance (normalize)
dat.aglo = tax_glom(hepdat, taxrank = "Phylum")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
prune.dat = prune_taxa(taxa_sums(dat.trans) > 0.02, dat.trans)
dat.dataframe = psmelt(prune.dat)
dat.agr = aggregate(Abundance~Sample_Label+Site_ID+Location_vs_WWTP+Phylum, data=dat.dataframe, FUN=mean)

#Edit labels for plot
dat.agr$Site_ID= factor(dat.agr$Site_ID, 
                        levels=c("REF1", "REF2", "REF3", "DSW1", "DSW2", "DSW3", "DSK1", "DSK2", "DSK3", "DSK4"))
dat.agr$Location_vs_WWTP= factor(dat.agr$Location_vs_WWTP, 
                                 levels=c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener"), 
                                 labels=c("Upstream", "Downstream Waterloo", "Downstream Kitchener"))

ggplot(dat.agr, aes(x=Site_ID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity") + 
  facet_grid(~Location_vs_WWTP, scales = "free_x", space = "free_x") + 
  scale_fill_manual(values=phylaColours1) +
  labs(y= "Relative Abundance", x= "Site") + theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
                                                   axis.text.y=element_text(size=12),
                                                   axis.title=element_text(size=12),
                                                   strip.text.x = element_text(size = 12),
                                                   legend.position = "none")

#Relative abundance proportions (Ex. Cyanobacteria)
with(dat.agr[dat.agr$Phylum == "Cyanobacteria", ], table(Abundance, Site_ID))

#For the normalized data (>2%), what is in the lower 2%? (Abundance column) - added to plot caption
dat.aglo = tax_glom(hepdat, taxrank = "Phylum")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.frame = psmelt(dat.trans)
dat.frame
dat.frame[200:300, 1:4]
tail(dat.frame, 10)

#########################################################
#Phylum - display legend category for <2% abundant taxa
dat.aglo = tax_glom(hepdat, taxrank = "Phylum")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~Sample_Label+Site_ID+Location_vs_WWTP+Phylum, data=dat.dataframe, FUN=mean)
dat.agr$Phylum <- as.character(dat.agr$Phylum)
#simple way to rename phyla with < 2% abundance
dat.agr$Phylum[dat.agr$Abundance < 0.02] <- "< 2% abund."

#Edit labels for plot
unique(dat.agr$Phylum) #what will show up in legend
dat.agr$Phylum <- factor(dat.agr$Phylum, levels = c("Bacteroidetes", "Cyanobacteria", "Firmicutes", "Proteobacteria", 
                                                    "RsaHf231", "Synergistetes", "Tenericutes", "< 2% abund."))
dat.agr$Site_ID= factor(dat.agr$Site_ID, 
                        levels=c("REF1", "REF2", "REF3", "DSW1", "DSW2", "DSW3", "DSK1", "DSK2", "DSK3", "DSK4"))
dat.agr$Location_vs_WWTP= factor(dat.agr$Location_vs_WWTP, 
                                 levels=c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener"), 
                                 labels=c("Upstream", "Downstream Waterloo", "Downstream Kitchener"))

ggplot(dat.agr, aes(x=Site_ID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity") + 
  facet_grid(~Location_vs_WWTP, scales = "free_x", space = "free_x") + 
  scale_fill_manual(values=phylaColours2) +
  labs(y= "Abundance", x= "Site")

##################################################### Cyanobacteria ############################################################

#Subset to Cyanobacteria
cyano = hepdat %>% subset_taxa(Phylum == "Cyanobacteria")
cyano

genusColours = c(Aphanizomenon_MDT14a = "#9e73f2", Calothrix_KVSF5 = "#bee148", `Calothrix_PCC-6303` = "#de62d9", 
                 Calothrix_UAM_374 = "#6bdb5e", `Chamaesiphon_PCC-7430` = "#ea54a4", 
                 `Chroococcidiopsis_PCC-6712` = "#6ce9a8", `Cyanobium_PCC-6307` = "#f34478", 
                 `Geminocystis_PCC-6308` = "#4ae8df", Gleocapsa = "#ee5f51", HAVOmat113 = "#248fe1", 
                 `JSC-12` = "#ebd530", Limnothrix = "#5686ed", `Microcystis_PCC-7914` = "#73ab37", 
                 Neosynechococcus = "#b486db", `Nodosilinea_PCC-7104` = "#c7b740", 
                 `Oscillatoria_SAG_1459-8` = "#8e9df0", Phormidesmis_ANT.LACV5.1 = "#da9835", 
                 Phormidium_CYN64 = "#62a6f0", `Planktothrix_NIVA-CYA_15` = "#e46f32", 
                 `Pleurocapsa_PCC-7319` = "#3ba7e5", Potamolinea_1PC = "#cce387", 
                 `Pseudanabaena_PCC-6802` = "#7a85d0", `Pseudanabaena_PCC-7429` = "#58b064", 
                 `Schizothrix_LEGE_07164` = "#de82c6", Snowella_0TU37S04 = "#3fba99", 
                 Synechococcus_MBIC10613 = "#e37589", `Synechocystis_CCALA_700` = "#4f8ed0", 
                 `Synechocystis_PCC-6803` = "#a49d47",`Tychonema_CCAP_1459-11B` = "#da925e")

#Genus
dat.aglo = tax_glom(cyano, taxrank = "Genus")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~Site_ID+Location_vs_WWTP+Genus, data=dat.dataframe, FUN=mean)

dat.agr$Site_ID= factor(dat.agr$Site_ID, 
                        levels=c("REF1", "REF2", "REF3", "DSW1", "DSW2", "DSW3", "DSK1", "DSK2", "DSK3", "DSK4"))
dat.agr$Location_vs_WWTP= factor(dat.agr$Location_vs_WWTP, 
                                 levels=c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener"), 
                                 labels=c("Upstream", "Downstream Waterloo", "Downstream Kitchener"))

ggplot(dat.agr, aes(x=Site_ID, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity") + 
  facet_grid(~Location_vs_WWTP, scales = "free_x", space = "free_x") + 
  scale_fill_manual(values=genusColours) +
  labs(y= "Relative Abundance", x= "Site") + theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
                                                   axis.text.y=element_text(size=12),
                                                   axis.title=element_text(size=12),
                                                   strip.text.x = element_text(size = 12),
                                                   legend.position = "none")

# Relative abundance matrices (Ex. Tychonema_CCAP_1459-11B)
with(dat.agr[dat.agr$Genus == "Tychonema_CCAP_1459-11B", ], table(Abundance, Site_ID))


####How many individuals, at which sites have cyano genera?
#Top 6 genera
top.taxa = sort(tapply(taxa_sums(cyano), tax_table(cyano)[, "Genus"], sum), TRUE)
top.taxa[1:6] #top 6 taxa by read count

cyanogen = c('Tychonema_CCAP_1459-11B', 'Cyanobium_PCC-6307', 'Chamaesiphon_PCC-7430',
             'Snowella_0TU37S04', 'Pleurocapsa_PCC-7319', 'Calothrix_KVSF5')
cygen = subset_taxa(hepdat, (Genus %in% cyanogen))
cygen = prune_samples(sample_sums(cygen) >0, cygen)
cglist = psmelt(cygen)
cglist = cglist[order(cglist$Genus, cglist$Site_ID, cglist$Sample),] #order by genus, site, sample
cglist = cglist[ , c('Sample','Genus','Abundance','Site_ID')] #select columns from metadata
cglist = cglist[!(cglist$Abundance=="0"),] #remove zero abundance entries
cglist


####################################################### Endosymbionts #############################################################

#Subset to Endosymbionts
endobac = c('Wolbachia', 'Arsenophonus', 'Spiroplasma', 'Rickettsia',
            'Candidatus_Cardinium', 'Rickettsiella')
endo = subset_taxa(hepdat, (Genus %in% endobac))
endo

#Top genera
top.taxa = sort(tapply(taxa_sums(endo), tax_table(endo)[, "Genus"], sum), TRUE)
top.taxa

endoColours = c(Arsenophonous = "#F8766D", Candidatus_Cardinium = "#B79F00", Rickettsia = "#00BA38", Rickettsiella = "#00BFC4", 
                Spiroplasma = "#619CFF", Wolbachia = "#F564E3")
endoColours = c(Arsenophonous = "#e07238", Candidatus_Cardinium = "#3a91fb", Rickettsia = "#cd62db", Rickettsiella = "#6ad165", 
                Spiroplasma = "#e7558d",Wolbachia = "#c6c643")

#Genus
dat.aglo = tax_glom(endo, taxrank = "Genus")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~Site_ID+Location_vs_WWTP+Genus, data=dat.dataframe, FUN=mean)

dat.agr$Site_ID= factor(dat.agr$Site_ID, 
                        levels=c("REF1", "REF2", "REF3", "DSW1", "DSW2", "DSW3", "DSK1", "DSK2", "DSK3", "DSK4"))
dat.agr$Location_vs_WWTP= factor(dat.agr$Location_vs_WWTP, 
                                 levels=c("Upstream", "Downstream_Waterloo", "Downstream_Kitchener"), 
                                 labels=c("Upstream", "Downstream Waterloo", "Downstream Kitchener"))

ggplot(dat.agr, aes(x=Site_ID, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity") + 
  facet_grid(~Location_vs_WWTP, scales = "free_x", space = "free_x") + 
  scale_fill_manual(values=endoColours) +
  labs(y= "Abundance", x= "Site") + theme(axis.text.x=element_text(angle=45, hjust=1))

#Relative abundance proportions (Ex. Rickettsiella)
with(dat.agr[dat.agr$Genus == "Rickettsiella", ], table(Abundance, Site_ID))


####################################################### DESeq2 for genera #####################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")
packageVersion("DESeq2")

DColours = c(Acidobacteria = "#F8766D", Actinobacteria = "#E18A00", Bacteroidetes = "#BE9C00", 
             Cyanobacteria = "#8CAB00", `Deinococcus-Thermus` = "#24B700", Firmicutes = "#00BE70", 
             Nitrospirae = "#00C1AB", Patescibacteria = "#00BBDA", Planctomycetes = "#00ACFC", 
             Proteobacteria = "#8B93FF", Synergistetes = "#D575FE", Tenericutes = "#F962DD", Verrucomicrobia = "#FF65AC")

#####Upstream and Downstream Waterloo 

REF_DSW = subset_samples(hepdat, Location_vs_WWTP=="Downstream_Waterloo" | Location_vs_WWTP=="Upstream")
REF_DSW <- prune_taxa(taxa_sums(REF_DSW) >0, REF_DSW)
REF_DSW #2002 ASVs, 60 samples

deseq2.REF_DSW <- phyloseq_to_deseq2(REF_DSW, ~Location_vs_WWTP)
deseq2.REF_DSW <- DESeq(deseq2.REF_DSW, test="Wald", fitType="parametric")
deseq2.REF_DSW

res = results(deseq2.REF_DSW, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(REF_DSW)[rownames(sigtab), ], "matrix"))
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
  theme(legend.title = element_text(size=13, color="black")) + ggtitle("Upstream vs Downstream Waterloo")+ 
  scale_colour_manual(values=DColours) + xlim(-33, 33)


#####Upstream and Downstream Kitchener

REF_DSK = subset_samples(hepdat, Location_vs_WWTP=="Downstream_Kitchener" | Location_vs_WWTP=="Upstream")
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
  scale_colour_manual(values=DColours) + xlim(-33, 33)


#####Downstream Waterloo and Kitchener

DSW_DSK = subset_samples(hepdat, Location_vs_WWTP=="Downstream_Waterloo" | Location_vs_WWTP=="Downstream_Kitchener")
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
  scale_colour_manual(values=DColours) + xlim(-33, 33)

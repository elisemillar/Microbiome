#
#title: "Water Microbiome"
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

#Subset to only ~Water~ samples
watdat <- subset_samples(data, Sample_Label=="Water" | Sample_Label=="W-WWTP" | Sample_Label=="K-WWTP")
watdat <- prune_taxa(taxa_sums(watdat) >0, watdat)
watdat #7740 ASVs, 42 samples

#Subset to only ~River~ samples
rivdat <- subset_samples(watdat, Site_ID != "WAT" & Site_ID != "KIT")
rivdat <- prune_taxa(taxa_sums(rivdat) >0, rivdat)
rivdat #6199 ASVs, 36 samples

#Subset to only ~Waterloo WWTP Effluent~ samples
wadat <- subset_samples(watdat, Site_ID =="WAT")
wadat <- prune_taxa(taxa_sums(wadat) >0, wadat)
wadat #1973 ASVs, 3 samples

#Subset to only ~Kitchener WWTP Effluent~ samples
kidat <- subset_samples(watdat, Site_ID =="KIT")
kidat <- prune_taxa(taxa_sums(kidat) >0, kidat)
kidat #262 ASVs, 3 samples

##################################################### Different Subsets #############################################################

#River samples containing WWTP Effluent Genera
#effluent = c('Ottowia', 'Fusibacter', 'Bifidobacterium',
             #'Actinomyces', 'Prevotellaceae_UCG-004', 'Cloacibacterium', 'Chiayiivirga', 
             #'Lelliottia', 'BD1-7_clade', 'Steroidobacter', 'Acidaminococcus', 
             #'Candidatus_Paenicardinium', 'Flavitalea', 'SWB02')
effluent = c('Ottowia', 'Fusibacter', 'Bifidobacterium',
             'Actinomyces', 'Prevotellaceae_UCG-004', 'Cloacibacterium', 'Chiayiivirga', 
             'Lelliottia', 'BD1-7_clade', 'Steroidobacter', 'Acidaminococcus', 
             'Candidatus_Paenicardinium', 'Flavitalea', 'SWB02', 'Acetoanaerobium', 'Aquimonas', 
             'AUTHM297', 'Bact-08', 'C1-B045', 'Candidatus_Cloacimonas', 'Candidatus_Protochlamydia', 
             'Desulfobacter', 'Leptotrichia', 'Mesotoga', 'Neochlamydia', 'Planctopirus', 'Proteiniclasticum', 
             'SC103', 'Succinivibrio', 'Thermovirga', 'Turneriella', 'U29-B03', 'XBB1006')
effgen = subset_taxa(rivdat, (Genus %in% effluent))
effgen <- prune_samples(sample_sums(effgen) >0, effgen)
efflist = psmelt(effgen)
efflist[1:81, c(2,3,5,33)]

#By Phylum (Ex. Cyanobacteria)
cyano = subset_taxa(watdat, Phylum=="Cyanobacteria")
cyano

#By Site (Ex. REF1)
REF1 <- subset_samples(watdat, Site_ID=="REF1")
REF1 <- prune_taxa(taxa_sums(REF1) >0, REF1)
REF1

#By Location (Ex. Upstream)
ups <- subset_samples(rivdat, Location_vs_WWTP=="Upstream")
ups <- prune_taxa(taxa_sums(ups) >0, ups)
ups

################################################ How many taxonomic levels ###############################################################
#Excluding NA taxa

get_taxa_unique(watdat, "Phylum") #40
get_taxa_unique(watdat, "Class") #80
get_taxa_unique(watdat, "Order") #215
get_taxa_unique(watdat, "Family") #313
get_taxa_unique(watdat, "Genus") #667

#################################################### Abundance Numbers ############################################################################

##Total number of sequence reads for dataset (including unassigned)
sum(sample_sums(watdat)) #906533
#Read sequence range for samples
sort(sample_sums(watdat)) #937-82058

#Number of sequence reads per rank (Ex. Phylum)
top.taxa = sort(tapply(taxa_sums(watdat), tax_table(watdat)[, "Phylum"], sum), TRUE)

top.taxa #highest amount of reads to lowest
top.taxa[1:10] #top 10 taxa by read count
sum(top.taxa) #total reads assigned to taxa = 905881

################
#Didn't need this code but could be useful 
#Reads per Species
taxa_sums(watdat)

##Number of clustered ASVs 
ntaxa(watdat)
tax_table(watdat)
top_taxa = table(tax_table(watdat)[, "Phylum"], exclude = NULL) #Ex. Phylum
top_taxa
sort(top_taxa, decreasing = TRUE)
head(sort(top_taxa, decreasing = TRUE))

#ASV counts per Species
asv_df <- t(otu_table(watdat))
colSums(asv_df != 0)
###############

###################################################### Alpha Diversity #########################################################################

#Sharokh's rarefy method
rare.depth <- min(sample_sums(watdat))
rare.depth #937
dat_r <- rarefy_even_depth(watdat, sample.size = rare.depth, rngseed=1414) #rngseed = random value
dat_r

#Plot by Site_ID
alpha_plot = plot_richness(dat_r, x="Site_ID", measures=c("Shannon", "Simpson"))
newlabelorder_Site_ID = c("REF1", "REF2", "REF3", "WAT", "DSW1", "DSW2", "G3", "DSW3", "KIT", "DSK1", 
                          "DSK2", "DSK3", "J1", "DSK4")
alpha_plot$data$Site_ID <- as.character(alpha_plot$data$Site_ID)
alpha_plot$data$Site_ID <- factor(alpha_plot$data$Site_ID, levels=newlabelorder_Site_ID)
alpha_plot + xlab("Site")

################################################# Alpha Diversity Metrics ###################################################################

alpha_diversity <- estimate_richness(dat_r, split=TRUE, measures=c("Shannon", "Simpson"))
alpha_div <- data.frame(alpha_diversity)
alpha_div

#Add Sample_ID column 
alpha_div["Sample_ID"] <- NA 
alpha_div$Sample_ID <- c(964:1005) #exact range for Water samples in metadata

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
SiteColours = c(REF1 = "#F8766D", REF2 = "#E38900", REF3 = "#C49A00", WAT = "#99A800", DSW1 = "#53B400", 
                DSW2 = "#00BC56", G3 = "#00C094", DSW3 = "#00BFC4", KIT = "#00B6EB", DSK1 = "#06A4FF", 
                DSK2 = "#A58AFF", J1 = "#DF70F8", DSK3 = "#FB61D7", DSK4 = "#FF66A8")

#Orders for boxplots
alpha_div_merge$Site_ID <- factor(alpha_div_merge$Site_ID, 
                                  levels=c("REF1", "REF2", "REF3", "WAT", "DSW1", "DSW2", "G3", "DSW3", "KIT", "DSK1", 
                                           "DSK2", "DSK3", "J1", "DSK4"))

alpha_div_merge$Location_vs_WWTP <- factor(alpha_div_merge$Location_vs_WWTP, 
                                           levels=c("Upstream", "Waterloo_WWTP", "Downstream_Waterloo", "Kitchener_WWTP", 
                                                    "Downstream_Kitchener"), 
                                           labels=c("Upstream", "W-WWTP", "Downstream Waterloo", "K-WWTP", 
                                                    "Downstream Kitchener"))

#Shannon Diversity boxplot
ggplot(data=alpha_div_merge, aes(x=Site_ID, y=Shannon)) + 
  geom_boxplot(stat="boxplot", fill="gray95", color=SiteColours) + 
  facet_grid(~Location_vs_WWTP, scale="free") + xlab("Site") + ylab("Shannon Diversity")
# + theme(strip.text = element_text(size = 11, color="black")) + theme(axis.text.x=element_text(angle=45, hjust=1, size=12, color="black"), axis.text=element_text(size=12, color="black"), axis.title=element_text(size=18))

#Simpson Diversity boxplot
ggplot(data=alpha_div_merge, aes(x=Site_ID, y=Simpson)) + 
  geom_boxplot(stat="boxplot", fill="gray95", color=SiteColours) + 
  facet_grid(~Location_vs_WWTP, scale="free") + xlab("Site") + ylab("Simpson Diversity")
# + theme(strip.text = element_text(size = 11, color="black")) + theme(axis.text.x=element_text(angle=45, hjust=1, size=12, color="black"), axis.text=element_text(size=12, color="black"), axis.title=element_text(size=18))


################################################### Beta Diversity ##################################################################

newlabelorder_Location = c("Upstream", "Waterloo_WWTP", "Downstream_Waterloo", "Kitchener_WWTP", "Downstream_Kitchener")
newlabelname_Location = c("Upstream", "W-WWTP", "Downstream Waterloo", "K-WWTP", "Downstream Kitchener")
newlabelorder_Site = c("REF1", "REF2", "REF3", "WAT", "DSW1", "DSW2", "G3", "DSW3", "KIT", "DSK1", 
                       "DSK2", "DSK3", "J1", "DSK4")

SiteColours = c(REF1 = "#F8766D", REF2 = "#E38900", REF3 = "#C49A00", WAT = "#99A800", DSW1 = "#53B400", 
                DSW2 = "#00BC56", G3 = "#00C094", DSW3 = "#00BFC4", KIT = "#00B6EB", DSK1 = "#06A4FF", 
                DSK2 = "#A58AFF", J1 = "#DF70F8", DSK3 = "#FB61D7", DSK4 = "#FF66A8")

# Transform data to proportions for Bray-Curtis distances
trans.dat <- transform_sample_counts(watdat, function(otu) otu/sum(otu))

# Data ordination
bray.ord = ordinate(trans.dat, method = "PCoA", distance = "bray")
bray.ord

# Plot Bray-Curtis
bray.plot = plot_ordination(trans.dat, bray.ord, color="Site_ID", title="Bray-Curtis Water", axes = c(1,2)) + 
  scale_color_manual(values = SiteColours, name="Site") 
bray.plot$data$Site_ID <- as.character(bray.plot$data$Site_ID)
bray.plot$data$Site_ID <- factor(bray.plot$data$Site_ID, levels=newlabelorder_Site)
bray.plot 

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

pairwise.adonis(bray.test, sample_data(watdat)$Location_vs_WWTP, perm = 99999) #by Location
pairwise.adonis(bray.test, sample_data(watdat)$Site_ID, perm = 99999) #by Site

################################################## Relative Abundance Plots ########################################################################

#colours
library(scales)
show_col(hue_pal()(29)) #ggplot colours, can adjust to amount of taxa. Didn't use but handy to have.

#Bar colours
phylaColours = c(Acidobacteria = "#c2c134", Actinobacteria = "#e84bd3", Bacteroidetes = "#e5e736", 
                  Chloroflexi = "#b16af2", Cyanobacteria = "#43e07b", Dependentiae = "#7488ec", Firmicutes = "#eab42c", 
                  Fusobacteria = "#5aac60", Chlamydiae = "#3ceabf", Epsilonbacteraeota = "#e8702d",
                  Nitrospirae = "#ea5da0", Patescibacteria = "#31abb2", Planctomycetes = "#ee603a", 
                  Proteobacteria = "#8fd9e5", Spirochaetes = "#d67ad6", Verrucomicrobia = "#c29937", 
                  `< 2% abund.` = "#A9A9A9")

########################################################
#Phylum - remove taxa in <2% abundance (normalize)
dat.aglo = tax_glom(watdat, taxrank = "Phylum")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
prune.dat = prune_taxa(taxa_sums(dat.trans) > 0.02, dat.trans)
dat.dataframe = psmelt(prune.dat)
dat.agr = aggregate(Abundance~Sample_Label+Site_ID+Location_vs_WWTP+Phylum, data=dat.dataframe, FUN=mean)

#Edit labels for plot
dat.agr$Site_ID= factor(dat.agr$Site_ID, 
                        levels=c("REF1", "REF2", "REF3", "WAT", "DSW1", "DSW2", "G3", "DSW3", "KIT", "DSK1", 
                                 "DSK2", "DSK3", "J1", "DSK4"))
dat.agr$Location_vs_WWTP= factor(dat.agr$Location_vs_WWTP, 
                                 levels=c("Upstream", "Waterloo_WWTP", "Downstream_Waterloo", "Kitchener_WWTP", 
                                          "Downstream_Kitchener"), 
                                 labels=c("Upstream", "W-WWTP", "Downstream Waterloo", "K-WWTP", 
                                          "Downstream Kitchener"))

ggplot(dat.agr, aes(x=Site_ID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity") + 
  facet_grid(~Location_vs_WWTP, scales = "free_x", space = "free_x") + 
  scale_fill_manual(values=phylaColours) +
  labs(y= "Abundance >2%", x= "Site")

#Relative abundance proportions (Ex. Cyanobacteria)
with(dat.agr[dat.agr$Phylum == "Cyanobacteria", ], table(Abundance, Site_ID))

#For the normalized data (>2%), what is in the lower 2%? (Abundance column) - added to plot caption
dat.aglo = tax_glom(watdat, taxrank = "Phylum")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.frame = psmelt(dat.trans)
dat.frame
dat.frame[200:300, 1:4]
tail(dat.frame, 10)

#########################################################
#Phylum - display legend category for <2% abundant taxa
dat.aglo = tax_glom(watdat, taxrank = "Phylum")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~Sample_Label+Site_ID+Location_vs_WWTP+Phylum, data=dat.dataframe, FUN=mean)
dat.agr$Phylum <- as.character(dat.agr$Phylum)
#simple way to rename phyla with < 2% abundance
dat.agr$Phylum[dat.agr$Abundance < 0.02] <- "< 2% abund."

#Edit labels for plot
unique(dat.agr$Phylum) #what will show up in legend
dat.agr$Phylum <- factor(dat.agr$Phylum, levels = c("Actinobacteria", "Bacteroidetes", "Cyanobacteria", "Dependentiae", 
                                                    "Epsilonbacteraeota", "Firmicutes", "Nitrospirae", "Patescibacteria",
                                                    "Proteobacteria", "< 2% abund."))
dat.agr$Site_ID= factor(dat.agr$Site_ID, 
                        levels=c("REF1", "REF2", "REF3", "WAT", "DSW1", "DSW2", "G3", "DSW3", "KIT", "DSK1", 
                                 "DSK2", "DSK3", "J1", "DSK4"))
dat.agr$Location_vs_WWTP= factor(dat.agr$Location_vs_WWTP, 
                                 levels=c("Upstream", "Waterloo_WWTP", "Downstream_Waterloo", "Kitchener_WWTP", 
                                          "Downstream_Kitchener"), 
                                 labels=c("Upstream", "W-WWTP", "Downstream Waterloo", "K-WWTP", 
                                          "Downstream Kitchener"))

ggplot(dat.agr, aes(x=Site_ID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity") + 
  facet_grid(~Location_vs_WWTP, scales = "free_x", space = "free_x") + 
  scale_fill_manual(values=phylaColours) +
  labs(y= "Abundance", x= "Site")

##################################################### Cyanobacteria ############################################################

#Subset to Cyanobacteria
cyano = watdat %>% subset_taxa(Phylum == "Cyanobacteria")
cyano

genusColours = c("#76d756","#e757c9","#d3da48","#6c81f4","#75af53","#b671e9","#64d79a","#4a9de3","#e96e4e","#44dcd1",
                 "#888ee2","#d69439","#d686cd","#c3b85d","#e9638b")

#Genus
dat.aglo = tax_glom(cyano, taxrank = "Genus")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~Site_ID+Location_vs_WWTP+Genus, data=dat.dataframe, FUN=mean)

dat.agr$Site_ID= factor(dat.agr$Site_ID, 
                        levels=c("REF1", "REF2", "REF3", "WAT", "DSW1", "DSW2", "G3", "DSW3", "KIT", "DSK1", 
                                 "DSK2", "DSK3", "J1", "DSK4"))
dat.agr$Location_vs_WWTP= factor(dat.agr$Location_vs_WWTP, 
                                 levels=c("Upstream", "Waterloo_WWTP", "Downstream_Waterloo", "Kitchener_WWTP", 
                                          "Downstream_Kitchener"), 
                                 labels=c("Upstream", "W-WWTP", "Downstream Waterloo", "K-WWTP", 
                                          "Downstream Kitchener"))

ggplot(dat.agr, aes(x=Site_ID, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity") + 
  facet_grid(~Location_vs_WWTP, scales = "free_x", space = "free_x") + 
  scale_fill_manual(values=genusColours) +
  labs(y="Abundance", x="Site")

# Relative abundance matrices (Ex. Tychonema_CCAP_1459-11B)
with(dat.agr[dat.agr$Genus == "Tychonema_CCAP_1459-11B", ], table(Abundance, Site_ID))

#Set working directory appropriately
setwd("~/Documents/PurdueLabStuff/Collaborations/Koziol/Koziol_mothuroutput/KoziolV4output/")

#Load libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(tidyr)


#Assign varibles for the paths of the data to import
sharedfile <- "Phylo.tx.1.pick.1.subsample.shared" #Phylo table
taxfile <- "Phylo.cons.taxonomy" #Phylo tax
metadata <- "Koziol_metadata.txt" #metadata
mdsfile <- "Tables/nmds.txt" #nmds for sorting

#first the otu table file
otu_subsample <- read.table(sharedfile, sep = "\t", header = T) #Read the OTU table into R
otu_subsample <- otu_subsample[,c(-1,-3)] #deleting unneeded columns
otu_subsample <- separate(data = otu_subsample, col = Group, into = c("sequencing_id", "Group"), sep = "_") #removing sequencing ID
otu_subsample <- otu_subsample[-c(1,2),-1] #deleting unneeded rows and columns
#setting rownames to sampleID and removing that column from array contents
row.names(otu_subsample) <- otu_subsample$Group
otu_subsample <- otu_subsample[,-1]


#then the taxonomy file
taxonomy <- read.table(taxfile, sep = "\t", header = T) #Read the taxonomy file into R
taxonomy <- separate(data = taxonomy, col = Taxonomy, into = c("kingdom", "phylum", "class", "family", "order", "genus", "species"), sep = ";")
str(taxonomy)

#then the metadata file
meta <- read.table(file = metadata, sep = '\t', header = TRUE)
str(meta)
meta$age <- factor(meta$age, levels = c("1", "2", "3", "4", "5+"))
str(meta)


#then the nmds file to merge with the meta data
nmds <- read.table(file = mdsfile, sep ='\t', header = TRUE)
metanmds <- merge(nmds, meta, by.x="group", by.y="group")


# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
   "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

#If you want different taxonomic level, find and replace the taxonomic level listed here
rm(otu.summary)
otu.summary <- prop.table(as.matrix(otu_subsample), 1) 
otu_abund <- colSums(otu.summary)
otu.summary <- rbind(otu_abund, otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]


#top 16 most abundant genera
num_genera <- 16 # enter the number of genera you want

melt_otu <- melt(otu.summary_sorted[,c(1:num_genera)])
colnames(melt_otu) <- c("Sample", "OTU", "Abundance")
tail(melt_otu)


#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu <- merge(metanmds, melt_otu, by.x = "group", by.y = "Sample")
meta_otu_tax <- merge(meta_otu, taxonomy)
str(meta_otu_tax)
summary(meta_otu_tax$group)
#sorting based on MDS1 from negative to positive (NMDS axis 1)
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$MDS1),]
#ordering samples based on NMDS axis 1
meta_otu_tax$group <- factor(meta_otu_tax$group, levels=unique(as.character(meta_otu_tax$group)) )


#MAKE A GRAPH! Plot individuals not group means
ggplot(meta_otu_tax, aes(x = group, y = Abundance, fill = genus)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_blank()) +
  ylab(paste0("Relative Abundance (top ", num_genera, " genera)")) +
  ggtitle("Genus Composition by Bull sorted by NMDS Axis 1") 
ggsave("GenusBarPlotNMDS1_Phylo_pwpt.png", width = 7, height = 3)




#making graph of next 16 most abundant genera
num_genera <- 16 # enter the number of genera you want

melt_otu32 <- melt(otu.summary_sorted[,c(17:32)])
colnames(melt_otu32) <- c("Sample", "OTU", "Abundance")
tail(melt_otu32)


#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu32 <- merge(metanmds, melt_otu32, by.x = "group", by.y = "Sample")
meta_otu32_tax <- merge(meta_otu32, taxonomy)
str(meta_otu32_tax)
#sorting based on MDS1 from negative to positive (NMDS axis 1)
meta_otu32_tax <- meta_otu32_tax[order(meta_otu32_tax$MDS1),]
#ordering samples based on NMDS axis 1
meta_otu32_tax$group <- factor(meta_otu32_tax$group, levels=unique(as.character(meta_otu32_tax$group)) )


ggplot(meta_otu32_tax, aes(x = group, y = Abundance, fill = genus)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=10)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_blank()) +
  ylab(paste0("Relative Abundance (next ", num_genera, " genera)"))+
  ggtitle("Genus Composition by Bull sorted by NMDS Axis 1") 
ggsave("GenusBarPlot_17to32_NMDS1_Phylo_pwpt.png", width = 8, height = 3)



#save plot to pdf for manuscript
pdf(file="Graphs/GenusBarPlots.pdf", height=4, width=10)
ggplot(meta_otu_tax, aes(x = group, y = Abundance, fill = genus)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab(paste0("Relative Abundance (top ", num_genera, " genera)")) +
  ggtitle("Genus Composition by Bull sorted by NMDS Axis 1") 
ggplot(meta_otu32_tax, aes(x = group, y = Abundance, fill = genus)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab(paste0("Relative Abundance (second ", num_genera, " top genera)")) +
  ggtitle("Genus Composition by Bull sorted by NMDS Axis 1") 
dev.off()



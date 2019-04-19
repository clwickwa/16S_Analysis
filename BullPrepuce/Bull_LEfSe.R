library(MASS)
library(tidyr)
library(ggplot2)
library(dplyr)


setwd("~/Documents/PurdueLabStuff/Collaborations/Koziol/Koziol_mothuroutput/KoziolV4output")
summary_lefse <- "Phylo.tx.1.pick.1.subsample.1.lefse_summary"
taxfile <- "Phylo.cons.taxonomy"

#read in lefse data
lefse_out <- read.table(file=summary_lefse, sep="\t", header=T)

#filter out unnwanted data - water samples, those without LDA, and p-values > 0.0001
lefse_noH2O <- lefse_out[lefse_out$Class != "water",]
lefse_NMDS <- lefse_noH2O[lefse_noH2O$Class != "-",]
lefse_NMDS_filt <- lefse_NMDS[which(lefse_NMDS$pValue <= 0.0001),]



#read in the taxonomy file and separate the taxa into columns
taxonomy<- read.table(file=taxfile, sep="\t", header=T)
taxonomy <- separate(data = taxonomy, col = Taxonomy, into = c("kingdom", "phylum", "class", "family", "order", "genus", "species"), sep = ";")

#merge the lefse and taxonomy by OTU
lefse_tax <- merge(lefse_NMDS_filt, taxonomy, by.x="OTU", by.y="OTU")

#order the data.frame by logLDA
lefse_tax$genus<- reorder(lefse_tax$genus, lefse_tax$LDA)
str(lefse_tax)
lefse_tax$Class <- factor(lefse_tax$Class, levels = c("LD", "HD"))
diversity <- c("Low Diversity", "High Diversity")
names(diversity) <- c("LD", "HD")


#plot and save for manuscript
#pdf(file="Bull_Lefse.pdf", width = 6, height = 3)
ggplot(lefse_tax, aes(x = genus, y = LDA, fill = Class)) +
  geom_bar(stat = "identity") +
  facet_grid(Class~., labeller = labeller(Class=diversity)) +
  #scale_fill_manual(values = my_colors) +
  # Remove x axis title
  #theme(axis.title.x = element_blank()) +
  ylab(paste0("Linear Discriminant Analysis Effect Size Score"))+
  #ylim(c(0,1)) +
  #guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5)) +
  #theme(legend.text=element_blank()) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 10)) +
  coord_flip()
ggsave(file="lefse_.png", width=6, height=6)
#dev.off()




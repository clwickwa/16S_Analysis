library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)

##### FUNCTIONS #####

pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method, permutations = 9999);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


###### DATA #####

setwd("~/Documents/PurdueLabStuff/Collaborations/Koziol/Koziol_mothuroutput/KoziolV4output/")
otu_table <-"OTU.opti_mcc.0.03.pick.0.03.subsample.shared" #rarefied OTU table
metadata <- "Koziol_metadata.txt" #metadata

#Read in OTU table
otu_subsample <- read.table(otu_table, header = TRUE)

#Remove sequence ID and leave just sample names
otu_subsample <- separate(data = otu_subsample, col = Group, into = c("sequencing_id", "Group"), sep = "_")
otu_subsample <- otu_subsample[-c(1,2),-2] #remove unneeded rows and columns

#Stores the sample name info as the rownames of the dataframe rather
rownames(otu_subsample) <- otu_subsample$Group

#Read in metadata
meta <- read.table(file = metadata, sep = '\t', header = TRUE)

#Makes sure that the meta table and the otu table have the same samples
meta <- meta[meta$group %in% rownames(otu_subsample),] 
otu_subsample <- otu_subsample[rownames(otu_subsample) %in% meta$group,]

# removes extra info that mothur includes in their OTU tables and outlier points
otu_subsample <- otu_subsample[,-c(1:3)]  

##################################################
##################################################

# this calculates the distance matrix using Bray-Curtis distances with vegan 
dist.matr.bray <- vegdist(otu_subsample, method = 'bray')

# this is vegan's function to make an NMDS ordination using k=2 dimensions

mds <- metaMDS(dist.matr.bray, k = 2,trymax = 1000, autotransform = FALSE)

#Calculation of the irdination stress

mds$stress

#merging MDS and metadata

nmds <-as.data.frame(mds$points)
nmds$group <- rownames(nmds)
metanmds <- merge(meta, nmds, by.x = 'group', by.y = 'group')
metanmds$age <- factor(metanmds$age)
str(metanmds)
write.csv(metanmds)

###### PLOTS #######

#General plots with basic facets

ggplot(metanmds, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=age))
ggplot(metanmds, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=breed))
ggplot(metanmds, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=diet))

ggplot(metanmds, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=breed, shape=age))


#Plot for manuscript

#pdf(file="NMDS.pdf",height = 5, width = 7)
ggplot(metanmds, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=breed, shape=age)) +
  labs(x='Axis 1', y= 'Axis 2', caption = paste('Ordination stress: ', round(mds$stress, digits = 2)))
ggsave(file="NMDS.png", height=4, width=5)

#dev.off()


#Save MDS data for later
write.table(nmds, file="Tables/nmds.txt")

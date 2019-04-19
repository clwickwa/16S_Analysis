# This is a demo for running the co-occurrence analysis much, much faster

#make sure you have these libraries
library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)
library(GGally)
library(intergraph)
library(tidyr)


setwd("~/Documents/PurdueLabStuff/Collaborations/Koziol/Koziol_mothuroutput/KoziolV4output/")
#Read in OTU table and meta data or design file
OTU <- read.table(file = "OTU.opti_mcc.0.03.pick.0.03.subsample.shared", sep = "\t", header = T)
meta <- read.table(file = "Koziol_metadata.txt", sep = "\t", header = T)

#load taxonomy file for later
taxonomy <- read.table (file = "OTU.cons.taxonomy", sep = "\t", header =T)
taxonomy <- separate(data = taxonomy, col = Taxonomy, into = c("kingdom", "phylum", "class", "family", "order", "genus", "species"), sep = ";")

#Split OTU table "Group" into two columns to match "group" in meta data
OTU <- separate(data = OTU, col = Group, into = c("sequencing_id", "Group"), sep = "_") 
#removing unnecessary columns and rows
OTU <- OTU[-c(1,2),-c(1,2,4)]
#merging OTUs with metadata table
dataset <- merge(meta, OTU, by.x = "group", by.y = "Group")


# we are going to create a network per treatment
datasetn<-dataset
datasetn[datasetn==0]<-NA
head(dataset[,1:12])
head(datasetn[,1:12])

occurrences <- colSums(datasetn[,-1:-11], na.rm=TRUE)
str(datasetn$Otu0001)

#NMDS_side<-as.vector(unique(dataset$X))
#NMDS_side
#NMDS_side[1]


rm(i)

#Can use loop if subsetting by treatment
final_results<-data.frame()
#for(i in 1:length(NMDS_side)){
	#subset the data for a particular treatment YOU MUST ENTER THE HEADER OF THE COLUMN THAT HAS THE DATA IN THIS CASE “NMDS_side”
	#temp<-subset(dataset, X==NMDS_side[i])
	#tempn<-subset(datasetn, X==NMDS_side[i])
	# making an object that has all the results in it (both rho and P values)
	results<-rcorr(as.matrix(dataset[,-c(1:9)]),type="spearman")
	resultsn<-rcorr(as.matrix(datasetn[,-c(1:9)]),type="spearman")
	#make two seperate objects for p-value and correlation coefficients
	rhos<-results$r
	ps<-results$P
	ns<-resultsn$n
	# going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
	ps_melt<-na.omit(melt(ps))
	#creating a qvalue based on FDR
	ps_melt$qval<-fdrtool(ps_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
	#making column names more relevant
	names(ps_melt)[3]<-"pval"
	# if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
	ps_sub<-subset(ps_melt, qval < 0.05)
	# now melting the rhos, note the similarity between ps_melt and rhos_melt
	rhos_melt<-na.omit(melt(rhos))
	names(rhos_melt)[3]<-"rho"
	ns_melt<-na.omit(melt(ns))
	names(ns_melt)[3]<-"ns"
	ns_sub<- subset(ns_melt, ns > nrow(dataset)/10)
	#merging together filter by rho value
	merged<-merge(ps_sub,subset(rhos_melt, abs(rho) > 0.5),by=c("Var1","Var2"))
	merged<-merge(merged, ns_sub, by=c("Var1","Var2"))
	#merged$X<-NMDS_side[i]
	final_results<-rbind(final_results, merged)
	#print(paste("finished ",NMDS_side[i],sep=""))
#}


Var1_tax <- merge(final_results, taxonomy, by.x = "Var1", by.y = "OTU")
Var1_tax[,c("Phylum1","Class1","Family1","Order1","Genus1")] <- Var1_tax[,c(9:13)]
Var1_tax<- Var1_tax[,-c(7:14)]
final_tax <- merge(Var1_tax, taxonomy, by.x = "Var2", by.y = "OTU")
final_tax[,c("Phylum2","Class2","Family2","Order2","Genus2")] <- final_tax[,c(14:18)]
final_tax<- final_tax[,-c(12:19)]

final_tax <- final_tax[,c(7,8,5,1:4,6)]
write.csv(final_tax, file = "Bull_correlation.csv")


fuso <- subset(final_tax, Genus1=="Fusobacterium(100)")
porph <- subset(final_tax, Genus1=="Porphyromonas(100)")
histo <- subset(final_tax, Genus1=="Histophilus(100)")
myco <- subset(final_tax, Genus1=="Mycoplasma(100)")
ecoli <- subset(final_tax, Genus1=="Escherichia(100)")
bacillus <- subset(final_tax, Genus1=="Bacillus(100)")
parvi <- subset(final_tax, Genus1=="Parvimonas(100)")
strepto <- subset(final_tax, Genus1=="Streptobacillus(100)")
brady <- subset(final_tax, Genus1=="Bradyrhizobium(100)")
rumino <- subset(final_tax, Genus1=="Ruminococcaceae_unclassified(100)")

#save tables for manuscript
write.csv(ecoli, file = "Ecoli_Correlation.csv")
write.csv(fuso, file = "Fusobacterium_Correlation.csv")
write.csv(porph, file = "Porphyromonas_Correlation.csv")
write.csv(histo, file = "Histophilus_Correlation.csv")
write.csv(myco, file = "Mycoplasma_Correlation.csv")
write.csv(bacillus, file = "Bacillus_Correlation.csv")
write.csv(parvi, file = "Parvimonas_Correlation.csv")
write.csv(strepto, file = "Streptobacillus_Correlation.csv")
write.csv(brady, file = "Bradyrhizobium_Correlation.csv")
write.csv(rumino, file = "Ruminococcaceae_Correlation.csv")

#run with subsetting

metanmds <- read.table("NMDS_sorted_summary.txt", sep='\t', header=T)

dataset2 <- merge(metanmds, OTU, by.x = "group", by.y = "Group")
dataset2 <- dataset2[,-c(2,4:19)] #getting rid of everything except for sample names, NMDS column, and OTU data

# we are going to create a network per treatment
datasetn2<-dataset2
datasetn2[datasetn2==0]<-NA
head(dataset2[,1:12])
head(datasetn2[,1:12])

occurrences2 <- colSums(datasetn2[,-c(1,2)], na.rm=TRUE)
str(datasetn$Otu0001)

NMDS_side<-as.vector(unique(dataset2$NMDS))
NMDS_side
NMDS_side[1]


rm(i)

#Can use loop if subsetting by treatment
final_results2<-data.frame()
for(i in 1:length(NMDS_side)){
#subset the data for a particular treatment YOU MUST ENTER THE HEADER OF THE COLUMN THAT HAS THE DATA IN THIS CASE “NMDS_side”
temp<-subset(dataset2, NMDS==NMDS_side[i])
tempn<-subset(datasetn2, NMDS==NMDS_side[i])
# making an object that has all the results in it (both rho and P values)
results<-rcorr(as.matrix(temp[,-c(1:2)]),type="spearman")
resultsn<-rcorr(as.matrix(tempn[,-c(1:2)]),type="spearman")
#make two seperate objects for p-value and correlation coefficients
rhos<-results$r
ps<-results$P
ns<-resultsn$n
# going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
ps_melt<-na.omit(melt(ps))
#creating a qvalue based on FDR
ps_melt$qval<-fdrtool(ps_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
#making column names more relevant
names(ps_melt)[3]<-"pval"
# if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
ps_sub<-subset(ps_melt, qval < 0.05)
# now melting the rhos, note the similarity between ps_melt and rhos_melt
rhos_melt<-na.omit(melt(rhos))
names(rhos_melt)[3]<-"rho"
ns_melt<-na.omit(melt(ns))
names(ns_melt)[3]<-"ns"
ns_sub<- subset(ns_melt, ns > nrow(dataset2)/10)
#merging together filter by rho value
merged<-merge(ps_sub,subset(rhos_melt, abs(rho) > 0.5),by=c("Var1","Var2"))
merged<-merge(merged, ns_sub, by=c("Var1","Var2"))
merged$NMDS<-NMDS_side[i]
final_results2<-rbind(final_results2, merged)
print(paste("finished ",NMDS_side[i],sep=""))
}



Var1_tax2 <- merge(final_results2, taxonomy, by.x = "Var1", by.y = "OTU")
Var1_tax2$Genus1 <- Var1_tax2$genus
Var1_tax2<- Var1_tax2[,-c(6,8:15)]
final_tax2 <- merge(Var1_tax2, taxonomy, by.x = "Var2", by.y = "OTU")
final_tax2$Genus2 <- final_tax2$genus
final_tax2<- final_tax2[,-c(8:15)]

final_tax2 <- final_tax2[,c(7,8,5,6,1:4)]
write.csv(final_tax2, file = "Bull_correlation_withNMDS.csv")


fuso2 <- subset(final_tax2, Genus1=="Fusobacterium(100)")
porph2 <- subset(final_tax2, Genus1=="Porphyromonas(100)")
histo2 <- subset(final_tax2, Genus1=="Histophilus(100)")
myco2 <- subset(final_tax2, Genus1=="Mycoplasma(100)")
ecoli2 <- subset(final_tax2, Genus1=="Escherichia(100)")
bacillus2 <- subset(final_tax2, Genus1=="Bacillus(100)")
parvi2 <- subset(final_tax2, Genus1=="Parvimonas(100)")
strepto2 <- subset(final_tax2, Genus1=="Streptobacillus(100)")
brady2 <- subset(final_tax2, Genus1=="Bradyrhizobium(100)")
rumino2 <- subset(final_tax2, Genus1=="Ruminococcaceae_unclassified(100)")











#put these in individual tables


library(ggplot2)
temp<-subset(dataset, X==NMDS_side[2])
qplot(Otu0002, Otu0006, data = temp)

merged_fortim<-merge(ps_sub,rhos_melt,by=c("Var1","Var2"))
absolute <- subset(rhos_melt, rho > abs(0.5))
final_tax$Genus1["Porphyromonas(100)"]

# # to make a network from the results
# # we need to pass the relevant relationships (pairs across columns Var1 and Var2) to graph.edgelist as a matrix, note I am subsetting for a particular treatment within the arguments
# head(final_results)
# g<-graph.edgelist(as.matrix(subset(final_results, diet=="W")[,1:2]),directed=F)
# 
# #are relationships different
# ggplot(final_results)+
#   geom_density(aes(rho,fill=diet),alpha=0.5)+
#   theme_bw(base_size=17)+
#   theme(aspect.ratio=1)+scale_fill_manual(name="Diet",values=c(1:7))
# 
# # now we can calculate stats for the network
# final_stats<-data.frame()
# for(i in 1:length(unique(final_results$diet))){
# 	temp<-subset(final_results, diet==as.vector(unique(final_results$diet))[i])
# 	temp.graph<-(graph.edgelist(as.matrix(temp[,c(1,2)]),directed=FALSE))
# 	E(temp.graph)$weight<-temp$rho
# 	temp.graph<-simplify(temp.graph)
# 	stats<-data.frame(row.names((as.matrix(igraph::degree(temp.graph,normalized=TRUE)))),(as.matrix(igraph::degree(temp.graph,normalized=TRUE))),(as.matrix(igraph::betweenness(temp.graph))))
# 	names(stats)<-c("otus","norm_degree","betweenness")
# 	stats$diet<-as.vector(unique(final_results$diet))[i]
# 	stats$clustering_coeff<-igraph::transitivity(temp.graph,type="global")
# 	stats$clustering_coeff_rand<-igraph::transitivity(igraph::erdos.renyi.game(length(V(temp.graph)),length(E(temp.graph)),type="gnm"))
# 	stats$cluster_ratio<-stats$clustering_coeff/stats$clustering_coeff_rand
# 	final_stats<-rbind(final_stats,stats)
# 	print(paste("finished ",as.vector(unique(final_results$diet))[i],sep=""))
# }
# 
# # is the structure changing?
# ddply(final_stats,.(diet),summarise,mean(cluster_ratio))
# 
# #WHICH TAXA HAVE HIGHEST DEGREE OR BETWEENNESS, ENTER THE TREATMENT TYPE
# head(arrange(subset(final_stats,diet=="non-foaming"), -norm_degree),10)
# OTU 166 d:Bacteria,p:Firmicutes,c:Negativicutes,o:Selenomonadales,f:Acidaminococcaceae
# OTU 2080 d:Bacteria,p:"Bacteroidetes"
# OTU 1864 d:Bacteria
# 
# #ENTER THE TREATMENT TYPE
# head(arrange(subset(final_stats,diet=="non-foaming"), -betweenness),10)
# OTU 1864 d:Bacteria
# OTU 1601  d:Bacteria,p:"Proteobacteria",c:Betaproteobacteria
# 
# # whats the relationship between these statistics?
# ggplot(final_stats)+geom_point(aes(x=norm_degree,y=betweenness,color=diet),alpha=0.5)+scale_y_log10()+theme_bw(base_size=17)+labs(x="Normalized Degree",y="Betweenness")+theme(aspect.ratio=1)+scale_colour_manual(name="Diet",values=c(1:7))
# 



# let's make plots with the most significant relationships
head(final_results)
dim(final_results)
hist(final_results$rho)
#YOU CAN CHOSE TO LOOK AT ONLY THE STRONGEST RELATIONSHIPS. YOU CAN CHANGE THE RHO USED
strong_results<-subset(final_results, rho >= 0.7)
dim(strong_results)
hist(strong_results$rho)

#ENTER THE CORRECT TREATMENT AND THEN RERUN WITH ALL OTHER TREATMENTS
temp.graph<-(graph.edgelist(as.matrix(subset(strong_results, X=="neg")[,c(1,2)]),directed=FALSE))
#ENTER THE CORRECT TREATMENT
E(temp.graph)$weight<-subset(strong_results, X=="neg")$rho
temp.graph<-simplify(temp.graph)
temp.graph

library(GGally)
library(intergraph)
library(network)
gnet<-asNetwork(temp.graph)
df<-asDF(gnet)
head(df$vertexes)

ggnet(gnet, size=0, method="kamadakawaii")+geom_point()

# lets color by phylum YOU HAVE TO HAVE A FILE THAT CATEGORIZES THE GENES OR OTUS BASED ON SOME COMMON FACTOR. I USED MGE, ARG, PHYLA. YOU COULD CATEGORIZE ON ARG TYPE OR WHATEVER YOU WANT.
vs<-df$vertexes
phyla<-read.table("~/Desktop/otu_phylum.txt",sep="\t",check.names=F,header=T)
head(phyla)
head(vs)
vs_phyla<-merge(vs, phyla, by.x="vertex.names",by.y="otus")
vs_phyla<-arrange(vs_phyla,intergraph_id)
head(vs_phyla)
##YOU HAVE TO HAVE THE CORRECT NUMBER OF COLORS INDICATED IN THE LAST SET OF PARENTHESIS OR THERE WILL BE AN ERROR.
ggnet(gnet, size=0, method="kamadakawaii")+geom_point(aes(colour=vs_phyla$phylum))+scale_colour_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","grey50"))
##YOU CAN ALSO GET LABELS INCLUDED WITH:
ggnet(gnet, label.nodes=T, size=0, method="kamadakawaii")+geom_point(aes(colour=vs_phyla$phylum))+scale_colour_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","grey50"))




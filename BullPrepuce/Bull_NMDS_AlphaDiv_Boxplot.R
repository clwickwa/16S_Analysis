library(ggplot2)

setwd("~/Documents/PurdueLabStuff/Collaborations/Koziol/Koziol_mothuroutput/KoziolV4output/")

#Read in alpha diversity table
V4summary <- read.table("OTU.opti_mcc.groups.summary", sep="\t", header=T)
V4summary <- separate(data = V4summary, col = group, into = c("sequencing_id", "Group"), sep = "_") #removing sequencing ID
V4summary <- V4summary[,c(-2)] #deleting unneeded columns

#Read in meta data table
meta <- read.table("Koziol_metadata.txt", sep="\t", header=T)

#Read in mds table
MDStable <- read.table("Tables/nmds.txt", sep="\t", header=T)

#merge tables together
V4summary_meta <- merge(meta, V4summary, by.x="group", by.y="Group")
V4_mds_sum_meta <- merge(MDStable, V4summary_meta, by.x = "group", by.y="group")
V4_mds_sum_meta <- V4_mds_sum_meta[order(V4_mds_sum_meta$MDS1),]
V4_mds_sum_meta$NMDS <- if_else(V4_mds_sum_meta$MDS1 >0, "positive","negative") 

#save graph for manuscript
#pdf(file="NMDS_shannon.pdf", width=4, height=4)
ggplot(V4_mds_sum_meta, aes(NMDS, shannon)) + 
  geom_boxplot(aes(color = NMDS)) +
  ylab(paste0("Shannon")) +
  theme(legend.position="none")
ggsave(file="NMDS_shannon.png", width = 3, height=3)

#dev.off()


#t-test of negative NMDS vs postive NMDS (left and right clusters) for Shannon
Shannon_ttest <- t.test(NMDS_sort_summary$shannon[1:42], NMDS_sort_summary$shannon[43:84], paired = TRUE)







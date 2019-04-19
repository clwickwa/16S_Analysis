library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("~/Documents/PurdueLabStuff/Collaborations/Koziol/Koziol_mothuroutput/KoziolV4output/")

#Rarefaction curve
rm(V4rfact)

meta <- read.table(file = "Koziol_metadata.txt", sep = "\t", header = TRUE)

#Before reading in the rarefaction file, delete unwanted samples. 
#I used Excel for this and used the search function for unwanted samples

V4rfact <- read_tsv(file = "OTU.opti_mcc.groups.rarefaction")%>%
  select(-contains("lci-"), -contains("hci-")) %>%
  gather(-numsampled, key=sample, value=coverage) %>%
  mutate(sample=str_replace_all(sample, pattern="0.03-", replacement="")) %>%
  drop_na()
V4rfact <- separate(data = V4rfact, col = sample, into = c("sequencing_id", "sample"), sep = "_") #removing sequencing ID
V4rfact <- V4rfact[,-2]

#merge rarefaction file with metadata
V4meta_rare <- meta%>%
  #sample_n(20) %>%
  merge(., V4rfact, by.x= "group", by.y = "sample")

#graph rarefaction plot with vertical line where subsampling cutoff is

#pdf("Rarefaction.pdf", width = 5, height = 4)
ggplot(V4meta_rare, aes(x=numsampled, y=coverage, group=group, color=age)) +
  geom_line()+
  geom_vline(xintercept=3000) +
  coord_cartesian(xlim=c(0,20000)) +
  labs(x="Number of Sequences Sampled per Subject",
       y="Number of OTUs per Subject") +
  theme_classic()
ggsave(file="Rarefaction_pwpt.png",width=4, height=3)
#dev.off()


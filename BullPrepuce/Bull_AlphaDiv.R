library(readr)
library(stringr)
library(ggplot2)
library(tidyr)

setwd("~/Documents/PurdueLabStuff/Collaborations/Koziol/Koziol_mothuroutput/KoziolV4output/")
alpha_div <- read.table(file = "OTU.opti_mcc.groups.summary", sep = "\t", header = T)
meta <- read.table(file = "Koziol_metadata.txt", sep = "\t", header = TRUE)

#removing sequencing ID and leaving just sample name
alpha_div <- separate(data = alpha_div, col = group, into = c("sequencing_id", "Group"), sep = "_")
alpha_div <- alpha_div[,c(-2)]

alpha_div_merge <- merge(meta, alpha_div, by.x = "group", by.y = "Group")
unique(alpha_div_merge$breed)
str(alpha_div_merge)
alpha_div_merge$age <- factor(alpha_div_merge$age, levels = c("1", "2", "3", "4", "5", "6+")) #setting age as a factor instead of int
str(alpha_div_merge)

#checking boxplots
qplot(age, chao, geom = "boxplot", colour = age, data = alpha_div_merge, size = I(0.3))
qplot(breed, chao, geom = "boxplot", colour = breed, data = alpha_div_merge, size = I(0.3))
qplot(age, shannon, geom = "boxplot", colour = age, data = alpha_div_merge, size = I(0.3))
qplot(breed, shannon, geom = "boxplot", colour = breed, data = alpha_div_merge, size = I(0.3))


#Get figures for manuscript.

#pdf(file="AlphaDiv.pdf", width=3, height=3)


ggplot(alpha_div_merge, aes(age, chao)) + 
  geom_boxplot(aes(color = age)) + 
  ylim(c(0,450)) +
  ylab(paste0("Chao1")) +
  theme(legend.position="none")
ggsave("AlphaDiv_chao.png", height = 3, width=3)
ggplot(alpha_div_merge, aes(age, shannon)) + 
  geom_boxplot(aes(color = age)) + 
  ylim(c(0,4)) +
  ylab(paste0("Shannon")) +
  theme(legend.position="none")
ggsave("AlphaDiv_shannon.png", height = 3, width=3)


#dev.off()


##Run the ANOVAs for statistics
breeds <- unique(alpha_div_merge$breed)
age_breeds <- unique(alpha_div_merge$age_breed)
ages <- unique(alpha_div_merge$age)


#Age
ad_metrics <- c("chao", "shannon")
for(b in breeds){
  print(b)
    for(m in ad_metrics){
      print(m)
      aov_temp <- aov(get(m) ~ age, data = subset(alpha_div_merge, breed == b))
      summary(aov_temp)
      anova_summary <- as.data.frame(summary(aov_temp)[[1]])
      write.table(anova_summary, file = paste0("Tables/anova_shannon", b, ".txt"), sep = "\t", quote = FALSE)
      if (summary(aov_temp)[[1]][["Pr(>F)"]][[1]]){
        tukey_out <- TukeyHSD(aov_temp)
        tukey_out_df <- as.data.frame(tukey_out$age)
        tukey_out_df$breed <- b
        tukey_out_df$ad_metric <- m
        if (exists("tukey_summary")) {
          tukey_summary <- rbind(tukey_summary, tukey_out_df)
        } else {
          tukey_summary <- tukey_out_df
        }
      }
    }
  }

tukey_summary$q.value <- p.adjust(tukey_summary$`p adj`, method = "BH")
write.table(tukey_summary, file = "Tables/tukey_summaryAge.txt", sep = "\t", quote = FALSE)


#Breed
ad_metrics <- c("chao", "shannon")
for(t in ages){
  print(t)
  for(m in ad_metrics){
    print(m)
    aov_temp <- aov(get(m) ~ breed, data = subset(alpha_div_merge, age == t))
    summary(aov_temp)
    anova_summary <- as.data.frame(summary(aov_temp)[[1]])
    write.table(anova_summary, file = paste0("Tables/anova_shannon", t, ".txt"), sep = "\t", quote = FALSE)
    if (summary(aov_temp)[[1]][["Pr(>F)"]][[1]]){
      tukey_out <- TukeyHSD(aov_temp)
      tukey_out_df <- as.data.frame(tukey_out$breed)
      tukey_out_df$age <- t
      tukey_out_df$ad_metric <- m
      if (exists("tukey_summary")) {
        tukey_summary <- rbind(tukey_summary, tukey_out_df)
      } else {
        tukey_summary <- tukey_out_df
      }
    }
  }
}

tukey_summary$q.value <- p.adjust(tukey_summary$`p adj`, method = "BH")
write.table(tukey_summary, file = "Tables/tukey_summaryBreed.txt", sep = "\t", quote = FALSE)




library(ggplot2)
library(dplyr)
library(reshape)

# TRAMIC diff abundance clr_genus_urine_only
# this file has the data_collapsed_genus_urine file as an input for clr in R




#make sure that taxa have numbers before the name for designated order; replace the sampleIDs with group types
data <-read.table("Bacteroides.csv", header=TRUE, sep=",")
melted_proportion <- melt(data)
write.csv(melted_proportion, file ="melted_Bacteroides.csv")
#adjust melted file: delete numbers in "variable", change labels in "id", "taxon", "variable", "value" 

data <- read.csv("melted_Bacteroides.csv", header=TRUE, sep=",", row.names="id")
data$variable <- factor(data$variable,levels = c("yes_adults","no_adults", "yes_elderly", "no_elderly", "yes_CEN", "no_CEN"))

#define group colors
group.colors <- c(yes_adults = "darkred", no_adults = "azure3", yes_elderly = "darkred", no_elderly = "azure3", yes_CEN = "darkred", no_CEN = "azure3")

genus_H <- ggplot(data, aes(x = taxon, y = (value), fill = variable)) +
  geom_boxplot() + 
  labs(fill = "group") + 
  geom_point(position=position_jitterdodge(),alpha=0.3, size=1) +
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  scale_fill_manual(values=group.colors)


genus_H +labs(x="", y="Bacteroides ", title="Differential abundance Bacteroides")


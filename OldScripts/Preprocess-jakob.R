# preproces data

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(stringr)
library(foreach)

###########################     Load Data     #################################

# Load dataframe with qPCR data. Note that Biorep1 i removed as it is an outlier as compared to the other:
data.qpcr <- as.data.frame(read.csv(file = "Data/Concat_QPCR.csv")[,-1])
data.qpcr1 <- data.qpcr[data.qpcr$Biorep != "biorep1",]

data.qpcr2 <- data.qpcr1 %>% group_by(Sample, Lineage, Culture) %>%
  summarise(k = mean(k), 
            t.gen = mean(t.gen), 
            endpoint = mean(Endpoint), 
            copynumb = mean(Copynumb), 
            ratio = mean(Ratio))


data.tsf <- read.csv(file = "Data/TosseForsøg.csv")[,-1]

str(data.cat)
data.cat <- merge(data.qpcr2, data.tsf[,c(4:8)], 
                  by.x = "Sample", by.y = "Sample", all.x = T)
data.cat <- data.cat %>% separate(Sample, c("ID", "Color"), sep = "-", remove = F)
#data.cat <- data.cat %>% unite(Special ,c("Culture", "Color"), sep = "_", remove = F)

data.omit_LL <- data.select

# Group and label lists
data.Culture <- data.cat$Culture
data.color <- data.cat$Color
data.special <- data.cat$Special
data.sample <- data.cat$Sample
data.linage <- data.cat$Lineage


##################   Outlier removal   ###########################
# names(data.cat)
# # Subset data.frame 
# data.select <- data.cat[,c(1,3:8,10,11,13,14)]
# 
# data.select$Cult.bin <- ifelse(data.select$Culture == "Coc",1,0) 
# data.select$Col.bin <- ifelse(data.select$Color == "l",1,0) 
# 
# # Outlier detection by cooks distance 
# # https://www.r-bloggers.com/outlier-detection-and-treatment-with-r/
# mod <- lm(ratio ~ ., data=data.select[,c(5:11)])
# cooksd <- cooks.distance(mod)
# 
# plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
# abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
# text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
# 
# # Outlier removal by Car test
# car::outlierTest(mod)
# 
# # Remove outlier
# data.select2 <- data.select[-c(16,22,40,70),]

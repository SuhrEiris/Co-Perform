# preproces data

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(stringr)
library(foreach)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut],
#     p = pmat[ut]
#   )
# }
# 

###########################     Load Data     #################################

# Load dataframe with qPCR data. Note that Biorep1 i removed as it is an outlier as compared to the other:
data.qpcr <- as.data.frame(read.csv(file = "Data/Concat_QPCR.csv")[,-1])
data.qpcr1 <- data.qpcr[data.qpcr$Biorep != "biorep1",]

# data.qpcr2 <- data.qpcr1 %>% group_by(Sample, Lineage, Culture) %>% 
#   summarise(k = mean(k), t.gen = mean(t.gen), endpoint = mean(Endpoint), copynumb = mean(Copynumb), ratio = mean(Ratio))


data.tsf <- read.csv(file = "Data/TosseForsøg.csv")[,-1]

str(data.cat)
data.cat <- merge(data.qpcr2, data.tsf[,c(4:8)], by.x = "Sample", by.y = "Sample", all.x = T)
data.cat <- data.cat %>% separate(Sample, c("ID", "Color"), sep = "-", remove = F)
#data.cat <- data.cat %>% unite(Special ,c("Culture", "Color"), sep = "_", remove = F)


###########################    PCA approach    ###############################

#library(ggbiplot)

# Group and label lists
data.Culture <- data.cat$Culture
data.color <- data.cat$Color
data.special <- data.cat$Special
data.sample <- data.cat$Sample
data.linage <- data.cat$Lineage


data.cat2 <- data.cat[5:13]
mtcars.pca <- prcomp(na.omit(data.cat2), center = TRUE,scale. = TRUE)
ggbiplot::ggbiplot(mtcars.pca)

# Overview plot 
ggbiplot::ggbiplot(mtcars.pca, ellipse = T, groups = data.special, labels = data.color)

data.cat3 <- data.cat2[,c(1,2,5,6,8)]
pca.reduc <- prcomp(data.cat3, center = TRUE,scale. = TRUE)

#
p <- ggbiplot(pca.reduc, ellipse = T, groups = data.Culture, labels = data.sample) + 
  geom_point(aes(shape = data.color))+
  theme_bw()
p

#
p <- ggbiplot(pca.reduc, ellipse = T, groups = data.special) + 
  theme_bw()
p

##################   Outlier removal   ###########################
names(data.cat)
# Subset data.frame 
data.select <- data.cat[,c(1,3:8,10,11,13,14)]

data.select$Cult.bin <- ifelse(data.select$Culture == "Coc",1,0) 
data.select$Col.bin <- ifelse(data.select$Color == "l",1,0) 

# Outlier detection by cooks distance 
# https://www.r-bloggers.com/outlier-detection-and-treatment-with-r/
mod <- lm(ratio ~ ., data=data.select[,c(5:11)])
cooksd <- cooks.distance(mod)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

# Outlier removal by Car test
car::outlierTest(mod)

# Remove outlier
data.select2 <- data.select[-c(16,22,40,70),]


mod <- lm(ratio ~ ., data=data.select2[,c(5:11)])
cooksd <- cooks.distance(mod)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels


### ----- Correlation of variables
# Simple Correlations on white variants
library("PerformanceAnalytics")

chart.Correlation(data.select2[,c(5:11)], histogram=TRUE, pch=19)

#library(car)
#pairs(~k+t.gen+Endpoint+Copynumb+Ratio,data=white[,c(5:9)], main="Simple Scatterplot Matrix")

### ----- Advanced correlations on white variants
# Make a correlation matrix with adjusted p-values, using library(igraph)
simple <- Hmisc::rcorr(as.matrix(data.select2[,c(5:11)]))
simple2 <- flattenCorrMatrix(simple$r, simple$P)
simple2$correct <- p.adjust(simple2$p, method = "fdr")
simple3 <- simple2[,c(1,2,3,5)]

g <- igraph::graph.data.frame(simple3[,c(1,2,3)], directed=FALSE)
a <- igraph::get.adjacency(g, attr="cor", sparse=FALSE)
table(a[,1]==simple$r[,1])  
a[a == 0] <- 1

g1 <- graph.data.frame(simple3[,c(1,2,4)], directed=FALSE)
b <- get.adjacency(g1, attr="correct", sparse=FALSE)
table(colnames(b) == colnames(a))   
b[b == 0] <- NA

corrplot(a, type = "upper", 
         tl.col = "black", tl.srt = 45, 
         p.mat = b, sig.level = c(.001, .01, .05), pch.cex = 2, insig = "label_sig", pch.col = "black",
         title = "Correlation matrix - white variants",
         mar=c(0,0,2,0))

######################      PCA cleaned data     ##########################

pca.clean <- prcomp(data.select[,c(5:9)], center = TRUE,scale. = TRUE)
ggbiplot(pca.clean)

# Overview plot 
ggbiplot(mtcars.pca, ellipse = T, groups = data.special, labels = data.color)

data.omit <- data.select
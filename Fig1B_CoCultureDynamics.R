#########################################################################################

#Packages


# Load library for data import
library(readxl)

# load library for data management
library(dplyr)
library(tidyr)
library(ggplot2)

#########################################################################################

#Data management

# Set working directory
#setwd("~/Dropbox/MME_Co-Adatpt/Experiments_2018/Evo_1/Evo1_CFU")
#setwd("C:/Users/ncz500/Dropbox/MME_Co-Perform/Experiments_2018/Evo_1/Evo1_CFU")
setwd("~/Dropbox/MME_Co-Adatpt/Experiments_2018/Evo_1/Evo1_CFU")

dir()

# Load data (create CSV file etc)
cfu <- read_excel("Evo1_CFU_JAHC.xlsx")

#Average technical replicates
cfu_data.1 <- cfu %>% group_by(Day, Species, Replicate, Count_of) %>% summarise(Count_tech= n(), Mean_cfu = mean(CFU_mL))

#df with xl and xs to be merged 
cfu_data.3 <-cfu_data.1[cfu_data.1$Count_of == "Xl",]
cfu_data.4 <-cfu_data.1[cfu_data.1$Count_of == "Xs",]
cfu_data.5 <- rbind(cfu_data.3,cfu_data.4)

#adding xl and xs 
cfu_data.6 <- cfu_data.5 %>% group_by(Day, Species, Replicate) %>% summarise(Count_tech=n(), Mean_cfu=sum(Mean_cfu))
cfu_data.6$Count_of <- "X"



#dataset without xl and xs
cfu_data.7 <- cfu_data.1[ which(cfu_data.1$Count_of !='Xs' & cfu_data.1$Count_of !='Xl'), ]

#rbind count of x with rest
vec.names <- colnames(cfu_data.7)
cfu_data.6 <- cfu_data.6[,vec.names]
cfu_data.8 <- rbind(cfu_data.6,cfu_data.7)

######## Counts day 1 and day 16 plot#######
data1 <- cfu_data.8[!grepl("Z", cfu_data.8$Species),]
data1$Species[data1$Species=="X"] <- "LaL (Single)"
data1$Species[data1$Species=="Y"] <- "LeM (Single)"
data1$Species[data1$Species=="XY"] <- "Co-culture"
data1$Count_of[data1$Count_of=="X"] <- "LaL"
data1$Count_of[data1$Count_of=="Y"] <- "LeM"

data16 <- data1[data1$Day == "16",]
data1 <-  data1[data1$Day== "1",]

dataDay <- rbind(data16,data1)

plot.day <- ggplot(data=dataDay[dataDay$Day == "16",], aes(x=Species, y=Mean_cfu, fill=Species))+
  geom_boxplot(width = 0.5)+
  facet_grid(.~Count_of, scales = "free_y")+
  scale_fill_brewer()+
  theme(axis.text.x = element_text(size=24))+
  theme(axis.title = element_text(size = 24))+
  theme(axis.text = element_text(size = 20, color = "Black"))+
  ylab ("CFU per mL")+
  theme(legend.position="bottom")+
  theme_bw()
plot.day

library("ggpubr")
#Test day16 CFU of LeM in co and single
#ay16.lem <- dataDay[grepl("16", dataDay$Day),]
#day16 <- day16.lem[grepl("LeM", day16.lem$Count_of),]
#t.test(day16$Mean_cfu ~ day16$Species)

#Test day16 CFU of LaL in co and single
Day16.lal <- dataDay[grepl("16", dataDay$Day),]
Day16.lal2 <- Day16.lal[grepl("LaL", Day16.lal$Count_of),]
t.test(Day16.lal2$Mean_cfu ~ Day16.lal2$Species)

########
dataDay2 <- dataDay %>% group_by(Day, Species, Count_of) %>% summarize(Count_bio= n(), CFU = mean(Mean_cfu), Std.dev_bio = sd(Mean_cfu))
dataco <- dataDay[grepl("Co-culture", dataDay$Species),]

#Average biological replicates
log_data = cfu_data.8
log_data$Mean_cfu = log10(log_data$Mean_cfu)

# rep A at day 13 and rep C at day 7 is zero (plating problem) taken out
cfu_data.8 = cfu_data.8[-38,]
cfu_data.8 = cfu_data.8[-65,]

log_data = cfu_data.8
log_data$Mean_cfu = log10(log_data$Mean_cfu)


cfu_data.9 <- log_data %>% group_by(Day, Species, Count_of) %>% summarize(Count_bio= n(), CFU = mean(Mean_cfu), Std.dev_bio = sd(Mean_cfu))
#cfu_data.9 <- cfu_data.8 %>% group_by(Day, Species, Count_of) %>% summarise(Count = n(), Mean = mean(Mean_cfu), test.sd = sd(Mean_cfu))
#write.table(cfu_data.9, file = "cfu_data_9.csv", row.names = TRUE, sep = ";")


#############   Plot Jakob  ################
pres <- cfu_data.9[!grepl("Z", cfu_data.9$Count_of),]
pres2 <- cfu_data.9[!grepl("Z", cfu_data.9$Species),]

pres2$Species[pres2$Species=="X"] <- "LaL (Single)"
pres2$Species[pres2$Species=="Y"] <- "LeM (Single)"
pres2$Species[pres2$Species=="XY"] <- "Co-culture"
pres2$Count_of[pres2$Count_of=="X"] <- "LaL"
pres2$Count_of[pres2$Count_of=="Y"] <- "LeM"

# SE for CI intervals

pres2$se = pres2$Std.dev_bio/sqrt(pres2$Count_bio)
pres2$CI = qnorm(0.975)*pres2$Std.dev_bio/sqrt(pres2$Count_bio)

#library("wesanderson")

pres.plot <- ggplot(data=pres2, aes(x=Day, y= CFU, fill = Count_of, color = Count_of))+
  geom_point(size=2)+
  geom_line(size=1)+
  facet_grid(~Species)+
  geom_ribbon(aes(ymax = CFU + CI, ymin = CFU - CI, alpha = 0.2, size = NA)) +
  labs(x="\nTime (days)",y="CFU per mL\n") + 
  scale_x_continuous(breaks = seq(1,18,by=3))+
  scale_y_continuous(breaks = seq(4,8, by =1))+
  #labs(shape = "Species") +
  #labs(color = "Species") +
  #scale_y_log10() +
  scale_color_manual(values= wes_palette("Moonrise2"))+
  scale_fill_manual(values = wes_palette("Moonrise2"))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_text(size = 20), 
        legend.text = element_text(size=20))+
  theme(axis.text = element_text(size = 18, color = "Black"), axis.title = element_text(size = 24),
        axis.text.x = element_text(size=20))
pres.plot

#################

# T.test different days 

# Make dataset nice and tidy

CFU.ttest = log_data[!grepl("Z", log_data$Count_of),]
CFU.ttest=  log_data[!grepl("Z", log_data$Species),]
CFU.ttest$Species[CFU.ttest$Species=="X"] <- "LaL (Single)"
CFU.ttest$Species[CFU.ttest$Species=="Y"] <- "LeM (Single)"
CFU.ttest$Species[CFU.ttest$Species=="XY"] <- "Co-culture"
CFU.ttest$Count_of[CFU.ttest$Count_of=="X"] <- "LaL"
CFU.ttest$Count_of[CFU.ttest$Count_of=="Y"] <- "LeM"

CFU.ttest.co = CFU.ttest[grepl("Co-culture", CFU.ttest$Species),]



# t.test of co-culture populations
d1.co = CFU.ttest[CFU.ttest$Day == "1",]
t.test(d1.co$Mean_cfu ~ d1.co$Count_of)

d4.co = CFU.ttest[CFU.ttest$Day == "4",]
t.test(d4.co$Mean_cfu ~ d4.co$Count_of)

d7.co = CFU.ttest[CFU.ttest$Day == "7",]
t.test(d7.co$Mean_cfu ~ d7.co$Count_of)

d10.co = CFU.ttest[CFU.ttest$Day == "10",]
t.test(d10.co$Mean_cfu ~ d10.co$Count_of)

d13.co = CFU.ttest[CFU.ttest$Day == "13",]
t.test(d13.co$Mean_cfu ~ d13.co$Count_of)

d16.co = CFU.ttest[CFU.ttest$Day == "16",]
t.test(d16.co$Mean_cfu ~ d16.co$Count_of)


# t.test of lal sin-culture populations
CFU.ttest.mono = CFU.ttest[!grepl("Co-culture", CFU.ttest$Species),]

d1.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "1",]
t.test(d1.mono$Mean_cfu ~ d1.mono$Count_of)

d4.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "4",]
t.test(d4.mono$Mean_cfu ~ d4.mono$Count_of)

d7.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "7",]
t.test(d7.mono$Mean_cfu ~ d7.mono$Count_of)

d10.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "10",]
t.test(d10.mono$Mean_cfu ~ d10.mono$Count_of)

d13.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "13",]
t.test(d13.mono$Mean_cfu ~ d13.mono$Count_of)

d16.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "16",]
t.test(d16.mono$Mean_cfu ~ d16.mono$Count_of)


# t.test of lal-culture populations co and single
CFU.ttest.lal = CFU.ttest[grepl("LaL", CFU.ttest$Count_of),]

d1.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "1",]
t.test(d1.lal$Mean_cfu ~ d1.lal$Species)

d4.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "4",]
t.test(d4.lal$Mean_cfu ~ d4.lal$Species)

d7.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "7",]
t.test(d7.lal$Mean_cfu ~ d7.lal$Species)

d10.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "10",]
t.test(d10.lal$Mean_cfu ~ d10.lal$Species) # not below

d13.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "13",]
t.test(d13.lal$Mean_cfu ~ d13.lal$Species)

d16.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "16",]
t.test(d16.lal$Mean_cfu ~ d16.lal$Species)



# t.test of lemthin-culture populations co and singler
CFU.ttest.lem = CFU.ttest[grepl("LeM", CFU.ttest$Count_of),]

d1.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "1",]
t.test(d1.lem$Mean_cfu ~ d1.lem$Species)

d4.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "4",] # not sig
t.test(d4.lem$Mean_cfu ~ d4.lem$Species)

d7.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "7",]
t.test(d7.lem$Mean_cfu ~ d7.lem$Species)

d10.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "10",]
t.test(d10.lem$Mean_cfu ~ d10.lem$Species)

d13.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "13",] # not sig
t.test(d13.lem$Mean_cfu ~ d13.lem$Species)

d16.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "16",]
t.test(d16.lem$Mean_cfu ~ d16.lem$Species)


#Morphotypes plot 

cfu_data.1.1 <- cfu_data.1[which(cfu_data.1$Species !="XZ" & cfu_data.1$Species != "Z"),]
cfu_data.x <- cfu_data.1.1[cfu_data.1.1$Species == "X",]
cfu_data.xy <- cfu_data.1.1[cfu_data.1.1$Species == "XY",]
cfu_data.xy <- cfu_data.xy[cfu_data.xy$Count_of != "Y",]

morphotype.plot.x <- ggplot(data=cfu_data.x, aes(x=Day, y= Mean_cfu, color = Count_of))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(~Replicate)+
  labs(x="Time (days)",y="CFU/mL") + 
  scale_x_continuous(breaks = seq(1,20,by=2))

morphotype.plot.x

morphotype.plot.xy <- ggplot(data=cfu_data.xy, aes(x=Day, y= Mean_cfu, fill = Count_of))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(~Replicate)+
  labs(x="Time (days)",y="CFU per mL") + 
  scale_x_continuous(breaks = seq(1,20,by=2))+
  scale_color_brewer()+
  theme_bw()+
  theme(legend.position = "bottom")
  

morphotype.plot.xy


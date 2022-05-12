### Growthcurves for wildtype La. lactis and Ln. mesenteroides ###

#Nathalie - MAC
setwd("~/Dropbox/MME_Co-Adatpt/Experiments_2018/Evo_1/Variant_data/GrowthCurves")

# Packs

#install.packages("growthcurver")
library("growthcurver")
# Install updated version of GrowthCurver
#install_github("Russel88/growthcurver", force = TRUE)

library("stringr")
library("broom")
library("multcomp")

# load data
biorep1 <- as.data.frame(read_excel("2018.09.26_VariantGrowthCurve_test1.xlsx", sheet = 3))
biorep4 <- as.data.frame(read_excel("2018.10.14_VariantGrowthCurve_test4.xlsx", sheet = 3))
biorep5 <- as.data.frame(read_excel("2018.10.17_VariantGrowthCurve_test5.xlsx", sheet = 3))

bio1 <- biorep1
bio4 <- biorep4
bio5 <- biorep5

bio1$Time <- as.numeric(bio1$Time)
bio4$Time <- as.numeric(bio4$Time)
bio5$Time <- as.numeric(bio5$Time)

bio1 <- filter(bio1, Time < 1140)
bio4 <- filter(bio4, Time < 1140)
bio5 <- filter(bio5, Time < 1140)

biorep1_gc <- SummarizeGrowthByPlate(bio1, plot_fit = TRUE, plot_file = "Biorep1.pdf")
biorep4_gc <- SummarizeGrowthByPlate(bio4, plot_fit = TRUE, plot_file = "Biorep4.pdf")
biorep5_gc <- SummarizeGrowthByPlate(bio5, plot_fit = TRUE, plot_file = "Biorep5.pdf")

biorep1_gc$biorep <- "biorep1"
biorep4_gc$biorep <- "biorep4"
biorep5_gc$biorep <- "biorep5"

data <- rbind(biorep1_gc, biorep4_gc, biorep5_gc)


data<- data[!grepl("blank", data$sample),]
datawild <- data[grepl("xw", data$sample),]
datawild$sample <- "LaL"

datawild2 <- data[grepl("yw", data$sample),]
datawild2$sample <- "LeM"

datawild2.2 <-rbind(datawild, datawild2)
datawild2.2$sample<-as.factor(datawild2.2$sample)

# Growth measures from growthcurver
data3 <- datawild2.2 %>% dplyr::group_by(sample) %>% dplyr::summarise(mean_CC = mean(k),
                                                                      mean_GT =mean(t_gen),
                                                                      sd_CC = sd(k),
                                                                      sd_GT = sd(t_gen))


# biorep 1 #
biorep1.new <- biorep1 %>% gather(Species, OD, -Time)
biorep1.new2 <- biorep1.new %>% separate(Species, into = c("Species", "Replicate"), sep = "_")
biorep1.new3 <- biorep1.new2 %>% group_by(Species, Time) %>% summarise(Count = n(), Average= mean(OD))
biorep1.new4 <- biorep1.new3

biorep1.xw <- biorep1.new4[grepl("xw", biorep1.new4$Species),]
biorep1.xw$Species <- str_sub(biorep1.xw$Species,1,2)
biorep1.xw2 <- biorep1.xw %>% group_by(Species, Time) %>% summarize(Count = n(), Average = mean(Average))
biorep1.new5 <-  rbind(biorep1.new4[!grepl("xw", biorep1.new4$Species),], biorep1.xw2)
biorep1.new5$biorep <- "biorep1"

# biorep4
biorep4.new <- biorep4 %>% gather(Species, OD, -Time)
biorep4.new2 <- biorep4.new %>% separate(Species, into = c("Species", "Replicate"), sep = "_")
biorep4.new3 <- biorep4.new2 %>% group_by(Species, Time) %>% summarise(Count = n(), Average= mean(OD))
biorep4.new4 <- biorep4.new3

biorep4.xw <- biorep4.new4[grepl("xw", biorep4.new4$Species),]
biorep4.xw$Species <- str_sub(biorep4.xw$Species,1,2)
biorep4.xw2 <- biorep4.xw %>% group_by(Species, Time) %>% summarize(Count = n(), Average = mean(Average))
biorep4.new5 <-  rbind(biorep4.new4[!grepl("xw", biorep4.new4$Species),], biorep4.xw2)
biorep4.new5$biorep <- "biorep4"

# biorep5
biorep5.new <- biorep5 %>% gather(Species, OD, -Time)
biorep5.new2 <- biorep5.new %>% separate(Species, into = c("Species", "Replicate"), sep = "_")
biorep5.new3 <- biorep5.new2 %>% group_by(Species, Time) %>% summarise(Count = n(), Average= mean(OD))
biorep5.new4 <- biorep5.new3

biorep5.xw <- biorep5.new4[grepl("xw", biorep5.new4$Species),]
biorep5.xw$Species <- str_sub(biorep5.xw$Species,1,2)
biorep5.xw2 <- biorep5.xw %>% group_by(Species, Time) %>% summarize(Count = n(), Average = mean(Average))
biorep5.new5 <-  rbind(biorep5.new4[!grepl("xw", biorep5.new4$Species),], biorep5.xw2)
biorep5.new5$biorep <- "Biorep5"

dat.comb <- rbind(biorep1.new5,
                  biorep4.new5,
                  biorep5.new5)

# lal and lem only 

data.wild <- dat.comb[grepl("w", dat.comb$Species),]
data.wild <- data.wild[!grepl("z", data.wild$Species),]
data.wild$Species[data.wild$Species=="xw"] <- "LaL"
data.wild$Species[data.wild$Species=="yw"] <- "LeM"

library("wesanderson")
library(Rmisc)

#Average Biorep (n = 3)
wildtypedata <- data.wild %>% 
  dplyr::group_by(Species, Time) %>% 
  dplyr::summarise(mean_OD = mean(Average),
                   sd_OD =sd(Average),
                   reps = n())
wildtypedata$CI =qnorm(0.975)*wildtypedata$sd_OD/sqrt(wildtypedata$reps) 

wildtypedata$Time <- as.numeric(wildtypedata$Time)
wildtypedata <- filter(wildtypedata, Time < 1140)
wildtypedata$hour <- wildtypedata$Time/60

wildtypedata$Species[wildtypedata$Species=="LaL"] <- "L. lactis"
wildtypedata$Species[wildtypedata$Species=="LeM"] <- "L. mesenteroides"


#FigS1B 

PS1B <- wildtypedata %>% 
  ggplot(aes(x = hour, y = mean_OD, color = Species))+
  geom_line(size = 1)+
  scale_fill_manual(values= wes_palette(name="Moonrise2"))+
  scale_color_manual(values= wes_palette(name="Moonrise2"))+
  geom_ribbon(aes(ymax = mean_OD + CI, ymin = mean_OD - CI, fill = Species), alpha = 0.2, color = NA) +
  theme_bw(base_size = 8)+
  scale_x_continuous(breaks = seq(0,16,by=4), limits=c(0, 16))+ 
  labs(x="\nTime(hour)", y= "OD600nm\n")+
  theme(legend.position = c(0.8,0.15),
        legend.title = element_blank(), 
        legend.text = element_text(face = "bold.italic"),
        legend.background = element_rect(colour = "Black"),
        axis.text = element_text(color = "Black", face = "bold"), 
        axis.title = element_text(color = "Black", face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
PS1B

#setwd("~/Desktop/PhD_nasuh/PhD_exp/Projects/Co-Perform")
ggsave("FigS1B.pdf",
       device = "pdf",
       plot = PS1B,
       units="mm",
       width=88,
       height=88,
       dpi = 300)



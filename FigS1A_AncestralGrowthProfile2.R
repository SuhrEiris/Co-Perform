### Growthcurves for wildtype La. lactis and Ln. mesenteroides ###

# Data
GrowthData <- as.data.frame(read_excel("Data/2022.05.11_growth_curves_LaL_LeM_anc_mono_co_cultures.xlsx", sheet = 3))

# Wrangle data

GrowthData2 = gather(GrowthData, Strain, OD, LaL_1:Cocu_4, factor_key = TRUE) %>%
  separate(Strain, c('Strain', 'BioRep'))

# Summarise
GrowthData3<- GrowthData2 %>% 
  dplyr::group_by(Time, Strain) %>% 
  dplyr::summarise(mean_OD = mean(OD),
                   sd_OD = sd(OD),
                   count = n())
GrowthData3$CI =qnorm(0.975)*GrowthData3$sd_OD/sqrt(GrowthData3$count) 
GrowthData3$Time = as.numeric(GrowthData3$Time)
GrowthData3<- filter(GrowthData3, Time < 15)
library("wesanderson")

GrowthData3$Strain[GrowthData3$Strain=="LaL"] <- "L. lactis"
GrowthData3$Strain[GrowthData3$Strain=="LeM"] <- "L. mesenteroides"
GrowthData3$Strain[GrowthData3$Strain=="Cocu"] <- "Co-culture"

#FigS1B 
GrowthData3$Strain <- factor(GrowthData3$Strain, 
                             levels = c("L. lactis", "L. mesenteroides", "Co-culture"))

PS1B <- GrowthData3 %>% 
  ggplot(aes(x = Time, y = mean_OD, color = Strain))+
  geom_line(size = 1.5)+
  scale_fill_manual(values= wes_palette(name="Moonrise2"))+
  scale_color_manual(values= wes_palette(name="Moonrise2"))+
  geom_ribbon(aes(ymax = mean_OD + CI, ymin = mean_OD - CI, fill = Strain), alpha = 0.2, color = NA) +
  theme_bw(base_size = 8)+
  #scale_x_continuous(breaks = seq(0,16,by=4), limits=c(0, 16))+ 
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



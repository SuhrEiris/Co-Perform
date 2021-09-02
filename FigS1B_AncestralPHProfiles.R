
ph_raw <- read.csv("Data/phprofiles.csv", header = T, sep = ";")

################################################################

## Data management ##

# Divide data
ph_species <- ph_raw[ph_raw$Species != "MRS",]
ph_MRS <- ph_raw[ph_raw$Species == "MRS",]

# Average blank
ph_MRS_timeblank <- ph_MRS %>% 
  group_by(Time) %>% 
  dplyr::summarize(Count_Blank= n(), 
            mean_blank_time = mean(pH), 
            sd_blank_time = sd(pH))

#Subtract mean_blank_time from time_blank (get difference per time per biorep)
ph_MRS$ph_difference <- (ph_MRS$pH - ph_MRS_timeblank$mean_blank_time)

#ad/subtract difference from pH per time per biorep
ph_species$ph_normalized <- (ph_species$pH + ph_MRS$ph_difference)

#####################################################################################################

library("wesanderson")
#Average all
ph_MRS1 <- ph_raw %>% group_by(Species, Time) %>% dplyr::summarize (Count = n(),
                                                    Mean_ph = mean(pH),                                                 
                                                    std.dev = sd(pH))

ph_MRS1$CI = qnorm(0.975)*ph_MRS1$std.dev/sqrt(ph_MRS1$Count)
ph_MRS1$Species <- as.character(ph_MRS1$Species)
ph_MRS1$Species[ph_MRS1$Species=="481"] <- "L. lactis"
ph_MRS1$Species[ph_MRS1$Species=="343"] <- "L. mesenteroides"
ph_MRS1$Species[ph_MRS1$Species=="Co"] <- "Co-culture"
ph_MRS1$Species[ph_MRS1$Species=="MRS"] <- "50%MRS"


# Plot profiles
ph_MRS1$Species <- factor(ph_MRS1$Species, levels = c("L. lactis", "L. mesenteroides", "Co-culture", "50%MRS"))

PS1A <- ggplot(ph_MRS1, aes(x = Time, y = Mean_ph, color = Species)) +
  geom_line (size =1.5) +
  geom_ribbon(aes(ymax = Mean_ph + CI, ymin = Mean_ph - CI, fill = Species), alpha = 0.2, colour = NA)+
  #geom_errorbar(aes(ymax = Mean_ph + std.dev, ymin = Mean_ph - std.dev), width = 0.3) + 
  ylab("pH\n")+
  xlab("\nTime (hour)") +
  theme_bw(base_size = 8)+
  scale_color_manual(values= c("#798E87", "#C27D38","#CCC591", "#29211F"))+
  scale_fill_manual(values= c("#798E87", "#C27D38","#CCC591", "#29211F"))+
  theme(axis.title = element_text(color = "Black", face = "bold"),
        axis.text = element_text(color = "Black", face = "bold"),
        legend.position = c(0.8, 0.55),
        legend.title = element_blank(), 
        legend.text = element_text(colour="Black", face = "bold"),
        legend.background = element_rect(colour = "Black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
PS1A 

ggsave("FigS1A.pdf",
       device = "pdf",
       plot = PS1A,
       units="mm",
       width=88,
       height=88,
       dpi = 300)

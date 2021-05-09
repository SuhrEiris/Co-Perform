
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
ph_MRS1$Species[ph_MRS1$Species=="481"] <- "LaL"
ph_MRS1$Species[ph_MRS1$Species=="343"] <- "LeM"
ph_MRS1$Species[ph_MRS1$Species=="co"] <- "Co-culture"


# Plot profiles
plot.profiles <- ggplot(ph_MRS1, aes(x = Time, y = Mean_ph, color = Species)) +
  geom_line (size =1) +
  geom_ribbon(aes(ymax = Mean_ph + CI, ymin = Mean_ph - CI, fill = Species), alpha = 0.2, colour = NA)+
  #geom_errorbar(aes(ymax = Mean_ph + std.dev, ymin = Mean_ph - std.dev), width = 0.3) + 
  ylab("pH\n")+
  xlab("\nTime (Hours)") +
  theme_bw()+
  #scale_fill_discrete(name = "Species", labels = c("LaL", "LeM", "Co-culture", "MRS"))+
  scale_color_manual(values= c("#CCC591", "#798E87", "#C27D38", "#29211F"))+
  scale_fill_manual(values= c("#CCC591", "#798E87", "#C27D38", "#29211F"))+
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 14, color = "Black"),
        legend.position = "right",legend.title = element_blank(), legend.text = element_text(colour="Black", size=14), legend.key.size = unit(1, "cm"))
  #guides(color=guide_legend(override.aes = list(size = 3)))
plot.profiles  

####################################################################################################


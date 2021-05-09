
evoplank <- read_excel("Data/df.plankexperimentCFU.xlsx")
evobead <- read_excel("Data/df.beadexperimentCFU.xlsx")


# Bead evolution dataframe ####

#Filter for Day 16
bead_D16 <- evobead[evobead$Day == "16",]
bead_D16_XY <- filter(bead_D16, Species == "XY")

#Remove Xl rows
bead_D16_XY <- bead_D16_XY %>%
  dplyr::select(Day, Species, Replicate, Count_of, CFU_mL) %>%
  dplyr::filter(Count_of != "Xl")

#df with Y
bead_xy_y <- bead_D16_XY %>%
  filter(Count_of == "Y") %>%
  tibble::rownames_to_column(var ="Sample") %>%
  dplyr::rename(LeM.cf = CFU_mL)

#df with X
bead_xy_x <- bead_D16_XY %>%
  filter(Count_of == "Xs") %>%
  tibble::rownames_to_column(var="Sample") %>%
  dplyr::rename( LaL.cf = CFU_mL)


#Merge x and y df
bead_ratio <- cbind(bead_xy_x, bead_xy_y)
bead_ratio.t = bead_ratio[,-(1:5)]
bead_ratio.t = bead_ratio.t[,-(6)]


bead_ratio.t$Total.cf <- bead_ratio.t$LaL.cf + bead_ratio.t$LeM.cf
bead_ratio.t$Relative.abundance <- bead_ratio.t$LaL.cf/bead_ratio.t$Total.cf
bead_ratio.t$Percentage <- bead_ratio.t$Relative.abundance*100

bead_rat <- bead_ratio.t %>%
  mutate(Experiment = "BEAD-D16")

bead_rat <- dplyr::select(bead_rat, -c(Sample))

## GROUP TECHNICAL
bead_rat2 <- bead_rat %>%
  group_by(Species, Replicate, Experiment) %>%
  dplyr::summarise(LaL.cfu = mean(LaL.cf),
            LeM.cfu = mean(LeM.cf),
            Total.cfu = mean(Total.cf),
            Relative.abund = mean(Relative.abundance),
            Percent = mean(Percentage))

bead_rat2 <-tibble::as_tibble(bead_rat2)


# Planktonic evolution dataframe####

# 30 % of the co-cultures does not contain LaL (6/20)

plank_co <- filter(evoplank, Plate == "Co")

plank_lal <- filter(plank_co, Count_of == "LaL_l")
plank_lal <- plank_lal %>%
  dplyr::rename (LaL = CFU_ml) %>%
  dplyr::select (Species, LaL)

plank_lem <- filter(plank_co, Count_of == "LeM")
plank_lem <- plank_lem %>%
  dplyr::rename (LeM = CFU_ml) %>%
  dplyr::select(Species, LeM)

plank_co_ratio <- cbind(plank_lem, plank_lal, by= "Species")
plank_co_ratio = plank_co_ratio[,2:4] 

plank_co_ratio$totalCFU <- plank_co_ratio$LaL + plank_co_ratio$LeM
plank_co_ratio$relative.abundance <- plank_co_ratio$LaL/plank_co_ratio$totalCFU
plank_co_ratio$percentage <- plank_co_ratio$relative.abundance*100

plank_co_ratio <- plank_co_ratio %>%
  dplyr::rename(LeM.cfu = LeM,
         LaL.cfu = LaL,
         Relative.abund = relative.abundance,
         Percent = percentage,
         Total.cfu = totalCFU,
         Replicate = Species) %>%
  mutate(Experiment = "PLANK-D16") %>%
  mutate(Species = "XY")

# Merge Evo1 and Evoplank ####
relative.abundance <- rbind(bead_rat2, plank_co_ratio)

#create percentage of LeM and Total
relative.abundance$Percentage.LeM <- (relative.abundance$LeM.cfu/relative.abundance$Total.cfu)*100
relative.abundance$Percentage.total <- (relative.abundance$Total.cfu/relative.abundance$Total.cfu)*100

# plot only LaL percentage  ####

# means

relative.abundance.2 <- relative.abundance %>%
  group_by(Experiment) %>%
  dplyr::summarise(Count = n(),
            Mean.Percentage = mean(Percent),
            Sd.Percentage = sd(Percent))

my_comparisons <- list( c("BEAD-D16", "PLANK-D16"))

plot.percentage <- ggplot(relative.abundance, aes(x = Experiment, y = Percent, fill = Experiment)) +
  geom_boxplot(alpha= 0.6) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1, binwidth = 0.4, fill = "#798E87", alpha = 0.5) +
  stat_summary(fun = mean, color = "black", size = 20, geom = "point", shape = 3)+
  theme_bw()+
  ylim(0,11.5)+
  stat_boxplot(geom ='errorbar', width = 0) +
  xlab("") +
  theme(axis.text.x = element_text(size=24))+
  theme(axis.title = element_text(size = 24))+
  theme(axis.text = element_text(size = 24, color = "Black"))+
  ylab ("Amount La. lactis (%)") +
  scale_fill_manual(values= c("#798E87", "#798E87"))+
  theme(legend.position="none")
plot.percentage


# Assumptions for unpaired two-sample t.test and wilcoxon ####
#LaL
#normal distribution
with(relative.abundance, shapiro.test(Percent[Experiment == "PLANK-D16"])) #p = 0.0004413
with(relative.abundance, shapiro.test(Percent[Experiment == "BEAD-D16"])) # p = 0.2

var.test(Percentage ~ Experiment, data = relative.abundance) #p = 0.006983

# t. test (P = 0.00006691)
t.test <- t.test(Percentage ~ Experiment, data = relative.abundance, var.equal = FALSE)
print(t.test)

# Wilcox.test ( P = 0.0056)
wilcox.lal <- wilcox.test(Percentage ~ Experiment, data = relative.abundance,
                   exact = FALSE)
print(wilcox.lal)

#LeM
wilcox.lem <- wilcox.test(Percentage.LeM ~ Experiment, data = relative.abundance,
                      exact = FALSE)
print(wilcox.lem)

with(relative.abundance, shapiro.test(Percentage.LeM[Experiment == "Planktonic"])) #p = 0.0004413
with(relative.abundance, shapiro.test(Percentage.LeM[Experiment == "Bead"])) # p = 0.2

var.test(Percentage.LeM ~ Experiment, data = relative.abundance) #p = 0.006983

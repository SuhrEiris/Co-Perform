### pH analysis - day14 variants ###


#load data
pH_raw <- as.data.frame(read_excel("Data/pH_measurement_R.xlsx", sheet = 1))
str(pH_raw)

#### Preprocessing data ####

#prepare data
pH_raw$type = as.factor(pH_raw$type)
pH_raw$biorep = as.factor(pH_raw$biorep)
pH_raw$color = as.factor(pH_raw$color)
pH_raw$lineage = as.factor(pH_raw$lineage)
str(pH_raw)

# remove lem and biorep 4-7

pH_clean = pH_raw[pH_raw$isolate != "lem" & pH_raw$biorep != 4 & pH_raw$biorep != 5 & pH_raw$biorep != 6 & pH_raw$biorep != 7,]

#plot

ggplot(pH_clean, aes(x = type, y = ph, color = type)) +
  geom_boxplot()

# Testing (wilcoxon test)

with(pH_clean, shapiro.test(ph[type == "co"]))
hist(pH_clean[pH_clean$type == "co",]$ph)

with(pH_clean, shapiro.test(ph[type == "wt"]))
hist(pH_clean[pH_clean$type == "wt",]$ph)

with(pH_clean, shapiro.test(ph[type == "mono"]))
hist(pH_clean[pH_clean$type == "mono",]$ph)

wilcox.test(pH_clean[pH_clean$type == "co",]$ph, pH_clean[pH_clean$type== "wt",]$ph)
wilcox.test(pH_clean[pH_clean$type == "mono",]$ph, pH_clean[pH_clean$type== "wt",]$ph)
wilcox.test(pH_clean[pH_clean$type == "mono",]$ph, pH_clean[pH_clean$type== "co",]$ph)

#divide into biorep 
biorep1 = pH_clean[pH_clean$biorep == "1",]
biorep2 = pH_clean[pH_clean$biorep == "2",]
biorep3 = pH_clean[pH_clean$biorep == "3",]

# pick out lal wt 
lal_b1 = biorep1[biorep1$isolate == "lal",] %>%
  summarise(ph = mean(ph))
lal_b2 = biorep2[biorep2$isolate == "lal",] %>%
  summarise(ph = mean(ph))
lal_b3 = biorep3[biorep3$isolate == "lal",] %>%
  summarise(ph = mean(ph))

# Logfold change
biorep1$logfc_ph = log2(biorep1$ph/lal_b1$ph)
biorep2$logfc_ph = log2(biorep2$ph/lal_b2$ph)
biorep3$logfc_ph = log2(biorep3$ph/lal_b3$ph)

# merge bioreps
combined_ph = rbind(biorep1, biorep2, biorep3)
combined_ph = combined_ph %>%  unite("ID1", isolate:type , na.rm = TRUE, remove = FALSE)
combined_ph = combined_ph %>%  unite("ID", ID1,color , na.rm = TRUE, remove = FALSE)
combined_ph = na.omit(combined_ph)

#Logfold change plot 
ggplot(combined_ph, aes(x = ID, y = logfc_ph,  fill = type)) +
  geom_bar(stat = "identity")+
  facet_grid(biorep~.) + 
  theme_bw()

# remove lal from raw data that has not been average across bioreps
combined_ph = combined_ph[combined_ph$isolate != "lal",]

# make binary for culture type and color
combined_ph$Col.bin <- ifelse(combined_ph$color == "l",1,0)
combined_ph$Cult.bin <- ifelse(combined_ph$type == "co",1,0) 

# Average bioreps
ave_combined_ph = combined_ph %>%
    group_by(isolate, type, color, lineage, ID, Cult.bin, Col.bin) %>%
    dplyr::summarise(mean_logfc_ph = mean(logfc_ph),
                     sd_logfc_ph = sd(logfc_ph),
                     counts = n())

# Remove outlier
mod <- lm(mean_logfc_ph ~ ., data=as.data.frame(ave_combined_ph[,8]))
cooksd <- cooks.distance(mod)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

# Outlier removal by Car test
car::outlierTest(mod)
# No outlier

set.seed(41)
lm <- lme(mean_logfc_ph ~ Cult.bin*Col.bin, random = ~ 1|lineage, data=ave_combined_ph)
lm.sum <- summary(lm)
lm.intervals <- intervals(lm,which = "fixed")

#library(multcomp)
contr <- rbind("COCvsInter"=c(1,1,0,0)) 
glht <- multcomp::glht(lm, linfct=contr)
glht.confint <- confint(glht)
glht.sum <- summary(glht)

# Make a summary dataframe for plots with single and co-culture 
dfx.ph <- data.frame(rbind(lm.intervals$fixed[1, ], glht.confint$confint[c(2,1,3)]))
df.ph.psing <- data.frame(lm.sum$tTable[1,5])
colnames(df.ph.psing) <- "pvalue"
df.ph.pcoc <- data.frame(glht.sum$test$pvalues)
colnames(df.ph.pcoc) <- "pvalue"
df.ph.pboth <- rbind(df.ph.psing,df.ph.pcoc)
dfx.ph2 <- cbind(dfx.ph,df.ph.pboth)
dfx.ph2$Var <- "Acidification (pH)"
dfx.ph2$Eff <- c("Mono-culture", "Co-culture")
dfx.ph2$plot <- "Acidification (pH)"

# Make a summary dataframe for plots with Culture and Color effect
df.ph <- data.frame(lm.intervals$fixed[-1, ])
df.ph.p <- data.frame(lm.sum$tTable[-1,5])
colnames(df.ph.p) <- "pvalue"
df.ph2 <- cbind(df.ph,df.ph.p)
df.ph2$Var <- "Acidification (pH)"
df.ph2$Eff <- rownames(df.ph)
df.ph2$plot <- "Acidification (variables)"
names(df.ph2)

# Unite the two data.frames 
df.ph.all <- rbind(dfx.ph2, df.ph2)

# FDR on all p-values
df.ph.all$p.correct <- p.adjust(df.ph.all$pvalue, method = "fdr")

#Change Eff names
rownames(df.ph.all) = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype")
df.ph.all$Eff = rownames(df.ph.all)
df.ph.all$Eff <- factor(df.ph.all$Eff, levels = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype"))

# change variables 
ave_combined_ph$type = as.character(ave_combined_ph$type)
ave_combined_ph$type[ave_combined_ph$type=="mono"] <- "Mono-culture"
ave_combined_ph$type[ave_combined_ph$type=="co"] <- "Co-culture"
ave_combined_ph$type = as.factor(ave_combined_ph$type)
str(ave_combined_ph)

#Figure
P2D= df.ph.all[df.ph.all$plot != "Acidification (variables)",] %>% ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_point(data = ave_combined_ph, aes(x = type, y = mean_logfc_ph), 
             alpha = 0.7, position = position_jitter(width = 0.1), color= "#798E87")+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0)+
  theme_bw(base_size = 8)+
  facet_grid(.~plot)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.6, alpha = 1)+
  labs(x="", y = "\n\n Log2(foldchange of pH)\n")+
  scale_y_continuous(limits = c(-0.1,0.1)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5),
        strip.text = element_text(color = "Black", face = "bold", size = 9))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = "white")) +
  annotate(geom = "text", label = "Padj = 0.393", x = 1, y = 0.07, family="sans", size = 2, fontface = 2)+
  annotate(geom = "text", label = "Padj = 0.147", x = 2, y = 0.07, family="sans", size = 2, fontface = 2)
 P2D

 #Plot Biofilm w other variables (FigS4C)
noave_2 = df.ph.all[df.ph.all$plot != "Acidification (pH)",] %>% ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0)+
  theme_bw(base_size = 8)+
  facet_grid(.~plot)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.6, alpha = 1)+
  labs(x="", y = "\n\nLog2(foldchange of pH)\n")+
  scale_y_continuous(limits = c(-0.1,0.1)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5),
        strip.text = element_text(color = "Black", face = "bold", size = 9))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white")) +
  #coord_flip()+
  annotate(geom = "text", label = "Padj = 0.147", x = 1, y = 0.095, family="sans", size = 2, fontface = 2) +
  annotate(geom = "text", label = "Padj= 0.882", x = 2, y = 0.046, family="sans", size = 2, fontface = 2) +
  annotate(geom = "text", label = "Padj= 0.294", x = 3, y = -0.083, family="sans", size = 2, fontface = 2)
noave_2

#### grid.arrange ####
gridExtra::grid.arrange(ave_1, ave_2, noave_1, noave_2, nrow = 2)

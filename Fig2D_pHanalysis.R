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

#### Average bioreps procedure ####
av_bioreps = combined_ph %>%
  group_by(isolate, type, color, lineage, ID) %>%
  dplyr::summarise(mean_logfc_ph = mean(logfc_ph),
                   sd_logfc_ph = sd(logfc_ph),
                   counts = n(),
                   mean_ph = mean(ph))

# T test of ph raw values
wilcox.test(av_bioreps[av_bioreps$type == "mono",]$mean_ph, av_bioreps[av_bioreps$type == "co",]$mean_ph)
#wilcox.test(av_bioreps[av_bioreps$Col.bin == "1",]$mean_cv, av_bioreps[av_bioreps$Col.bin == "0",]$mean_cv)

# Plot check
ggplot(av_bioreps, aes(x = color, y = mean_ph,  fill = color)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8))+
  facet_grid(type~.) +
  theme_bw()

# remove lal
av_bioreps = av_bioreps[av_bioreps$isolate != "lal",]

# make binary for culture type and color
av_bioreps$Col.bin <- ifelse(av_bioreps$color == "l",1,0)
av_bioreps$Cult.bin <- ifelse(av_bioreps$type == "co",1,0) 

lm <- lme(mean_logfc_ph ~ Cult.bin*Col.bin, random = ~ 1|lineage, data=av_bioreps)
#lm2 <- lm(mean_logfc_ph ~ Cult.bin*Col.bin, data=av_bioreps)
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
dfx.ph2$Var <- "ph"
dfx.ph2$Eff <- c("Mono-culture", "Co-culture")
dfx.ph2$plot <- "Plot1"

# Make a summary dataframe for plots with Culture and Color effect
df.ph <- data.frame(lm.intervals$fixed[-1, ])
df.ph.p <- data.frame(lm.sum$tTable[-1,5])
colnames(df.ph.p) <- "pvalue"
df.ph2 <- cbind(df.ph,df.ph.p)
df.ph2$Var <- "ph"
df.ph2$Eff <- rownames(df.ph)
df.ph2$plot <- "Plot2"
names(df.ph2)

# Dataframe for lineage 
lm.fit.ph.null <- lmer(mean_logfc_ph ~ Cult.bin*Col.bin + (1|lineage), data=av_bioreps, REML=FALSE)
lm.fit.ph.model <- lm(mean_logfc_ph ~ Cult.bin*Col.bin , data=av_bioreps)
linage <- anova(lm.fit.ph.null, lm.fit.ph.model)
linage2 <- as.data.frame(rbind(c(NA,NA,NA,linage$`Pr(>Chisq)`[2])))
colnames(linage2) <- c("lower",  "est.",   "upper",  "pvalue")
linage2$Var <- "ph"
linage2$Eff <- "Linage"
linage2$plot <- "Plot2"

# Unite the two data.frames 
df.ph.all <- rbind(dfx.ph2, df.ph2,linage2)

# FDR on all p-values
df.ph.all$p.correct <- p.adjust(df.ph.all$pvalue, method = "fdr")

#Change Eff names
rownames(df.ph.all) = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype", "Lineage")
df.ph.all$Eff = rownames(df.ph.all)
df.ph.all$Eff <- factor(df.ph.all$Eff, levels = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype", "Lineage"))
df.ph.all2 = df.ph.all[df.ph.all$Eff != "Lineage",]


#Figure
ave_1 = df.ph.all2[df.ph.all2$plot != "Plot2",] %>% ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.3)+
  theme_bw(base_size = 8)+
  facet_grid(plot~.)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "Black", size = 0.6, alpha = 1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of pH)", 
       title = "average of bioreps (N = 12)\n")+
  scale_y_continuous(limits = c(-0.1,0.1)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate(geom = "text", label = "Padj = 0.378", x = 1, y = 0.1, family="sans", size = 2, fontface = 2)+
  annotate(geom = "text", label = "Padj = 0.117", x = 2, y = 0.1, family="sans", size = 2, fontface = 2)


#Plot Biofilm w other variables (FigS4C)

ave_2 = df.ph.all2[df.ph.all2$plot != "Plot1",] %>% ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.3)+
  theme_bw(base_size = 8)+
  facet_grid(plot~.)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "Black", size = 0.6, alpha = 1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of pH)", 
       title = "average of bioreps (N = 12)\n")+
  scale_y_continuous(limits = c(-.15,.15)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate(geom = "text", label = "Padj= 0.117", x = 1, y = 0.15, family="sans", size = 2, fontface = 2) +
  annotate(geom = "text", label = "Padj= 0.883", x = 2, y = 0.15, family="sans", size = 2, fontface = 2) +
  annotate(geom = "text", label = "Padj= 0.265", x = 3, y = 0.15, family="sans", size = 2, fontface = 2)


#### No average bioreps ####

# remove lal from raw data that has not been average across bioreps
combined_ph = combined_ph[combined_ph$isolate != "lal",]

# make binary for culture type and color
combined_ph$Col.bin <- ifelse(combined_ph$color == "l",1,0)
combined_ph$Cult.bin <- ifelse(combined_ph$type == "co",1,0) 

lm <- lme(logfc_ph ~ Cult.bin*Col.bin, random = ~ 1|lineage, data=combined_ph)
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
dfx.ph2$Var <- "ph"
dfx.ph2$Eff <- c("Mono-culture", "Co-culture")
dfx.ph2$plot <- "Plot1"

# Make a summary dataframe for plots with Culture and Color effect
df.ph <- data.frame(lm.intervals$fixed[-1, ])
df.ph.p <- data.frame(lm.sum$tTable[-1,5])
colnames(df.ph.p) <- "pvalue"
df.ph2 <- cbind(df.ph,df.ph.p)
df.ph2$Var <- "ph"
df.ph2$Eff <- rownames(df.ph)
df.ph2$plot <- "Plot2"
names(df.ph2)

# Dataframe for lineage 
lm.fit.ph.null <- lmer(logfc_ph ~ Cult.bin*Col.bin + (1|lineage), data=combined_ph, REML=FALSE)
lm.fit.ph.model <- lm(logfc_ph ~ Cult.bin*Col.bin , data=combined_ph)
linage <- anova(lm.fit.ph.null, lm.fit.ph.model)
linage2 <- as.data.frame(rbind(c(NA,NA,NA,linage$`Pr(>Chisq)`[2])))
colnames(linage2) <- c("lower",  "est.",   "upper",  "pvalue")
linage2$Var <- "ph"
linage2$Eff <- "Linage"
linage2$plot <- "Plot2"

# Unite the two data.frames 
df.ph.all <- rbind(dfx.ph2, df.ph2,linage2)

# FDR on all p-values
df.ph.all$p.correct <- p.adjust(df.ph.all$pvalue, method = "fdr")

#Change Eff names
rownames(df.ph.all) = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype", "Lineage")
df.ph.all$Eff = rownames(df.ph.all)
df.ph.all$Eff <- factor(df.ph.all$Eff, levels = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype", "Lineage"))
df.ph.all2 = df.ph.all[df.ph.all$Eff != "Lineage",]


#Figure
 noave_1 = df.ph.all2[df.ph.all2$plot != "Plot2",] %>% ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.3)+
  theme_bw(base_size = 8)+
  facet_grid(plot~.)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "Black", size = 0.6, alpha = 1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of pH)", 
       title = "NO average of bioreps (n = 36)\n")+
  scale_y_continuous(limits = c(-0.1,0.1)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate(geom = "text", label = "Padj = 0.240", x = 1, y = 0.1, family="sans", size = 2, fontface = 2)+
  annotate(geom = "text", label = "Padj = 0.041", x = 2, y = 0.1, family="sans", size = 2, fontface = 2)

#Plot Biofilm w other variables (FigS4C)

noave_2 = df.ph.all2[df.ph.all2$plot != "Plot1",] %>% ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.3)+
  theme_bw(base_size = 8)+
  facet_grid(plot~.)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "Black", size = 0.6, alpha = 1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of pH)", 
       title = "NO average of bioreps (N = 36)\n")+
  scale_y_continuous(limits = c(-.15,.15)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate(geom = "text", label = "Padj= 0.005", x = 1, y = 0.15, family="sans", size = 2, fontface = 2) +
  annotate(geom = "text", label = "Padj= 0.895", x = 2, y = 0.15, family="sans", size = 2, fontface = 2) +
  annotate(geom = "text", label = "Padj= 0.056", x = 3, y = 0.15, family="sans", size = 2, fontface = 2)

#### grid.arrange ####
gridExtra::grid.arrange(ave_1, ave_2, noave_1, noave_2, nrow = 2)
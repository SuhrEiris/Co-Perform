### CV analysis - day14 variants ###


#load data
cv_raw <- as.data.frame(read_excel("Data/CV.xlsx", sheet = 2))
str(cv_raw)

###

#prepare data
cv_raw$type = as.factor(cv_raw$type)
cv_raw$biorep = as.factor(cv_raw$biorep)
cv_raw$color = as.factor(cv_raw$color)
cv_raw$lineage = as.factor(cv_raw$lineage)
str(cv_raw)

# remove medium, unknown and lem 

cv_clean = cv_raw[cv_raw$type !="unknown" & cv_raw$isolate != "medium" & cv_raw$isolate != "lem",]

#plot

ggplot(cv_clean, aes(x = type, y = bgcorrect_cv_adjposvalue, color = type)) +
  geom_boxplot()

# Testing (wilcoxon test)

with(cv_clean, shapiro.test(bgcorrect_cv_adjposvalue[type == "co"]))
hist(cv_clean[cv_clean$type == "co",]$bgcorrect_cv_adjposvalue)

with(cv_clean, shapiro.test(bgcorrect_cv_adjposvalue[type == "wt"]))
hist(cv_clean[cv_clean$type == "wt",]$bgcorrect_cv_adjposvalue)

with(cv_clean, shapiro.test(bgcorrect_cv_adjposvalue[type == "mono"]))
hist(cv_clean[cv_clean$type == "mono",]$bgcorrect_cv_adjposvalue)

wilcox.test(cv_clean[cv_clean$type == "co",]$bgcorrect_cv_adjposvalue, cv_clean[cv_clean$type== "wt",]$bgcorrect_cv_adjposvalue)
wilcox.test(cv_clean[cv_clean$type == "mono",]$bgcorrect_cv_adjposvalue, cv_clean[cv_clean$type== "wt",]$bgcorrect_cv_adjposvalue)
wilcox.test(cv_clean[cv_clean$type == "mono",]$bgcorrect_cv_adjposvalue, cv_clean[cv_clean$type== "co",]$bgcorrect_cv_adjposvalue)


#divide into biorep 
biorep1 = cv_clean[cv_clean$biorep == "1",]
biorep2 = cv_clean[cv_clean$biorep == "2",]
biorep3 = cv_clean[cv_clean$biorep == "3",]

# pick out lal wt 
lal_b1 = biorep1[biorep1$isolate == "lal",]
lal_b2 = biorep2[biorep2$isolate == "lal",]
lal_b3 = biorep3[biorep3$isolate == "lal",]

# # log2 foldchange (using the OD per CV value without background correction n_cv/n_cv_wildtype)
# biorep1$logfc_cv = log2(biorep1$n_cv/lal_b1$n_cv)
# biorep2$logfc_cv = log2(biorep2$n_cv/lal_b2$n_cv)
# biorep3$logfc_cv = log2(biorep3$n_cv/lal_b3$n_cv)

#Background corrected logfc (all bgcorrected CV values that are negative are set to zero) a pseudocount of one is added 
# biorep1$logfc_cv = log2((biorep1$n_bgcorrect_cv_adjneg+1)/(lal_b1$n_bgcorrect_cv_adjneg+1))
# biorep2$logfc_cv = log2((biorep2$n_bgcorrect_cv_adjneg+1)/(lal_b2$n_bgcorrect_cv_adjneg+1))
# biorep3$logfc_cv = log2((biorep3$n_bgcorrect_cv_adjneg+1)/(lal_b3$n_bgcorrect_cv_adjneg+1))

# Logfold change with only CV values (That normalized to OD), negative values are set to the lowest CV value measured in the dataset
biorep1$logfc_cv = log2(biorep1$bgcorrect_cv_adjposvalue/lal_b1$bgcorrect_cv_adjposvalue)
biorep2$logfc_cv = log2(biorep2$bgcorrect_cv_adjposvalue/lal_b2$bgcorrect_cv_adjposvalue)
biorep3$logfc_cv = log2(biorep3$bgcorrect_cv_adjposvalue/lal_b3$bgcorrect_cv_adjposvalue)

# merge bioreps
combined_cv = rbind(biorep1, biorep2, biorep3)
combined_cv = combined_cv %>%  unite("ID1", isolate:type , na.rm = TRUE, remove = FALSE)
combined_cv = combined_cv %>%  unite("ID", ID1,color , na.rm = TRUE, remove = FALSE)
combined_cv = na.omit(combined_cv)

#Logfold change plot 
ggplot(combined_cv, aes(x = ID, y = logfc_cv,  fill = type)) +
  geom_bar(stat = "identity")+
  facet_grid(biorep~.) + 
  theme_bw()

# T test of cv raw values
wilcox.test(combined_cv[combined_cv$type == "mono",]$bgcorrect_cv_adjposvalue, combined_cv[combined_cv$type == "co",]$bgcorrect_cv_adjposvalue)
#wilcox.test(combined_cv[combined_cv$Col.bin == "1",]$bgcorrect_cv_adjposvalue, combined_cv[combined_cv$Col.bin == "1",]$bgcorrect_cv_adjposvalue)

# Plot check
ggplot(combined_cv, aes(x = color, y = bgcorrect_cv_adjposvalue,  fill = color)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8))+
  facet_grid(type~.) +
  theme_bw()

#Average bioreps
ave_combined_cv = combined_cv %>%
  group_by(isolate, type, color, lineage, ID) %>%
  dplyr::summarise(mean_logfc_cv = mean(logfc_cv),
                   sd_logfc_cv = sd(logfc_cv),
                   counts = n(),
                   mean_cv = mean(bgcorrect_cv_adjposvalue))

# remove lal
ave_combined_cv = ave_combined_cv[ave_combined_cv$isolate != "lal",]

# make binary for culture type and color
ave_combined_cv$Col.bin <- ifelse(ave_combined_cv$color == "l",1,0)
ave_combined_cv$Cult.bin <- ifelse(ave_combined_cv$type == "co",1,0)

# Remove outlier
mod <- lm(mean_logfc_cv ~ ., data=as.data.frame(ave_combined_cv[,6]))
cooksd <- cooks.distance(mod)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

# Outlier removal by Car test
car::outlierTest(mod)
# No outliers

set.seed(41)
lm <- lme(mean_logfc_cv ~ Cult.bin*Col.bin, random = ~ 1|lineage, data=ave_combined_cv)
lm2 <- lm(mean_logfc_cv ~ Cult.bin*Col.bin, data=ave_combined_cv)
lm.sum <- summary(lm)
lm.intervals <- intervals(lm,which = "fixed")

#library(multcomp)
contr <- rbind("COCvsInter"=c(1,1,0,0)) 
glht <- multcomp::glht(lm, linfct=contr)
glht.confint <- confint(glht)
glht.sum <- summary(glht)

# Make a summary dataframe for plots with single and co-culture 
dfx.cv <- data.frame(rbind(lm.intervals$fixed[1, ], glht.confint$confint[c(2,1,3)]))
df.cv.psing <- data.frame(lm.sum$tTable[1,5])
colnames(df.cv.psing) <- "pvalue"
df.cv.pcoc <- data.frame(glht.sum$test$pvalues)
colnames(df.cv.pcoc) <- "pvalue"
df.cv.pboth <- rbind(df.cv.psing,df.cv.pcoc)
dfx.cv2 <- cbind(dfx.cv,df.cv.pboth)
dfx.cv2$Var <- "cv"
dfx.cv2$Eff <- c("Mono-culture", "Co-culture")
dfx.cv2$plot <- "Biofilm formation"

# Make a summary dataframe for plots with Culture and Color effect
df.cv <- data.frame(lm.intervals$fixed[-1, ])
df.cv.p <- data.frame(lm.sum$tTable[-1,5])
colnames(df.cv.p) <- "pvalue"
df.cv2 <- cbind(df.cv,df.cv.p)
df.cv2$Var <- "cv"
df.cv2$Eff <- rownames(df.cv)
df.cv2$plot <- "Biofilm formation (variables)"
names(df.cv2)

# Unite the two data.frames 
df.cv.all <- rbind(dfx.cv2, df.cv2)

# FDR on all p-values
df.cv.all$p.correct <- p.adjust(df.cv.all$pvalue, method = "fdr")

#Change Eff names
rownames(df.cv.all) = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype")
df.cv.all$Eff = rownames(df.cv.all)
df.cv.all$Eff <- factor(df.cv.all$Eff, levels = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype"))

#Make jitter ready
ave_combined_cv$type = as.character(ave_combined_cv$type)
ave_combined_cv$type[ave_combined_cv$type=="mono"] <- "Mono-culture"
ave_combined_cv$type[ave_combined_cv$type=="co"] <- "Co-culture"
ave_combined_cv$type = as.factor(ave_combined_cv$type)
str(ave_combined_cv)


#Figure
P2C = df.cv.all[df.cv.all$plot != "Biofilm formation (variables)",] %>% ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_point(data = ave_combined_cv, aes(x = type, y = mean_logfc_cv), 
             alpha = 0.7, position = position_jitter(width = 0.1), color= "#798E87")+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0)+
  theme_bw(base_size = 8)+
  facet_grid(.~plot)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.6, alpha = 1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of biofilm formation)\n")+
  scale_y_continuous(limits = c(-7,7), breaks = seq(-7,7)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5),
        strip.text = element_text(color = "Black", face = "bold", size = 9))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"))+ 
  annotate(geom = "text", label = "Padj < 0.0001", x = 1, y = 7, family="sans", size = 2, fontface = 2) + 
  annotate(geom = "text", label = "Padj < 0.0001", x = 2, y = 7, family="sans", size = 2, fontface = 2) 
P2C

#Plot Biofilm w other variables (FigS4C)

PS4C = df.cv.all[df.cv.all$plot != "Biofilm formation",] %>% ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0)+
  theme_bw(base_size = 8)+
  facet_grid(.~plot)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.6, alpha = 1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of biofilm formation)\n")+
  scale_y_continuous(limits = c(-4,4), breaks = seq(-4,4)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5),
        strip.text = element_text(color = "Black", face = "bold", size = 9))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = "white")) +
  annotate(geom = "text", label = "Padj = 0.002", x = 1, y = 3.5, family="sans", size = 2, fontface = 2) + 
  annotate(geom = "text", label = "Padj < 0.0001", x = 2, y = 3.5, family="sans", size = 2, fontface = 2) +
  annotate(geom = "text", label = "Padj = 0.007", x = 3, y = -4, family="sans", size = 2, fontface = 2)

PS4C

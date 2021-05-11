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

# Average bioreps

av_bioreps = combined_cv %>%
  group_by(isolate, type, color, lineage, ID) %>%
  dplyr::summarise(mean_logfc_cv = mean(logfc_cv),
                   sd_logfc_cv = sd(logfc_cv),
                   counts = n(),
                   mean_cv = mean(bgcorrect_cv_adjposvalue))

# T test of cv raw values
wilcox.test(av_bioreps[av_bioreps$type == "mono",]$mean_cv, av_bioreps[av_bioreps$type == "co",]$mean_cv)

# Plot check
ggplot(av_bioreps, aes(x = type, y = mean_cv,  fill = type)) +
  geom_boxplot()+
  theme_bw()

# remove lal
av_bioreps = av_bioreps[av_bioreps$isolate != "lal",]

# make binary for culture type and color
av_bioreps$Col.bin <- ifelse(av_bioreps$color == "l",1,0)
av_bioreps$Cult.bin <- ifelse(av_bioreps$type == "co",1,0) 

lm <- lme(mean_logfc_cv ~ Cult.bin*Col.bin, random = ~ 1|lineage, data=av_bioreps)
lm2 <- lm(mean_logfc_cv ~ Cult.bin*Col.bin, data=av_bioreps)
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
dfx.cv2$Eff <- c("Mono", "CoCulture")
dfx.cv2$plot <- "Plot1"

# Make a summary dataframe for plots with Culture and Color effect
df.cv <- data.frame(lm.intervals$fixed[-1, ])
df.cv.p <- data.frame(lm.sum$tTable[-1,5])
colnames(df.cv.p) <- "pvalue"
df.cv2 <- cbind(df.cv,df.cv.p)
df.cv2$Var <- "cv"
df.cv2$Eff <- rownames(df.cv)
df.cv2$plot <- "Plot2"
names(df.cv2)

# Dataframe for lineage 
lm.fit.cv.null <- lmer(mean_logfc_cv ~ Cult.bin*Col.bin + (1|lineage), data=av_bioreps, REML=FALSE)
lm.fit.cv.model <- lm(mean_logfc_cv ~ Cult.bin*Col.bin , data=av_bioreps)
linage <- anova(lm.fit.cv.null, lm.fit.cv.model)
linage2 <- as.data.frame(rbind(c(NA,NA,NA,linage$`Pr(>Chisq)`[2])))
colnames(linage2) <- c("lower",  "est.",   "upper",  "pvalue")
linage2$Var <- "cv"
linage2$Eff <- "Linage"
linage2$plot <- "Plot2"

# Unite the two data.frames 
df.cv.all <- rbind(dfx.cv2, df.cv2,linage2)

# FDR on all p-values
df.cv.all$p.correct <- p.adjust(df.cv.all$pvalue, method = "fdr")

#Order the variables for plotting
df.cv.all$Eff <- factor(df.cv.all$Eff, levels = c("Mono", "CoCulture", "Cult.bin","Col.bin", "Cult.bin:Col.bin", "Linage"))

#Figure
df.cv.all[df.cv.all$Eff != "Linage",] %>% ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.3)+
  theme_bw()+
  geom_hline(yintercept = 0, color = "dodgerblue4", size = 0.8)+
  labs(title = "CV analysis")+
  facet_grid(~plot, scale= "free_x", space = "free_x")



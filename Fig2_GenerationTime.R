
#Load dataframe with culture yield values
data.omit  <- readRDS("Data/df.phenotypicdata.rds")

############  Data from generation time ######
lm.fit.t.gen <- lme(t.gen ~ Cult.bin*Col.bin, random = ~ 1|Lineage, data=data.omit)
jakob.2 <- summary(lm.fit.t.gen)
jakob <- intervals(lm.fit.t.gen)

library(multcomp)
contr <- rbind("COCvsInter"=c(1,1,0,0)) 
test <- multcomp::glht(lm.fit.t.gen, linfct=contr)
jakob1 <- confint(test)
jakob3 <- summary(test)

# Make a summary dataframe for plots with single and co-culture 
dfx.t.gen <- data.frame(rbind(jakob$fixed[1, ], jakob1$confint[c(2,1,3)]))
df.t.gen.psing <- data.frame(jakob.2$tTable[1,5])
colnames(df.t.gen.psing) <- "pvalue"
df.t.gen.pcoc <- data.frame(jakob3$test$pvalues)
colnames(df.t.gen.pcoc) <- "pvalue"
df.t.gen.pboth <- rbind(df.t.gen.psing,df.t.gen.pcoc)
dfx.t.gen2 <- cbind(dfx.t.gen,df.t.gen.pboth)
dfx.t.gen2$Var <- "t.gen"
dfx.t.gen2$Eff <- c("Single", "CoCulture")
dfx.t.gen2$plot <- "Plot1"

# Make a summary dataframe for plots with Culture and Color effect
df.t.gen <- data.frame(jakob$fixed[-1, ])
df.t.gen.p <- data.frame(jakob.2$tTable[-1,5])
colnames(df.t.gen.p) <- "pvalue"
df.t.gen2 <- cbind(df.t.gen,df.t.gen.p)
df.t.gen2$Var <- "t.gen"
df.t.gen2$Eff <- rownames(df.t.gen)
df.t.gen2$plot <- "Plot2"
names(df.t.gen2)

# Dataframe for linage:
lm.fit.t.gen.null <- lmer(t.gen ~ Cult.bin*Col.bin + (1|Lineage), data=data.select, REML=FALSE)
lm.fit.t.gen.model <- lm(t.gen ~ Culture*Color , data=data.select)
linage <- anova(lm.fit.t.gen.null, lm.fit.t.gen.model)
linage2 <- as.data.frame(rbind(c(NA,NA,NA,linage$`Pr(>Chisq)`[2])))
colnames(linage2) <- c("lower",  "est.",   "upper",  "pvalue")
linage2$Var <- "t.gen"
linage2$Eff <- "Linage"
linage2$plot <- "Plot2"

# Unite the two data.frames 
df.t.gen.all <- rbind(dfx.t.gen2, df.t.gen2,linage2)
df.t.gen.all$p.correct <- p.adjust(df.t.gen.all$pvalue, method = "fdr")

df.t.gen.all$Eff <- factor(df.t.gen.all$Eff, levels = c("Single", "CoCulture", "Cult.bin","Col.bin", "Cult.bin:Col.bin", "Linage"))


p <- df.t.gen.all[df.t.gen.all$Eff != "Linage",] %>% ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.3)+
  theme_bw()+
  geom_hline(yintercept = 0, color = "dodgerblue4", size = 0.8)+
  labs(title = "GenerationTime")+
  scale_y_continuous(limits = c(-0.6,0.5))+
  facet_grid(~plot, scale= "free_x", space = "free_x")
p

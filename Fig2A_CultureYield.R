
#Load dataframe with culture yield values
data.omit  <- readRDS("Data/df.phenotypicdata.rds")

lm.fit.k <- lme(k ~ Cult.bin*Col.bin, random = ~ 1|Lineage, data=data.omit)
jakob.2 <- summary(lm.fit.k)
jakob <- intervals(lm.fit.k)

contr <- rbind("COCvsInter"=c(1,1,0,0)) 
test <- multcomp::glht(lm.fit.k, linfct=contr)
jakob1 <- confint(test)
jakob3 <- summary(test)

# Make a summary dataframe for plots with single and co-culture 
dfx.k <- data.frame(rbind(jakob$fixed[1, ], jakob1$confint[c(2,1,3)]))
df.k.psing <- data.frame(jakob.2$tTable[1,5])
colnames(df.k.psing) <- "pvalue"
df.k.pcoc <- data.frame(jakob3$test$pvalues)
colnames(df.k.pcoc) <- "pvalue"
df.k.pboth <- rbind(df.k.psing,df.k.pcoc)
dfx.k2 <- cbind(dfx.k,df.k.pboth)
dfx.k2$Var <- "k"
dfx.k2$Eff <- c("Mono-culture", "Co-culture")
dfx.k2$plot <- "Plot1"

# Make a summary dataframe for plots with Culture and Color effect
df.k <- data.frame(jakob$fixed[-1, ])
df.k.p <- data.frame(jakob.2$tTable[-1,5])
colnames(df.k.p) <- "pvalue"
df.k2 <- cbind(df.k,df.k.p)
df.k2$Var <- "k"
df.k2$Eff <- rownames(df.k)
df.k2$plot <- "Plot2"
names(df.k2)

# Dataframe for linage:
lm.fit.k.null <- lmer(k ~ Cult.bin*Col.bin + (1|Lineage), data=data.omit, REML=FALSE)
lm.fit.k.model <- lm(k ~ Cult.bin*Col.bin , data=data.omit)
linage <- anova(lm.fit.k.null, lm.fit.k.model)
linage2 <- as.data.frame(rbind(c(NA,NA,NA,linage$`Pr(>Chisq)`[2])))
colnames(linage2) <- c("lower",  "est.",   "upper",  "pvalue")
linage2$Var <- "k"
linage2$Eff <- "Linage"
linage2$plot <- "Plot2"

# Unite the two data.frames 
df.k.all <- rbind(dfx.k2, df.k2,linage2)
df.k.all$p.correct <- p.adjust(df.k.all$pvalue, method = "fdr")

#Change Eff names
rownames(df.k.all) = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype", "Lineage")
df.k.all$Eff = rownames(df.k.all)
df.k.all$Eff <- factor(df.k.all$Eff, levels = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype", "Lineage"))
df.k.all2 = df.k.all[df.k.all$Eff != "Lineage",]

P2A <- df.k.all2[df.k.all2$plot != "Plot2",] %>% 
  ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.2)+
  theme_bw(base_size = 8)+
  facet_grid(plot~.)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "black", size = 0.6, alpha =1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of culture yield)", 
       title = "Culture yield\n")+
  scale_y_continuous(limits = c(-0.4,0.6)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+ 
  annotate(geom = "text", label = "Padj= 0.676", x = 1, y = 0.3, family="sans", size = 2, fontface = 2) + 
  annotate(geom = "text", label = "Padj=0.008", x = 2, y = 0.1, family="sans", size = 2, fontface = 2)
P2A

# Plot culture yield other variables (figS4A)
PS4A <- df.k.all2[df.k.all2$plot != "Plot1",] %>% 
  ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.2)+
  theme_bw(base_size = 8)+
  facet_grid(plot~.)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "black", size = 0.6, alpha =1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of culture yield)", 
       title = "Culture yield\n")+
  scale_y_continuous(limits = c(-0.6,0.6)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+ 
  annotate(geom = "text", label = "Padj= 0.013", x = 1, y = 0.6, family="sans", size = 2, fontface = 2) + 
  annotate(geom = "text", label = "Padj=0.676", x = 2, y = 0.22, family="sans", size = 2, fontface = 2) +
  annotate(geom = "text", label = "Padj=0.676", x = 3, y = 0.3, family="sans", size = 2, fontface = 2)
PS4A

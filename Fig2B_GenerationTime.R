
#Load dataframe with culture yield values
data.omit  <- readRDS("Data/df.phenotypicdata.rds")
data.omit<-data.omit[!(data.omit$Sample=="XB.3-m"),]

############  Data from generation time ######
set.seed(41)
lm.fit.t.gen <- lme(t.gen ~ Cult.bin*Col.bin, random = ~ 1|Lineage, data=data.omit)
jakob.2 <- summary(lm.fit.t.gen)
jakob <- intervals(lm.fit.t.gen)

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
dfx.t.gen2$Eff <- c("Mono-culture", "Co-culture")
dfx.t.gen2$plot <- "Generation time"

# Make a summary dataframe for plots with Culture and Color effect
df.t.gen <- data.frame(jakob$fixed[-1, ])
df.t.gen.p <- data.frame(jakob.2$tTable[-1,5])
colnames(df.t.gen.p) <- "pvalue"
df.t.gen2 <- cbind(df.t.gen,df.t.gen.p)
df.t.gen2$Var <- "t.gen"
df.t.gen2$Eff <- rownames(df.t.gen)
df.t.gen2$plot <- "Generation time (variables)"
names(df.t.gen2)

# Dataframe for linage:
lm.fit.t.gen.null <- lmer(t.gen ~ Cult.bin*Col.bin + (1|Lineage), data=data.omit, REML=FALSE)
lm.fit.t.gen.model <- lm(t.gen ~ Culture*Color , data=data.omit)
linage <- anova(lm.fit.t.gen.null, lm.fit.t.gen.model)
linage2 <- as.data.frame(rbind(c(NA,NA,NA,linage$`Pr(>Chisq)`[2])))
colnames(linage2) <- c("lower",  "est.",   "upper",  "pvalue")
linage2$Var <- "t.gen"
linage2$Eff <- "Linage"
linage2$plot <- "Generation time (variables"

# Unite the two data.frames 
df.t.gen.all <- rbind(dfx.t.gen2, df.t.gen2,linage2)
df.t.gen.all$p.correct <- p.adjust(df.t.gen.all$pvalue, method = "fdr")

#Change names of Eff
rownames(df.t.gen.all) = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype", "Lineage")
df.t.gen.all$Eff = rownames(df.t.gen.all)
df.t.gen.all$Eff <- factor(df.t.gen.all$Eff, levels = c("Mono-culture", "Co-culture","Culture","Morphotype", "Culture:Morphotype", "Lineage"))
df.t.gen.all2 = df.t.gen.all[df.t.gen.all$Eff != "Lineage",]

data.omit$Culture = as.character(data.omit$Culture)
data.omit$Culture[data.omit$Culture=="Single"] <- "Mono-culture"
data.omit$Culture[data.omit$Culture=="Coc"] <- "Co-culture"
data.omit$Culture = as.factor(data.omit$Culture)
str(data.omit)
data.omit$Culture <- factor(data.omit$Culture, 
                            levels = c("Mono-culture", "Co-culture"))


#Generation time plot with main variables (Fig2B)
P2B <- df.t.gen.all2[df.t.gen.all2$plot != "Generation time (variables)",] %>% 
  ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_point(data = data.omit, aes(x = Culture, y = t.gen),
             alpha = 0.7, position = position_jitter(width = 0.1), color= "#798E87")+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0)+
  theme_bw(base_size = 8)+
  facet_grid(.~plot)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.6, alpha =1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of generation time)\n")+
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.6, alpha =1)+
  scale_y_continuous(limits = c(-1,1)) +
  theme_bw(base_size = 8)+
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5),
        strip.text = element_text(color = "Black", face = "bold", size = 9))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"))+
  annotate(geom = "text", label = "Padj = 0.3", x = 1, y = 0.7, family="sans", size = 2, fontface = 2) + 
  annotate(geom = "text", label = "Padj = 0.0004", x = 2, y = 0.7, family="sans", size = 2, fontface = 2)
P2B

#Generation time plot with other variables (FigS4B)
PS4B <- df.t.gen.all2[df.t.gen.all2$plot != "Generation time",] %>% 
  ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
  theme_bw()+
  facet_grid(.~plot)+
  labs(x="", 
       y = "\n\nLog2(foldchange of generation time)\n")+
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.6, alpha =1)+
  scale_y_continuous(limits = c(-0.6,0.6)) +
  theme_bw(base_size = 8)+
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5),
        strip.text = element_text(color = "Black", face = "bold", size = 9))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = "white")) +
  annotate(geom = "text", label = "Padj < 0.0001", x = 1, y = 0.5, family="sans", size = 2, fontface = 2) + 
  annotate(geom = "text", label = "Padj = 0.0014", x = 2, y = 0.5, family="sans", size = 2, fontface = 2) +
  annotate(geom = "text", label = "Padj < 0.0001", x = 3, y = -0.6, family="sans", size = 2, fontface = 2)
PS4B

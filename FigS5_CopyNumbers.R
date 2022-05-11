
#Load dataframe with culture yield values
data.omit  <- readRDS("Data/df.phenotypicdata.rds")

#remove XB.3-m due to contamination
data.omit<-data.omit[!(data.omit$Sample=="XB.3-m"),]

lm.fit.copy <- lme(copynumb ~ Cult.bin*Col.bin, random = ~ 1|Lineage, data=data.omit)
lm.fit.copy2 <- lm(copynumb ~ Cult.bin*Col.bin, data=data.omit)
jakob.2 <- summary(lm.fit.copy)
jakob <- intervals(lm.fit.copy,which = "fixed")

#library(multcomp)
contr <- rbind("COCvsInter"=c(1,1,0,0)) 
test <- multcomp::glht(lm.fit.copy, linfct=contr)
jakob1 <- confint(test)
jakob3 <- summary(test)

# Make a summary dataframe for plots with single and co-culture 
dfx.copy <- data.frame(rbind(jakob$fixed[1, ], jakob1$confint[c(2,1,3)]))
df.copy.psing <- data.frame(jakob.2$tTable[1,5])
colnames(df.copy.psing) <- "pvalue"
df.copy.pcoc <- data.frame(jakob3$test$pvalues)
colnames(df.copy.pcoc) <- "pvalue"
df.copy.pboth <- rbind(df.copy.psing,df.copy.pcoc)
dfx.copy2 <- cbind(dfx.copy,df.copy.pboth)
dfx.copy2$Var <- "copy"
dfx.copy2$Eff <- c("Mono-culture", "Co-culture")
dfx.copy2$plot <- "Plot1"

# Make a summary dataframe for plots with Culture and Color effect
df.copy <- data.frame(jakob$fixed[-1, ])
df.copy.p <- data.frame(jakob.2$tTable[-1,5])
colnames(df.copy.p) <- "pvalue"
df.copy2 <- cbind(df.copy,df.copy.p)
df.copy2$Var <- "copy"
df.copy2$Eff <- rownames(df.copy)
df.copy2$plot <- "Plot2"
names(df.copy2)

# Dataframe for linage:
lm.fit.copy.null <- lmer(copynumb ~ Cult.bin*Col.bin + (1|Lineage), data=data.omit, REML=FALSE)
lm.fit.copy.model <- lm(copynumb ~ Cult.bin*Col.bin , data=data.omit)
linage <- anova(lm.fit.copy.null, lm.fit.copy.model)
linage2 <- as.data.frame(rbind(c(NA,NA,NA,linage$`Pr(>Chisq)`[2])))
colnames(linage2) <- c("lower",  "est.",   "upper",  "pvalue")
linage2$Var <- "copy"
linage2$Eff <- "Linage"
linage2$plot <- "Plot2"

# Unite the two data.frames 
df.copy.all <- rbind(dfx.copy2, df.copy2,linage2)
df.copy.all$p.correct <- p.adjust(df.copy.all$pvalue, method = "fdr")

df.copy.all$Eff <- factor(df.copy.all$Eff, levels = c("Single", "CoCulture", "Cult.bin","Col.bin", "Cult.bin:Col.bin", "Linage"))

# 
# p <- df.copy.all[df.copy.all$Eff != "Linage",] %>% ggplot(aes(x= Eff, y =est.))+
#   geom_point(size = 2)+
#   geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.3)+
#   theme_bw()+
#   geom_hline(yintercept = 0, color = "dodgerblue4")+
#   labs(title = "copy")+
#   #scale_y_continuous(limits = c(-0.55,1))+
#   facet_grid(~plot, scale= "free_x", space = "free_x")
# p

#Figure S5

data.omit$Culture = as.character(data.omit$Culture)
data.omit$Culture[data.omit$Culture=="Single"] <- "Mono-culture"
data.omit$Culture[data.omit$Culture=="Coc"] <- "Co-culture"
data.omit$Culture = as.factor(data.omit$Culture)
str(data.omit)
data.omit$Culture <- factor(data.omit$Culture, 
                            levels = c("Mono-culture", "Co-culture"))

dfx.copy2$Eff <- factor(dfx.copy2$Eff, levels = c("Mono-culture", "Co-culture"))

PS5 = dfx.copy2 %>% 
  ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_point(data = data.omit, aes(x = Culture, y =copynumb ),
             alpha = 0.7, position = position_jitter(width = 0.1), color= "#798E87")+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0)+
  theme_bw(base_size = 8)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.6, alpha = 1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of copynumbers)\n")+
  scale_y_continuous(limits = c(-2,2)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5),
        strip.text = element_text(color = "Black", face = "bold", size = 9))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"))+ 
  annotate(geom = "text", label = "Padj = 0.193", x = 1, y = 1.5, family="sans", size = 2, fontface = 2) + 
  annotate(geom = "text", label = "Padj < 0.0001", x = 2, y = 1.5, family="sans", size = 2, fontface = 2)
PS5

ggsave("FigS5.pdf",
       device = "pdf",
       plot = PS5,
       units="mm",
       width=88,
       height=100,
       dpi = 300)

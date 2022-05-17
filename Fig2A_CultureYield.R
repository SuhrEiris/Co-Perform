
#Load dataframe with culture yield values
data.omit  <- readRDS("Data/df.phenotypicdata.rds")

# Outlier removal
test2 = data.omit[,c(1,5)]
mod <- lm(k ~ ., data=test2)
cooksd <- cooks.distance(mod)

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

# Outlier removal by Car test
car::outlierTest(mod)

#Remove outlier
data.omit2 <- data.omit[-c(16,22,40,70),] # large cook distance

#MELM
set.seed(41)
lm.fit.k <- lme(k ~ Cult.bin*Col.bin, random = ~ 1|Lineage, data=data.omit2)
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
dfx.k2$plot <- "Culture yield"

# Make a summary dataframe for plots with Culture and Color effect
df.k <- data.frame(jakob$fixed[-1, ])
df.k.p <- data.frame(jakob.2$tTable[-1,5])
colnames(df.k.p) <- "pvalue"
df.k2 <- cbind(df.k,df.k.p)
df.k2$Var <- "k"
df.k2$Eff <- rownames(df.k)
df.k2$plot <- "Culture yield (variables)"
names(df.k2)

# Unite the two data.frames 
df.k.all <- rbind(dfx.k2, df.k2)
df.k.all$p.correct <- p.adjust(df.k.all$pvalue, method = "fdr")

#Change Eff names
rownames(df.k.all) = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype")
df.k.all$Eff = rownames(df.k.all)
df.k.all$Eff <- factor(df.k.all$Eff, 
                       levels = c("Mono-culture", "Co-culture", "Culture", "Morphotype", "Culture:Morphotype"))

data.omit2$Culture = as.character(data.omit2$Culture)
data.omit2$Culture[data.omit2$Culture=="Single"] <- "Mono-culture"
data.omit2$Culture[data.omit2$Culture=="Coc"] <- "Co-culture"
data.omit2$Culture = as.factor(data.omit2$Culture)
str(data.omit2)
data.omit2$Culture <- factor(data.omit2$Culture, 
                            levels = c("Mono-culture", "Co-culture"))

#Plot

P2A <- df.k.all[df.k.all$plot != "Culture yield (variables)",] %>% 
  ggplot(aes(x= Eff, y =est.))+
  geom_point(data = data.omit2, aes(x = Culture, y = k),
             alpha = 0.7, position = position_jitter(width = 0.1), color= "#798E87")+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0)+
  theme_bw(base_size = 8)+
  facet_grid(.~plot)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.6, alpha =1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of culture yield)\n")+
  scale_y_continuous(limits = c(-1.3,1.3)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5),
        strip.text = element_text(color = "Black", face = "bold", size = 9))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"))+ 
  annotate(geom = "text", label = "Padj = 0.175", x = 1, y = 1, family="sans", size = 2, fontface = 2) + 
  annotate(geom = "text", label = "Padj < 0.0001", x = 2, y = 1, family="sans", size = 2, fontface = 2)
P2A

# Plot culture yield other variables (figS4A)
PS4A <- df.k.all[df.k.all$plot != "Culture yield",] %>% 
  ggplot(aes(x= Eff, y =est.))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0)+
  theme_bw(base_size = 8)+
  facet_grid(.~plot)+
  theme(strip.text.y = element_text(color="black", face="bold"))+
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.6, alpha =1)+
  labs(x="", 
       y = "\n\nLog2(foldchange of culture yield)\n")+
  scale_y_continuous(limits = c(-0.5,0.5)) +
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"),
        plot.title = element_text(color = "Black", face ="bold", hjust = 0.5),
        strip.text = element_text(color = "Black", face = "bold", size = 9))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = "white")) +
  annotate(geom = "text", label = "Padj = 0.007", x = 1, y = 0.45, family="sans", size = 2, fontface = 2) + 
  annotate(geom = "text", label = "Padj = 0.168", x = 2, y = 0.1, family="sans", size = 2, fontface = 2) +
  annotate(geom = "text", label = "Padj = 0.928", x = 3, y = 0.3, family="sans", size = 2, fontface = 2)
PS4A

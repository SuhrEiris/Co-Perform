
# Load data
cfu <- read_excel("Data/df.beadexperimentCFU.xlsx")

#Average technical replicates
cfu_data.1 <- cfu %>% 
  group_by(Day, Species, Replicate, Count_of) %>% 
  dplyr::summarise(Count_tech= n(), 
            Mean_cfu = mean(CFU_mL))

#df with xl and xs to be merged 
cfu_data.3 <-cfu_data.1[cfu_data.1$Count_of == "Xl",]
cfu_data.4 <-cfu_data.1[cfu_data.1$Count_of == "Xs",]
cfu_data.5 <- rbind(cfu_data.3,cfu_data.4)

#adding xl and xs 
cfu_data.6 <- cfu_data.5 %>% 
  group_by(Day, Species, Replicate) %>% 
  dplyr::summarise(Count_tech=n(), 
            Mean_cfu=sum(Mean_cfu))

cfu_data.6$Count_of <- "X"

#dataset without xl and xs
cfu_data.7 <- cfu_data.1[ which(cfu_data.1$Count_of !='Xs' & cfu_data.1$Count_of !='Xl'), ]

#rbind count of x with rest
vec.names <- colnames(cfu_data.7)
cfu_data.6 <- cfu_data.6[,vec.names]
cfu_data.8 <- rbind(cfu_data.6,cfu_data.7)

# rep A at day 13 and rep C at day 7 is zero (plating problem) taken out
cfu_data.8 = cfu_data.8[-38,]
cfu_data.8 = cfu_data.8[-65,]

log_data = cfu_data.8
log_data$Mean_cfu = log10(log_data$Mean_cfu)

#Average biological replicates
cfu_data.9 <- log_data %>% 
  group_by(Day, Species, Count_of) %>% 
  dplyr::summarize(Count_bio= n(), 
            CFU = mean(Mean_cfu), 
            Std.dev_bio = sd(Mean_cfu))

pres <- cfu_data.9[!grepl("Z", cfu_data.9$Count_of),]
pres2 <- cfu_data.9[!grepl("Z", cfu_data.9$Species),]

pres2$Species[pres2$Species=="X"] <- "Mono-culture (L. lactis)"
pres2$Species[pres2$Species=="Y"] <- "Mono-culture (L. mesenteroides)"
pres2$Species[pres2$Species=="XY"] <- "Co-culture"
pres2$Count_of[pres2$Count_of=="X"] <- "L. lactis"
pres2$Count_of[pres2$Count_of=="Y"] <- "L. mesenteroides"

# SE for CI intervals

pres2$se = pres2$Std.dev_bio/sqrt(pres2$Count_bio)
pres2$CI = qnorm(0.975)*pres2$Std.dev_bio/sqrt(pres2$Count_bio)

#Labels for factgrid

pres.plot <- ggplot(data=pres2, aes(x=Day, y= CFU, fill = Count_of, color = Count_of))+
  geom_point(size=1)+
  geom_line(size=0.5)+
  facet_grid(~Species)+
  geom_ribbon(aes(ymax = CFU + CI, ymin = CFU - CI, alpha = 0.1, size = NA)) +
  labs(x="\nTime (days)",y="CFU/mL\n") + 
  scale_x_continuous(breaks = seq(1,18,by=3))+
  scale_y_continuous(breaks = seq(4,8, by =1)) +
  scale_color_manual(values= wes_palette("Moonrise2"))+
  scale_fill_manual(values = wes_palette("Moonrise2"))+
  scale_alpha(guide=FALSE)+
  theme_bw(base_size = 8)+
  theme(legend.position = c(0.91, 0.15),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold.italic"),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"))+
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold")) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
pres.plot

ggsave("Fig1B.pdf",
       plot = pres.plot,
       units="cm",
       width=18,
       height=10,
       dpi = 300)

# T.test different days 

# Make dataset nice and tidy
CFU.ttest = log_data[!grepl("Z", log_data$Count_of),]
CFU.ttest=  log_data[!grepl("Z", log_data$Species),]
CFU.ttest$Species[CFU.ttest$Species=="X"] <- "LaL (Single)"
CFU.ttest$Species[CFU.ttest$Species=="Y"] <- "LeM (Single)"
CFU.ttest$Species[CFU.ttest$Species=="XY"] <- "Co-culture"
CFU.ttest$Count_of[CFU.ttest$Count_of=="X"] <- "LaL"
CFU.ttest$Count_of[CFU.ttest$Count_of=="Y"] <- "LeM"

CFU.ttest.co = CFU.ttest[grepl("Co-culture", CFU.ttest$Species),]

# t.test of co-culture populations
d1.co = CFU.ttest[CFU.ttest$Day == "1",]
t.test(d1.co$Mean_cfu ~ d1.co$Count_of)

d4.co = CFU.ttest[CFU.ttest$Day == "4",]
t.test(d4.co$Mean_cfu ~ d4.co$Count_of)

d7.co = CFU.ttest[CFU.ttest$Day == "7",]
t.test(d7.co$Mean_cfu ~ d7.co$Count_of)

d10.co = CFU.ttest[CFU.ttest$Day == "10",]
t.test(d10.co$Mean_cfu ~ d10.co$Count_of)

d13.co = CFU.ttest[CFU.ttest$Day == "13",]
t.test(d13.co$Mean_cfu ~ d13.co$Count_of)

d16.co = CFU.ttest[CFU.ttest$Day == "16",]
t.test(d16.co$Mean_cfu ~ d16.co$Count_of)


# t.test of lal sin-culture populations
CFU.ttest.mono = CFU.ttest[!grepl("Co-culture", CFU.ttest$Species),]

d1.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "1",]
t.test(d1.mono$Mean_cfu ~ d1.mono$Count_of)

d4.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "4",]
t.test(d4.mono$Mean_cfu ~ d4.mono$Count_of)

d7.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "7",]
t.test(d7.mono$Mean_cfu ~ d7.mono$Count_of)

d10.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "10",]
t.test(d10.mono$Mean_cfu ~ d10.mono$Count_of)

d13.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "13",]
t.test(d13.mono$Mean_cfu ~ d13.mono$Count_of)

d16.mono = CFU.ttest.mono[CFU.ttest.mono$Day == "16",]
t.test(d16.mono$Mean_cfu ~ d16.mono$Count_of)

# t.test of lal-culture populations co and single
CFU.ttest.lal = CFU.ttest[grepl("LaL", CFU.ttest$Count_of),]

d1.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "1",]
t.test(d1.lal$Mean_cfu ~ d1.lal$Species)

d4.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "4",]
t.test(d4.lal$Mean_cfu ~ d4.lal$Species)

d7.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "7",]
t.test(d7.lal$Mean_cfu ~ d7.lal$Species)

d10.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "10",]
t.test(d10.lal$Mean_cfu ~ d10.lal$Species) # not below

d13.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "13",]
t.test(d13.lal$Mean_cfu ~ d13.lal$Species)

d16.lal = CFU.ttest.lal[CFU.ttest.lal$Day == "16",]
t.test(d16.lal$Mean_cfu ~ d16.lal$Species)

# t.test of lemthin-culture populations co and singler
CFU.ttest.lem = CFU.ttest[grepl("LeM", CFU.ttest$Count_of),]

d1.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "1",]
t.test(d1.lem$Mean_cfu ~ d1.lem$Species)

d4.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "4",] # not sig
t.test(d4.lem$Mean_cfu ~ d4.lem$Species)

d7.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "7",]
t.test(d7.lem$Mean_cfu ~ d7.lem$Species)

d10.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "10",]
t.test(d10.lem$Mean_cfu ~ d10.lem$Species)

d13.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "13",] # not sig
t.test(d13.lem$Mean_cfu ~ d13.lem$Species)

d16.lem = CFU.ttest.lem[CFU.ttest.lem$Day == "16",]
t.test(d16.lem$Mean_cfu ~ d16.lem$Species)
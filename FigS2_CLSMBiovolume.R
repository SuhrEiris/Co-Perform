
###############################     General Settings    ###################

## The package
#devtools::install_github("Russel88/RCon3D")
library(RCon3D)

## Other packages
library(foreach)
library(reshape2) # melt
library(rgl) # plot3d
library(COEF) # fancy_scientific (devtools::install_github("Russel88/COEF"))


######################################### Load and prepare images ###########################################
# Specify paths of images
#path <- "C:/Users/ncz500/Dropbox/MME_Co-Perform/Experiments_2019/Confocal/ConvertData/"
path <- "~/Desktop/Github/Co-Perform/Confocal/ConvertData"

setwd(gsub("/ConvertData","",path))

# Convert tiff to R-friendly RDS files (first time only)
# This might give warnings about TIFF image tags that cannot be read. This should not be a problem
# This creates the RDS files, and when we run the script again we can just use "findIMG" to get the paths for the RDS files
# img contains the paths for the RDS files, and this is the input for all subsequent analyses
img <- loadIMG(path, channels = c("C1"))
# img <- findIMG(path)

# We test the dimensions of the loaded images to ensure they are as expected
lapply(1:length(img), function(x) dim(readRDS(img[x])))

# Thresh-holding
imgT <- threshIMG(img, method = "Otsu")

# See the out-put of thresh-holding
mmand::display(readRDS(imgT[1])[,,5])

# Image details in microns (check these!)
pwidth <- 1.78 * 2
zstep <- 2

######################################## Biomass #########################################
# naming makes variables by looking in the names of the image files. 
biomass <- quant(imgT, channels = c("C1"), cores = 1,
                     naming = list(Type = c("LAL","COC","LEM"),
                                   Date = c("0409","0509","2908"),
                                   Rep=c("RepA","RepB","RepC","RepD","RepE","RepF","RepG")))

biomass2 <- biomass %>% group_by(Date, Type, Rep) %>% dplyr::summarise(Observations = n(),
                                                           biomass = sum(Count))

biomass3 <- biomass2 %>% group_by(Date, Type) %>% dplyr::summarise(Observations = n(),
                                                            biomass = mean(biomass))


# The biomass in cubic micron
biomass3$Microns <- biomass3$biomass * pwidth^2 * zstep
biomass3$biomass_pr_cubicmicrons = biomass3$Microns/(1.225*10^6)
biomass3$log_microns = log10(biomass3$Microns)

biomass4 = biomass3 %>%
  group_by(Type) %>%
  dplyr::summarise(observations = n(),
                   mean_microns = mean(Microns),
                   log_mean_microns = mean(log_microns),
                   sd_microns = sd(Microns),
                   sd_log_microns = sd(log_microns))

biomass4$CI =  qnorm(0.975)*biomass4$sd_microns/sqrt(biomass4$observations)
biomass4$log_CI = qnorm(0.975)*biomass4$sd_log_microns/sqrt(biomass4$observations)

biomass3$Type[biomass3$Type=="LAL"] <- "Mono-culture"
biomass3$Type[biomass3$Type=="COC"] <- "Co-culture"
# Plot
PS2 <- biomass3[biomass3$Type != "LEM",] %>% 
  ggplot(aes(x = Type, y = log_microns, fill = Type, color = Type)) +
  geom_point(alpha = 0.9, position = position_jitter(width = 0.1), size = 2)+
  stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), geom = "pointrange", color = "black", 
               size = 0.5, alpha = 1)+
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#CCC591", "#798E87"))+
  scale_color_manual(values = c("#CCC591", "#798E87"))+
  ylim(5,7)+
  xlab("")+
  theme(legend.position = "none",
        axis.text = element_text(color = "Black", face = "bold"), 
        axis.title = element_text(color = "Black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  coord_flip()+
  ylab(expression(paste("Biovolume (",mu,m^3,")")))+
  theme(axis.text.y = element_text(color = "Black", face = "bold"))+
  geom_segment(aes(x = 1, y = 7, xend = 2, yend = 7), color = "black", size = 0.35)+
  geom_segment(aes(x = 1, y = 6.95, xend = 1, yend = 7), color = "black", size = 0.35)+
  geom_segment(aes(x = 2, y = 6.95, xend = 2, yend = 7), color = "black", size = 0.35)+
  annotate(geom = "text", label = "Padj = 0.023", x = 1.5, y = 6.7, family="sans", size = 2, fontface = 2)
PS2

ggsave("FigS2.pdf",
       device = "pdf",
       plot = PS2,
       units="mm",
       width=88,
       height=88,
       dpi = 300)

#T-test
stat = biomass3 %>%
  group_by(Type) %>%
  summarise(count = n(),
            mean_Microns = mean(Microns),
            std_Microns = sd(Microns))

### Regular stats
# Test approach - paired with equal variance
L = t.test(log(biomass3[biomass3$Type == "Mono-culture",]$Microns), log(biomass3[biomass3$Type == "Co-culture",]$Microns), 
       paired = T, var.equal = TRUE)

p.adjust(L$p.value, method = "fdr")

# # Linear model approach 
# biomass3$Type <- as.factor(as.character(biomass3$Type))
# fit.lm <- lm(Microns ~ Type, data = biomass3[biomass3$Type != "LEM",])
# summary(fit.lm)
# 
# ### Foldchange approach:
# fold <- cbind(biomass3[biomass3$Type == "LAL",], biomass3[biomass3$Type == "COC",])
# fold2 <- fold %>% dplyr::mutate(Ratio = log2(Microns / Microns1)) %>% select(c(1,11))
# 
# # linear model approach
# fit.ratio <- lm(Ratio ~ 1, data = fold2, offset = rep(0, nrow(fold2)))
# summary(fit.ratio)
# 
# # t-test approach
# t.test(fold2$Ratio, mu=0)



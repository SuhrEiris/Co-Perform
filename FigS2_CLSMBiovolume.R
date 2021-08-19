
###############################     General Settings    ###################

# Set working directory
setwd("C:/Users/ncz500/Dropbox/MME_Co-Perform/Experiments_2019/Confocal")

# Set working directory (Nathalie)
setwd("C:/Users/nasuh/Dropbox/MME_Co-Adatpt/Experiments_2019/Confocal")

## The package
#devtools::install_github("Russel88/RCon3D")
library(RCon3D)

## Other packages
library(foreach)
library(ggplot2) 
library(reshape2) # melt
library(rgl) # plot3d
library(COEF) # fancy_scientific (devtools::install_github("Russel88/COEF"))
library(dplyr)
library(tidyr)


######################################### Load and prepare images ###########################################
# Specify paths of images
#path <- "C:/Users/ncz500/Dropbox/MME_Co-Perform/Experiments_2019/Confocal/ConvertData/"
path <- "C:/Users/nasuh/Dropbox/MME_Co-Adatpt/Experiments_2019/Confocal/ConvertData/"

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

biomass3$log_microns = log(biomass3$Microns)

biomass4 = biomass3 %>%
  group_by(Type) %>%
  dplyr::summarise(observations = n(),
                   mean_microns = mean(Microns),
                   log_mean_microns = mean(log_microns),
                   sd_microns = sd(Microns),
                   sd_log_microns = sd(log_microns))

biomass4$CI =  qnorm(0.975)*biomass4$sd_microns/sqrt(biomass4$observations)
biomass4$log_CI = qnorm(0.975)*biomass4$sd_log_microns/sqrt(biomass4$observations)


# Plot
p <- biomass3[biomass3$Type != "LEM",] %>% ggplot(aes(Type, Microns, fill = Type)) +
  geom_boxplot(alpha= 1) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.15, fill = "black", binwidth = 0.4, alpha = 1) +
  stat_summary(fun = mean, color = "black", shape = 3)+
  #stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange")+
  theme_bw() +
  scale_fill_manual(values = c("#CCC591", "#798E87"))+
  scale_y_log10()+
  xlab("")+
  theme(legend.title = element_text(size = 20), 
        legend.text = element_text(size=20),
        axis.text = element_text(size = 18, color = "Black"), 
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size=20))+
  ylab(expression(paste("Biovolume (",mu,m^3,")")))
p

#T-test
stat = biomass3 %>%
  group_by(Type) %>%
  summarise(count = n(),
            mean_Microns = mean(Microns),
            std_Microns = sd(Microns))

### Regular stats
# Test approach - paired with equal variance
t.test(log(biomass3[biomass3$Type == "LAL",]$Microns), log(biomass3[biomass3$Type == "COC",]$Microns), 
       paired = T, var.equal = TRUE)

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



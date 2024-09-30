# TL QUAD Data analysis 

# Written by T. Lindsay 
# 9.25.2024

# Libraries & Data  -------------------------------------------------------

library(tidyverse)
library(dplyr)
library(car)
library(performance) # for check_model()
library(lmtest) # for breusht-pagan test
library(nlme) # for stepping gls models 
library(MASS) # for stepAIC()

library(ggpubr)# for R2 and pvaules 
library(ggplot2)
library('gridExtra')  #to arrange plots 
library('grid')       #arranging plots 
library(ggpmisc)

rm(list=ls())
setwd("/Users/tayrlindsay/Desktop/GITHUB/Pub-Astrangia-Gradient/")

# Data prep  ---------------------------------------------------------

# read and edit data 
data <- read.csv('TL_Quad_raw_data.csv')  %>%
  dplyr::select(algae, corrected_depth_m, mean_sym, mean_apo, light) %>%
  mutate(Sym = mean_sym*4) %>%
  mutate(Apo = mean_apo*4) %>%
  mutate(all_colonies = Sym + Apo) %>%
  mutate(bins = cut(corrected_depth_m, 
                    breaks = c(0,2,4,6,8,10,12,14,16,18,20,22,24))) 

# remove outliers 
data1 <- data %>%
  filter(!is.na(mean_sym )) %>%
  mutate(zscore = (.$mean_sym - mean(.$mean_sym ))/sd(.$mean_sym )) %>%
  filter(abs(zscore)<3) %>%
  filter(!is.na(mean_apo)) %>%
  mutate(zscore2 = (.$mean_apo- mean(.$mean_apo))/sd(.$mean_apo)) %>%
  filter(abs(zscore2)<3)

### CALCULATE PERCENT COVER 
# Percent cover based on sizes 
# colonies were 0.053629137 m^2 
data1 <- data1 %>%
  mutate(percent_sym = mean_sym*0.0053629137*100) %>%
  mutate(percent_apo = mean_apo*0.0053629137*100) %>%
  mutate(percent_AP = percent_apo + percent_sym)


# Model Light  ------------------------------------------------------------

# This analysis was based on an exponential decay curve 
# https://douglas-watson.github.io/post/2018-09_exponential_curve_fitting/ 

# fit a non linear curve using nls
# use the SSasymp function to guess where to start with the basic variables 
fit <- nls(light ~ SSasymp(corrected_depth_m, yf, y0, log_alpha), data = data)
summary(fit)

# formula   =   y(t)∼yf+(y0−yf)e−αt
# yf        =   -248.0086               # this value is the asymptotic value when depth approaches infinity 
# y0        =   10237.8397              # initial value when depth = 0
# log_alpha =   -1.3191                 # logarithm decay rate (natural log of the decay constant)
# alpha     =   0.2673758
# equation  =   light(depth) = -248.0086 + (10237.8397 - -248.0086)e^(-0.2673758*depth)
# simple    =   light(depth) = -248.0086 + (10485.8483)e^(-0.2673758*depth)

# create this as a function 
simplified_function <- function(depth) {
  -248.0086 + (10237.8397 - (-248.0086))*exp(-0.2673758*depth)
}

# Calculate the values for each depth based on the equation
data1 <- data1 %>%
  mutate(modeled_light = 
           -248.0086 + (10237.8397 - (-248.0086))*exp(-0.2673758*corrected_depth_m))

# Calculate mean abundance  -----------------------------------------------

### MEAN ABUNDANCE OF APO AND SYM 
# pivot longer
data2 <- data1 %>%
  pivot_longer(cols = c(Apo, Sym), names_to = "colony_type", values_to = "abundance") %>%
  dplyr::select(bins, corrected_depth_m, colony_type, abundance, modeled_light, algae)

# Calculate mean abundance by depth bins 
summary_merged <- data2 %>% 
  group_by(bins,colony_type) %>% 
  summarise_at(vars(abundance), 
               list(av = ~mean(., na.rm = TRUE), 
                    sd = ~sd(., na.rm = TRUE)
               ))

# create new df of numbers for break points 
break_points <- c(0,0,2,2,4,4,6,6,8,8,10,10,12,12,14,14,16,16,18,18,20,20,22,22)
#add to summary df 
summary_merged$breaks = break_points

### MEAN TOTAL ABUNDANCE
# Calculate means by bins 
# summarize by bins 
summary_merged_all <- data1 %>% 
  group_by(bins) %>% 
  summarise_at(vars(all_colonies), 
               list(av = ~mean(., na.rm = TRUE), 
                    sd = ~sd(., na.rm = TRUE)))

# create new df of numbers for break points 
break_points <- c(0,2,4,6,8,10,12,14,16,18,20, 22)
#add to summary df 
summary_merged_all$breaks = break_points

### MEAN PERCENT COVER 
# pivot longer
percent_cover_long <- data1 %>%
  pivot_longer(cols = c(percent_apo, percent_sym, algae), names_to = "cover_type", values_to = "percent_cover") %>%
  dplyr::select(bins, corrected_depth_m, cover_type, percent_cover)

# Calculate mean abundance by bins 
# summarize by bins 
percent_cover_summary <- percent_cover_long  %>% 
  group_by(bins,cover_type) %>% 
  summarise_at(vars(percent_cover), 
               list(av = ~mean(., na.rm = TRUE), 
                    sd = ~sd(., na.rm = TRUE)))

# create new df of numbers for break points 
break_points3 <- c(0,0,0,2,2,2,4,4,4,6,6,6,8,8,8,10,10,10,12,12,12,14,14,14,16,16,16,18,18,18,20,20,20,22,22,22)
#add to summary df 
percent_cover_summary$breaks = break_points3

# Figure 2. Light curve  --------------------------------------------------

#### Create light figure 

# PLOT Light x Depth 
light_modeled <- ggplot(data1) +
  # DATA
  geom_point(aes(x=corrected_depth_m,y=light), color = "black", alpha = 0.6, size=1, shape=1) + 
  stat_function(fun = simplified_function, color = "orange") + 
  # AESTHETICS 
  labs( y = expression(paste("Light Intensity (lum ", m^{-2}, ")")), x = "Depth (m below MLLW)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25000)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  theme_bw() +
  theme(text = element_text(size=15)) 
light_modeled

ggsave("TLAP_Quad_fig2_light.jpg", plot = light_modeled, path = 'Figures', width =4, height = 6)

#### Calculating R2 value 

# calculate the residuals 
residuals <- residuals(fit)
plot(data$corrected_depth_m, residuals)

# calculate sum of squares of residuals 
RSS <- sum(residuals^2)
# Calculate total sum of squares 
TSS <- sum((data$light - mean(data$light))^2)

# Calculate R-squared
r_squared <- 1 - RSS / TSS
r_squared

# Figure 3. Abundance  ------------------------------------------------------------------

# total abundance 
total_line <- ggplot(summary_merged_all, aes(x=breaks,y=av)) + 
  #DATA
  geom_point(size=3) + 
  geom_line(linetype="dashed") + 
  geom_errorbar(aes(ymin = av - sd, ymax = av + sd), width = 0.1) +
  #AESTHETICS
  labs(x = "Depth (m)", y="Mean colony density" ~(m^2)) + 
  theme_bw() + 
  theme(text = element_text( size=15),  legend.position = "bottom") 
total_line
#ggsave("TLAP_Quad_lines.jpg", plot = line_by_ecotype , path = 'Figures', width =10, height = 5)

# abundance x ecotype 
line_by_ecotype <- ggplot(summary_merged, aes(x=breaks,y=av, color=colony_type)) + 
  #DATA
  geom_point(size=5) + 
  geom_line(linetype="dashed", size=2) + 
  geom_vline(xintercept = 10.5, linetype="dashed", linewidth=0.8) + 
  geom_vline(xintercept = 13, linetype="dashed", linewidth=0.8) + 
  geom_errorbar(aes(ymin = av - sd, ymax = av + sd, color=colony_type), width = 0.1) +
  #AESTHETICS
  scale_color_manual(
    values = c("Apo" = "#bf9e72", "Sym" = "#7F1734"),
    labels=c('Aposymbiotic', 'Symbiotic'), 
    name = "Ecotype") + 
  labs(x = "Depth (m)", y="Mean colony density" ~(m^2)) + 
  theme_bw() + 
  theme(text = element_text(size=15),  legend.position = c(0.8, 0.8)) 
line_by_ecotype 

ggsave("TLAP_Quad_fig3_ecotype.jpg", plot = line_by_ecotype , path = 'Figures', width =10, height = 5)

# Figure 4. Percent Cover  ----------------------------------------------------------------

percent_cover <- ggplot(percent_cover_summary, aes(x=breaks,y=av, fill=cover_type)) + 
  #DATA
  geom_area(color="black") + 
  geom_vline(xintercept = 10.5, linetype="dashed", linewidth=0.8) + 
  geom_vline(xintercept = 13, linetype="dashed", linewidth=0.8) + 
  #geom_line(linetype="dashed") + 
  #geom_errorbar(aes(ymin = av - sd, ymax = av + sd), width = 0.1) +
  #AESTHETICS
  scale_fill_manual(
    values = c("percent_apo" = "#bf9e72", "percent_sym" = "#7F1734", algae="#1B6B22"),
    labels = c("Macroalgae","Aposymbiotic", "Symbiotic"), 
    name = "Ecotype") + 
  labs(x = "Depth (m)", y="Benthic cover (%)") + 
  theme_bw() + 
  theme(text = element_text( size=15),  
        legend.position = c(0.8, 0.8)) 
percent_cover

ggsave("TLAP_Quad_fig4_percent_cover.jpg", plot = percent_cover, path = 'Figures', width =10, height = 5)


# Figure 5. Linear Regressions ----------------------------------------------------------------

# set formula to be a line 
formula <- y ~ x   

# filter data into the three depth categories 
data_s <- data2 %>% filter(corrected_depth_m < 10.5)
data_m <- data2 %>% 
  filter(corrected_depth_m > 10.5) %>% 
  filter(corrected_depth_m < 13)
data_d <- data2 %>% filter(corrected_depth_m > 13) %>%
  filter(colony_type == "Apo")

# generate legend  
for_the_legend <- ggplot(data_s, aes(x=corrected_depth_m, y = abundance, color = colony_type)) + 
  # DATA
  geom_point(size = 2) + 
  geom_smooth(method="lm", size=2) +
  # AESTHETICS 
  theme_bw() + 
  scale_color_manual(
    values = c("Apo" = "#bf9e72", "Sym" = "#7F1734"),
    labels = c("Aposymbiotic", "Symbiotic"), 
    name = "Ecotype") + 
  theme(text = element_text(size=15),
        #theme(legend.position = c(0.6, 0.9),
        legend.title=element_text(size=25, face= "bold"),
        legend.text=element_text(size=25),
        legend.key.width = unit(3,"cm"))

get_only_legend <- function(plot) { 
  # get tabular interpretation of plot 
  plot_table <- ggplot_gtable(ggplot_build(plot))  
  #  Mark only legend in plot 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")  
  # extract legend 
  legend <- plot_table$grobs[[legend_plot]] 
  # return legend 
  return(legend)  
}

led <- get_only_legend(for_the_legend)

# MACROALGAE ZONE x DEPTH  
plot_shallow_depth <- ggplot(data_s, aes(x=corrected_depth_m, y = abundance, color = colony_type)) + 
  # DATA
  geom_point(size = 2) + 
  geom_smooth(method="lm", size=2) + 
  # AESTHETICS 
  theme_bw() + 
  scale_color_manual(
    values = c("Apo" = "#bf9e72", "Sym" = "#7F1734"),
    labels = c("Aposymbiotic", "Symbiotic"), 
    name = "Ecotype") + 
  ggtitle("Macroalgae Zone") + 
  labs(x = "Depth (m)", y="") + 
  theme(text = element_text(size=25),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))

# CORAL ZONE x DEPTH  
plot_middle_depth <- ggplot(data_m, aes(x=corrected_depth_m, y = abundance, color = colony_type)) + 
  # DATA
  geom_point(size = 2) + 
  geom_smooth(method="lm", size=2) + 
  # AESTHETICS 
  theme_bw() + 
  scale_color_manual(
    values = c("Apo" = "#bf9e72", "Sym" = "#7F1734"),
    labels = c("Aposymbiotic", "Symbiotic"), 
    name = "Ecotype") + 
  ggtitle("Coral Zone") + 
  labs(x = "Depth (m)", y="")+  #title = "Middle Zone (10.5-13m)"
  theme(text = element_text(size=25),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))

# DEEP ZONE x DEPTH  
plot_deep_depth <- ggplot(data_d, aes(x=corrected_depth_m, y = abundance, color = colony_type)) + 
  # DATA
  geom_point(size = 2) + 
  geom_smooth(method="lm", size=2) + 
  # AESTHETICS 
  theme_bw() + 
  scale_color_manual(
    values = c("Apo" = "#bf9e72", "Sym" = "#7F1734"),
    labels = c("Aposymbiotic", "Symbiotic"), 
    name = "Ecotype") + 
  ggtitle("Deep Zone") + 
  labs(x = "Depth (m)", y="") + #title = "Deep Zone (13-24m)", 
  theme(text = element_text(size=25),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))

plot_shallow_depth
plot_middle_depth
plot_deep_depth

# MACROALGAE ZONE x LIGHT
plot_shallow_light <- ggplot(data_s, aes(x=modeled_light, y = abundance, color = colony_type)) + 
  # DATA
  geom_point(size = 2) + 
  geom_smooth(method="lm", size=2) + 
  # AESTHETICS 
  theme_bw() + 
  scale_color_manual(
    values = c("Apo" = "#bf9e72", "Sym" = "#7F1734"),
    labels = c("Aposymbiotic", "Symbiotic"), 
    name = "Ecotype") + 
  labs( x = expression(paste("Light Intensity (lum ", m^{-2}, ")")), y="") + #Mean Abundance" ~(m^2))
  theme(text = element_text(size=25),
        legend.position = "none")

# CORAL ZONE x LIGHT
plot_middle_light <-ggplot(data_m, aes(x=modeled_light, y = abundance, color = colony_type)) + 
  # DATA
  geom_point(size = 2) + 
  geom_smooth(method="lm", size=2) + 
  # AESTHETICS 
  theme_bw() + 
  scale_color_manual(
    values = c("Apo" = "#bf9e72", "Sym" = "#7F1734"),
    labels = c("Aposymbiotic", "Symbiotic"), 
    name = "Ecotype") + 
  labs( x = expression(paste("Light Intensity (lum ", m^{-2}, ")")), y="") +
  theme(text = element_text(size=25),
        legend.position = "none")

# DEEP ZONE x LIGHT  
plot_deep_light <-ggplot(data_d, aes(x=modeled_light, y = abundance, color = colony_type)) + 
  # DATA
  geom_point(size = 2) + 
  geom_smooth(method="lm", size=2) +
  # AESTHETICS 
  theme_bw() + 
  scale_color_manual(
    values = c("Apo" = "#bf9e72", "Sym" = "#7F1734"),
    labels = c("Aposymbiotic", "Symbiotic"), 
    name = "Ecotype") + 
  labs(x = expression(paste("Light Intensity (lum ", m^{-2}, ")")), y="") +
  theme(text = element_text(size=25),
        legend.position = "none")

plot_shallow_light
plot_middle_light
plot_deep_light

# MACROALGAE ZONE x MACROAGALE
plot_shallow_algae <-ggplot(data_s, aes(x=algae, y = abundance, color = colony_type)) + 
  # DATA
  geom_point(size = 2) + 
  geom_smooth(method="lm", size=2) + 
  # AESTHETICS 
  theme_bw() + 
  scale_color_manual(
    values = c("Apo" = "#bf9e72", "Sym" = "#7F1734"),
    labels = c("Aposymbiotic", "Symbiotic"), 
    name = "Ecotype") + 
  labs(x ="Macroalgae Cover (%)", y="")+  # Mean Abundance" ~(m^2))
  theme(text = element_text(size=25), 
        legend.position = "none")

plot_shallow_algae

# Arrange plots into one 
#top <- textGrob("Macroalgae Zone           Coral Zone           Deep Zone",gp = gpar(fontsize = 25))
yleft <- textGrob(expression(bold("Mean Abundance" ~(colonies/m^2))), rot = 90, gp = gpar(fontsize = 25))
tp_arrange <- grid.arrange(plot_shallow_depth, plot_middle_depth, plot_deep_depth, plot_shallow_light, plot_middle_light, plot_deep_light, plot_shallow_algae,led, nrow=3,  left=yleft)

tp_arrange

ggsave("TLAP_Quad_fig5_linear_models.jpg", plot = tp_arrange, path = 'Figures', width =20, height = 15)

# normality testing -------------------------------------------------------

##### SYM 

#model 
model_sym <- lm(mean_sym ~ algae *corrected_depth_m* modeled_light, data = data1)
vif(model_sym)

# is data normal? 
ggdensity(data1$mean_sym, xlab = "mean sym")
#Shapiro-wilks for 
shapiro.test(data1$mean_sym) #p= <0.000001, data is not normal 

# is data homoscedastic? 
ggqqplot(data1$mean_sym)
#Breusch-Pagan test
bptest(model_sym) #p = <0.000001, data is heteroscetastic (variance of residuals is not constant)

##### APO

#model 
model_apo <- lm(mean_apo ~ algae + corrected_depth_m + modeled_light, data = data1)
vif(model_apo)

# is data normal? 
ggdensity(data1$mean_apo, xlab = "mean apo")
#Shapiro-wilks for 
shapiro.test(data1$mean_apo) #p= <0.000001, data is not normal 

# is data homoscedastic? 
ggqqplot(data1$mean_apo)
#Breusch-Pagan test
bptest(model_apo) #p = 0.005, data is heteroscetastic (variance of residuals is not constant)

##### ALL COLONIES 

#model 
model_all <- lm(all_colonies ~ algae + corrected_depth_m + modeled_light, data = data1)
vif(model_all)

# is data normal? 
ggdensity(data1$all_colonies, xlab = "mean apo")
#Shapiro-wilks for 
shapiro.test(data1$all_colonies) #p= <0.000001, data is not normal 

# is data homoscedastic? 
ggqqplot(data1$all_colonies)
#Breusch-Pagan test
bptest(model_all) #p = 0.001, data is heteroscetastic (variance of residuals is not constant)


# Generalized Linear Models -----------------------------------------------

# I am using a GLm because it can handle non-normal data
# add 0.0001 to all data so that we can use gamma distribution 
data3 <- data1 %>%
  mutate(mean_sym=mean_sym+0.0001) %>%
  mutate(mean_apo=mean_apo+0.0001) %>%
  mutate(all_colonies=all_colonies+0.0001)

#all
glm_all <- glm(all_colonies ~ corrected_depth_m * modeled_light * algae, data = data3, family = Gamma())
summary(glm_all)
glm_all_best <- stepAIC(glm_all, direction = "both") # 
summary(glm_all_best) # summary of best model 
glm_all_r <- resid(glm_all_best) # pull residuals of best model 
ggqqplot(glm_all_r) # plot residuals in q-q plot 
vif(glm_all_best) # test for multicolinearity 

#apo
glm_apo <- glm(mean_apo ~ corrected_depth_m * modeled_light*algae, data = data3, family = Gamma())
summary(glm_apo)
glm_apo_best <- stepAIC(glm_apo, direction = "both") # 
summary(glm_apo_best) # summary of best model 
glm_apo_r <- resid(glm_apo_best) # pull residuals of best model 
ggqqplot(glm_apo_r) # plot residuals in q-q plot 
vif(glm_apo_best) # test for multicolinearity 

#sym
glm_sym <- glm(mean_sym ~ corrected_depth_m * modeled_light*algae, data = data3, family = Gamma())
summary(glm_sym)
glm_sym_best <- stepAIC(glm_sym, direction = "both") # 
summary(glm_sym_best) # summary of best model 
glm_sym_r <- resid(glm_sym_best) # pull residuals of best model 
ggqqplot(glm_sym_r) # plot residuals in q-q plot 
vif(glm_sym_best) # test for multicolinearity 

# Supplemental ------------------------------------------------------------

### A priori symbiont density 

#AP_Apriori
ap_raw <- read.csv('AP_Sym_Apriori.csv')

# boxplot
a_priori_plot <- ggplot(ap_raw, aes(x=ecotype, y=Cells.cm2, color=ecotype)) +
  # DATA 
  geom_boxplot() +
  # AESTHETICS 
  theme_bw()+
  labs(x= "Ecotype", y= expression(paste("Symbiont cells per ", cm^{-2})))+
  scale_color_manual(
    values = c("Apo" = "#bf9e72", "Sym" = "#7F1734"),
    labels = c("Aposymbiotic", "Symbiotic")) + 
  theme(text = element_text(size=25),
        legend.position = "none")
a_priori_plot        

# t-test 
a_priori_plot2 = a_priori_plot + stat_compare_means(method = "t.test", size = 5)
a_priori_plot2

ggsave("TL_Quad_S1_apriori.jpg", plot = a_priori_plot2, path = 'Figures', height = 8, width = 5)


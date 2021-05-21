library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggridges)
library(lme4)
library(car)

setwd("~/Google Drive/Rudolf_MRE/results/ch4/cusper/")

## load data
s1_1 = read.csv("./cusper_liqculture/exp1/S1DATA-ADJ-FIXED.csv", sep = ",", header = TRUE)
s1_2 = read.csv("./cusper_liqculture/exp1/S2DATA-ADJ-FIXED.csv", sep = ",", header = TRUE)
s1_3 = read.csv("./cusper_liqculture/exp1/S3DATA-ADJ-FIXED.csv", sep = ",", header = TRUE)
s2_1 = read.csv("./cusper_liqculture/exp2/S1Data-Adjusted.csv", sep = ",", header = TRUE)
s2_2 = read.csv("./cusper_liqculture/exp2/S2Data-Adjusted.csv", sep = ",", header = TRUE)
s2_3 = read.csv("./cusper_liqculture/exp2/S3Data-Adjusted.csv", sep = ",", header = TRUE)

data = rbind(s1_1, s1_2, s1_3, s2_1, s2_2, s2_3) %>% 
      mutate(timepoint = case_when(
            time == "T0" ~ 0,
            time == "T1" ~ 1,
            time == "T2" ~ 2,
            time == "T3" ~ 3,
            time == "T4" ~ 4))

summ = data %>% 
      filter(time == "T0") %>% 
      group_by(exp, sample, rep) %>% 
      summarise(mean_t0 = mean(intensity),
                median_t0 = median(intensity))
head(summ)

## reproductive success
cusper <-  data %>% 
      filter(intensity > 0) %>% 
      inner_join(., summ, by = c("exp", "sample", "rep")) %>% 
      mutate(success = log2(mean_t0/intensity),
             success2 = log2(median_t0/intensity)) %>% 
      dplyr::select(exp, sample, rep, time, timepoint, success, success2)
cusper

cusper_mean = ggplot(cusper[cusper$exp=="exp2",], aes(x=success, y=time))+
      geom_density_ridges(alpha = 0.3, color = "dark grey", scale = 1.5)+
      scale_x_continuous(limit = c(-2, 7), breaks = seq(0,6,1))+
      scale_y_discrete(expand = expansion(mult = c(0, 0.3))) +
      theme_ridges(grid = TRUE, center_axis_labels = TRUE)

cusper_median = ggplot(cusper[cusper$exp=="exp2",], aes(x=success2, y=time))+
      geom_density_ridges(alpha = 0.3, color = "dark grey", scale = 1.5)+
      scale_x_continuous(limit = c(-2, 7), breaks = seq(0,6,1))+
      scale_y_discrete(expand = expansion(mult = c(0, 0.3))) +
      theme_ridges(grid = TRUE, center_axis_labels = TRUE)

cusper_mean + cusper_median


## model
lm = lmer(success ~ timepoint + (1|exp/sample), data = cusper)
summary(lm)
qqnorm(resid(lm))

Anova(lm, type = "III")

## save plots
pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_liq.pdf", height = 10, width = 12, useDingbats = F)
cusper_mean
dev.off()

library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggridges)
library(lme4)
library(car)
library(gamlss)

setwd("~/Google Drive/Rudolf_MRE/results/ch4/cusper/")

###         LIQUID CULTURE

## load data
s1_1 = read.csv("./cusper_liqculture/exp1/S1DATA-ADJ-FIXED.csv", sep = ",", header = TRUE)
s1_2 = read.csv("./cusper_liqculture/exp1/S2DATA-ADJ-FIXED.csv", sep = ",", header = TRUE)
s1_3 = read.csv("./cusper_liqculture/exp1/S3DATA-ADJ-FIXED.csv", sep = ",", header = TRUE)
s2_1 = read.csv("./cusper_liqculture/exp2/S1Data-Adjusted.csv", sep = ",", header = TRUE)
s2_2 = read.csv("./cusper_liqculture/exp2/S2Data-Adjusted.csv", sep = ",", header = TRUE)
s2_3 = read.csv("./cusper_liqculture/exp2/S3Data-Adjusted.csv", sep = ",", header = TRUE)

liq = rbind(s1_1, s1_2, s1_3, s2_1, s2_2, s2_3) %>% 
      mutate(timepoint = case_when(
            time == "T0" ~ 0,
            time == "T1" ~ 1,
            time == "T2" ~ 2,
            time == "T3" ~ 3,
            time == "T4" ~ 4))

summ_liq = liq %>% 
      filter(time == "T0") %>% 
      group_by(exp, sample, rep) %>% 
      summarise(mean_t0 = mean(intensity),
                median_t0 = median(intensity))
head(summ_liq)

## reproductive success
cusper_liq <-  liq %>% 
      filter(intensity > 0) %>% 
      inner_join(., summ_liq, by = c("exp", "sample", "rep")) %>% 
      mutate(success = log2(mean_t0/intensity),
             success2 = log2(median_t0/intensity)) %>% 
      dplyr::select(exp, sample, rep, time, timepoint, success, success2)
head(cusper_liq)


liq_lod = liq %>% 
      filter(exp == "exp2") %>% 
      filter(time == "T4") %>% 
      group_by(sample, rep) %>% 
      summarise(fl = min(mean),
            bg = mean(bakground),
            int = abs(fl - bg)) %>% 
      mutate(lod = log2(mean(summ_liq$mean_t0)/int)) %>% 
      ungroup() %>% 
      summarise(mean_lod = mean(lod),
                sd_lod = sd(lod),
                up = mean_lod + sd_lod,
                low = mean_lod - sd_lod)
liq_lod

cusper_liq_n = cusper_liq %>% filter(exp == "exp2") %>% group_by(sample, time)
cusper_liq_n = sample_n(cusper_liq_n, 100, replace = TRUE)

cusper_liq_plot = ggplot(cusper_liq[cusper_liq$exp=="exp2",], aes(x=success, y=time))+
      geom_density_ridges(alpha = 0.5, color = "dark grey", scale = 1.5)+
      scale_x_continuous(limit = c(-2, 10), breaks = seq(0,10,1))+
      scale_y_discrete(expand = expansion(mult = c(0, 0.3))) +
      theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
      geom_vline(data = liq_lod, aes(xintercept = mean_lod))+
      geom_vline(data = liq_lod, aes(xintercept = up), linetype = "dashed")+
      geom_vline(data = liq_lod, aes(xintercept = low), linetype = "dashed")
cusper_liq_plot

qq_liq = ggplot(cusper_liq_n)+
      geom_qq(aes(sample = success, color = time), alpha = 0.6, size = 3)+
      theme_bw()+
      scale_y_continuous(name = "Reproductive success (# generations)",
                         limits = c(-2,10), breaks = seq(0,10,2), expand = c(0,0))+
      scale_x_continuous(name = "Theoretical quantiles",
                         limits = c(-4,4), expand = c(0,0))+
      scale_color_brewer(palette = "Set2")+
      geom_hline(data = liq_lod, aes(yintercept = mean_lod))+
      geom_hline(data = liq_lod, aes(yintercept = up), linetype = "dashed")+
      geom_hline(data = liq_lod, aes(yintercept = low), linetype = "dashed")+
      coord_fixed(ratio = 8/12)
qq_liq

###         PLANT

## load data
ath = read.csv("./cusper_arabidopsis/plant_cusper.csv", sep = ",", header = TRUE)

ath = ath %>% 
      mutate(timepoint = case_when(
            time == "T0" ~ 0,
            time == "T1" ~ 16,
            time == "T2" ~ 24,
            time == "T3" ~ 48))

summ_ath = ath %>% 
      filter(time == "T0") %>% 
      group_by(sample, rep) %>% 
      summarise(mean_t0 = mean(intensity),
                median_t0 = median(intensity))
head(summ_ath)

## reproductive success
cusper_ath <-  ath %>% 
      filter(intensity > 0) %>% 
      inner_join(., summ, by = c("sample", "rep")) %>% 
      mutate(success = log2(mean_t0/intensity),
             success2 = log2(median_t0/intensity)) %>% 
      dplyr::select(sample, rep, time, timepoint, success, success2)
head(cusper_ath)

ath_lod = ath %>% filter(intensity>10) %>% slice_min(intensity) %>% 
      mutate(lod = log2(mean(summ_ath$mean_t0)/(intensity/2)))
head(ath_lod)

ath_lod = ath %>% 
      filter(time == "T3") %>% 
      filter(rep != "3") %>% 
      group_by(rep) %>% 
      summarise(fl = min(mean),
                bg = mean(bakground),
                int = abs(fl - bg)) %>% 
      mutate(lod = log2(mean(summ_liq$mean_t0)/int)) %>% 
      ungroup() %>% 
      summarise(mean_lod = mean(lod),
                sd_lod = sd(lod),
                up = mean_lod + sd_lod,
                low = mean_lod - sd_lod)
ath_lod

cusper_ath_plot = ggplot(cusper_ath, aes(x=success, y=time))+
      geom_density_ridges(alpha = 0.5, color = "dark grey", scale = 1.5)+
      scale_x_continuous(limit = c(-2, 10), breaks = seq(0,10,1))+
      scale_y_discrete(expand = expansion(mult = c(0, 0.3))) +
      theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
      geom_vline(data = ath_lod, aes(xintercept = mean_lod))+
      geom_vline(data = ath_lod, aes(xintercept = up), linetype = "dashed")+
      geom_vline(data = ath_lod, aes(xintercept = low), linetype = "dashed")
cusper_ath_plot

qq_ath = ggplot(cusper_ath, aes(sample = success, color = time))+
      geom_qq(alpha = 0.6, size = 3)+
      theme_bw()+
      scale_y_continuous(name = "Reproductive success (# generations)",
                         limits = c(-2,10), breaks = seq(0,10,2), expand = c(0,0))+
      scale_x_continuous(name = "Theoretical quantiles",
                         limits = c(-4,4), expand = c(0,0))+
      scale_color_brewer(palette = "Set2")+
      coord_fixed(ratio = 8/12)+
      geom_hline(data = ath_lod, aes(yintercept = mean_lod))+
      geom_hline(data = ath_lod, aes(yintercept = up), linetype = "dashed")+
      geom_hline(data = ath_lod, aes(yintercept = low), linetype = "dashed")
qq_ath


##    MODELS
# fitting probability distributions for liq_cusper
hist(cusper_liq$success)
liq_hist1 = histDist(cusper_liq$success, "NO", density = TRUE) 
liq_hist2 = histDist(cusper_liq$success, "ST1", density = TRUE) 
liq_hist3 = histDist(cusper_liq$success, "ST2", density = TRUE) 
liq_hist4 = histDist(cusper_liq$success, "TF", density = TRUE)

GAIC(liq_hist1, liq_hist2, liq_hist3, liq_hist4)

#model liquid culture
cusper_liq$exp = factor(cusper_liq$exp)
cusper_liq$sample = factor(cusper_liq$sample)
cusper_liq$rep = factor(cusper_liq$rep)

model_liq = gamlss(success ~ timepoint + re(random= ~1|exp/sample/rep), data = cusper_liq, family = ST1)
summary(model_liq)
qqnorm(resid(model_liq))



# fitting probability distributions for ath_cusper
hist(cusper_ath$success)
ath_hist1 = histDist(cusper_ath$success, "NO", density = TRUE) 
ath_hist2 = histDist(cusper_ath$success, "ST1", density = TRUE) 
ath_hist3 = histDist(cusper_ath$success, "ST2", density = TRUE) 
ath_hist4 = histDist(cusper_ath$success, "TF", density = TRUE)

GAIC(ath_hist1, ath_hist2, ath_hist3, ath_hist4)

#model arabidopsis
cusper_ath$sample = factor(cusper_ath$sample)
cusper_ath$rep = factor(cusper_ath$rep)

model_ath = gamlss(success ~ timepoint + re(random= ~1|rep), data = cusper_ath, family = ST1)
summary(model_ath)
qqnorm(resid(model_ath))


## save plots
pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper.pdf", height = 12, width = 8, useDingbats = F)
cusper_liq_plot/cusper_ath_plot
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_qqplots.pdf", height = 10, width = 6, useDingbats = F)
qq_liq/qq_ath
dev.off()

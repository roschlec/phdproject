# Pe299R competition in planta

library(tidyverse)
library(patchwork)
library(emmeans)
library(gamlss)


setwd("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/")

df <- read.csv("comp_ath.csv", header = T)[,-6] %>% 
      mutate(comp = log10(pe299r)) %>% 
      data.frame(stringsAsFactors = FALSE)
head(df)

mono <- df %>% 
      filter(trt == "cus_mono") %>% 
      dplyr::select(time, rep, mono = comp)
wt <- df %>% 
      filter(trt == "cus_wt") %>% 
      dplyr::select(time, rep, wt = comp)
df2 <- df[,-5] %>% 
      filter(trt != "cus_mono" & trt != "cus_wt") %>% 
      left_join(., mono, by = c("time", "rep")) %>% 
      left_join(., wt, by = c("time", "rep")) %>% 
      pivot_longer(cols = c(comp, mono, wt), names_to = "type", values_to = "cfu") %>% 
      arrange(trt, type)
df2

sum_df <- df2 %>% 
      group_by(trt, time, hpi, type) %>% 
      summarise(mean = mean(cfu),
                sd = sd(cfu))
sum_df


plot = ggplot(df2, aes(x=hpi, y = cfu, color = type))+
      facet_wrap(~ trt, ncol = 3)+
      geom_point(alpha = 0.8, size = 1)+
      geom_line(data = sum_df, aes(y=mean), size = 1, alpha = 0.8)+
      theme_bw()+
      scale_y_continuous(limits = c(3,9), breaks = seq(4,8,2))+
      scale_x_continuous(limits = c(0,50), breaks = seq(0,48,24))+
      coord_fixed(ratio = 50/6)

df_aov = aov(comp ~ trt * time, data = df)
summary(df_aov)
qqnorm(df$comp)

emm_df = emmeans(df_aov, ~ trt | time, adjust = "bonferroni")
contrast(emm_df, "trt.vs.ctrl", ref = 3)
contrast(emm_df, "trt.vs.ctrl", ref = 8)

##  SAVE FILES
pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/cfu.pdf", width = 8, height = 8, useDingbats=FALSE)
plot
dev.off()

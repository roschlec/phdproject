### AREA UNDER THE CURVE FROM FLUORESCENCE CURVES

library(tidyverse)
library(MESS)
library(RColorBrewer)
library(circlize)
library(patchwork)
library(lme4)
library(MuMIn)
library(car)
library(gamlss)
library(reshape)

### FUNCTIONS
coalesce_all_columns <- function(df) {
  tibble(
        #tag = df$tag[1],
        exp = df$exp[1],
        #channel = df$channel[1],
        resource = df$resource[1],
        mix = df$mix[1],
        chg_sp1 = na.omit(df$chg_sp1),
        sp1 = na.omit(df$sp1),
        chg_sp2 = na.omit(df$chg_sp2),
        sp2 = na.omit(df$sp2)
  )
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))} # rescaling function (min-max rescaling)


##### COMPETITION FOR CARBON SOURCES #####

#setwd("~/Rudolf drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/comp_assay_1/combined/")    # windows
setwd("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/")      # mac
all_data <- read.csv("all_reps.csv", header = TRUE, sep = ",")
all_data = all_data %>% 
      mutate(tag = channel) %>% 
      separate(tag, into =c("tag", "exp"), sep = "_")

all_data$plateID <- factor(all_data$plateID)
all_data$channel <- factor(all_data$channel)
all_data$mix <- factor(all_data$mix)
all_data$strain <- factor(all_data$strain)
all_data$resource <- factor(all_data$resource)
all_data$tag <- factor(all_data$tag)
all_data$exp <- factor(all_data$exp)


pair = read.csv("./pairs.csv", header = TRUE) %>% 
      dplyr::select(mix, sp1, sp2, pair) %>% 
      pivot_longer(cols=c(sp1,sp2), names_to = "sp", values_to = "strain") %>% 
      mutate(tag = case_when(
            sp == "sp1" ~ "red",
            sp == "sp2" ~ "yellow"))
head(pair)
pair$sp = factor(pair$sp)
pair$tag = factor(pair$tag)

all_data = all_data %>% inner_join(., pair, by=c("mix", "strain", "tag"))
head(all_data)

## background correction
# create a column to identify wild type strain for filtering
all_data = all_data %>% 
  mutate(
    wt = case_when(
      str_detect(mix, "_wt") ~ "yes",
      TRUE ~ "no"))

# correct fluorescence measurements by subtracting the mean of wt fluorescence at each timepoint for each species and each condition 
data_wt = all_data %>% 
      group_by(channel) %>% 
      filter(wt == "yes") %>% 
      group_by(tag, exp, channel, strain, resource, time_days) %>% 
      summarise(mean_wt = mean(rfu)) # obtain the mean fluorescence of each wild type strain in each condition per time point

data2 = merge(all_data, data_wt, by = c("tag", "exp", "channel", "strain", "resource", "time_days")) %>% 
      group_by(tag, exp, channel, strain, resource, time_days) %>% 
      mutate(rfu_corr = rfu - mean_wt) # correct fluorescence by subtracting autofluorescence from fluorescence measurement
head(data2)

      # data inspection > plotting corrected fluorescence data
data2 %>% filter(tag == "red" & strain == "smfr1") %>% # select channel and strain to check
  ggplot(., aes(x=time, y=rfu_corr, color = exp))+
  facet_wrap(resource ~ mix, ncol = 8)+
  geom_point()+
  theme_bw()+
  labs(title = "Background-corrected fluorescence")

##  AUC
# Area under the curve of each replicate
auc_data <- na.omit(data2) %>% 
      group_by(tag, exp, channel, mix, strain, resource, pair, sp, rep) %>%
      summarise(AUC=auc(time, rfu_corr)) %>% ungroup()

      # data inspection > quick visualisation plot
auc_data %>% filter(tag == "red") %>% 
      ggplot(., aes(x=resource, y=AUC, fill=strain))+
      facet_wrap( ~ mix)+
      geom_point(pch=21, size=3, alpha=0.8, position=position_dodge(width=0.78))+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90),
            text = element_text(size=20))

# Area under the curve from no-carbon source growth condition
auc_noC <- auc_data %>%
      filter(resource == "noC") %>%
      group_by(tag, exp, mix, strain, channel, pair, sp) %>%
      summarise(AUC_noC = mean(AUC),
                qth = quantile(AUC, probs = .9),
                stdv = sd(AUC),
                mean_2sd = AUC_noC + 1.97*stdv)

# Filtering AUC by limit of detection (DL), which is the AUC of each combination under no carbon source media
auc_data2 = auc_data %>% 
      filter(resource != "noC") %>% 
      group_by(tag, exp, channel, mix, strain, resource, pair, sp) %>% 
      summarise(AUC = mean(AUC)) %>% 
      merge(., auc_noC, by=c("tag", "exp", "mix", "strain", "channel", "pair", "sp")) %>%
      mutate(
            AUC_DL = case_when(
                  AUC <= qth ~ qth,
                  TRUE ~ AUC)) %>% 
      dplyr::select(tag, exp, mix, strain, channel, pair, sp, resource, AUC_DL)
head(auc_data2)

# Relative change to INTRASPECIFIC COMPETITION
auc_intra = auc_data2 %>% 
      filter(!str_detect(mix, "mc")) %>% 
      separate(mix, into=c("sp1", "sp2"), sep = "_") %>% 
      filter(sp1 == sp2) %>% 
      group_by(tag, exp, channel, sp1, sp2, strain, resource, pair, sp) %>%
      summarise(AUC_intra = AUC_DL)

auc_all = auc_data2 %>% 
      filter(!str_detect(mix, "mc")) %>% 
      inner_join(., auc_intra[,c(1:3,6:7,10)], by=c("tag", "exp", "channel", "strain", "resource")) %>% 
      mutate(relchg = log2(AUC_DL/AUC_intra))

      
      # plot: relative change
ggplot(auc_all, aes(x=resource, y=relchg, color = tag))+
      facet_wrap(strain ~ mix, ncol=5, labeller = label_both)+
      geom_hline(yintercept=0, linetype= 2, alpha = 0.8)+
      geom_point(size = 3, alpha = 0.5)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90),
            text = element_text(size=20))+
      ylim(-3,3)+
      labs(title = "Relative change to monoculture",
           y = "Relative change (log2 scale)",
           x = "Resource")

# Geometry analysis -- Cartesian coordinates
all = auc_all %>% 
      dplyr::select(-AUC_DL, -AUC_intra, -tag, -channel) %>% 
      group_by(pair) %>% 
      pivot_wider(names_from = sp, values_from = c(relchg, strain)) %>% 
      mutate(chg_sp1 = case_when(
                  pair == "P11" ~ relchg_sp2,
                  pair == "P12" ~ relchg_sp2,
                  pair == "P13" ~ relchg_sp2,
                  pair == "P14" ~ relchg_sp2,
                  pair == "P15" ~ relchg_sp2,
                  pair == "P06" ~ relchg_sp2, #16 to 06
                  pair == "P07" ~ relchg_sp2, #17 to 07
                  pair == "P18" ~ relchg_sp2,
                  pair == "P19" ~ relchg_sp2,
                  pair == "P20" ~ relchg_sp2,
                  TRUE ~ relchg_sp1),
             chg_sp2 = case_when(
                   pair == "P11" ~ relchg_sp1,
                   pair == "P12" ~ relchg_sp1,
                   pair == "P13" ~ relchg_sp1,
                   pair == "P14" ~ relchg_sp1,
                   pair == "P15" ~ relchg_sp1,
                   pair == "P06" ~ relchg_sp1, #16 to 06
                   pair == "P07" ~ relchg_sp1, #17 to 07
                   pair == "P18" ~ relchg_sp1,
                   pair == "P19" ~ relchg_sp1,
                   pair == "P20" ~ relchg_sp1,
                   TRUE ~ relchg_sp2),
             sp1 = case_when(
                  pair == "P11" ~ strain_sp2,
                  pair == "P12" ~ strain_sp2,
                  pair == "P13" ~ strain_sp2,
                  pair == "P14" ~ strain_sp2,
                  pair == "P15" ~ strain_sp2,
                  pair == "P06" ~ strain_sp2, #16 to 06
                  pair == "P07" ~ strain_sp2, #17 to 07
                  pair == "P18" ~ strain_sp2,
                  pair == "P19" ~ strain_sp2,
                  pair == "P20" ~ strain_sp2,
                  TRUE ~ strain_sp1),
             sp2  = case_when(
                   pair == "P11" ~ strain_sp1,
                   pair == "P12" ~ strain_sp1,
                   pair == "P13" ~ strain_sp1,
                   pair == "P14" ~ strain_sp1,
                   pair == "P15" ~ strain_sp1,
                   pair == "P06" ~ strain_sp1, #16 to 06
                   pair == "P07" ~ strain_sp1, #17 to 07
                   pair == "P18" ~ strain_sp1,
                   pair == "P19" ~ strain_sp1,
                   pair == "P20" ~ strain_sp1,
                   TRUE ~ strain_sp2)) %>% 
      dplyr::select(exp, mix, pair, resource, chg_sp1, chg_sp2, sp1, sp2) %>% 
      ungroup() 
head(all)

all_cart = all %>%
      mutate(magnitude = sqrt(chg_sp1^2 + chg_sp2^2),
             radian = atan2(chg_sp2, chg_sp1),
             rad2 = case_when(
                   radian < 0 ~ radian + 2*pi,
                   TRUE ~ radian),
             degree = rad2 * 180/pi,
             coord_sp = case_when(
                   0 <= degree & degree <= 45 ~ paste(sp2,sp1),
                   degree > 45 & degree <= 90 ~ paste(sp1,sp2),
                   degree > 90 & degree <= 180 ~ paste(sp1,sp2),
                   degree > 180 & degree <= 225 ~ paste(sp1,sp2),
                   degree > 225 & degree <= 270 ~ paste(sp2,sp1),
                   degree > 270 & degree <= 360 ~ paste(sp2,sp1)),
             coord = case_when(
                   0 <= degree & degree <= 45 ~ paste(chg_sp2,chg_sp1),
                   degree > 45 & degree <= 90 ~ paste(chg_sp1,chg_sp2),
                   degree > 90 & degree <= 180 ~ paste(chg_sp1,chg_sp2),
                   degree > 180 & degree <= 225 ~ paste(chg_sp1,chg_sp2),
                   degree > 225 & degree <= 270 ~ paste(chg_sp2,chg_sp1),
                   degree > 270 & degree <= 360 ~ paste(chg_sp2,chg_sp1)),
             relation = case_when(
                   sp1 == "meth85" & sp2 == "meth92" ~ "Intra-Methylobacteriaceae",
                   sp1 == "meth85" & sp2 == "mr01" ~ "Intra-Methylobacteriaceae",
                   sp1 == "meth92" & sp2 == "mr01" ~ "Intra-Methylobacteriaceae",
                   sp1 == "smfr1" & sp2 == "spfa2" ~ "Intra-Sphingomonadaceae",
                   sp2 == "smfr1" & sp1 == "spfa2" ~ "Intra-Sphingomonadaceae",
                   TRUE ~ "Inter-Family"),
             rel2 = case_when(
                   relation == "Intra-Methylobacteriaceae" ~ "Intra-Family",
                   relation == "Intra-Sphingomonadaceae" ~ "Intra-Family",
                   TRUE ~ "Inter-Family"),
             met_group = case_when(
                   sp1 == "meth85" & sp2 == "meth92" ~ "Inter-Clade",
                   sp2 == "meth85" & sp1 == "meth92" ~ "Inter-Clade",
                   sp1 == "meth85" & sp2 == "mr01" ~ "Inter-Clade",
                   sp2 == "meth85" & sp1 == "mr01" ~ "Inter-Clade",
                   sp1 == "smfr1" & sp2 == "meth92" ~ "Inter-Clade",
                   sp2 == "smfr1" & sp1 == "meth92" ~ "Inter-Clade",
                   sp1 == "spfa2" & sp2 == "meth92" ~ "Inter-Clade",
                   sp2 == "spfa2" & sp1 == "meth92" ~ "Inter-Clade",
                   sp1 == "spfa2" & sp2 == "mr01" ~ "Inter-Clade",
                   sp2 == "spfa2" & sp1 == "mr01" ~ "Inter-Clade",
                   sp1 == "smfr1" & sp2 == "mr01" ~ "Inter-Clade",
                   sp2 == "smfr1" & sp1 == "mr01" ~ "Inter-Clade",
                   sp1 == "meth92" & sp2 == "mr01" ~ "Clade2",
                   sp2 == "meth92" & sp1 == "mr01" ~ "Clade2",
                   sp1 == "smfr1" & sp2 == "spfa2" ~ "Clade1",
                   sp2 == "smfr1" & sp1 == "spfa2" ~ "Clade1",
                   sp1 == "smfr1" & sp2 == "meth85" ~ "Clade1",
                   sp2 == "smfr1" & sp1 == "meth85" ~ "Clade1",
                   sp1 == "meth85" & sp2 == "spfa2" ~ "Clade1",
                   sp2 == "meth85" & sp1 == "spfa2" ~ "Clade1",
                   TRUE ~ "Conspecific"))

all_cart$coord_sp2 = all_cart$coord_sp
all_cart = separate(all_cart, coord_sp, into=c("s1","s2"), sep=" ") %>% 
      separate(., coord, into = c("a21", "a12"), sep = " ")
all_cart$a21 = as.numeric(all_cart$a21)
all_cart$a12 = as.numeric(all_cart$a12)

all_cart = all_cart %>% mutate(radian2 = atan2(a12, a21),
                               rad3 = case_when(
                                     radian2 < 0 ~ radian2 + 2*pi,
                                     TRUE ~ radian2),
                               degree_corr = rad3 * 180/pi)

int_str = all_cart %>% 
      dplyr::select(mix, resource, a21, a12, degree = degree_corr, magnitude, relation) %>% 
      filter(degree > 0) %>% 
      mutate(relation = case_when(
            relation == "Intra-Methylobacteriaceae" ~ "Intra-Family",
            relation == "Intra-Sphingomonadaceae" ~ "Intra-Family",
            TRUE ~ relation),
            int = (degree - 45)/(225 - 45))
int_str$mix = factor(int_str$mix)

  # summary
allSummary = all %>% 
      group_by(mix, pair, resource, sp1, sp2) %>% 
      summarise(median1 = median(chg_sp1),
                median2 = median(chg_sp2),
                q25_1 = quantile(chg_sp1, probs=0.25),
                q75_1 = quantile(chg_sp1, probs=0.75),
                q25_2 = quantile(chg_sp2, probs=0.25),
                q75_2 = quantile(chg_sp2, probs=0.75),
                mean1 = mean(chg_sp1),
                mean2 = mean(chg_sp2),
                sd1 = sd(chg_sp1),
                sd2 = sd(chg_sp2)) %>% 
      separate(mix, into=c("strain1", "strain2"), sep="_")
all_cart_summ_res = all_cart %>% 
      group_by(resource) %>% 
      summarise(med_a21 = median(a21),
                med_a12 = median(a12),
                q25_a21 = quantile(a21, probs=0.25),
                q75_a21 = quantile(a21, probs=0.75),
                q25_a12 = quantile(a12, probs=0.25),
                q75_a12 = quantile(a12, probs=0.75),
                magnitude = sqrt(med_a21^2+med_a12^2),
                radian = atan2(med_a12, med_a21),
                degree = radian * 180/pi,
                mean_a21 = mean(a21),
                mean_a12 = mean(a12),
                sd_a21 = sd(a21),
                sd_a12 = sd(a12))
all_cart_summ = all_cart %>% 
      group_by(rel2) %>% 
      summarise(med_a21 = median(a21),
                med_a12 = median(a12),
                q25_a21 = quantile(a21, probs=0.25),
                q75_a21 = quantile(a21, probs=0.75),
                q25_a12 = quantile(a12, probs=0.25),
                q75_a12 = quantile(a12, probs=0.75),
                magnitude = sqrt(med_a21^2+med_a12^2),
                radian = atan2(med_a12, med_a21),
                degree = radian * 180/pi,
                mean_a21 = mean(a21),
                mean_a12 = mean(a12),
                sd_a21 = sd(a21),
                sd_a12 = sd(a12))

all_cart_summ_intra2 = all_cart %>% filter(!grepl('CON', pair) & !grepl('MC', pair)) %>% 
      ungroup() %>% 
      summarise(mean_a21 = mean(a21),
                mean_a12 = mean(a12),
                sd_a21 = sd(a21),
                sd_a12 = sd(a12),
                magnitude = sqrt(mean_a21^2+mean_a12^2),
                radian = atan2(mean_a12, mean_a21),
                degree = radian * 180/pi,
                int = (degree - 45)/(225 - 45))
t(all_cart_summ)
all$mix = as.character(all$mix)


# Plots
      #     interaction plot by resource
p_interaction_resource = ggplot()+
      #facet_wrap(~met_group)+
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)+
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 1)+
      geom_point(data=all_cart, aes(x=a21, y=a12, color=resource), alpha=0.7, size=2)+
      geom_point(data=all_cart_summ_res, aes(x=mean_a21, y=mean_a12))+ #mean
      geom_errorbar(data=all_cart_summ_res, 
                    aes(x=mean_a21, 
                        y=mean_a12, 
                        ymin=mean_a12-sd_a12, 
                        ymax=mean_a12+sd_a12,
                        group = resource), 
                    width = 0.05)+ # error bar mean
      geom_errorbarh(data=all_cart_summ_res, 
                     aes(x=mean_a21, 
                         y=mean_a12, 
                         xmin=mean_a21-sd_a21, 
                         xmax=mean_a21+sd_a21,
                         group = resource), 
                     height = 0.05)+ # error bar mean
      geom_point(data=all_cart_summ_res, aes(x=mean_a21, y=mean_a12, fill = resource), pch = 21, size = 2.5)+ #mean
      theme_bw()+
      guides(fill = "none")+
      labs(title = 'Species interactions\nhomogeneous environment\nFold change intraspecific competition',
           x = "Fold change sp. 1",
           y = "Fold change sp. 2")+
      coord_fixed(ratio = 1)+
      theme(text = element_text(size = 18))+
      scale_color_discrete(name = "Resource", labels = c("Fru", "Fuma", "Gluc", "Glu", "Meth"))+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-3.1,3.1))+
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-3.1,3.1))
p_interaction_resource

      #     interaction plot by phylogenetic relationship
p_interaction_phyl = ggplot()+
      #facet_wrap(~met_group)+
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)+
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 1)+
      geom_point(data=all_cart, aes(x=a21, y=a12, color=rel2), alpha=0.7, size=2.5)+
      geom_errorbar(data=all_cart_summ, 
                    aes(x=mean_a21, 
                        y=mean_a12, 
                        ymin=mean_a12-sd_a12, 
                        ymax=mean_a12+sd_a12, group = rel2), width = 0.2)+ # vertical error bars
      geom_errorbarh(data=all_cart_summ, 
                     aes(x = mean_a21, 
                         y = mean_a12, 
                         xmin = mean_a21-sd_a21, 
                         xmax = mean_a21+sd_a21, group = rel2), height = 0.2)+ # horizontal error bars
      geom_point(data=all_cart_summ, aes(x=mean_a21, y=mean_a12, fill = rel2), pch = 21, size = 2.5)+ #mean
      geom_point(data=all_cart_summ_intra2, aes(x=mean_a21, y=mean_a12), color = "black", shape = 18, size = 5, alpha = 0.8)+
      geom_errorbar(data=all_cart_summ_intra2, 
                    aes(x=mean_a21, 
                        y=mean_a12, 
                        ymin=mean_a12-sd_a12, 
                        ymax=mean_a12+sd_a12), width = 0.2)+
      geom_errorbarh(data=all_cart_summ_intra2, 
                     aes(x = mean_a21, 
                         y = mean_a12, 
                         xmin = mean_a21-sd_a21, 
                         xmax = mean_a21+sd_a21), height = 0.2)+
      theme_bw()+
      guides(fill = "none")+
      labs(title = 'Species interactions\nhomogeneous environment\nFold change intraspecific competition',
           x = "Fold change sp. 1",
           y = "Fold change sp. 2")+
      coord_fixed(ratio = 1)+
      theme(text = element_text(size = 18))+
      scale_color_discrete(name = "Relation", labels = c("Inter", "Intra"))+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-3.1,3.1))+
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(-3.1,3.1))
p_interaction_phyl

      #     plot pairwise comparison
p_int_sp = ggplot(data=allSummary)+
      facet_grid(rows=vars(sp2), cols=vars(sp1), switch = "both")+
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "dashed")+
      geom_vline(xintercept = 0, alpha = 0.5, linetype = "dashed")+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5)+
      stat_ellipse(data=all, aes(x = chg_sp1, y = chg_sp2), geom = "polygon", alpha = 0.1, na.rm = T)+
      geom_point(data=all, aes(x=chg_sp1, y=chg_sp2, fill=resource), pch=21, alpha=0.7, size=1.5)+
      #geom_point(aes(x=mean1, y=mean2, fill=resource), pch=21, size=1.5, alpha=0.8)+
      #geom_errorbar(aes(y = mean2, x = mean1, ymin = mean2-sd2, ymax = mean2+sd2), 
      #              width=.2, alpha = 0.5, size=0.4)+
      #geom_errorbarh(aes(x = mean1, y=mean2, xmin = mean1-sd1, xmax = mean1+sd1), 
      #               height=.2, alpha = 0.5, size=0.4)+
      theme_bw()+
      labs(title = 'Species interactions\nhomogeneous environment\nFold change intraspecific competition',
           x = "Fold change sp. 1",
           y = "Fold change sp. 2")+
      scale_x_continuous(limits = c(-4,4), breaks = c(-3, 0, 3))+
      scale_y_continuous(limits = c(-4,4), breaks = c(-3, 0, 3))+
      coord_fixed(ratio = 1)+
      theme(text = element_text(size = 18))+
      scale_colour_discrete(name = "Resource", 
                          labels = c("Fructose", "Fumarate", "Glucose", "Glutamate", "Methanol"))
p_int_sp

      #     TYPES OF INTERACTIONS
int_str_intra1 = ggplot(int_str, aes(x = relation, y = int, fill = relation))+
      geom_hline(yintercept = 0.25, linetype = "dashed", alpha = 0.5)+
      geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5)+
      geom_hline(yintercept = 0.75, linetype = "dashed", alpha = 0.5)+
      geom_boxplot(alpha = 0.5)+
      geom_point(alpha = 0.5)+
      ylim(0,1)+
      theme_bw()+
      labs(y = "Interaction type",
           x = "Relation")+
      coord_fixed(ratio = 2)+
      scale_fill_discrete(name = "Relation")+
      theme(text = element_text(size = 15))
int_str_intra1

int_str_intra2 = ggplot(int_str, aes(x = relation, y = int, fill = relation))+
      facet_wrap(~ resource, ncol = 5)+
      geom_hline(yintercept = 0.25, linetype = "dashed", alpha = 0.5)+
      geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5)+
      geom_hline(yintercept = 0.75, linetype = "dashed", alpha = 0.5)+
      geom_boxplot(alpha = 0.5)+
      geom_point(alpha = 0.5)+
      ylim(0,1)+
      theme_bw()+
      labs(y = "Interaction type",
           x = "Relation")+
      coord_fixed(ratio = 2)+
      scale_fill_discrete(name = "Relation")+
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 45, hjust = 1))
int_str_intra2

int_str_intra3 = ggplot(int_str, aes(x = resource, y = int))+
      geom_hline(yintercept = 0.25, linetype = "dashed", alpha = 0.5)+
      geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5)+
      geom_hline(yintercept = 0.75, linetype = "dashed", alpha = 0.5)+
      geom_boxplot(alpha = 0.8)+
      geom_point(alpha = 0.5)+
      ylim(0,1)+
      theme_bw()+
      labs(y = "Interaction type",
           x = "Resource")+
      coord_fixed(ratio = 5)+
      theme(text = element_text(size = 15))
int_str_intra3

int_density = ggplot(int_str)+
      #facet_wrap(~relation)+
      geom_density(aes(x=int, color = relation), size=1.5)+
      theme_bw()+
      theme(text = element_text(size = 18))+
      labs(y = "Kernel density",
           x = "Interaction type")+
      guides(colour = "none")+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,3))+
      scale_x_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1.0) )
int_density

      #     STRENGTH OF INTERACTIONS
str_intra1 = ggplot(int_str, aes(x = relation, y = magnitude, fill = relation))+
      geom_boxplot(alpha = 0.5)+
      geom_point(alpha = 0.5)+
      theme_bw()+
      labs(y = "Interaction strength",
           x = "Relation")+
      coord_fixed(ratio = 1/3)+
      scale_fill_discrete(name = "Relation")+
      theme(text = element_text(size = 15))+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01), limits = c(0,6))
str_intra1

str_intra2 = ggplot(int_str, aes(x = relation, y = magnitude, fill = relation))+
      facet_wrap(~ resource, ncol = 5)+
      geom_boxplot(alpha = 0.5)+
      geom_point(alpha = 0.5)+
      theme_bw()+
      labs(y = "Interaction strength",
           x = "Relation")+
      coord_fixed(ratio = 1/3)+
      scale_fill_discrete(name = "Relation")+
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01), limits = c(0,6))
str_intra2

str_intra3 = ggplot(int_str, aes(x = resource, y = magnitude))+
      geom_boxplot(alpha = 0.8)+
      geom_point(alpha = 0.5)+
      theme_bw()+
      labs(y = "Interaction strength",
           x = "Resource")+
      coord_fixed(ratio = 5/6)+
      theme(text = element_text(size = 15))+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01), limits = c(0,6))
str_intra3

strength_density = ggplot(int_str)+
      #facet_wrap(~relation)+
      geom_density(aes(x=magnitude, color = relation), size=1.5)+
      theme_bw()+
      theme(text = element_text(size = 18))+
      labs(y = "Kernel density",
           x = "Interaction strength")+
      guides(colour = "none")+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,1.2))+
      scale_x_continuous(limits = c(0,6))
strength_density



#  MODELS
#     data
int = read.csv("matrix_int_str.csv", header = T)

# phylogenetic distances
phyl = read.csv("~/Google Drive/Rudolf_MRE/results/phylogeny/kbase/amphora_genes/phyl_distance.csv", header = T) %>% 
      .[,c(2,3)] %>% 
      mutate(community = case_when(
            community == "smfr1_meth92" ~ "meth92_smfr1",
            TRUE ~ as.character(community)
      ))
phyl$community = factor(phyl$community)
colnames(phyl) = c("mix", "phyl")

# metabolic distances 
metab = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/csource_growth/5resources/distmatrix_biomass_final.csv", header = T, check.names = T)

#     Re-labeling
sp.label = c("1" = "meth85",
             "2" = "meth92",
             "3" = "mr01",
             "4" = "smfr1",
             "5" = "spfa2") 
sp.label2 = data.frame(id1 = c(1:5),
                       id2 = c(1:5),
                       label = sp.label)
colnames(metab) = c(1:5)

#     Convert matrices into data frames for dataset 1 and 2
df1 <- melt(as.matrix(metab), varnames = c("id1", "id2"))

#     Re-label dataframes
df_met = inner_join(df1, sp.label2[,c(1,3)], by = "id1") %>% 
      inner_join(., sp.label2[,c(2,3)], by = "id2") %>% 
      dplyr::select(met_dist = value, sp1 = label.x, sp2 = label.y) %>% 
      unite("mix", sp1:sp2, sep="_")

phyl_met = inner_join(phyl, df_met, by = "mix")

int2 = inner_join(int, phyl_met, by="mix") %>% 
      unite("exp_group", exp:group, sep = "_")
int3 = left_join(int2, unique(int_str[,c(1,7)]), by = "mix") %>% 
      mutate(dummy_relation = case_when(
                  relation == "Intra-Family" ~ 0,
                  relation == "Inter-Family" ~ 1),
            dummy_resource = case_when(
                  resource == "fructose" ~ 0,
                  resource == "fumarate" ~ 1,
                  resource == "glucose" ~ 2,
                  resource == "glutamate" ~ 3,
                  resource == "methanol" ~ 4),
            scale_int = scale(int),
            scale_mag = scale(magnitude),
            scale_phyl = scale(phyl),
            scale_met = scale(met_dist))
head(int3)

int4 = int3 %>% 
      filter(relation == "Inter-Family") %>% 
      mutate(methylo = case_when(
            sp1 == "smfr1" ~ as.character(sp2),
            TRUE ~ as.character(sp1)),
            sphingo = case_when(
                  sp2 == "meth92" ~ as.character(sp1),
                  sp2 == "mr01" ~ as.character(sp1),
                  TRUE ~ as.character(sp2)))

#factors
int3$exp_group = factor(int3$exp_group)
int3$mix = factor(int3$mix)
int3$pair = factor(int3$pair)
int4$methylo = factor(int4$methylo)
int4$sphingo = factor(int4$sphingo)

# logistic regression (beta distribution) : interaction type
be0 = gamlss(int ~ 1 + re(random = ~1|exp_group/mix), data = int3, family = BE)
summary(be0)
exp(coef(be0))
Rsq(be0)
confint(be0)
be0$df.fit
qqnorm(resid(be0))

be1 = gamlss(int ~ dummy_relation + re(random = ~1|exp_group/mix), data = int3, family = BE)
summary(be1)
exp(coef(be1))
Rsq(be1)
confint(be1)
be1$df.fit
qqnorm(resid(be1))

be2 = gamlss(int ~ dummy_resource + re(random = ~1|exp_group/mix), data = int3, family = BE)
summary(be2)
exp(coef(be2))
Rsq(be2)
confint(be2)
be2$df.fit
qqnorm(resid(be2))

be3 = gamlss(int ~ dummy_relation + dummy_resource + re(random = ~1|exp_group/mix), data = int3, family = BE)
summary(be3)
qqnorm(resid(be3))

be4 = gamlss(int ~ dummy_relation:dummy_resource + re(random = ~1|exp_group/mix), data = int3, family = BE)
summary(be4)
qqnorm(resid(be4))


beta_fru = gamlss(int ~ dummy_relation + re(random = ~1|exp_group), data = int3[int3$resource=="fructose",], family=BE)
beta_fum = gamlss(int ~ dummy_relation + re(random = ~1|exp_group), data = int3[int3$resource=="fumarate",], family=BE)
beta_glc = gamlss(int ~ dummy_relation + re(random = ~1|exp_group), data = int3[int3$resource=="glucose",], family=BE)
beta_glu = gamlss(int ~ dummy_relation + re(random = ~1|exp_group), data = int3[int3$resource=="glutamate",], family=BE)
beta_met = gamlss(int ~ dummy_relation + re(random = ~1|exp_group), data = int3[int3$resource=="methanol",], family=BE)

summary(beta_fru)
summary(beta_fum)
summary(beta_glc)
summary(beta_glu)
summary(beta_met)

      # model comparison
AICc(be0, be1, be2, be3, be4)

# linear models : interaction strength 
ga0 = gamlss(magnitude ~ 1 + re(random = ~1|exp_group/mix), data = int3, family = GA)
summary(ga0)
exp(coef(ga0))
Rsq(ga0)
confint(ga0)
ga0$df.fit
plot(ga0)

ga1 = gamlss(magnitude ~ dummy_relation + re(random = ~1|exp_group/mix), data = int3, family = GA)
summary(ga1)
exp(coef(ga1))
Rsq(ga1)
confint(ga1)
ga1$df.fit
plot(ga1)

ga2 = gamlss(magnitude ~ dummy_resource + re(random = ~1|exp_group), data = int3, family = GA)
summary(ga2)
exp(coef(ga2))
Rsq(ga2)
confint(ga2)
ga2$df.fit
plot(ga2)

ga3 = gamlss(magnitude ~ dummy_relation + dummy_resource + re(random = ~1|exp_group), data = int3, family = GA)
summary(ga3)
ga4 = gamlss(magnitude ~ dummy_relation*dummy_resource + re(random = ~1|exp_group), data = int3, family = GA)
summary(ga4)

      # model comparison
AICc(ga0, ga1, ga2, ga3, ga4)



      # INTERACTION TYPES VS PHYLOGENETIC + METABOLIC DISTANCE
bpm1 = gamlss(int ~ phyl + re(random = ~1|exp_group/resource), data = int3, family = BE)
summary(bpm1)

bpm2 = gamlss(int ~ met_dist + re(random = ~1|exp_group/resource), data = int3, family = BE)
summary(bpm2)

bpm3 = gamlss(int ~ phyl + met_dist + re(random = ~1|exp_group/resource), data = int3, family = BE)
summary(bpm3)
exp(coef(bpm3))
confint(bpm3)
Rsq(bpm3)
plot(resid(bpm3))
plot(bpm3)

bpm4 = gamlss(int ~ phyl * met_dist + re(random = ~1|exp_group/resource), data = int3, family = BE)
summary(bpm4)
exp(coef(bpm4))
confint(bpm4)
Rsq(bpm4)
plot(resid(bpm4))
plot(bpm4)

AICc(bpm1, bpm2, bpm3, bpm4)

ggplot(int3, aes(x=phyl, y = int, color = met_dist))+
      facet_wrap(~resource, ncol = 5)+
      geom_point()


int3$pred = predict(bpm3, int3, what = "mu")

int3 = int3 %>% 
      mutate(met = case_when(
            met_dist < 0.8 ~ "met1",
            met_dist >= 0.8 ~ "met2"
      ))

pred_int1 = ggplot(int3, aes(x=phyl, y=pred))+
      geom_point(alpha = 0.8, color = "grey")+
      geom_smooth(aes(y=pred), method = "lm", se = F, fullrange = T, size = 1.5, color = "black")+
      scale_y_continuous(name = "Log (odds)", limits = c(-2, 2))+
      scale_x_continuous(name = "Phylogenetic distance", limits = c(0, 1.5))+
      scale_color_brewer(name = "Resource\noverlap", palette = "Set2", label = c("High", "Low"))+
      theme_bw()+
      coord_fixed(ratio = 1.5/4)+
      theme(text = element_text(size = 15))
      
pred_int2 = ggplot(int3, aes(x=met_dist, y=pred, color = relation))+
      geom_point()+
      geom_smooth(method = "lm", se = F, fullrange = T, size = 1.5)+
      scale_y_continuous(name = "Log (odds)", limits = c(-2, 2))+
      scale_x_continuous(name = "Metabolic distance", limits = c(0, 1.5))+
      scale_color_discrete(name = "Relation")+
      theme_bw()+
      theme(text = element_text(size = 15))

pred_int1 + pred_int2

# INTERACTION STRENGTHS VS PHYLOGENETIC + METABOLIC DISTANCE
sbpm1 = gamlss(magnitude ~ phyl + re(random = ~1|exp_group/resource), data = int3)
summary(sbpm1)

sbpm2 = gamlss(magnitude ~ met_dist + re(random = ~1|exp_group/resource), data = int3)
summary(sbpm2)

sbpm3 = gamlss(magnitude ~ phyl + met_dist + re(random = ~1|exp_group/resource), data = int3)
summary(sbpm3)

sbpm4 = gamlss(magnitude ~ phyl * met_dist + re(random = ~1|exp_group/resource), data = int3)
summary(sbpm4)
exp(coef(sbpm4))
confint(sbpm4)
Rsq(sbpm4)

AICc(sbpm1, sbpm2, sbpm3, sbpm4)


# INTERACTION TYPES SPHINGOMONADS WITH METHYLOS
ggplot(int4, aes(x=sphingo, y=int))+
      facet_grid(methylo ~ resource) + 
      geom_boxplot() +
      geom_point(aes(color = methylo))

sphi0 = gamlss(int ~ 1 + re(random = ~1|exp_group/resource/methylo), data = int4, family = BE)
summary(sphi0)

sphi1 = gamlss(int ~ sphingo + re(random = ~1|exp_group/resource/methylo), data = int4, family = BE)
summary(sphi1)
exp(coef(sphi1))
confint(sphi1)
Rsq(sphi1)
plot(resid(sphi1))
qqnorm(resid(sphi1))

#### EXPORT DATA
write.csv(all_cart, "~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/foldchange_intra.csv")
write.csv(int_str, "~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/interaction_intra.csv")

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_resource.pdf", width = 6, height = 6, useDingbats=FALSE)
p_interaction_resource
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_phyl.pdf", width = 6, height = 6, useDingbats=FALSE)
p_interaction_phyl
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_pairwise.pdf", width = 12, height = 12, useDingbats=FALSE)
p_int_sp
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_int1.pdf", width = 5, height = 5, useDingbats=FALSE)
int_str_intra1
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_int2.pdf", width = 10, height = 5, useDingbats=FALSE)
int_str_intra2
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_int3.pdf", width = 10, height = 5, useDingbats=FALSE)
int_str_intra3
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_density.pdf", width = 5, height = 5, useDingbats=FALSE)
int_density
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_str1.pdf", width = 5, height = 5, useDingbats=FALSE)
str_intra1
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_str2.pdf", width = 10, height = 5, useDingbats=FALSE)
str_intra2
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_str3.pdf", width = 10, height = 5, useDingbats=FALSE)
str_intra3
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_intra_str_density.pdf", width = 5, height = 5, useDingbats=FALSE)
strength_density
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/plot_int_pred.pdf", width = 6, height = 4, useDingbats=FALSE)
pred_int1
dev.off()

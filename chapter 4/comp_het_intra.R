####### CFU ANALYSIS 2

##    Libraries
library(tidyverse) #ggplot2, tibble, tidyr, readr, purrr, dplyr, stringr, forcats
library(lme4)
library(lmtest)
library(emmeans)
library(robustHD)
library(glmmTMB)
library(gamlss)

### FUNCTIONS
#FOR INTRA
coalesce_all_columns <- function(df) {
      tibble(
            community = df$community[1],
            species = df$species[1],
            dpi = df$dpi[1:2],
            strain_1 = na.omit(df$strain_1),
            strain_2 =na.omit(df$strain_2),
            chg_sp_1 = na.omit(df$chg_sp_1),
            chg_sp_2 = na.omit(df$chg_sp_2))}

range01 <- function(x){(x-min(x))/(max(x)-min(x))} # rescaling function (min-max rescaling)


setwd("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/")
cfu <- read.csv("sp1sp2.csv", header = T) %>% 
      mutate(s1 = sp1,
             s2 = sp2,
             inoc1log = log10(inoc1),
             inoc2log = log10(inoc2),
             cfu1log = log10(cfu1),
             cfu2log = log10(cfu2)) %>% 
      separate(s1, into = c("s1", "tag1"), sep = "_") %>% 
      separate(s2, into = c("s2", "tag2"), sep = "_") %>% 
      mutate(s1 = case_when(
                  community == "smfr1_meth92" ~ "meth92",
                  community == "smfr1_spfa2" ~ "smfr1",
                  TRUE ~ s1),
            s2 = case_when(
                  community == "smfr1_meth92" ~ "smfr1",
                  community == "smfr1_spfa2" ~ "spfa2",
                  TRUE ~ s2)) %>% 
      dplyr::select(community, type, exp_date, dpi, s1, s2, inoc1log, inoc2log, cfu1, cfu2, cfu1log, cfu2log)
cfu$community = factor(cfu$community)
cfu$exp_date = factor(cfu$exp_date)
cfu$dpi = factor(cfu$dpi)
head(cfu)

plot_cfu = ggplot(cfu, aes(x=cfu1log, y=cfu2log, fill=dpi))+
      facet_grid(s2 ~ s1, switch = "both")+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.8)+
      geom_point(alpha=0.5, size=2, pch=21)+
      theme_bw()+
      labs(title = 'Species interactions in\nheterogeneous environment',
           x = "Colony counts\nsp. 1 (CFU/gFW)",
           y = "Colony counts\nsp. 2 (CFU/gFW)")+
      scale_x_continuous(limits = c(4,10), breaks = c(5, 7, 9))+
      scale_y_continuous(limits = c(4,10), breaks = c(5, 7, 9))+
      coord_fixed(ratio = 1)+
      theme(text = element_text(size = 18))+
      scale_fill_discrete(name = "Days post\ninoculation",
                          labels = c("7 dpi", "14 dpi"))


###   COMPARED TO INTRASPECIFIC
cfu2 <- read.csv("cfu_bincomm.csv", header = T) %>% 
      mutate(community = case_when(
            community == "smfr1_meth92" & sp =="meth92_1" ~ "smfr1_meth92_1",
            community == "smfr1_meth92" & sp =="smfr1_2" ~ "smfr1_meth92_1",
            community == "smfr1_meth92" & sp =="meth92_2" ~ "smfr1_meth92_2",
            community == "smfr1_meth92" & sp =="smfr1_1" ~ "smfr1_meth92_2",
            TRUE ~ as.character(community))) %>% 
      filter(type != "con") %>% 
      group_by(community, sp, dpi) %>% 
      summarise(mean_cfu = mean(cfu),
                log_cfu = log10(mean_cfu))

cfu_intra <- read.csv("cfu_bincomm.csv", header = T) %>% 
      filter(type == "con") %>% 
      group_by(sp, dpi) %>% 
      summarise(mean_intra = mean(cfu),
                log_intra = log10(mean_intra))

intra = inner_join(cfu2, cfu_intra, by = c("sp", "dpi")) %>% 
      mutate(species = sp,
             chg_sp = log2(mean_cfu/mean_intra)) %>% 
      separate(sp, into = c("strain", "tag"), sep = "_")
intra$tag = factor(intra$tag)

comp_intra = pivot_wider(intra, id_cols = c(community, tag, dpi, species), 
                         names_from = tag, values_from = c(strain, chg_sp)) %>% 
      do(coalesce_all_columns(.)) %>% 
      dplyr::select(community:dpi, s1 = strain_1, s2 = strain_2, chg_sp1 = chg_sp_1, chg_sp2 = chg_sp_2)
head(comp_intra)

all_cart_intra = comp_intra %>%
      mutate(magnitude = sqrt(chg_sp1^2+chg_sp2^2),
             radian = atan2(chg_sp2, chg_sp1),
             rad2 = case_when(
                   radian < 0 ~ radian + 2*pi,
                   TRUE ~ radian),
             degree = rad2 * 180/pi,
             coord_sp = case_when(
                   0 <= degree & degree <= 45 ~ paste(s2,s1),
                   degree > 45 & degree <= 90 ~ paste(s1,s2),
                   degree > 90 & degree <= 180 ~ paste(s1,s2),
                   degree > 180 & degree <= 225 ~ paste(s1,s2),
                   degree > 225 & degree <= 270 ~ paste(s2,s1),
                   degree > 270 & degree <= 360 ~ paste(s2,s1)),
             coord = case_when(
                   0 <= degree & degree <= 45 ~ paste(chg_sp2,chg_sp1),
                   degree > 45 & degree <= 90 ~ paste(chg_sp1,chg_sp2),
                   degree > 90 & degree <= 180 ~ paste(chg_sp1,chg_sp2),
                   degree > 180 & degree <= 225 ~ paste(chg_sp1,chg_sp2),
                   degree > 225 & degree <= 270 ~ paste(chg_sp2,chg_sp1),
                   degree > 270 & degree <= 360 ~ paste(chg_sp2,chg_sp1)),
             relation = case_when(
                   s1 == "meth85" & s2 == "meth92" ~ "Intra-Methylobacteriaceae",
                   s1 == "meth85" & s2 == "mr01" ~ "Intra-Methylobacteriaceae",
                   s1 == "meth92" & s2 == "mr01" ~ "Intra-Methylobacteriaceae",
                   s1 == "smfr1" & s2 == "spfa2" ~ "Intra-Sphingomonadaceae",
                   s1 == "spfa2" & s2 == "smfr1" ~ "Intra-Sphingomonadaceae",
                   TRUE ~ "Inter-Family"),
             rel2 = case_when(
                   relation == "Intra-Methylobacteriaceae" ~ "Intra-Family",
                   relation == "Intra-Sphingomonadaceae" ~ "Intra-Family",
                   TRUE ~ "Inter-Family"))
all_cart_intra$coord_sp2 = all_cart_intra$coord_sp
all_cart_intra = separate(all_cart_intra, coord_sp, into=c("s1","s2"), sep=" ") %>% 
      separate(., coord, into = c("a21", "a12"), sep = " ")
all_cart_intra$a21 = as.numeric(all_cart_intra$a21)
all_cart_intra$a12 = as.numeric(all_cart_intra$a12)
all_cart_intra = all_cart_intra %>% mutate(radian2 = atan2(a12, a21),
                               rad3 = case_when(
                                     radian2 < 0 ~ radian2 + 2*pi,
                                     TRUE ~ radian2),
                               degree_corr = rad3 * 180/pi)

min(all_cart_intra$degree_corr)
max(all_cart_intra$degree_corr)

all_cart_summ_intra = all_cart_intra %>% ungroup() %>% 
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

all_cart_summ_intra2 = all_cart_intra %>% 
      ungroup() %>% 
      summarise(mean_a21 = mean(a21),
                mean_a12 = mean(a12),
                sd_a21 = sd(a21),
                sd_a12 = sd(a12),
                magnitude = sqrt(mean_a21^2+mean_a12^2),
                radian = atan2(mean_a12, mean_a21),
                degree = radian * 180/pi,
                int = (degree - 45)/(225 - 45))


int_str_intra = all_cart_intra %>% 
      ungroup() %>% 
      dplyr::select(community, dpi, a21, a12, degree = degree_corr, magnitude, relation) %>% 
      mutate(relation = case_when(
                  relation == "Intra-Methylobacteriaceae" ~ "Intra-Family",
                  relation == "Intra-Sphingomonadaceae" ~ "Intra-Family",
                  TRUE ~ relation),
            community = case_when(
                  community == "smfr1_meth92_1" ~ "smfr1_meth92",
                  community == "smfr1_meth92_2" ~ "smfr1_meth92",
                  TRUE ~ community),
            int = (degree - 45)/(225 - 45))
int_str_intra$community = factor(int_str_intra$community)
int_str_intra$dpi = factor(int_str_intra$dpi)


## plots
plot_cfu_intra = ggplot()+
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)+
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 1)+
      geom_point(data=all_cart_intra, aes(x=a21, y=a12, color=rel2), alpha=0.7, size=2.5)+
      geom_errorbar(data=all_cart_summ_intra, 
                    aes(x=mean_a21, 
                        y=mean_a12, 
                        ymin=mean_a12-sd_a12, 
                        ymax=mean_a12+sd_a12, group = rel2), width = 0.2)+ # vertical error bars
      geom_errorbarh(data=all_cart_summ_intra, 
                     aes(x = mean_a21, 
                         y = mean_a12, 
                         xmin = mean_a21-sd_a21, 
                         xmax = mean_a21+sd_a21, group = rel2), height = 0.2)+ # horizontal error bars
      geom_point(data=all_cart_summ_intra, aes(x=mean_a21, y=mean_a12, fill = rel2), pch = 21, size = 2.5)+ #mean
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
      labs(title = 'Species interactions\nheterogeneous environment\nFold change intraspecific interaction',
           x = "Relative change sp. 1",
           y = "Relative change sp. 2")+
      scale_x_continuous(limits = c(-5.5, 5.5), breaks = c(-5, -2.5, 0, 2.5, 5))+
      scale_y_continuous(limits = c(-5.5, 5.5), breaks = c(-5, -2.5, 0, 2.5, 5))+
      coord_fixed(ratio = 1)+
      theme(text = element_text(size = 18))+
      scale_color_discrete(name = "Relation", labels = c("Inter", "Intra"))
plot_cfu_intra

plot_int_intra = ggplot(int_str_intra, aes(x = relation, y = int, fill = relation))+
      facet_wrap(~dpi)+
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
plot_int_intra

plot_str_intra = ggplot(int_str_intra, aes(x = relation, y = magnitude, fill = relation))+
      facet_wrap(~dpi)+
      geom_boxplot(alpha = 0.5)+
      geom_point(alpha = 0.5)+
      theme_bw()+
      labs(y = "Interaction strength",
           x = "Relation")+
      coord_fixed(ratio = 2/6)+
      scale_fill_discrete(name = "Relation")+
      theme(text = element_text(size = 15))+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01), limits = c(0,6))
plot_str_intra

plot_int_intra/plot_str_intra + plot_layout(guides = "collect")

plot_density_intra = ggplot(int_str_intra)+
      #facet_wrap(~relation)+
      geom_density(aes(x=int, color = relation), size=1.5)+
      theme_bw()+
      theme(text = element_text(size = 18))+
      labs(y = "Kernel density",
           x = "Interaction type")+
      guides(colour = "none")+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,3.01))
plot_density_intra

plot_str_density_intra = ggplot(int_str_intra)+
      #facet_wrap(~relation)+
      geom_density(aes(x=magnitude, color = relation), size=1.5)+
      theme_bw()+
      theme(text = element_text(size = 18))+
      labs(y = "Kernel density",
           x = "Interaction strength")+
      guides(colour = "none")+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,1.2))+
      scale_x_continuous(limits = c(0,6))
plot_str_density_intra


#  MODELS
#     data
# phylogenetic distances
phyl = read.csv("~/Google Drive/Rudolf_MRE/results/phylogeny/kbase/amphora_genes/phyl_distance.csv", header = T) %>% 
      .[,c(2,3)]
phyl$community = factor(phyl$community)
colnames(phyl) = c("community", "phyl")

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
      unite("community", sp1:sp2, sep="_")

phyl_met = inner_join(phyl, df_met, by = "community")

int3 = inner_join(int_str_intra, phyl_met, by="community") %>% 
      mutate(dummy_relation = case_when(
            relation == "Intra-Family" ~ 0,
            relation == "Inter-Family" ~ 1),
            scale_int = scale(int),
            scale_mag = scale(magnitude),
            scale_phyl = scale(phyl),
            scale_met = scale(met_dist))
head(int3)
str(int3)

#factors
int3$community = factor(int3$community)

# logistic regression (beta distribution) : interaction type
be0 = gamlss(int ~ 1, data = int3, family = BE)
summary(be0)
exp(coef(be0))
Rsq(be0)
confint(be0)
be0$df.fit
qqnorm(resid(be0))

be1 = gamlss(int ~ dummy_relation, data = int3, family = BE)
summary(be1)
exp(coef(be1))
Rsq(be1)
confint(be1)
be1$df.fit
qqnorm(resid(be1))

be2 = gamlss(int ~ dpi, data = int3, family = BE)
summary(be2)
exp(coef(be2))
Rsq(be2)
confint(be2)
be2$df.fit
qqnorm(resid(be2))

be3 = gamlss(int ~ dummy_relation * dpi, data = int3, family = BE)
summary(be3)
exp(coef(be3))
Rsq(be3)
confint(be3)
be3$df.fit
qqnorm(resid(be3))

# model comparison
AICc(be0, be1, be2, be3)

# gamma linear models : interaction strength 
ga0 = gamlss(magnitude ~ 1, data = int3, family = GA)
summary(ga0)
exp(coef(ga0))
Rsq(ga0)
confint(ga0)
ga0$df.fit
plot(ga0)

ga1 = gamlss(magnitude ~ dummy_relation, data = int3, family = GA)
summary(ga1)
exp(coef(ga1))
Rsq(ga1)
confint(ga1)
ga1$df.fit
plot(ga1)

ga2 = gamlss(magnitude ~ dpi, data = int3, family = GA)
summary(ga2)
exp(coef(ga2))
Rsq(ga2)
confint(ga2)
ga2$df.fit
plot(ga2)

ga3 = gamlss(magnitude ~ dummy_relation * dpi, data = int3, family = GA)
summary(ga3)
exp(coef(ga3))
Rsq(ga3)
confint(ga3)
ga3$df.fit
plot(ga3)

# model comparison
AICc(ga0, ga1, ga2, ga3)


# INTERACTION TYPES VS PHYLOGENETIC + METABOLIC DISTANCE
bpm1 = gamlss(int ~ phyl + dpi, data = int3, family = BE)
summary(bpm1)

bpm2 = gamlss(int ~ met_dist + dpi, data = int3, family = BE)
summary(bpm2)

bpm3 = gamlss(int ~ phyl + met_dist + dpi, data = int3, family = BE)
summary(bpm3)

bpm4 = gamlss(int ~ phyl : met_dist * dpi, data = int3, family = BE)
summary(bpm4)
exp(coef(bpm4))
plot(resid(bpm4))


#### EXPORT DATA
## CFU
write.csv(cfu, "~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/cfu.csv")

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/plot_cfu_mono.pdf", 
    width = 6, height = 6, useDingbats=FALSE)
plot_cfu
dev.off()


## FOLD CHANGE TO INTRASPECIFIC COMPETITION
write.csv(intra, "~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/foldchange_intra_long.csv")
write.csv(all_cart_intra, "~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/foldchange_intra.csv")
write.csv(int_str_intra, "~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/interaction_intra.csv")

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/plot_intra_cfu.pdf", 
    width = 6, height = 6, useDingbats=FALSE)
plot_cfu_intra
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/plot_intra_int.pdf", 
    width = 5, height = 5, useDingbats=FALSE)
plot_int_intra
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/plot_intra_density.pdf", 
    width = 5, height = 5, useDingbats=FALSE)
plot_density_intra
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/plot_intra_str.pdf", 
    width = 5, height = 5, useDingbats=FALSE)
plot_str_intra
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/plot_intra_density_str.pdf", 
    width = 5, height = 5, useDingbats=FALSE)
plot_str_density_intra
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/plot_intra_int_str.pdf", 
    width = 5, height = 5, useDingbats=FALSE)
plot_int_intra/plot_str_intra + plot_layout(guides = "collect")
dev.off()


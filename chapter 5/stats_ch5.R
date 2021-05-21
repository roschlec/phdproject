### STATS CHAPTER 4 ### PE299R COMPETITION 

# Library
library(tidyverse)
library(lme4)
library(mgcv)
library(car)
library(MuMIn)
library(ggrepel)
library(patchwork)
library(vegan)
library(boot)
library(mgcv.helper)
library(gamlss)
library(gamlss.tr)

setwd("~/Google Drive/Rudolf_MRE/results/ch4")

strain.label = c("acid84" = "AcidoL84", "aero245" = "AeromL245", "agre335" = "AgreiL335", "arth145" = "ArthrL145",
                 "brad396" = "BradyL396", "meth85" = "MethyL85", "meth92" = "MethyL92", "micr320" = "MicroL320", 
                 "pa299r" = "Pe299R", "pkore" = "PkP19E3", "pssb728a" = "PssB728a", "rhod225" = "RhodoL225",
                 "smelo" = "SmFR1", "sphi17" = "SphinL17", "sphi357" = "SphinL357", "sphyl" = "SpFA2")
names = data.frame(strain = strain.label)
names$sp1 = rownames(names)

comp <- read.csv("./competition/comp_chg.csv", header = T)[,-1:-2] %>% 
      group_by(strain) %>% 
      summarise(comp = mean(chg)) %>% 
      dplyr::select(sp1 = strain, comp)
comp2 <- read.csv("./competition/comp_chg.csv", header = T)[,-1:-2] %>% 
      group_by(strain) %>% 
      dplyr::select(sp1 = strain, comp = chg) %>% 
      filter(sp1 != "Monoculture")
met <- read.csv("./growth_assay/Cprof/all/met_dist.csv", header = T)[,-1] %>% 
      dplyr::select(sp1, sp2, met_dist)
fitness <- read.csv("./growth_assay/5xC/all/fitness5C.csv", header = T)[,-1] %>% 
      dplyr::select(sp1, fitness = value)
cfu <- read.csv("~/Google Drive/Rudolf_MRE/results/ch2/in planta growth/r_K_df.csv", header = T)[,-1] %>% 
      inner_join(., names, by = "sp1") %>% 
      dplyr::select(sp1 = strain, r, K)
phyl <- read.csv("~/Google Drive/Rudolf_MRE/THESIS/ch5:cusper/Phylogeny/phyl_distance.csv", header = T)[,-1]
rs <- read.csv("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/dist_rs.csv", header = T)[,-1]

df = inner_join(met, phyl, by = c("sp1", "sp2")) %>% 
      filter(sp2 == "Pe299R" & sp1 != "Pe299R") %>% 
      left_join(., df.color[-1,], by = "sp1") %>% 
      droplevels()
df2 = inner_join(comp2, met[met$sp2 == "Pe299R",], by = "sp1") %>% 
      left_join(., fitness, by = "sp1") %>% 
      left_join(., cfu, by = "sp1") %>%
      left_join(., phyl[phyl$sp2 == "Pe299R",], by = c("sp1", "sp2")) %>%
      left_join(., rs[rs$sp2 == "Pe299R",], by = c("sp1", "sp2")) %>%
      left_join(., df.color[-1,], by = "sp1") %>% 
      droplevels()

sp1 = data.frame(sp1 = levels(df2$sp1))
sp1 = data.frame(sp1 = levels(df$sp1))
palette = as.character(na.omit(left_join(df.color[-1,], sp1, by = "sp1")[,2]))
palette2 = as.character(na.omit(left_join(df.color[-1,], sp1, by = "sp1")[,2]))

palette3 = df2 %>% dplyr::select(sp1, palette) %>% unique() %>% data.frame()
palette3 = as.character(palette3[c(7,8,6,3,4,1,10,2,9,5),2])
levels(df2$sp1)

#     Data exploration
qqnorm(scale(df2$comp))
qqnorm(unique(df2$met_dist))
qqnorm(unique(df2$fitness))
qqnorm(unique(df2$r))
qqnorm(unique(df2$K))

#     Models
# METABOLIC DISTANCE
lm_met = lmer(comp ~ met_dist + (1|sp1), data = df2)
summary(lm_met)
qqnorm(resid(lm_met))
plot(resid(lm_met))
r.squaredGLMM(lm_met)
r.squaredLR(lm_met)[1]
confint(lm_met)
Anova(lm_met, type = "III")

poly3_met = lmer(comp ~ poly(met_dist,3) + (1|sp1), data = df2)
summary(poly3_met)
qqnorm(resid(poly3_met))
plot(resid(poly3_met))
r.squaredGLMM(poly3_met)
r.squaredLR(poly3_met)[1]
Anova(poly3_met, type = "III")

poly4_met = lmer(comp ~ poly(met_dist,4) + (1|sp1), data = df2)
summary(poly4_met)
qqnorm(resid(poly4_met))
plot(resid(poly4_met))
r.squaredGLMM(poly4_met)
r.squaredLR(poly4_met)[1]
Anova(poly4_met, type = "III")

gam_met = gam(comp ~ s(met_dist, bs="cr", k = 4) + s(sp1, bs = "re"), data = df2)
summary(gam_met)
qqnorm(resid(gam_met))
plot(resid(gam_met))
r.squaredGLMM(gam_met)
r.squaredLR(gam_met)[1]
confint.gam(gam_met)


AIC(lm_met, poly3_met, poly4_met, gam_met)

plot1 = ggplot(df2, aes(x = met_dist, y = comp))+
      #geom_smooth(method = "lm", formula = y ~ x, color = "grey", linetype = "dashed", se = FALSE)+
      geom_hline(aes(yintercept=0), linetype = "dashed")+
      geom_smooth(method = "lm", formula = y ~ x, color = "dark grey", se = FALSE)+
      geom_point(aes(fill = sp1), pch = 21, alpha = 0.8, size = 2)+
      theme_bw()+
      theme(text = element_text(size = 15))+
      geom_text_repel(data = df, aes(label = sp1), 
                      max.overlaps = Inf, 
                      box.padding = 1.5,
                      segment.size = 0.2)+
      scale_y_continuous(name = "Fold change\nPe299R::mSc (log2-scale)",
                         limits = c(-1, 1))+
      scale_x_continuous(name = "Metabolic distance",
                         limits = c(0, 75), breaks = seq(0,75,15))+
      guides(fill = "none")+
      coord_fixed(ratio = 75/2)
      #scale_fill_manual(values = palette)


# FITNESS 5XC
lm_fit = lmer(comp ~ fitness + (1|sp1), data = df2)
summary(lm_fit)
qqnorm(resid(lm_fit))
plot(resid(lm_fit))
r.squaredGLMM(lm_fit)
r.squaredLR(lm_fit)[1]
Anova(lm_fit, type = "III")

gam_fit = gam(comp ~ s(fitness, bs="cr", k = 4) + s(sp1, bs = "re"), data = df2)
summary(gam_fit)
qqnorm(resid(gam_fit))
plot(resid(gam_fit))
r.squaredGLMM(gam_fit)
r.squaredLR(gam_fit)[1]

AIC(lm_fit, gam_fit)

      
plot2 = ggplot(df2, aes(x = fitness, y = comp))+
      #geom_smooth(method = "lm", formula = y ~ x, color = "grey", linetype = "dashed", se = FALSE)+
      geom_hline(aes(yintercept=0), linetype = "dashed")+
      geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 4), color = "dark grey", se = FALSE)+
      geom_point(aes(fill = sp1), pch = 21, alpha = 0.8, size = 2)+
      theme_bw()+
      theme(text = element_text(size = 15))+
      geom_text_repel(data = df, aes(label = sp1), 
                      max.overlaps = Inf, 
                      box.padding = 1.5,
                      segment.size = 0.2)+
      scale_y_continuous(name = "Fold change\nPe299R::mSc (log2-scale)",
                         limits = c(-1, 1))+
      scale_x_continuous(name = "Growth differences",
                         limits = c(0, 15))+
      guides(fill = "none")+
      coord_fixed(ratio = 15/2)+
      scale_fill_manual(values = palette)

# PHYLOGENETIC DISTANCE
lm_phyl = lmer(comp ~ phyl_dist + (1|sp1), data = df2)
summary(lm_phyl)
qqnorm(resid(lm_phyl))
plot(resid(lm_phyl))
r.squaredGLMM(lm_phyl)
r.squaredLR(lm_phyl)[1]
confint(lm_phyl)
Anova(lm_phyl, type = "III")

poly3_phyl = lmer(comp ~ poly(phyl_dist,3) + (1|sp1), data = df2)
summary(poly3_phyl)
qqnorm(resid(poly3_phyl))
plot(resid(poly3_phyl))
r.squaredGLMM(poly3_phyl)
r.squaredLR(poly3_phyl)[1]
confint(poly3_phyl)
Anova(poly3_phyl, type = "III")

poly4_phyl = lmer(comp ~ poly(phyl_dist,4) + (1|sp1), data = df2)
summary(poly4_phyl)
qqnorm(resid(poly4_phyl))
plot(resid(poly4_phyl))
r.squaredGLMM(poly4_phyl)
r.squaredLR(poly4_phyl)[1]
Anova(poly4_phyl, type = "III")

gam_phyl = gam(comp ~ s(phyl_dist, bs="cr", k = 4) + s(sp1, bs = "re"), data = df2)
summary(gam_phyl)
qqnorm(resid(gam_phyl))
plot(resid(gam_phyl))
r.squaredGLMM(gam_phyl)
r.squaredLR(gam_phyl)[1]
confint.gam(gam_phyl)

AIC(lm_phyl, poly3_phyl, poly4_phyl, gam_phyl)

plot3 = ggplot(df2, aes(x = phyl_dist, y = comp))+
      #geom_smooth(method = "lm", formula = y ~ x, color = "grey", linetype = "dashed", se = FALSE)+
      geom_hline(aes(yintercept=0), linetype = "dashed")+
      geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 4), color = "dark grey", se = FALSE)+
      geom_point(aes(fill = sp1), pch = 21, alpha = 0.8, size = 2)+
      theme_bw()+
      theme(text = element_text(size = 15))+
      geom_text_repel(aes(label = sp1), 
                      max.overlaps = Inf, 
                      box.padding = 1.5,
                      segment.size = 0.2)+
      scale_y_continuous(name = "Fold change\nPe299R::mSc (log2-scale)",
                         limits = c(-1, 1))+
      scale_x_continuous(name = "Phylogenetic distance",
                         limits = c(0, 1.5), breaks = seq(0,1.5,.5))+
      guides(fill = "none")+
      coord_fixed(ratio = 1.5/2)+
      scale_fill_manual(values = palette)


plot1 + plot3 + plot2

## METABOLISM AND PHYLOGENY IN COMPETITION ASSAY
gam_met_phyl = gam(comp ~ s(phyl_dist, bs="cr", k = 4) + s(met_dist, bs="cr", k = 4) + s(sp1, bs = "re"), data = df2)
summary(gam_met_phyl)
qqnorm(resid(gam_met_phyl))
plot(resid(gam_met_phyl))
r.squaredGLMM(gam_met_phyl)
r.squaredLR(gam_met_phyl)[1]
confint.gam(gam_met_phyl)

lm_met_phyl = lm(comp ~ phyl_dist * met_dist, data = df2)
summary(lm_met_phyl)
Anova(lm_met_phyl, type = "III")

df3 = df2 %>% 
      mutate(relation = case_when(
                   sp1 == "Pe299R" ~ "Gammaproteobacteria",
                   sp1 == "PssB728a" ~ "Gammaproteobacteria",
                   sp1 == "PkP19E3" ~ "Gammaproteobacteria",
                   sp1 == "MethyL85" ~ "Alphaproteobacteria",
                   sp1 == "MethyL92" ~ "Alphaproteobacteria",
                   sp1 == "SmFR1" ~ "Alphaproteobacteria",
                   sp1 == "BradyL396" ~ "Alphaproteobacteria",
                   TRUE ~ "Actinobacteria"))
df3$pred = predict(lm_met_phyl, df3)


plot_met_phyl = ggplot(df3, aes(x = met_dist))+
      geom_hline(aes(yintercept=0), linetype = "dashed")+
      geom_smooth(aes(y = pred, color = relation), method = "lm", fullrange = T,
                  se = F, formula = y ~ x, alpha = 0.2, size = 1)+
      geom_point(aes(y = comp, fill = sp1), pch = 21, size = 2, color = "black")+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_y_continuous(name = "Fold change\nPe299R::mSc (log2-scale)",
                         limits = c(-1, 1))+
      scale_x_continuous(name = "Metabolic distance",
                         limits = c(0, 75), breaks = seq(0,75,15))+
      coord_fixed(ratio = 75/2)+
      scale_color_brewer(palette = "Set2")+
      scale_fill_manual(values = palette3)
plot_met_phyl

## METABOLISM AND PHYLOGENY IN REPSUC
df3 = inner_join(met, rs, by = c("sp1", "sp2")) %>% 
      inner_join(., phyl, by = c("sp1", "sp2")) %>% 
      filter(sp2 == "Pe299R") %>% 
      filter(sp1 != "Pe299R") %>% 
      mutate(met = 1/met_dist,
             relation = case_when(
                         phyl_dist < 1 ~ "close",
                         phyl_dist >= 1 ~ "dist"))

ggplot(df3, aes(x = met_dist, y = rs))+
      geom_point()
ggplot(df3, aes(x = met, y = rs))+
      geom_point()+
      geom_smooth(method = "lm", se = F)
ggplot(df3, aes(x = phyl_dist, y = rs))+
      geom_point()+
      geom_smooth(method = "lm", se = F)

gamma1 = gamlss(rs ~ met_dist, data = df3, family = GA)
gamma2 = gamlss(rs ~ phyl_dist, data = df3, family = GA)
gamma3 = gamlss(rs ~ met_dist + phyl_dist, data = df3, family = GA)
gamma4 = gamlss(rs ~ met_dist * phyl_dist, data = df3, family = GA)
summary(gamma1)
summary(gamma2)
summary(gamma3)
summary(gamma4)

AIC(gamma1, gamma2, gamma3, gamma4)

plot(resid(lm3))
cor(df3$rs, df3$met, method = "pearson")
cor(df3$rs, df3$phyl_dist, method = "pearson")
cor(df3$rs, df3$met_dist, method = "pearson")

df3$pred = predict(gamma4, df3, what = "mu", type = "response", se.fit = TRUE)$fit
df3$se = predict(gamma4, df3, what = "mu", type = "response", se.fit = TRUE)$se.fit


pred1 = data.frame(met_dist = seq(1, 100, .1),
                   phyl_dist = seq(0.02, 2, 0.002))
pred1$pred = predict(glm(rs ~ met_dist, data = df3, family=Gamma), type = "response", pred1)
pred1$pred2 = predict(glm(rs ~ phyl_dist, data = df3, family=Gamma), type = "response", pred1)
pred1 = pred1 %>% 
      mutate(relation = case_when(
            phyl_dist < 1 ~ "close",
            phyl_dist >= 1 ~ "dist"
      ))


rs_met = ggplot(pred1, aes(x = met_dist, y = pred))+
      geom_line(size = 1)+
      geom_point(data = df3, aes(x = met_dist, y = rs), color = "black")+
      theme_bw()+
      theme(
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 15)
      )+
      scale_y_continuous(name = "Reproductive success\nsimilarity",
                         limits = c(0, 0.3))+
      scale_x_continuous(name = "Metabolic distance",
                         limits = c(0, 100))+
      coord_fixed(ratio = 100/.3)
rs_phyl = ggplot(pred1, aes(x = phyl_dist, y = pred2))+
      geom_line(size = 1)+
      geom_point(data = df3, aes(x = phyl_dist, y = rs), color = "black")+
      theme_bw()+
      theme(
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 15)
      )+
      scale_y_continuous(name = "",
                         limits = c(0,0.3))+
      scale_x_continuous(name = "Phylogenetic distance",
                         limits = c(0,1.5))+
      coord_fixed(ratio = 1.5/.3)

### METABOLIC VS PHYLOGENETIC DISTANCE
cor(df$met_dist, df$phyl_dist, method = "pearson")
lm2 = lm(met_dist ~ phyl_dist, data = df)
summary(lm2)
plot(resid(lm2))

dfmet = met %>% arrange(sp2, sp1) %>% 
      pivot_wider(id_cols = sp1, names_from = sp2, values_from = met_dist)
matmet = as.matrix(dfmet[,-1])
rownames(matmet) = dfmet$sp1
dmet = dist(matmet, diag = F)
mdsmet = monoMDS(dmet)

dfphyl = phyl %>% arrange(sp2, sp1) %>% 
      pivot_wider(id_cols = sp1, names_from = sp2, values_from = phyl_dist)
matphyl = as.matrix(dfphyl[,-1])
rownames(matphyl) = dfphyl$sp1
dphyl = dist(matphyl, diag = F)
mdsphyl = monoMDS(dphyl)

## Mantel test
mantel(dmet, dphyl, permutations = 9999)

## Procrustes test
pro1 = protest(mdsmet, mdsphyl,permutations = 999)

pro <- vector(length = 1000)
for (i in 1:length(pro)){
      test <- protest(mdsmet, mdsphyl, permutations = 999)
      pro[i] <- test$signif
}
pro
quantile(pro, c(0.025, 0.975))
median(pro)
mean(pro)

## Plot
metphyl = inner_join(met, phyl, by = c("sp1", "sp2"))
met_phyl1 = ggplot(data = df, aes(x = phyl_dist, y = met_dist)) +
      geom_point(alpha = 0.5)+
      geom_text_repel(data = df, aes(label = sp1))+
      theme_bw()+
      scale_y_continuous(name = "Metabolic distance",
                         limits = c(0, 100), breaks = seq(0,100,25))+
      scale_x_continuous(name = "Phylogenetic distance",
                         limits = c(0, 1.5), breaks = seq(0,1.5,.5))+
      guides(fill = "none")+
      coord_fixed(ratio = 1.5/100)
met_phyl1

met_phyl2 = ggplot(metphyl, aes(x = phyl_dist, y = met_dist)) +
      geom_point(alpha = 0.5)+
      geom_point(data = df, aes(x = phyl_dist, y = met_dist), 
                 pch = 21, alpha = 0.8, size = 2.5, fill = "magenta")+
      geom_text_repel(data = df, aes(label = sp1), 
                      max.overlaps = Inf, 
                      box.padding = 1.5,
                      segment.size = 0.2)+
      theme_bw()+
      scale_y_continuous(name = "Metabolic distance",
                         limits = c(0, 100), breaks = seq(0,100,25))+
      scale_x_continuous(name = "Phylogenetic distance",
                         limits = c(0, 1.5), breaks = seq(0,1.5,.5))+
      guides(fill = "none")+
      coord_fixed(ratio = 1.5/100)

## SAVE PLOT
pdf("~/Google Drive/Rudolf_MRE/results/ch4/model1.pdf", width = 5, height = 5, useDingbats=FALSE)
plot1
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/model2.pdf", width = 5, height = 5, useDingbats=FALSE)
plot2
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/model3.pdf", width = 5, height = 5, useDingbats=FALSE)
plot3
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/met_phyl1.pdf", width = 5, height = 5, useDingbats=FALSE)
met_phyl1
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/met_phyl.pdf", width = 5, height = 5, useDingbats=FALSE)
met_phyl2
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/rs_met_phyl.pdf", width = 10, height = 5, useDingbats=FALSE)
rs_met + rs_phyl
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/plot_met_phyl.pdf", width = 6, height = 5, useDingbats=FALSE)
plot_met_phyl
dev.off()

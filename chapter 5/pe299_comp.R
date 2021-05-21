# Pe299R competition

library(tidyverse)
library(MESS)
library(RColorBrewer)
library(circlize)
library(patchwork)
library(ggdendro)
library(ComplexHeatmap)
library(emmeans)
library(gamlss)


setwd("~/Google Drive/Rudolf_MRE/results/ch4/competition/")

strain.labelA = c(x01 = "Monoculture", 
                 x02 = "Pe299R",
                 x03 = "PssB728a",
                 x04 = "MethyL85",
                 x05 = "MethyL92",
                 x06 = "SmFR1",
                 x07 = "BradyL396",
                 x08 = "AgreiL335",
                 x09 = "ArthrL145", 
                 x10 = "RhodoL225",
                 x11 = "PkP19E3")
strain.labelB = c(x01 = "Pe299R + Pe299R", 
                  x02 = "AgreiL335 + Pe299R",
                  x03 = "AgreiL335 + PssB728a",
                  x04 = "AgreiL335 + MethyL85",
                  x05 = "AgreiL335 + MethyL92",
                  x06 = "AgreiL335 + SmFR1",
                  x07 = "AgreiL335 + BradyL396",
                  x08 = "AgreiL335 + AgreiL335",
                  x09 = "AgreiL335 + ArthL145", 
                  x10 = "AgreiL335 + RhodL225",
                  x11 = "AgreiL335 + PkP19E3")
cond.label = c(A = "1 Competitor",
               B = "2 Competitors",
               C = "Control")

data = read.csv("comp1.csv", header = T) %>% 
      pivot_longer(cols = -1, names_to = c("cond", "sample", "rep"), values_to = "rfu", names_sep = "_") %>%
      na.omit()
head(data)

data$cond = factor(data$cond)
data$sample = factor(data$sample)
data$rep = factor(data$rep)

autofl = data %>% filter(sample == "wt") %>% 
      group_by(time) %>% 
      summarise(autofl = mean(rfu))
head(autofl)

data_corr = data %>% 
      filter(cond != "C") %>% 
      filter(cond != "B" | sample != "x11" | rep != "4") %>% 
      inner_join(., autofl, by = "time") %>% 
      mutate(rfu_corr = rfu - autofl)

plot1 = ggplot(data_corr, aes(x = time, y=rfu_corr, color = rep))+
      facet_grid(cond ~ sample)+
      geom_point()

data_auc = data_corr %>% 
      group_by(cond, sample, rep) %>% 
      summarise(auc = auc(time, rfu_corr))

plot2 = ggplot(data_auc, aes(x=sample, y=auc))+
      facet_wrap(~ cond)+
      geom_boxplot()
plot2

### one competitor
mean_intraA = mean(data_auc$auc[data_auc$cond == "A" & data_auc$sample == "x02"])

aucA = data_auc %>% 
      filter(cond == "A") %>% 
      ungroup() %>% 
      mutate(chg = log2(auc/mean_intraA),
             dummy = case_when(
                   sample == "x01" ~ 1,
                   sample == "x02" ~ 0,
                   sample == "x03" ~ 2,
                   sample == "x04" ~ 3,
                   sample == "x05" ~ 4,
                   sample == "x06" ~ 5,
                   sample == "x07" ~ 6,
                   sample == "x08" ~ 7,
                   sample == "x09" ~ 8,
                   sample == "x10" ~ 9,
                   sample == "x11" ~ 10))

names = data.frame(strain = strain.labelA)
names$sample = rownames(names)

aucA = inner_join(aucA, names, by = "sample")

aucA_summ = aucA %>% 
      group_by(sample) %>% 
      summarise(mean_fold = mean(chg),
                sd_fold = sd(chg)) %>% 
      arrange(mean_fold)

new.order = aucA_summ$sample
new.order[1] = aucA_summ$sample[10]
new.order[2:10] = aucA_summ$sample[1:9]
aucA$sample = factor(aucA$sample, levels = new.order)

palette = vector(length = 11)
palette[2:11] = brewer.pal(10, "Set3")
palette[1] = "black"

strain = data.frame(sp1 = strain.labelA)
strain$strain = rownames(strain)
df.color = data.frame(strain = new.order, palette = palette)
df.color = left_join(df.color, strain, by = "strain")

plotA = ggplot(aucA, aes(x = sample, y = chg))+
      geom_point(aes(fill=sample), pch=21, size = 2, alpha = 0.8)+
      #geom_errorbar(data = aucA_summ, aes(y = mean_fold, ymin = mean_fold - sd_fold, ymax = mean_fold + sd_fold), 
      #              width = 0.3, alpha = 0.6)+
      #geom_point(data = aucA_summ, aes(y = mean_fold), shape = 3, size = 2, alpha = 0.6)+
      theme_bw()+
      theme(text = element_text(size = 15),
            plot.title =  element_text(size = 10),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
      coord_fixed(ratio = 11/2)+
      scale_y_continuous(name = "Fold change\nPe299R::mSc", limits = c(-1, 1))+
      scale_x_discrete(labels = strain.labelA, name ="")+
      labs(title = "Two-species competition")+
      scale_fill_manual(values = palette)


      # models
hist(aucA$chg)

nullA = lm(chg ~ 1, data = aucA)
summary(nullA)

lmA = lm(chg ~ sample, data = aucA)
summary(lmA)
plot(resid(lmA))
qqnorm(resid(lmA))

summary(aov(lmA))

Anova(lmA, nullA, type = "III")

emA = emmeans(lmA, specs = trt.vs.ctrl ~ sample, ref = 2, adjust = "bonferroni")
emA
effA = eff_size(emmeans(lmA, "sample"), sigma = sigma(lmA), edf = df.residual(lmA))
effA


### two competitors
mean_intraB = mean(data_auc$auc[data_auc$cond == "B" & data_auc$sample == "x01"])

aucB = data_auc %>% 
      filter(cond == "B") %>% 
      group_by(sample) %>% 
      mutate(chg = log2(auc/mean_intraB))

aucB_summ = aucB %>% 
      summarise(mean_fold = mean(chg),
                sd_fold = sd(chg))

plotB = ggplot(aucB, aes(x = sample, y = chg))+
      geom_point(size = 1.5, alpha = 0.8)+
      geom_errorbar(data = aucB_summ, aes(y = mean_fold, ymin = mean_fold - sd_fold, ymax = mean_fold + sd_fold), 
                    width = 0.3)+
      geom_point(data = aucB_summ, aes(y = mean_fold), shape = 23, size = 2, fill = "red", alpha = 0.6)+
      theme_bw()+
      theme(text = element_text(size = 15),
            plot.title =  element_text(size = 10),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
      coord_fixed(ratio = 11/2)+
      scale_y_continuous(name = "Fold change\nPe299R::mSc", limits = c(-1, 1))+
      scale_x_discrete(labels = strain.labelB, name ="")+
      labs(title = "Three-species competition")


      # models
lmB = lm(chg ~ sample, data = aucB)
summary(lmB)

emB = emmeans(lmB, specs = trt.vs.ctrl ~ sample, ref = 1, adjust = "bonferroni")

effB = eff_size(emmeans(lmB, "sample"), sigma = sigma(lmB), edf = df.residual(lmB))



## matrices
# Palette
my_palette <- colorRamp2(c(-1, 0, 1), c("#004488", "#ffffff", "#40004B"))

matA = matrix(aucA_summ$mean_fold)
row.names(matA) = strain.labelA
clusA = as.dendrogram(hclust(dist(matA)))

matB = matrix(aucB_summ$mean_fold)
row.names(matB) = strain.labelB
clusB = as.dendrogram(hclust(dist(matB)))


# Ward Hierarchical Clustering
dA <- dist(matA, method = "euclidean") %>% 
      hclust(., method="ward.D") 
dB <- dist(matB, method = "euclidean") %>% 
      hclust(., method="ward.D") 


# HM
hA = Heatmap(matA,
             rect_gp = gpar(col= "#7e7e7e"),
             col = my_palette,
             border = TRUE,
             show_heatmap_legend = FALSE,
             cluster_rows = dA,
             cluster_columns = TRUE,
             row_dend_width  = unit(5, "cm"),
             cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f", matA[i, j]), x, y, gp = gpar(fontsize = 10))
             })
hB = Heatmap(matB,
             rect_gp = gpar(col= "#7e7e7e"),
             col = my_palette,
             border = TRUE,
             show_heatmap_legend = FALSE,
             cluster_rows = dB,
             cluster_columns = TRUE,
             row_dend_width  = unit(5, "cm"),
             cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f", matB[i, j]), x, y, gp = gpar(fontsize = 10))
             })

hA+hB

#legends
#lgd1 main legend
lgd = Legend(title = "Fold\nchange",
             title_gap = unit(5, "mm"),
             title_gp = gpar(fontsize=15),
             col=my_palette,
             labels_gp = gpar(fontsize = 13),
             at = c(-1, 0, 1),
             legend_height = unit(2, "cm"),
             border = TRUE)
dev.off()
draw(hA, padding = unit(c(2, 30, 2, 2), "mm"))
draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))


##  SAVE FILES
pdf("~/Google Drive/Rudolf_MRE/results/ch4/competition/compA.pdf", width = 4, height = 4, useDingbats=FALSE)
plotA
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/competition/compA_hm.pdf", width = 5, height = 4, useDingbats=FALSE)
draw(hA, padding = unit(c(1, 40, 1, 1), "mm"))
draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/competition/compB.pdf", width = 5, height = 5, useDingbats=FALSE)
plotB
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/competition/compB_hm.pdf", width = 6, height = 4, useDingbats=FALSE)
draw(hB, padding = unit(c(1, 45, 1, 1), "mm"))
draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
dev.off()

write.csv(aucA, "~/Google Drive/Rudolf_MRE/results/ch4/competition/comp_chg.csv")

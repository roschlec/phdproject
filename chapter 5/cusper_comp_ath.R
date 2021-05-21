###   Resource competition and reproductive success  ###

##    libraries
library(tidyverse)
library(patchwork)
library(ggridges)
library(ggbeeswarm)
library(gamlss)
library(reshape2)
library(fpc)
library(ggdendro)
library(gridExtra)

##    labels
sp.lab = c(mono = "Monoculture",
           pk = "PkP19E3",
           pss = "PssB728a",
           arthr = "ArthrL145",
           wt = "Pe299R",
           meth = "MethyL85",
           smfr = "SmFR1",
           rhod = "RhodoL225")

pal = inner_join(data.frame(sp1 = sp.lab), df.color, by = "sp1")

##    load data
setwd("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/")
int = read.csv("intensity.csv", header = T)
bg = read.csv("background.csv", header = T)
df = inner_join(int, bg, by = c("sample", "time", "rep", "image")) %>% 
      dplyr::select(sample, time, rep, image, number, img, int, background) %>% 
      mutate(int_corr = int - background)
      #mutate(int_corr = int)
df$sample = factor(df$sample, levels = c("mono", "pk", "pss", "arthr", "wt", "meth", "smfr", "rhod"), labels = sp.lab)

##    initial fluorescence
t0 = df %>% 
      filter(time == "t0") %>% 
      group_by(sample, rep) %>% 
      summarise(mean_t0 = mean(int_corr))
t0

##    CUSPER calculation
cusper = df %>% 
      filter(int_corr > 0) %>% 
      inner_join(., t0, by = c("sample", "rep")) %>% 
      mutate(gfp = int_corr/mean_t0,
            success = log2(mean_t0/int_corr)) %>% 
      dplyr::select(sample, time, rep, image, number, img, gfp, success)
cusper

##    limit of detection
lod = df %>% 
      group_by(time, sample, rep) %>% 
      summarise(bg = mean(background),
                sd = sd(background),
                loq = case_when(
                      sd != NA ~ (bg+1.96*1),
                      TRUE ~ (bg+1.96*sd))) %>% 
      inner_join(df, ., by = c("time", "sample", "rep")) %>% 
      mutate(int_lod = if_else(loq - background < 0, 1, loq - background),
            lod = log2(mean(t0$mean_t0)/int_lod)) %>% 
      filter(int_lod != 1) %>% 
      ungroup() %>% 
      summarise(mean_lod = mean(lod),
                sd_lod = sd(lod),
                up = mean_lod + sd_lod,
                low = mean_lod - sd_lod)
lod

cusper_n = cusper %>% group_by(sample, time)
cusper_n = sample_n(cusper_n, 200, replace = TRUE)


##    Plotting
qq1 = ggplot(cusper_n, aes(sample = success, color = sample, shape = time))+
      geom_qq(size = 2, alpha = 0.6)+
      theme_bw()+
      scale_y_continuous(name = "Reproductive success (# generations)", limits = c(-2,10), breaks = seq(0,10,2), expand = c(0,0))+
      scale_x_continuous(name = "Theoretical quantiles", limits = c(-4,4), expand = c(0,0))+
      scale_color_manual(values = as.character(pal[,3]), name = "Competitor")+
      coord_fixed(ratio = 8/12)+
      geom_hline(data = lod, aes(yintercept = mean_lod))
qq1

qq2 = ggplot(cusper, aes(sample = success, fill = sample, group = interaction(sample,time)))+
      geom_qq(size = 3, alpha = 0.6, shape = 21, color = "dark grey")+
      theme_bw()+
      scale_y_continuous(name = "Reproductive success (# generations)", limits = c(-2,4), breaks = seq(0,10,2), expand = c(0,0))+
      scale_x_continuous(name = "Theoretical quantiles", limits = c(-2.5,0), expand = c(0,0))+
      scale_fill_manual(values = as.character(pal[,3]), name = "Competitor")+
      coord_fixed(ratio = 8/12)
qq2

plot = ggplot(cusper, aes(x = success, y = reorder(sample, desc(sample)), fill = time))+
      geom_density_ridges(alpha = 0.2, scale = 5, panel_scaling = T)+
      scale_x_continuous(limit = c(-2, 8), breaks = seq(0,8,1))+
      scale_y_discrete(expand = expansion(mult = c(0, 0.3))) +
      theme_ridges(grid = TRUE, center_axis_labels = TRUE)+
      geom_vline(data = lod, aes(xintercept = mean_lod))
plot

viol = ggplot(cusper, aes(x = sample, y = success))+
      facet_wrap(~time, ncol=1)+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray", alpha = 0.2)+
      geom_point(alpha = 0.1, size = 1)+
      theme_bw()+
      scale_x_discrete(name = "Competitor")+
      scale_y_continuous(name = "Reproductive\nsuccess")

##    Models
cus1 = cusper %>% filter(time == "t1" & sample != "Monoculture" & success < 3 & success < 6)
cus1$sample = factor(cus1$sample, levels = c("Pe299R", "PkP19E3", "PssB728a", "ArthrL145", "MethyL85", "SmFR1", "RhodoL225"))
cus1$rep = factor(cus1$rep)
cus1$image = factor(cus1$image)
cus2 = cusper %>% filter(time == "t1" & sample != "Monoculture" & success >= 3 & success < 6)
cus2$sample = factor(cus2$sample, levels = c("Pe299R", "PkP19E3", "PssB728a", "ArthrL145", "MethyL85", "SmFR1", "RhodoL225"))
cus2$rep = factor(cus2$rep)
cus2$image = factor(cus2$image)
hist(cus1$success)
hist(cus2$success)

dens1 = ggplot(cus1, aes(x = success))+
      facet_wrap( ~ sample, ncol = 7)+
      geom_density(size = 1)+
      theme_bw()+
      scale_x_continuous(name = "Reproductive success", limits = c(-2,7), breaks = seq(0,6,2))+
      scale_y_continuous(name = "Kernel density")

dens2 = ggplot(cus2, aes(x = success))+
      facet_wrap( ~ sample, ncol = 7)+
      geom_density(size = 1)+
      theme_bw()+
      scale_x_continuous(name = "Reproductive success", limits = c(-2,7), breaks = seq(0,6,2))+
      scale_y_continuous(name = "Kernel density")

viol1 = ggplot(cus1, aes(x = sample, y = success))+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray", alpha = 0.2)+
      geom_point(alpha = 0.1, size = 1)+
      theme_bw()+
      scale_x_discrete(name = "Competitor")+
      scale_y_continuous(name = "Reproductive\nsuccess")

viol2 = ggplot(cus2, aes(x = sample, y = success))+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray", alpha = 0.2)+
      geom_point(alpha = 0.1, size = 1)+
      theme_bw()+
      scale_x_discrete(name = "Competitor")+
      scale_y_continuous(name = "Reproductive\nsuccess")

#     model competition
model1 = gamlss(success ~ sample + re(random= ~1|rep/image), data = cus1, family = ST1)
summary(model1)
qqnorm(resid(model1))
plot(resid(model1))

model2 = gamlss(success ~ sample + re(random= ~1|rep/image), data = cus2, family = ST1)
summary(model2)
qqnorm(resid(model2))
plot(resid(model2))

### Fractions
##GFP/sum(GFP)
head(cusper)

cus3 = cusper %>% 
      group_by(sample, time, rep) %>% 
      summarise(rf = gfp/sum(gfp),
                og = log2(1/gfp)) %>% 
      mutate(rs = case_when(
            og < 0.5 ~ 0,
            og >= 0.5 & og < 1.5 ~ 1,
            og >= 1.5 & og < 2.5 ~ 2,
            og >= 2.5 & og < 3.5 ~ 3,
            og >= 3.5 & og < 4.5 ~ 4,
            og >= 4.5 & og < 5.5 ~ 5,
            og >= 5.5 & og < 6.1 ~ 6,
            og >= 6.1 ~ 7))

sample = unique(cus3$sample)
replicate = seq(1,4)
rs = seq(0,7)
n = 0

tb = data.frame(sample = sort(rep(sample, times = length(replicate*rs))),
           rep = sort(rep(replicate, times = length(rs))),
           rs = rs,
           n = 0)

cus4 = cus3 %>% 
      filter(time == "t1") %>% 
      group_by(sample, rep, rs) %>% 
      tally() %>% 
      data.frame()

cus5 = anti_join(tb, cus4, by = c("sample", "rep", "rs")) %>% 
      rbind(., cus4) %>% 
      arrange(sample, rep, rs) %>% 
      group_by(sample, rep) %>% 
      mutate(fraction = n/sum(n))

cus6 = cus5 %>% 
      group_by(sample, rs) %>% 
      summarise(mean = mean(fraction),
                sd = sd(fraction))
cus6$sample = factor(cus6$sample, levels = c("Pe299R", "Monoculture", "PkP19E3", "PssB728a", "ArthrL145", "MethyL85", "SmFR1", "RhodoL225"))

cus5 %>% group_by(sample, rs) %>% tally() %>% print(n = Inf)

fraction_plot = ggplot(cus6, aes(x = as.factor(rs), y = mean))+
      facet_wrap(~ sample, ncol = 4)+
      geom_col(alpha = 0.8)+
      geom_point(data = cus5, aes(x=as.factor(rs), y=fraction))+
      scale_y_continuous(limits = c(0,1), name = "Relative fraction of\nindividuals at t = 0")+
      scale_x_discrete(name = "Reproductive success", labels = c("0", "1", "2", "3", "4", "5", "6", ">6"))+
      theme_bw()+
      theme(axis.title = element_text(size = 18),
            axis.text = element_text(size = 15))

# model
cus5$sample = factor(cus5$sample, levels = c("Pe299R", "Monoculture", "PkP19E3", "PssB728a", "ArthrL145", "MethyL85", "SmFR1", "RhodoL225"))
cus5$sample = factor(cus5$sample, levels = c("Monoculture","Pe299R", "PkP19E3", "PssB728a", "ArthrL145", "MethyL85", "SmFR1", "RhodoL225"))
glm1 = gamlss(n ~ sample, data = cus5, family = ZAP)
summary(glm1)

## Clustering
cus7 = cus6 %>% 
      pivot_wider(id_cols = sample, names_from = rs, values_from = mean) %>% 
      as.data.frame()
colnames(cus7) = c("sample", "RS0", "RS1", "RS2", "RS3", "RS4", "RS5", "RS6", "RS7")
rownames(cus7) = cus7$sample
matcus = as.matrix(cus7[,-1])

#   Create the distance matrix.
d <- dist(matcus, method="euclidean") 

#   Do the clustering. 
pfit <- hclust(d, method = "ward.D")   
dfit <- as.dendrogram(pfit)

#   Plot the dendrogram.
plot(pfit)

# set the desired number of clusters                               
kbest.p<-3

#   Run clusterboot() with hclust 
#   ('clustermethod=hclustCBI') using Ward's method 
#   ('method="ward"') and kbest.p clusters 
#   ('k=kbest.p'). Return the results in an object 
#   called cboot.hclust.
cboot.hclust <- clusterboot(matcus, clustermethod = hclustCBI,
                            method = "ward.D", k = 3, B = 1000)

cboot.hclust$bootmean 
cboot.hclust$bootbrd

ddata <- dendro_data(dfit, type = "rectangle")
p1 = ggplot(segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      coord_flip() + 
      scale_y_reverse(expand = c(.5, .2))+
      theme_dendro()+
      geom_text(data = label(ddata), 
                aes(x = x, y = y, label = label, hjust = 0), 
                size = 6)

dendro.order = c("PkP19E3","ArthrL145", "Pe299R", "PssB728a", "MethyL85", "Monoculture", "SmFR1", "RhodoL225")
cus6$sample = factor(cus6$sample, 
                     levels = dendro.order)
p2 = ggplot(cus6, aes(x = as.factor(rs), y = sample, size = mean, color = mean))+
      geom_point()+
      scale_color_viridis_c(direction = -1) +
      theme_minimal() +
      scale_x_discrete(name = "Reproductive success", 
                       labels = c("0", "1", "2", "3", "4", "5", "6", ">6"),
                                  breaks = seq(0, 7, 1)) +
      theme(axis.text.y = element_blank()) +
      ylab(NULL)

grid.arrange(p1, p2, ncol = 2, widths = 3:2)


dist_rs = melt(as.matrix(dist(matcus, method = "euclidean")), varnames = c("sp1", "sp2"), value.name = "rs")
dist_rs[dist_rs$sp2 == "Pe299R",]

#     Save plots
pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/qq1.pdf", height = 6, width = 6, useDingbats = F)
qq1
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/qq2.pdf", height = 6, width = 6, useDingbats = F)
qq2
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/plot.pdf", height = 8, width = 12, useDingbats = F)
plot
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/viol.pdf", height = 10, width = 12, useDingbats = F)
viol
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/density.pdf", height = 10, width = 12, useDingbats = F)
dens1/dens2
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/violin.pdf", height = 10, width = 12, useDingbats = F)
viol1/viol2
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/dist_rs.pdf", height = 6, width = 12, useDingbats = F)
fraction_plot
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/cluster_rs.pdf", height = 4, width = 10, useDingbats = F)
grid.arrange(p1, p2, ncol = 2, widths = 3:2)
dev.off()

write.csv(dist_rs, "~/Google Drive/Rudolf_MRE/results/ch4/cusper/cusper_comp_ath/dist_rs.csv")

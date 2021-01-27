### STATS CHAPTER 4 ### FBA

# Library
library(tidyverse)
library(ggrepel)
library(patchwork)
library(vegan)
library(ROCR)

setwd("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final")

label = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/labels.csv", header = T)

phyl = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_phyl.csv", header = T)[,-1] %>% 
      filter(sp1 != "Pseudomonas_koreensis_P19E3" & sp2 != "Pseudomonas_koreensis_P19E3") %>% 
      arrange(sp2, sp1)

fba = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_fba.csv", header = T)[,-1] %>% 
      arrange(sp2, sp1)

no = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/pianka.csv", header = T)[,-1] %>% 
      arrange(sp2, sp1)

upt = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_uptake.csv", header = T)[,-1] %>% 
      arrange(sp2, sp1)

exc = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_excretion.csv", header = T)[,-1] %>% 
      arrange(sp2, sp1)

#RENAMING

phyl$sp1 = str_replace_all(as.character(phyl$sp1), setNames(as.character(label$nameB), as.character(label$nameC)))
phyl$sp2 = str_replace_all(as.character(phyl$sp2), setNames(as.character(label$nameB), as.character(label$nameC)))
fba$sp1 = str_replace_all(as.character(fba$sp1), setNames(as.character(label$nameB), as.character(label$name)))
fba$sp2 = str_replace_all(as.character(fba$sp2), setNames(as.character(label$nameB), as.character(label$name)))
no$sp1 = str_replace_all(as.character(no$sp1), setNames(as.character(label$nameB), as.character(label$name)))
no$sp2 = str_replace_all(as.character(no$sp2), setNames(as.character(label$nameB), as.character(label$name)))
upt$sp1 = str_replace_all(as.character(upt$sp1), setNames(as.character(label$nameB), as.character(label$name)))
upt$sp2 = str_replace_all(as.character(upt$sp2), setNames(as.character(label$nameB), as.character(label$name)))
exc$sp1 = str_replace_all(as.character(exc$sp1), setNames(as.character(label$nameB), as.character(label$name)))
exc$sp2 = str_replace_all(as.character(exc$sp2), setNames(as.character(label$nameB), as.character(label$name)))

## data frames
df = inner_join(phyl, fba, by = c("sp1", "sp2")) %>% 
      left_join(., no, by = c("sp1", "sp2")) %>% 
      left_join(., upt, by = c("sp1", "sp2")) %>%
      left_join(., exc, by = c("sp1", "sp2")) %>%
      filter(sp1 != sp2) %>% 
      droplevels()

df_phyl = pivot_wider(phyl, id_cols = sp1, names_from = sp2, values_from = value)[,-1]
df_fba = pivot_wider(fba, id_cols = sp1, names_from = sp2, values_from = fba)[,-1]
df_upt = pivot_wider(upt, id_cols = sp1, names_from = sp2, values_from = uptake)[,-1]
df_exc = pivot_wider(exc, id_cols = sp1, names_from = sp2, values_from = excreted)[,-1]
df_no = pivot_wider(no, id_cols = sp1, names_from = sp2, values_from = no)[,-1]

## matrices
mat_phyl = as.dist(df_phyl)
mat_fba = as.dist(df_fba)
mat_upt = as.dist(df_upt)
mat_exc = as.dist(df_exc)
mat_no = as.dist(df_no)
mat_no2 = dist(df_no, method = "euclidean")

## NMDS
mds_phyl = monoMDS(mat_phyl)
mds_fba = monoMDS(mat_fba)
mds_upt = monoMDS(mat_upt)
mds_exc = monoMDS(mat_exc)
mds_no = monoMDS(mat_no2)


## MANTEL
#     PHYL VS
mantel(mat_phyl, mat_fba, permutations = 999)
mantel(mat_phyl, mat_upt, permutations = 999)
mantel(mat_phyl, mat_exc, permutations = 999)
mantel(mat_phyl, mat_no2, permutations = 999)

#     FBA VS
mantel(mat_fba, mat_upt, permutations = 999)
mantel(mat_fba, mat_exc, permutations = 999)
mantel(mat_fba, mat_no2, permutations = 999)

#     UPT VS
mantel(mat_upt, mat_exc, permutations = 999)
mantel(mat_upt, mat_no2, permutations = 999)

#     EXC VS
mantel(mat_exc, mat_no2, permutations = 999)

## PROCRUSTES
##Figure 6
protest(mds_phyl, mds_fba, permutations = 999, scores = "sites")
protest(mds_phyl, mds_no, permutations = 999, scores = "sites")
protest(mds_phyl, mds_upt, permutations = 999, scores = "sites")
protest(mds_phyl, mds_exc, permutations = 999, scores = "sites")

##Figure7
protest(mds_no, mds_fba, permutations = 999, scores = "sites")
protest(mds_no, mds_upt, permutations = 999, scores = "sites")
protest(mds_no, mds_exc, permutations = 999, scores = "sites")


#plots
#phyl vs
p1 = ggplot(df, aes(x = value, y = fba))+
      geom_point(alpha = 0.2)+
      theme_bw()+
      scale_x_continuous(name = "Phylogenetic distance",
                         limits = c(0,1.5))+
      scale_y_continuous(name = "Metabolic distance",
                         limits = c(0,3))+
      coord_fixed(ratio = .5)

p2 = ggplot(df, aes(x = value, y = no))+
      geom_point(alpha = 0.2)+
      theme_bw()+
      scale_x_continuous(name = "Phylogenetic distance",
                         limits = c(0,1.5))+
      scale_y_continuous(name = "Resource overlap",
                         limits = c(0,1))+
      coord_fixed(ratio = 1.5)

p3 = ggplot(df, aes(x = value, y = uptake))+
      geom_point(alpha = 0.2)+
      theme_bw()+
      scale_x_continuous(name = "Phylogenetic distance",
                         limits = c(0,1.5))+
      scale_y_continuous(name = "Resource uptake dissimilarity",
                         limits = c(0,1.2))+
      coord_fixed(ratio = 1.25)

p4 = ggplot(df, aes(x = value, y = excreted))+
      geom_point(alpha = 0.2)+
      theme_bw()+
      scale_x_continuous(name = "Phylogenety distance",
                         limits = c(0,1.5))+
      scale_y_continuous(name = "Byproduct excretion dissimilarity",
                         limits = c(0,1.5))+
      coord_fixed(ratio = 1)


#metabolism vs resource overlap
p5 = ggplot(df, aes(x = fba, y = no))+
      geom_point(alpha = 0.2)+
      theme_bw()+
      scale_x_continuous(name = "Metabolic distance",
                         limits = c(0,3))+
      scale_y_continuous(name = "Resource overlap",
                         limits = c(0,1))+
      coord_fixed(ratio = 3)

p6 = ggplot(df, aes(x = uptake, y = no))+
      geom_point(alpha = 0.2)+
      theme_bw()+
      scale_x_continuous(name = "Resource uptake dissimilarity",
                         limits = c(0,1.2),
                         breaks = seq(0, 1.2, 0.4))+
      scale_y_continuous(name = "Resource overlap",
                         limits = c(0,1))+
      coord_fixed(ratio = 1.2)

p7 = ggplot(df, aes(x = excreted, y = no))+
      geom_point(alpha = 0.2)+
      theme_bw()+
      scale_x_continuous(name = "Byproduct excretion dissimilarity",
                         limits = c(0,1.5))+
      scale_y_continuous(name = "Resource overlap",
                         limits = c(0,1))+
      coord_fixed(ratio = 1.5)

## save plots!
pdf("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/plot_phylvs.pdf", height = 8, width = 8, useDingbats = F)
(p1+p2)/(p3+p4)
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/plot_met.pdf", height = 4, width = 8, useDingbats = F)
p5+p6+p7
dev.off()




### FBA vs EMPIRICAL

cprof1 = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_cprof.csv", header = T, check.names=FALSE)[,-1] %>% 
      arrange(sp2, sp1)
cprof2 = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_cprof2.csv", header = T, check.names=FALSE)[,-1] %>% 
      arrange(sp1)

fba1 = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_fba_short.csv", header = T)[,-1] %>% 
      arrange(sp2, sp1)
fba2 = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_fba2_short.csv", header = T)[,-1] %>% 
      arrange(sp1)

# RENAMING

fba1$sp1 = str_replace_all(as.character(fba1$sp1), setNames(as.character(label$nameB), as.character(label$name)))
fba1$sp2 = str_replace_all(as.character(fba1$sp2), setNames(as.character(label$nameB), as.character(label$name)))
fba2$sp1 = str_replace_all(as.character(fba2$sp1), setNames(as.character(label$nameB), as.character(label$name)))

cprof2 = cprof2 %>% arrange(sp1)
fba2 = fba2 %>% arrange(sp1)

## data
empi = pivot_wider(cprof2, id_cols = resource, names_from = sp1, values_from = cprof)[,-1]
pred = pivot_wider(fba2, id_cols = resource, names_from = sp1, values_from = pred.cprof)[,-1]

## matrices
mat_empi = dist(empi, "euclidean")
mat_pred = dist(pred, "euclidean")

## NMDS
mds_empi = monoMDS(mat_empi)
mds_pred = monoMDS(mat_pred)


## DISTANCE MATRIX FBA VS AUC
mantel(mat_empi, mat_pred, permutations = 1000)
protest(mds_empi, mds_pred, permutations = 9999, scores = "sites")
pro = procrustes(mds_empi, mds_pred)
plot(pro)



##ROC curve analysis
fba_bin = fba2 %>% 
      mutate(bin_fba = case_when(
            pred.cprof >= 0.1 ~ 1,
            TRUE ~ 0))
cprof_bin = cprof2 %>% 
      mutate(bin_auc = case_when(
            cprof >= 0.1 ~ 1,
            TRUE ~ 0))
binary <- inner_join(fba_bin, cprof_bin, by = c("resource", "sp1"))

pred <- prediction(binary$pred.cprof, binary$bin_auc)
pred

perf <- performance(pred, "tpr", "fpr")
perf

performance(pred, measure="tpr", x.measure = "fpr") %>% plot()
performance(pred, measure="fnr") %>% plot()
performance(pred, measure="acc") %>% plot()

auc.perf <- performance(pred, measure="auc")
auc.perf@y.values

acc.perf <- performance(pred, measure="acc")

opt.cut = function(perf, pred){
      cut.ind = mapply(FUN=function(x, y, p){
            d = (x - 0)^2 + (y-1)^2
            ind = which(d == min(d))
            c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
              cutoff = p[[ind]])
      }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(perf, pred))

cost.perf = performance(pred, "cost")
pred@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]

ind = which.max(slot(acc.perf, "y.values")[[1]] )
acc = slot(acc.perf, "y.values")[[1]][ind]
cutoff = slot(acc.perf, "x.values")[[1]][ind]
accu = data.frame(accuracy= acc, cutoff = cutoff) %>% format(digits = 2)
m <- data.frame(x = perf@x.values, y = perf@y.values)
colnames(m) <- c("x", "y")

rocplot <- ggplot(m, aes(x,y))+
      geom_line(size=1.5)+
      geom_text(data = accu, size = 5, hjust=0, color = "dark blue",
                mapping = aes(x=0.02, y=0.95, label = paste('Accuracy ==', accuracy)), parse = T)+
      theme_bw()+
      coord_fixed(ratio = 1)+
      scale_x_continuous(name = "1-Specificity (False Positive Rate)",
                         limits=c(0,1), expand = expansion(0))+
      scale_y_continuous(name = "Sensitivity (True Positive Rate)",
                         limits=c(0,1), expand = expansion(0))
rocplot


pdf("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/plot_roc.pdf", height = 4, width = 4, useDingbats = F)
rocplot
dev.off()


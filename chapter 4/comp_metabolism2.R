####        CFU BINARY COMMUNITIES        ####
####        RELATIONSHIP BETWEEN SP INTERACTION AND METABOLIC SIMILARITY ####

##    Libraries
library(tidyverse) #ggplot2, tibble, tidyr, readr, purrr, dplyr, stringr, forcats
library(lme4)
library(nlme)
library(lmtest)
library(emmeans)
library(car)
library(reshape)
library(MuMIn)
library(mgcv)
library(ecodist)
library(ciTools)
library(ade4)
library(patchwork)

##    Set working directory and load data
#    File contains the composition of the binary community, the date of the experiment, and cfu counts of each bacterial strain over time (0, 4, 7, and 14 dpi)
#    Strain IDs:
#          RS4-60	      smfr1_1	::Tn5::mSc::GmR
#          RS5-17	      smfr1_2	::Tn5::mTq2::TcR
#          RS3-59	      spfa2_1	::Tn5::mSc::GmR
#          RS4-1	      spfa2_2	::Tn5::mTq2::TcR
#          RS4-9	      meth85_1	::Tn5::mSc::GmR
#          RS4-14	      meth85_2	::Tn5::sYFP2::KmR
#          RS5-3	      meth92_1	::Tn5::mSc::GmR
#          RS4-21	      meth92_2	::Tn5::mTq2::KmR
#          RS5-1	      mr01_1	::Tn5::mSc::GmR
#          RS4-24	      mr01_2	::Tn5::mTq2::KmR

## COMUNITY LABELS

comm.label = c("meth85_meth85", "meth85_meth92", "mr01_meth85",
               "meth92_meth92", "mr01_meth92", "mr01_mr01",
               "meth85_smfr1", "smfr1_meth92", "mr01_smfr1",
               "smfr1_smfr1", "smfr1_spfa2", "meth85_spfa2",
               "meth92_spfa2", "mr01_spfa2", "spfa2_spfa2")
comm.label2 = c("meth85_meth92", "mr01_meth85", "mr01_meth92", "meth85_smfr1", "smfr1_meth92", 
                "mr01_smfr1", "smfr1_spfa2", "meth85_spfa2", "meth92_spfa2", "mr01_spfa2")
comm.label3 = c("meth85_meth92", "meth85_smfr1", "meth85_spfa2", "smfr1_meth92", "meth92_spfa2",
                "mr01_meth85", "mr01_meth92", "mr01_smfr1", "mr01_spfa2", "smfr1_spfa2")

      
#     CFU data INTRASPECIFIC
setwd("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/")
het_comp = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/interaction_intra.csv",
                    header = T, check.names = F)[,-1]
het_comp$dpi = factor(het_comp$dpi)
het_comp$community = factor(het_comp$community, labels = comm.label2)  
het_comp = het_comp[het_comp$community %in% comm.label2, ]
head(het_comp)

#     In vitro data
homo_comp = read.csv("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/interaction_intra.csv", 
                   header = T, check.names = F)[,-1]
colnames(homo_comp)[1] = "community"
homo_comp$community = factor(homo_comp$community, levels = c("meth85_meth92", "mr01_meth85",
                                                             "mr01_meth92", "meth85_smfr1",
                                                             "meth92_smfr1", "mr01_smfr1", 
                                                             "smfr1_spfa2", "meth85_spfa2",
                                                             "meth92_spfa2", "mr01_spfa2"))  
homo_comp$community = factor(homo_comp$community, labels = comm.label2)  
homo_comp = homo_comp[homo_comp$community %in% comm.label2, ]


##    Heterogeneous vs Homogeneous environment
levels(het_comp$community)
levels(homo_comp$community)

comp = inner_join(het_comp, homo_comp, by = c("community", "relation"), suffix = c("_het", "_hom"))
comp

comp_intrafam <- comp %>% filter(relation == "Intra-Family")
comp_intrafam_7day <- comp_intrafam %>% filter(dpi == "7")
comp_intrafam_14day <- comp_intrafam %>% filter(dpi == "14")

comp_interfam <- comp %>% filter(relation == "Inter-Family")
comp_interfam_7day <- comp_interfam %>% filter(dpi == "7")
comp_interfam_14day <- comp_interfam %>% filter(dpi == "14")


## interaction types
qqnorm(comp$int_het)
qqnorm(comp$int_hom)

model_intrafam_7 = glmmTMB(int_het ~ int_hom + (1|resource/community), data = comp_intrafam_7day, family=beta_family(link = "logit"))
summary(model_intrafam_7)
qqnorm(resid(model_intrafam_7))
plot(resid(model_intrafam_7))

model_intrafam_14 = glmmTMB(int_het ~ int_hom + (1|resource/community), data = comp_intrafam_14day, family=beta_family(link = "logit"))
summary(model_intrafam_14)
qqnorm(resid(model_intrafam_14))
plot(resid(model_intrafam_14))

model_interfam_7 = glmmTMB(int_het ~ int_hom + (1|resource/community), data = comp_interfam_7day, family=beta_family(link = "logit"))
summary(model_interfam_7)
qqnorm(resid(model_interfam_7))
plot(resid(model_interfam_7))

model_interfam_14 = glmmTMB(int_het ~ int_hom + (1|resource/community), data = comp_interfam_14day, family=beta_family(link = "logit"))
summary(model_interfam_14)
qqnorm(resid(model_interfam_14))
plot(resid(model_interfam_14))



plot_intrafam_7day = ggplot(comp_intrafam_7day, aes(x=int_hom, y = int_het))+
      geom_point(aes(fill = resource), size = 2, alpha = 0.7, pch = 21)+
      annotate(geom="text", x = 0.8, y = 1,
               label = paste("p-value =", round(summary(model_intrafam_7)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15),
            title = element_text(size = 12))+
      scale_fill_discrete(name = "Resource", labels = c("Fructose", "Fumarate", "Glucose", "Glutamate", "Methanol"))+
      xlim(0,1)+
      ylim(0,1)+
      labs(title = "Intra-Family (7 dpi)",
           x = "Interaction index\n(homogeneous environment)",
           y = "Interaction index\n(heterogeneous environment)")

plot_intrafam_14day = ggplot(comp_intrafam_14day, aes(x=int_hom, y = int_het))+
      geom_point(aes(fill = resource), size = 2, alpha = 0.7, pch = 21)+
      annotate(geom="text", x = 0.8, y = 1,
               label = paste("p-value =", round(summary(model_intrafam_14)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15),
            title = element_text(size = 12))+
      scale_fill_discrete(name = "Resource", labels = c("Fructose", "Fumarate", "Glucose", "Glutamate", "Methanol"))+
      xlim(0,1)+
      ylim(0,1)+
      labs(title = "Intra-Family (14 dpi)",
           x = "Interaction index\n(homogeneous environment)",
           y = "Interaction index\n(heterogeneous environment)")

plot_interfam_7day = ggplot(comp_interfam_7day, aes(x=int_hom, y = int_het))+
      geom_point(aes(fill = resource), size = 2, alpha = 0.7, pch = 21)+
      annotate(geom="text", x = 0.8, y = 1,
               label = paste("p-value =", round(summary(model_interfam_7)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15),
            title = element_text(size = 12))+
      scale_fill_discrete(name = "Resource", labels = c("Fructose", "Fumarate", "Glucose", "Glutamate", "Methanol"))+
      xlim(0,1)+
      ylim(0,1)+
      labs(title = "Inter-Family (7 dpi)",
           x = "Interaction index\n(homogeneous environment)",
           y = "Interaction index\n(heterogeneous environment)")

plot_interfam_14day = ggplot(comp_interfam_14day, aes(x=int_hom, y = int_het))+
      geom_point(aes(fill = resource), size = 2, alpha = 0.7, pch = 21)+
      annotate(geom="text", x = 0.8, y = 1,
               label = paste("p-value =", round(summary(model_interfam_14)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15),
            title = element_text(size = 12))+
      scale_fill_discrete(name = "Resource", labels = c("Fructose", "Fumarate", "Glucose", "Glutamate", "Methanol"))+
      xlim(0,1)+
      ylim(0,1)+
      labs(title = "Inter-Family (14 dpi)",
           x = "Interaction index\n(homogeneous environment)",
           y = "Interaction index\n(heterogeneous environment)")

(plot_intrafam_7day+plot_intrafam_14day) / (plot_interfam_7day+plot_interfam_14day) +
      plot_layout(guides = "collect")


## interaction strengths
qqnorm(comp$magnitude_het)
qqnorm(comp$magnitude_hom)

str_intrafam_7 = lmer(magnitude_het ~ magnitude_hom + (1|resource/community), data = comp_intrafam_7day)
summary(str_intrafam_7)
qqnorm(resid(str_intrafam_7))
plot(resid(str_intrafam_7))

str_intrafam_14 = lmer(magnitude_het ~ magnitude_hom + (1|resource/community), data = comp_intrafam_14day)
summary(str_intrafam_14)
qqnorm(resid(str_intrafam_14))
plot(resid(str_intrafam_14))

str_interfam_7 = lmer(magnitude_het ~ magnitude_hom + (1|resource/community), data = comp_interfam_7day)
summary(str_interfam_7)
qqnorm(resid(str_interfam_7))
plot(resid(str_interfam_7))

str_interfam_14 = lmer(magnitude_het ~ magnitude_hom + (1|resource/community), data = comp_interfam_14day)
summary(str_interfam_14)
qqnorm(resid(str_interfam_14))
plot(resid(str_interfam_14))






comp2 = comp %>% dplyr::select(community, dpi, relation, resource, magnitude_het, magnitude_hom) %>% 
      pivot_longer(cols = c(magnitude_het, magnitude_hom), names_to = "env", values_to = "value" ) %>% 
      arrange(env) %>% 
      mutate(resource = case_when(
                  env == "magnitude_het" ~ "complex",
                  TRUE ~ as.character(resource)))
comp2 = unique(comp2)
comp2$resource = factor(comp2$resource)

str_env = lme(value ~ env*relation, random = ~1|community/resource, data = comp2)
summary(str_env)
plot(resid(str_env))

emmeans(str_env, specs = pairwise ~ env|relation)

plot_str_env = ggplot(comp2, aes(x=env, y=value))+
      geom_boxplot(fill = "grey", alpha = 0.5)+
      geom_point(aes(color = relation), alpha = 0.5, size = 3)+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_fill_discrete(name = "Resource", labels = c("Complex", "Fructose", "Fumarate", "Glucose", "Glutamate", "Methanol"))+
      scale_x_discrete(name = "Environment", labels = c("Heterogeneous", "Homogeneous"))+
      scale_y_continuous(name = "Interaction strength", limits = c(0,6))
      

comp3 = comp %>% dplyr::select(community, dpi, relation, resource, int_het, int_hom) %>% 
      pivot_longer(cols = c(int_het, int_hom), names_to = "env", values_to = "value" ) %>% 
      arrange(env) %>% 
      mutate(resource = case_when(
            env == "int_het" ~ "complex",
            TRUE ~ as.character(resource)))
comp3 = unique(comp3)
comp3$resource = factor(comp3$resource)

int_env = glmmTMB(value ~ env*relation + (1|community/resource), data = comp3, family=beta_family(link = "logit"))
summary(int_env)
plot(resid(int_env))

plot_int_env = ggplot(comp3, aes(x=env, y=value))+
      geom_boxplot(fill = "grey", alpha = 0.5)+
      geom_point(aes(color = relation), alpha = 0.5, size = 3)+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_fill_discrete(name = "Resource", labels = c("Complex", "Fructose", "Fumarate", "Glucose", "Glutamate", "Methanol"))+
      scale_x_discrete(name = "Environment", labels = c("Heterogeneous", "Homogeneous"))+
      scale_y_continuous(name = "Interaction type", limits = c(0,1))


plot_int_env + plot_str_env + plot_layout(guides = "collect")

######     METABOLIC SIMILARITIES
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

#     Filters according to community label
met_dist = df_met[df_met$mix %in% comm.label2,] %>% .[,c(2,1)]
met_dist$mix = factor(met_dist$mix)

######     PHYLOGENETIC SIMILARITIES
phyl = read.csv("~/Google Drive/Rudolf_MRE/results/phylogeny/kbase/amphora_genes/phyl_distance.csv", header = T) %>% 
      .[,c(2,3)]
phyl$community = factor(phyl$community)

## phylo vs met
phyl_met = inner_join(met_dist, phyl, by = "community") %>% 
      filter(phyl_dist != 0)
phyl_met$community = factor(phyl_met$community, labels = comm.label2)

####  mantel test

# Biomass yield metabolic similarity output
metab
colnames(metab) = sp.label
rownames(metab) = colnames(metab)
dist_met = as.dist(metab, diag = T)
dist_met

# Phylogeny
phyl_mat = read.csv("~/Google Drive/Rudolf_MRE/results/phylogeny/kbase/amphora_genes/phyl_matrix.csv", header = T) %>% 
      .[c(4,5,3,2,1),c(5,6,4,3,2)]
colnames(phyl_mat) = sp.label
rownames(phyl_mat) = colnames(phyl_mat)
dist_phyl = as.dist(phyl_mat, diag = T)
dist_phyl

# statistical test
set.seed(041120)
mantel <- mantel.rtest(dist_met, dist_phyl, nrepet = 9999)
mantel # p = 0.239

mantel(dist_met ~ dist_phyl, nperm=10000, nboot=500, mrank = F)

model = lm(dist_met~dist_phyl)
summary(model)
qqnorm(resid(model))
confint(model)

phylmet_plot = ggplot(phyl_met, aes(x=phyl_dist, y=met_dist1))+
      geom_point()+
      geom_smooth(method = "lm", se = F, color = "black", fullrange=TRUE)+
      theme_bw()+
      theme(text = element_text(size=15))+
      ylim(0,1.5)+
      xlim(0,1.5)+
      labs(y = "Metabolic distance",
           x = "Phylogenetic distance")+
      coord_fixed(ratio=1)

pdf("./phylmet_plot.pdf", height = 4, width = 4, useDingbats = F)
phylmet_plot
dev.off()


##    ALL
all <- inner_join(comp, phyl_met, by="community") %>% 
      mutate(comm = community) %>% 
      separate(comm, into = c("s1", "s2"), sep = "_")
head(all)

all_het <- inner_join(het_comp, phyl_met, by="community") %>% 
      mutate(comm = community) %>% 
      separate(comm, into = c("s1", "s2"), sep = "_")
head(all_het)

all_hom <- inner_join(homo_comp, phyl_met, by="community") %>% 
      mutate(comm = community) %>% 
      separate(comm, into = c("s1", "s2"), sep = "_")
head(all_hom)

#     models
qqnorm(unique(all$int_het))
qqnorm(unique(all$int_hom))
qqnorm(unique(all$met_dist1))
qqnorm(unique(all$phyl_dist))

      # homogeneous env data - interaction types
            # null models
null_hom = glmmTMB(int ~ 1, data = all_hom, family=beta_family(link = "logit"))
null_hom_res = glmmTMB(int ~ 1 + (1|resource), data = all_hom, family=beta_family(link = "logit"))
null_hom_comm = glmmTMB(int ~ 1 + (1|community), data = all_hom, family=beta_family(link = "logit"))
null_hom_res_comm = glmmTMB(int ~ 1 + (1|resource) + (1|community), data = all_hom, family=beta_family(link = "logit"))
null_hom_res_comm2 = glmmTMB(int ~ 1 + (1|resource/community), data = all_hom, family=beta_family(link = "logit"))
            # models
res_hom = glmmTMB(int ~ resource + (1|community), data = all_hom, family=beta_family(link = "logit"))
rel_hom = glmmTMB(int ~ relation + (1|resource) + (1|community), data = all_hom, family=beta_family(link = "logit"))
rel_res_hom = glmmTMB(int ~ relation*resource + (1|community), data = all_hom, family=beta_family(link = "logit"))
phyl_hom = glmmTMB(int ~ phyl_dist + (1|resource) + (1|community), data = all_hom, family=beta_family(link = "logit"))
phyl_res_hom = glmmTMB(int ~ phyl_dist*resource + (1|community), data = all_hom, family=beta_family(link = "logit"))
met_hom = glmmTMB(int ~ met_dist1 + (1|resource) + (1|community), data = all_hom, family=beta_family(link = "logit"))
met_res_hom = glmmTMB(int ~ met_dist1*resource + (1|community), data = all_hom, family=beta_family(link = "logit"))
phyl_met_res = glmmTMB(int ~ met_dist1*phyl_dist*resource + (1|community), data = all_hom, family=beta_family(link = "logit"))

aic_hom = AICc(null_hom, null_hom_res, null_hom_comm, null_hom_res_comm, null_hom_res_comm2, #null models
               res_hom, rel_hom, rel_res_hom, phyl_hom, phyl_res_hom, met_hom, met_res_hom,
               phyl_met_res) %>% 
      mutate(model = rownames(.),
            dAIC = AICc - max(AICc)) %>% 
      arrange(-dAIC) %>% 
      dplyr::select(model, df, AICc, dAIC)
aic_hom

anova(null_hom_res_comm, rel_hom)
Anova(rel_hom, type = "III")
anova(null_hom_res_comm, met_hom)
Anova(met_hom, type = "III")
anova(null_hom_res_comm, phyl_hom)
Anova(phyl_hom, type = "III")

summary(null_hom_res_comm)
qqnorm(resid(null_hom_res_comm))
summary(rel_hom)
qqnorm(resid(rel_hom))
summary(met_hom)
qqnorm(resid(met_hom))
summary(phyl_hom)
qqnorm(resid(phyl_hom))

fixed.effects(rel_hom)
random.effects(rel_hom)


      # heterogeneous env data - interaction types
            # null model
null_het = glmmTMB(int ~ 1, data = all_het, family=beta_family(link = "logit"))
null_het_dpi = glmmTMB(int ~ 1 + (1|dpi), data = all_het, family=beta_family(link = "logit"))
null_het_comm = glmmTMB(int ~ 1 + (1|community), data = all_het, family=beta_family(link = "logit"))
null_het_dpi_comm = glmmTMB(int ~ 1 + (1|dpi) + (1|community), data = all_het, family=beta_family(link = "logit"))
null_het_dpi_comm2 = glmmTMB(int ~ 1 + (1|dpi/community), data = all_het, family=beta_family(link = "logit"))

rel_het = glmmTMB(int ~ relation, data = all_het, family=beta_family(link = "logit"))
phyl_het = glmmTMB(int ~ phyl_dist, data = all_het, family=beta_family(link = "logit"))
met_het = glmmTMB(int ~ met_dist1, data = all_het, family=beta_family(link = "logit"))
phyl_met_het = glmmTMB(int ~ met_dist1*phyl_dist, data = all_het, family=beta_family(link = "logit"))

rel_comm_het = glmmTMB(int ~ relation + (1|community), data = all_het, family=beta_family(link = "logit"))
phyl_comm_het = glmmTMB(int ~ phyl_dist + (1|community), data = all_het, family=beta_family(link = "logit"))
met_comm_het = glmmTMB(int ~ met_dist1 + (1|community), data = all_het, family=beta_family(link = "logit"))
phyl_met_comm_het = glmmTMB(int ~ met_dist1:phyl_dist + (1|community), data = all_het, family=beta_family(link = "logit"))

aic_het = AICc(null_het,null_het_dpi, null_het_comm, null_het_dpi_comm, null_het_dpi_comm2,
               rel_het, phyl_het, met_het, phyl_met_het,
               rel_comm_het, phyl_comm_het, met_comm_het, phyl_met_comm_het) %>% 
      mutate(model = rownames(.),
             dAIC = AICc - max(AICc)) %>% 
      arrange(-dAIC) %>% 
      dplyr::select(model, df, AICc, dAIC)
aic_het

summary(null_het)
qqnorm(resid(null_het))
summary(met_het)
qqnorm(resid(met_het))
summary(phyl_het)
qqnorm(resid(phyl_het))
summary(rel_het)
qqnorm(resid(rel_het))

#     mean comparisons
emm_rel_hom = emmeans(rel_hom, specs = pairwise ~ relation, adjust = "mvt", type = "response")
emm_rel_hom
emm_res_hom = emmeans(res_hom, specs = pairwise ~ resource, adjust = "mvt", type = "response")
emm_res_hom
emm_rel_res_hom = emmeans(rel_res_hom, specs = pairwise ~ relation|resource, adjust = "mvt", type = "response")
emm_rel_res_hom
emm_rel_het = emmeans(rel_het, specs = pairwise ~ relation, adjust = "mvt", type = "response")
emm_rel_het

#     confidence intervals
label_models = c("rel_hom", "phyl_hom", "met_hom", "rel_het", "phyl_het", "met_het")
df = data.frame()

for(i in 1:length(label_models)){
      data <- data.frame(confint(get(label_models[i])))[2,]
      data[, ncol(data) + 1] <- label_models[i]
      df <- rbind(df, data)
}

colnames(df) = c("2.5% CI", "97.5% CI", "Estimate", "Model")
df$Variable = rownames(df)
rownames(df) = NULL
df = df[, c(4, 5, 3, 1 ,2)]

df

#     PLOTS INTERACTIONS
      #     by relations/life strategy   
ggplot(all, aes(x=relation, y = int_het))+
      geom_point()+
      theme_bw()

ggplot(all, aes(x=relation, y = int_hom))+
      geom_point()+
      theme_bw()

      #     by phylogeny 
plot_phyl_het = ggplot(all_het, aes(x=phyl_dist, y=int))+
      geom_point(pch = 21, alpha = 0.8)+
      annotate(geom="text", x = 0.55, y = 1,
               label = paste("p-value =", round(summary(phyl_het)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_x_continuous(limits = c(0, 1.2), breaks = c(0.2, 0.4, 0.6,0.8, 1.0))+
      scale_y_continuous(limits = c(0, 1))+
      labs(title = "Heterogeneous\nenvironment",
           x = "Phylogenetic distance",
           y = "Interaction type")

plot_phyl_hom = ggplot(all_hom, aes(x=phyl_dist, y=int))+
      geom_point(pch = 21, alpha = 0.8)+
      annotate(geom="text", x = 0.55, y = 1,
               label = paste("p-value =", round(summary(phyl_hom)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_x_continuous(limits = c(0, 1.2), breaks = c(0.2, 0.4, 0.6,0.8, 1.0))+
      scale_y_continuous(limits = c(0, 1))+
      labs(title = "Homogeneous\nenvironment",
           x = "Phylogenetic distance",
           y = "Interaction type")

      #     by metabolism
plot_met_het = ggplot(all_het, aes(x=met_dist1, y=int))+
      geom_point(pch = 21, alpha = 0.8)+
      annotate(geom="text", x = 1.1, y = 1,
               label = paste("p-value =", round(summary(met_het)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_x_continuous(limits = c(0.4, 1.2), breaks = c(0.4, 0.6, 0.8, 1.0))+
      scale_y_continuous(limits = c(0, 1))+
      labs(title = "Heterogeneous\nenvironment",
           x = "Metabolic distance",
           y = "Interaction type")

plot_met_hom = ggplot(all_hom, aes(x=met_dist1, y=int))+
      geom_point(pch = 21, alpha = 0.8)+
      annotate(geom="text", x = 1.1, y = 1,
               label = paste("p-value =", round(summary(met_hom)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_x_continuous(limits = c(0.4, 1.2), breaks = c(0.4, 0.6, 0.8, 1.0))+
      scale_y_continuous(limits = c(0, 1))+
      labs(title = "Heterogeneous\nenvironment",
           x = "Metabolic distance",
           y = "Interaction type")








#     MODELS INTERACTION STRENGTH
qqnorm(unique(all$magnitude_het))
qqnorm(unique(all$magnitude_hom))
qqnorm(unique(all$met_dist1))
qqnorm(unique(all$phyl_dist))

# homogeneous env data - interaction strength
# null models
null_hom = lm(magnitude ~ 1, data = all_hom)
null_hom_res = lmer(magnitude ~ 1 + (1|resource), data = all_hom)
null_hom_comm = lmer(magnitude ~ 1 + (1|community), data = all_hom)
null_hom_res_comm = lmer(magnitude ~ 1 + (1|resource) + (1|community), data = all_hom)
null_hom_res_comm2 = lmer(magnitude ~ 1 + (1|resource/community), data = all_hom)
# models
res_hom = lmer(magnitude ~ resource + (1|community), data = all_hom)
rel_hom = lmer(magnitude ~ relation + (1|resource/community), data = all_hom)
rel_res_hom = lmer(magnitude ~ relation*resource + (1|community), data = all_hom)
phyl_hom = lmer(magnitude ~ phyl_dist + (1|resource/community), data = all_hom)
phyl_res_hom = lmer(magnitude ~ phyl_dist*resource + (1|community), data = all_hom)
met_hom = lmer(magnitude ~ met_dist1 + (1|resource/community), data = all_hom)
met_res_hom = lmer(magnitude ~ met_dist1*resource + (1|community), data = all_hom)
phyl_met_res = lmer(magnitude ~ met_dist1*phyl_dist*resource + (1|community), data = all_hom)

aic_hom = AICc(null_hom, null_hom_res, null_hom_comm, null_hom_res_comm, null_hom_res_comm2, #null models
               res_hom, rel_hom, rel_res_hom, phyl_hom, phyl_res_hom, met_hom, met_res_hom,
               phyl_met_res) %>% 
      mutate(model = rownames(.),
             dAIC = AICc - min(AICc)) %>% 
      arrange(dAIC) %>% 
      dplyr::select(model, df, AICc, dAIC)
aic_hom

anova(null_hom_res_comm2, rel_hom)
Anova(rel_hom, type = "III")
anova(null_hom_res_comm2, met_hom)
Anova(met_hom, type = "III")
anova(null_hom_res_comm2, phyl_hom)
Anova(phyl_hom, type = "III")

summary(null_hom_res_comm2)
qqnorm(resid(null_hom_res_comm2))
summary(rel_hom)
qqnorm(resid(rel_hom))
summary(met_hom)
qqnorm(resid(met_hom))
summary(phyl_hom)
qqnorm(resid(phyl_hom))

# heterogeneous env data - interaction strength
# null model
null_het = lm(magnitude ~ 1, data = all_het)
null_het_dpi = lmer(magnitude ~ 1 + (1|dpi), data = all_het)
null_het_comm = lmer(magnitude ~ 1 + (1|community), data = all_het)
null_het_dpi_comm = lmer(magnitude ~ 1 + (1|dpi) + (1|community), data = all_het)
null_het_dpi_comm2 = lmer(magnitude ~ 1 + (1|dpi/community), data = all_het)

rel_het = lm(magnitude ~ relation, data = all_het)
phyl_het = lm(magnitude ~ phyl_dist, data = all_het)
met_het = lm(magnitude ~ met_dist1, data = all_het)
phyl_met_het = lm(magnitude ~ met_dist1*phyl_dist, data = all_het)

rel_comm_het = lmer(magnitude ~ relation + (1|community), data = all_het)
phyl_comm_het = lmer(magnitude ~ phyl_dist + (1|community), data = all_het)
met_comm_het = lmer(magnitude ~ met_dist1 + (1|community), data = all_het)
phyl_met_comm_het = lmer(magnitude ~ met_dist1:phyl_dist + (1|community), data = all_het)

aic_het = AICc(null_het,null_het_dpi, null_het_comm, null_het_dpi_comm, null_het_dpi_comm2,
               rel_het, phyl_het, met_het, phyl_met_het,
               rel_comm_het, phyl_comm_het, met_comm_het, phyl_met_comm_het) %>% 
      mutate(model = rownames(.),
             dAIC = AICc - min(AICc)) %>% 
      arrange(dAIC) %>% 
      dplyr::select(model, df, AICc, dAIC)
aic_het

summary(null_het)
qqnorm(resid(null_het))
summary(rel_het)
qqnorm(resid(rel_het))
summary(met_het)
qqnorm(resid(met_het))
summary(phyl_het)
qqnorm(resid(phyl_het))

#     mean comparisons
emm_rel_hom = emmeans(rel_hom, specs = pairwise ~ relation, adjust = "mvt", type = "response")
emm_rel_hom
emm_res_hom = emmeans(res_hom, specs = pairwise ~ resource, adjust = "mvt", type = "response")
emm_res_hom
emm_rel_res_hom = emmeans(rel_res_hom, specs = pairwise ~ relation|resource, adjust = "mvt", type = "response")
emm_rel_res_hom

emm_rel_het = emmeans(rel_het, specs = pairwise ~ relation, adjust = "mvt")
emm_rel_het
emm_met_het = emmeans(met_het, specs = ~ met_dist1)
emm_met_het
emm_phyl_het = emmeans(phyl_het, specs = ~ phyl_dist, adjust = "mvt")
emm_phyl_het

#     confidence intervals
label_models = c("rel_hom", "phyl_hom", "met_hom", "rel_het", "phyl_het", "met_het")
df = data.frame()

for(i in 1:length(label_models)){
      data <- data.frame(confint(get(label_models[i])))[2,]
      data[, ncol(data) + 1] <- label_models[i]
      df <- rbind(df, data)
}

colnames(df) = c("2.5% CI", "97.5% CI", "Estimate", "Model")
df$Variable = rownames(df)
rownames(df) = NULL
df = df[, c(4, 5, 3, 1 ,2)]

df


#     PLOTS STRENGTH
#     by relations/life strategy   
ggplot(all, aes(x=relation, y = magnitude_het))+
      geom_point()+
      theme_bw()

ggplot(all, aes(x=relation, y = magnitude_hom))+
      geom_point()+
      theme_bw()

#     by phylogeny 
plot_str_phyl_het = ggplot(all_het, aes(x=phyl_dist, y=magnitude))+
      geom_point(pch = 21, alpha = 0.8)+
      #annotate(geom="text", x = 0.55, y = 1,
      #         label = paste("p-value =", round(summary(phyl_het)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_x_continuous(limits = c(0, 1.2), breaks = c(0.2, 0.4, 0.6,0.8, 1.0))+
      scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6))+
      labs(title = "Heterogeneous\nenvironment",
           x = "Phylogenetic distance",
           y = "Interaction strength")

plot_str_phyl_hom = ggplot(all_hom, aes(x=phyl_dist, y=magnitude))+
      geom_point(pch = 21, alpha = 0.8)+
      #annotate(geom="text", x = 0.55, y = 1,
      #         label = paste("p-value =", round(summary(phyl_hom)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_x_continuous(limits = c(0, 1.2), breaks = c(0.2, 0.4, 0.6,0.8, 1.0))+
      scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6))+
      labs(title = "Homogeneous\nenvironment",
           x = "Phylogenetic distance",
           y = "Interaction strength")

#     by metabolism
plot_str_met_het = ggplot(all_het, aes(x=met_dist1, y=magnitude))+
      geom_point(pch = 21, alpha = 0.8)+
      #annotate(geom="text", x = 1.1, y = 1,
      #         label = paste("p-value =", round(summary(met_het)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_x_continuous(limits = c(0.4, 1.2), breaks = c(0.4, 0.6, 0.8, 1.0))+
      scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6))+
      labs(title = "Heterogeneous\nenvironment",
           x = "Metabolic distance",
           y = "Interaction strength")

plot_str_met_hom = ggplot(all_hom, aes(x=met_dist1, y=magnitude))+
      geom_point(pch = 21, alpha = 0.8)+
      #annotate(geom="text", x = 1.1, y = 1,
      #         label = paste("p-value =", round(summary(met_hom)$coefficients$cond[2,4], digits = 2)))+
      theme_bw()+
      theme(text = element_text(size = 15))+
      scale_x_continuous(limits = c(0.4, 1.2), breaks = c(0.4, 0.6, 0.8, 1.0))+
      scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6))+
      labs(title = "Heterogeneous\nenvironment",
           x = "Metabolic distance",
           y = "Interaction strength")


# EXPORT PLOTS
pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/int_hom_het.pdf", width = 10, height = 5, useDingbats=FALSE)
plot_interfam + plot_intrafam + plot_layout(guides = 'collect')
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/hom_phyl.pdf", width = 10, height = 5, useDingbats=FALSE)
plot_phyl_hom + plot_str_phyl_hom
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/hom_met.pdf", width = 10, height = 5, useDingbats=FALSE)
plot_met_hom + plot_str_met_hom
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/het_phyl.pdf", width = 10, height = 5, useDingbats=FALSE)
plot_phyl_het + plot_str_phyl_het
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/het_met.pdf", width = 10, height = 5, useDingbats=FALSE)
plot_met_het + plot_str_met_het
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/competition in planta/analysis/int_str_hom_het.pdf", width = 10, height = 5, useDingbats=FALSE)
plot_int_env + plot_str_env + plot_layout(guides = "collect")
dev.off()

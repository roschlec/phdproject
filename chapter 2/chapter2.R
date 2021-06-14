###   CHAPTER 2
###   Chromatic bacteria — a broad host-range plasmid and chromosomal 
###   insertion toolbox for fluorescent protein expression in bacteria

##    DATA ANALYSIS AND REPRESENTATION
#     Libraries
library(tidyverse)
library(patchwork)
library(egg)
library(grid)


#     Figure 2.2. Normalised absorption and emission spectra
spectra = read.csv("~/Google Drive/Rudolf_MRE/THESIS/ch2:chromaticBac/Figures/spectra.csv", header = T) %>% 
      pivot_longer(cols = -1, names_to = "var", values_to = "val") %>% 
      na.omit() %>% 
      separate(col = var, into = c("type", "fp"))

col.pal = c("#1c0d82", "#23abc7", "#269a31", "#13491e", "#dede0b", "#f57f12", "#f70c0c", "#7f097d")
spectra$fp = factor(spectra$fp,
                    levels = c("mTagBFP2", "mTq2", "sGFP2", "mCl3", "sYFP2", "mO2", "mSc", "mCar"))
lab.fp = c("mTB2", "mTq2", "sGFP2", "mCl3", "sYFP2", "mO2", "mSc", "mCar")
lab.type = c("Abs" = "Absorption", "Em" = "Emission")

      #     Plots
spectra_plot_a = ggplot(spectra, aes(x = nm, y = val, color = fp))+
      facet_wrap(~ type, ncol = 1, labeller = labeller(type = lab.type))+
      geom_line(size = 1, alpha = 0.8)+
      scale_x_continuous(name = "Wavelength (nm)", limits = c(350,700), expand = c(0,0))+
      scale_y_continuous(name = "Normalised fluorescence intensity (%)")+
      theme_bw()+
      theme(legend.title = element_blank(),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 10))+
      scale_color_manual(values = col.pal, labels = lab.fp)

spectra_plot_b = ggplot(spectra, aes(x = nm, y = val, color = type))+
      facet_wrap(~ fp, ncol = 4)+
      geom_line(size = 1, alpha = 0.8)+
      scale_x_continuous(name = "Wavelength (nm)", limits = c(350,700), expand = c(0,0))+
      scale_y_continuous(name = "Normalised fluorescence intensity (%)")+
      theme_bw()+
      theme(legend.title = element_blank(),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 10))+
      scale_color_brewer(palette = "Set2", labels = c("Absorption", "Emission"))


#     FIGURE 2.3: Fluorescence intensity of E. coli DH5α cultures expressing mScarlet-I from different pMRE plasmid series
rfu = read.csv("~/Google Drive/Rudolf_MRE/THESIS/ch2:chromaticBac/Figures/fig2-3.csv", header = T)      
      #     dataset 1
rfu_1 = rfu %>% 
      filter(strain != "tn5-145" &
                   strain != "tn7-145")
rfu_1$strain = factor(rfu_1$strain, 
                      levels = c("wt", "pmre135", "pmre145", "pmre155", "pmre165"))
strain.lab1 = c("E. coli DH5a", "E. coli DH5a (pMRE135)", "E. coli DH5a (pMRE145)", "E. coli DH5a (pMRE155)", "E. coli DH5a (pMRE165)")

      #     dataset 2
rfu_2 = rfu %>% 
      filter(strain == "tn5-145" |
                   strain == "tn7-145" |
                   strain == "wt" |
                   strain == "pmre145") %>% 
      mutate(facet = case_when(
            rfu > 20000 ~ "B",
            TRUE ~ "A"
      ))
rfu_2$strain = factor(rfu_2$strain, 
                      levels = c("wt", "pmre145", "tn5-145", "tn7-145"))
strain.lab2 = c("E. coli DH5a", "E. coli DH5a (pMRE145)", "E. coli DH5a::Tn5-145", "E. coli DH5a::Tn7-145")

      #     Plots
plot2.3.a = ggplot(rfu_1, aes(x = strain, y = rfu))+
      geom_boxplot()+
      geom_point(alpha = 0.5)+
      theme_bw()+
      scale_x_discrete(labels = strain.lab1, name = "")+
      scale_y_continuous( name = "Fluorescence\nintensity (a.u.)",
                          limits = c(-100,50000))+
      theme(legend.title = element_blank(),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

#Transform the data onto the display scale
trans <- function(x){pmin(x,4000) + 0.1*pmax(x-4000,0)}
rfu_2$rfu2 <- trans(rfu_2$rfu)

yticks <- c(0, 2000, 4000, 30000, 50000)

plot2.3.b = ggplot(rfu_2, aes(x = strain, y = rfu2))+
      geom_boxplot()+
      geom_point(alpha = 0.5)+
      geom_rect(aes(xmin=1, xmax=3, ymin=4000, ymax=5000), fill="white") +
      theme_bw()+
      scale_x_discrete(labels = strain.lab2, name = "")+
      scale_y_continuous(name = "Fluorescence\nintensity (a.u.)", limits = c(-100,9000), 
                         breaks=trans(yticks), labels=yticks)+
      theme(legend.title = element_blank(),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

#     FIGURE 2.5: Effect of antibiotic pressure in growth and fluorescence of Pantoea eucalypti 299R (Pe299R)
#     Figure 2.5.a
c <- read.csv("~/Google Drive/Rudolf_MRE/results/ch1/chromatic bacteria paper/plasmid stability/cfu_pa299r.csv", header = T)
c$o.n.culture <- factor(c$o.n.culture,levels = c("NB", "Cm"))

sum <- c %>% 
      group_by(strain, o.n.culture) %>% 
      summarise(
            ave = mean(log),
            stdev = sd(log))

sum$strain <- factor(sum$strain,levels = c("Pa299R", "Pa299R(pMRE135)", "Pa299R::tn5-145", "Pa299R::tn7-145"))
sum$o.n.culture <- factor(sum$o.n.culture,levels = c("NB", "Cm"))
name.labs = c("Pe299R", "Pe299R\n(pMRE)", "Pe299R::\nMRE-Tn5", "Pe299R::\nMRE-Tn7")

cfu = ggplot(c, aes(x=strain, y=log, fill=o.n.culture))+
      geom_boxplot(data = sum, aes(x=strain, y = ave, fill = o.n.culture), alpha = 0.5)+
      geom_point(size = 2, pch = 21, position = position_dodge(width=0.78))+
      theme_bw()+
      theme(axis.title = element_text(size=10),
            axis.text = element_text(size=8),
            axis.text.x = element_text(size=10))+
      scale_fill_manual(name="Growth\nCondition", labels = c("NB", "NB+Cm"), values = c("#757575", "#e2e2e2"))+
      scale_y_continuous(limits = c(4,10), name = "Bacterial density\n(log10 CFU/mL)")+
      scale_x_discrete(labels = name.labs, name = "")

fix_cfu <- set_panel_size(cfu,
                          width  = unit(4, "in"),
                          height = unit(1, "in"))
grid.newpage()
grid.draw(fix_cfu)

#     Figure 2.b.a
d <- read.csv("~/Google Drive/Rudolf_MRE/results/ch1/chromatic bacteria paper/plasmid stability/data_pa299r.csv", header = T)
d$log <- log10(d$rfu_exp)

lab.name = c("Pe299R\n(pMRE)", "Pe299R::\nMRE-Tn5", "Pe299R::\nMRE-Tn7")
d$strain <- factor(d$strain,levels = c("Pa299R(pMRE)", "Pa299R::Tn5", "Pa299R::Tn7"))
q <- ggplot(d, aes(x=strain, y=log, fill=growthcondition))+
      geom_violin(alpha=0.5,
                  draw_quantiles = 0.5,
                  size=0.3)+
      geom_point(size = 0.5, alpha = 0.25, position = position_dodge(width=0.9))+
      theme_bw()+
      theme(axis.title = element_text(size=10),
            axis.text.x = element_text(size=10),
            axis.text = element_text(size=8))+
      scale_fill_manual(name="Growth\nCondition", labels = c("NB", "NB+Cm"), values = c("#757575", "#e2e2e2"))+
      labs(y = "log10 Normalised\nfluorescence\nintensity (a.u.)", x = "")+
      scale_x_discrete(labels = lab.name)
q

fix_q <- set_panel_size(q,
                        width  = unit(4, "in"),
                        height = unit(1, "in"))
grid.newpage()
grid.draw(fix_q)

      #     Stats 2.5.b
shapiro.test(d$log)
bartlett.test(log~interaction(strain,growthcondition), data=d)

m <- aov(log~strain*growthcondition, data=d)
summary(m)

pairwise.t.test(strain*growthcondition, log, p.adj = "bonferroni")
TukeyHSD(m)

kruskal.test(rfu_exp~growthcondition, data=subset(d, d$strain=="Pa299R(pMRE)"))
wilcox.test(rfu_exp~growthcondition, data=subset(d, d$strain=="Pa299R::Tn7"))


d2 <- subset(d, d$strain=="Pa299R::Tn7")
pairwise.wilcox.test(d2$rfu_exp, d2$growthcondition, p.adjust.method = "BH")

d2 <- subset(d, d$growthcondition=="NB")
pairwise.wilcox.test(d2$rfu_exp, d2$strain, p.adjust.method = "BH")


##    Figure 2.8: 
dat <- read.csv("~/Google Drive/Rudolf_MRE/results/ch1/chromatic bacteria paper/fluorescence single cell/data.csv", header = T)
dat$species = factor(dat$species, levels = unique(dat$species))

sp.lab = c("Ec (pMRE145)",
           "Ec::Tn5-145",
           "Ec::Tn7-145",
           "BradyL396::Tn5-145",
           "EaCFBP1430 (pMRE135)",
           "EaCFBP1430::Tn5-145",
           "EaCFBP1430::Tn7-145",
           "MethyL92::Tn5-145",
           "Pe299R (pMRE145)",
           "Pe299R::Tn5-145",
           "Pe299R::Tn7-145",
           "PcP3B5 (pMRE145)",
           "PcB3B5::Tn5-145",
           "PcB3B5::Tn7-145",
           "PssB728a (pMRE145)",
           "PssB728a::Tn5-145",
           "PssB728a::Tn7-145",
           "SmFR1::Tn5-145",
           "SmFR1::Tn7-145",
           "SpFA2 (pMRE135)",
           "SpFA2::Tn5-145")

singlecellplot = ggplot(dat, aes(x = species, y=milifl.exp))+
      geom_violin(alpha=0.3,
                  draw_quantiles = 0.5,
                  size=0.2)+
      geom_point(size = 0.5, alpha = 0.25, position = position_dodge(width=0.9))+
      theme_bw()+
      theme(axis.title = element_text(size=10),
            axis.text.x = element_text(size=10, angle = 45, hjust = 1),
            axis.text = element_text(size=8))+
      scale_y_continuous(name = "Normalised fluorescence\nintensity (milli a.u.)")+
      scale_x_discrete(labels = sp.lab, name = "")
fix_singlecellplot <- set_panel_size(singlecellplot,
                                     width  = unit(5, "in"),
                                     height = unit(1.5, "in"))
grid.newpage()
grid.draw(fix_singlecellplot)


###   SAVING PLOTS
#Figure 2.2.a
pdf(file = "~/Google Drive/Rudolf_MRE/THESIS/figure/illustrator_files/2-2a.pdf", width = 6, height = 3, useDingbats = FALSE)
spectra_plot_a
dev.off()
#Figure 2.2.b
pdf(file = "~/Google Drive/Rudolf_MRE/THESIS/figure/illustrator_files/2-2b.pdf", width = 6, height = 3, useDingbats = FALSE)
spectra_plot_b
dev.off()
#Figure2.3.a
pdf(file = "~/Google Drive/Rudolf_MRE/THESIS/figure/illustrator_files/2-3a.pdf", width = 3, height = 3, useDingbats = FALSE)
plot2.3.a
dev.off()
#Figure2.3.b
pdf(file = "~/Google Drive/Rudolf_MRE/THESIS/figure/illustrator_files/2-3b.pdf", width = 3, height = 3, useDingbats = FALSE)
plot2.3.b
dev.off()
#Figure2.5.a
pdf(file = "~/Google Drive/Rudolf_MRE/THESIS/figure/illustrator_files/2-5a.pdf", width = 6, height = 3, useDingbats = FALSE)
grid.newpage()
grid.draw(fix_cfu)
dev.off()
#Figure2.5.b
pdf(file = "~/Google Drive/Rudolf_MRE/THESIS/figure/illustrator_files/2-5b.pdf", width = 6, height = 3, useDingbats = FALSE)
grid.newpage()
grid.draw(fix_q)
dev.off()
#Figure2.8
pdf(file = "~/Google Drive/Rudolf_MRE/THESIS/figure/illustrator_files/2-8.pdf", width = 6, height = 3, useDingbats = FALSE)
grid.newpage()
grid.draw(fix_singlecellplot)
dev.off()

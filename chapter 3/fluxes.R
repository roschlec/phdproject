#fluxes FBA

#     libraries

library(g2f)
library(Matrix)
library(lattice)
library(glpkAPI)
library(sybil)
library(sybilSBML)
library(sybilDynFBA)
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#     functions
fx <- function(x){
      flx_frac = x/sum(x)}

#     set directory
setwd("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/fba/sbml/test")

#     load files
min <- read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/fba/media/min.csv", header = T, sep = ",") %>% 
      dplyr::select(cpd = compound, react_id, name, formula, charge) #read minimal medium formulation
leaf <- read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/fba/media/csources3.csv", header = T, sep = ",") %>% 
      dplyr::select(cpd, react_id, name, formula, charge)#read leaf-derived c-sources list
resources <- rbind(min, leaf)

cpd <- read_tsv("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/dev/Biochemistry/compounds.tsv")
cpd <- as.data.frame(cpd[,1:3])
cpd$cpd <- paste("EX", cpd$id, "e0", sep="_")
head(cpd)

spname <- read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/fba/spname.csv", header = F, sep = ",")
spname <- as.character(spname$V1)


#     FBA models
temp = list.files(pattern="*.sbml")
myfiles = lapply(temp, readSBMLmod)
myfiles = lapply(myfiles, changeObjFunc, "biomass2")
myex = lapply(myfiles, findExchReact)
myupt = lapply(myex, uptReact)

for (i in 1:length(myfiles)){
      myfiles[[i]] = changeBounds(myfiles[[i]], myex[[i]][myupt[[i]]], lb = 0)
      myfiles[[i]] = addExchReact(myfiles[[i]], met = c("cpd00012[e0]"), lb=c(-100))
      lowbnd(myfiles[[i]])[react_id(myfiles[[i]]) %in% resources$react_id] = -1000
}

mynewex = lapply(myfiles, findExchReact)
mynewupt = lapply(mynewex, uptReact)

fbas = lapply(myfiles, optimizeProb, algorithm = "fba")

df = data.frame(NULL)

for (i in 1:length(myfiles)){
      flx <- getFluxDist(fbas[[i]], mynewex[[i]])
      dframe <- data.frame(model = myfiles[[i]]@mod_name, 
                           flx)
      dframe$cpd = rownames(dframe)
      rownames(dframe) = NULL
      df = rbind(df, dframe)
      
}

head(df)

df$model = gsub(pattern = "*_SEED_model", replacement = "", x = df$model)

fluxes = inner_join(df, cpd, by = "cpd") %>% 
      mutate(dir = case_when(
                  flx > 0 ~ "exc",
                  flx == 0 ~ "bal",
                  flx < 0 ~ "upt"
      )) %>%
      filter(dir != "bal") %>% 
      filter(name != "H2O" & name != "H+") %>% 
      data.frame()
fluxes$model = factor(fluxes$model)
fluxes$flx = as.numeric(fluxes$flx)
head(fluxes)

uptflx <- fluxes %>% 
      filter(dir == "upt") %>% 
      mutate(flx = as.numeric(abs(flx))) %>% 
      group_by(model) %>% 
      mutate(fraction = fx(flx),
            fraction = case_when(
            fraction < 0.01 ~ 0,
            TRUE ~ fraction
      )) %>% 
      dplyr::select(model, name, fraction) %>%  
      pivot_wider(id_cols = name, names_from = model, values_from = fraction, values_fill = 0)

excflx <- fluxes %>% 
      filter(dir == "exc") %>% 
      mutate(flx = as.numeric(abs(flx))) %>% 
      group_by(model) %>% 
      mutate(fraction = fx(flx),
             fraction = case_when(
                   fraction < 0.01 ~ 0,
                   TRUE ~ fraction
             )) %>% 
      dplyr::select(model, name, fraction) %>%  
      pivot_wider(id_cols = name, names_from = model, values_from = fraction, values_fill = 0)



#Clustering
my_palette <- colorRamp2(c(0, .5, 1), c("#ffffff", "#40004B", "#000000"))


## UPTAKE FLUXES from main carbon sources
uptmatrix <- as.matrix(uptflx[, 2:ncol(uptflx)])
row.names(uptmatrix) = uptflx$name

uptmatrix = uptmatrix[rowSums(uptmatrix[, -1]) > 0,]

# Ward Hierarchical Clustering
uptd <- dist(t(uptmatrix), method = "euclidean") # distance matrix
uptfit <- hclust(uptd, method="ward.D") 
plot(uptfit) # display dendogram

hm_upt <- Heatmap(uptmatrix,
              col = my_palette,
              rect_gp = gpar(col= "#f2f2f2"),
              heatmap_legend_param = list(at = c(0, 0.25, 0.5, 0.75, 1),
                                          title = "Fraction of\nUptake Resources",
                                          legend_height = unit(20, "mm")),
              cluster_rows = T,
              row_dend_reorder = T,
              show_row_dend = F,
              cluster_columns = uptfit,
              column_dend_height=unit(20,"mm"))

hm_upt

## EXCRETION FLUXES from main carbon sources
excmatrix <- as.matrix(excflx[, 2:ncol(excflx)])
row.names(excmatrix) = excflx$name

excmatrix = excmatrix[rowSums(excmatrix[, -1]) > 0,]

# Ward Hierarchical Clustering
excd <- dist(t(excmatrix), method = "euclidean") # distance matrix
excfit <- hclust(excd, method="ward.D") 
plot(excfit) # display dendogram

hm_exc <- Heatmap(excmatrix,
                  col = my_palette,
                  rect_gp = gpar(col= "#f2f2f2"),
                  heatmap_legend_param = list(at = c(0, 0.25, 0.5, 0.75, 1),
                                              title = "Fraction of\nBy-products",
                                              legend_height = unit(20, "mm")),
                  cluster_rows = T,
                  row_dend_reorder = F,
                  show_row_dend = F,
                  column_dend_reorder = T,
                  cluster_columns = excfit,
                  column_dend_height=unit(20,"mm"))

hm_exc


## SAVE FILES
### plots
pdf("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/heatmap_upt.pdf", height = 6, width = 8, useDingbats = F)
draw(hm_upt)
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/heatmap_exc.pdf", height = 8, width = 8, useDingbats = F)
draw(hm_exc)
dev.off()

## save files
write.csv(fluxes, "~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/fluxes.csv")
write.csv(uptflx, "~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/flx_uptake.csv")
write.csv(excflx, "~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/flx_excretion.csv")
write.csv(melt(as.matrix(uptd), varnames = c("sp1", "sp2"), value.name = "uptake"),
          "~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_uptake.csv")
write.csv(melt(as.matrix(excd), varnames = c("sp1", "sp2"), value.name = "excreted"),
          "~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_excretion.csv")



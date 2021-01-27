#FBA

#libraries

library(g2f)
library(Matrix)
library(lattice)
library(glpkAPI)
library(sybil)
library(sybilSBML)
library(sybilDynFBA)
library(gplots)
library(dplyr)
library(tidyr)
library(plyr)
library(magrittr)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)

#set directory
setwd("~/Rudolf drive/Rudolf_MRE/results/ch2/metabolism/fba/")
setwd("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/fba/")


#     load files
min <- read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/fba/media/min.csv", header = T, sep = ",") %>% 
      dplyr::select(cpd = compound, react_id, name, formula, charge) #read minimal medium formulation
leaf <- read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/fba/media/csources.csv", header = T, sep = ",") %>% 
      dplyr::select(cpd, react_id, name, formula, charge, type)#read leaf-derived c-sources list
resources <- rbind(min, leaf)
leafc <- as.character(leaf$react_id)

cpd <- read_tsv("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/dev/Biochemistry/compounds.tsv")
cpd <- as.data.frame(cpd[,1:3])
cpd$react_id <- paste("EX", cpd$id, "e0", sep="_")
colnames(cpd)[3] = "cpd"
head(cpd)

spname <- read.csv("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/fba/spname.csv", header = F, sep = ",")
spname <- as.character(spname$V1)


#     FBA models
setwd("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/fba/sbml/test")
temp = list.files(pattern="*.sbml")
myfiles = lapply(temp, readSBMLmod)
myfiles = lapply(myfiles, changeObjFunc, "biomass2")
myex = lapply(myfiles, findExchReact)
myupt = lapply(myex, uptReact)

#min
for (i in 1:length(myfiles)){
      myfiles[[i]] = changeBounds(myfiles[[i]], myex[[i]][myupt[[i]]], lb = 0)
      myfiles[[i]] = addExchReact(myfiles[[i]], met = c("cpd00012[e0]"), lb=-100)
      lowbnd(myfiles[[i]])[react_id(myfiles[[i]]) %in% min$react_id] = -1000
}

mynewex = lapply(myfiles, findExchReact)
mynewupt = lapply(mynewex, uptReact)

fbas = lapply(myfiles, optimizeProb, algorithm = "fba")


#resources
maxgrowth <- data.frame(matrix(NA, ncol=length(myfiles), nrow=2))
row.names(maxgrowth) = c("all", "min")

cprof <- data.frame(matrix(NA, ncol=length(myfiles), nrow=length(leafc)), row.names = leaf$react_id)

for (i in 1:length(myfiles)){
      for (j in 1:length(leafc)){
            myfiles[[i]] = changeBounds(myfiles[[i]], myex[[i]][myupt[[i]]], lb = 0)
            myfiles[[i]] = addExchReact(myfiles[[i]], met = c("cpd00012[e0]"), lb=-100)
            lowbnd(myfiles[[i]])[react_id(myfiles[[i]]) %in% min$react_id] = -Inf
            lowbnd(myfiles[[i]])[react_id(myfiles[[i]]) == leaf$react_id[j]] <- -Inf
            opti <- optimizeProb(myfiles[[i]], algorithm = "fba", retOptSol=F)
            cprof[j,i] <- opti$obj
            colnames(cprof)[i] = myfiles[[i]]@mod_name
      }
      
}

for (i in 1:length(myfiles)){
            myfiles[[i]] = changeBounds(myfiles[[i]], myex[[i]][myupt[[i]]], lb = 0)
            myfiles[[i]] = addExchReact(myfiles[[i]], met = c("cpd00012[e0]"), lb=-100)
            lowbnd(myfiles[[i]])[react_id(myfiles[[i]]) %in% min$react_id] = -Inf
            optinull <- optimizeProb(myfiles[[i]], algorithm = "fba", retOptSol=F)
            lowbnd(myfiles[[i]])[react_id(myfiles[[i]]) %in% resources$react_id] <- -Inf
            opti <- optimizeProb(myfiles[[i]], algorithm = "fba", retOptSol=F)
            maxgrowth[1,i] <- opti$obj
            maxgrowth[2,i] <- optinull$obj
            colnames(maxgrowth)[i] = myfiles[[i]]@mod_name
            }

maxgrowth2 = maxgrowth %>% 
      mutate(resources = row.names(.)) %>% 
      pivot_longer(cols = -resources) %>% 
      pivot_wider(id_cols = name, names_from = resources, values_from = value)

biomass = cprof %>%
      mutate(react_id = row.names(.)) %>% 
      pivot_longer(cols = -react_id) %>% 
      inner_join(., maxgrowth2, by = "name") %>% 
      arrange(name) %>% 
      mutate(yield = (value - min)/(all - min)) %>% 
      inner_join(., cpd, by = "react_id")
biomass$name = gsub(pattern = "*_SEED_model", replacement = "", x = biomass$name)


#Clustering

# palette
my_palette <- colorRamp2(c(0, 1), c("#ffffff", "#40004B"))

# matrix dataset
mat <- biomass %>% 
      pivot_wider(id_cols = cpd, names_from = name, values_from = yield) %>% 
      data.frame()
row.names(mat) <- mat$cpd
mat <- as.matrix(mat[,2:ncol(mat)])
mat <- mat[rowSums(mat[, -1]) > 0,]

# Ward Hierarchical Clustering
matd <- dist(t(mat), method = "euclidean") # distance matrix
fit <- hclust(matd, method="ward.D") 
plot(fit) # display dendogram

hm <- Heatmap(mat,
              col = my_palette,
              rect_gp = gpar(col= "#f2f2f2"),
              heatmap_legend_param = list(at = c(0, 0.25, 0.5, 0.75, 1),
                                          title = "Biomass Yield",
                                          legend_height = unit(20, "mm")),
              cluster_rows = T,
              row_dend_reorder = T,
              show_row_dend = F,
              cluster_columns = fit,
              column_dend_height=unit(20,"mm"))
      
hm

# Save files!
write.csv(biomass,"~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/biomass_fba.csv")
write.csv(melt(as.matrix(matd), varnames = c("sp1", "sp2"), value.name = "fba"),
          "~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/dist_fba.csv")
pdf("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/heatmap_fba.pdf", height = 8, width = 8, useDingbats = F)
draw(hm)
dev.off()


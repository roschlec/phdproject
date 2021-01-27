#nutrient overlaps

#nutrient overlap
library(EcoSimR)
library(spaa)
library(tidyverse)
library(reshape2)


setwd("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final")
no <- read.csv("flx_uptake.csv", header = T, sep = ",")           #read file (col=species, row=nutrient)
rownames(no) <- no$name                                           #add row names to data frame
no <- data.matrix(no[,3:ncol(no)])                                #create numeric matrix 
no1 <- t(no)                                                      #transpose matrix

n <- as.matrix(niche.overlap(no, method = c("pianka")))           #calculate pianka nutrient overlap
n1 <- melt(as.matrix(n))                                          #tidy up output, ordered as sp1-sp2-value

#plotting
ggplot(n1, aes(x=Var1, y=Var2, fill=value))+
      geom_tile()

# Get lower triangle of the correlation matrix
get_lower_tri<-function(n){
      n[upper.tri(n)] <- NA
      return(n)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(n){
      n[lower.tri(n)]<- NA
      return(n)
}

upper_tri <- get_upper_tri(n)
upper_tri

#melt matrix
melted_n <- melt(upper_tri, na.rm=T)

reorder_n <- function(n){
      # Use correlation between variables as distance
      dd <- as.dist((1-n)/2)
      hc <- hclust(dd)
      n <-n[hc$order, hc$order]
}

# Reorder the  matrix
n <- reorder_n(n)
upper_tri <- get_upper_tri(n)

# Melt the  matrix
melted_n <- melt(upper_tri, na.rm = TRUE)

#heatmap
my_palette1 <- colorRamp2(c(0, 1), c("#ffffff", "#40004B"))

p <- ggplot(data=melted_n, aes(Var2, Var1, fill=value))+
      geom_tile(color="white")+
      scale_fill_gradient2(low = "#ffffff", high = "#40004B")+
      coord_fixed()+
      labs(y = NULL, x = NULL)+
      #geom_text(aes(Var2, Var1, label=sprintf("%0.1f", round(value, digits=2))), color="white", size=3)+
      scale_y_discrete(position = "right")+
      theme(
            axis.text.x = element_text(angle=90,size=8, hjust=1, vjust=0.5, color="black"),
            axis.text.y = element_text(size=8, color="black"),
            axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            legend.text = element_text(size=8),
            legend.title = element_text(size=10),
            legend.justification = c(0.5,0.5))+
      guides(fill=guide_colorbar(barwidth = 1, 
                                 barheight = 5,
                                 title.position = "top", 
                                 direction = "vertical"))
p

#null model
myModel <- niche_null_model(speciesData = t(no),  #dataframe
                            algo = "ra2",         #algorithm
                            metric = "pianka",    #nutrient overlap index
                            suppressProg = F,
                            saveSeed = T,
                            nReps = 10000)         #bootstrap
summary(myModel)
plot(myModel)
nicheplot <- plot(myModel, type="niche")

sim <- myModel$Sim

max.len <- max(length(sim), length(n1$value))
sim = c(sim, rep(NA, max.len - length(sim)))
obs = c(n1$value, rep(NA, max.len - length(n1$value)))
data <- data.frame(cbind(sim, obs))
data2 <- melt(data, na.rm = T)
head(data2)

data2$variable <- factor(data2$variable,levels = c("obs", "sim"))
pnm <- ggplot(data2, aes(x=variable, y=value))+
      geom_boxplot()+
      geom_point(alpha = 0.05)+
      annotate("text", x = 2.2, y = 1, label = "italic(p) == 0.004", 
               parse = T, size = 4, color = "black")+
      theme_bw()+
      theme(text = element_text(size=15))+
      scale_x_discrete(labels = c("Sample", "Random"))+
      labs(y = expression(paste("Pianka's Index")),
           x = NULL)
pnm


## SAVE PLOTS
pdf("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/nutrientoverlap.pdf", height = 12, width = 12, useDingbats = F)
p
dev.off()

pdf("~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/nullmodel_NO.pdf", height = 5, width = 3, useDingbats = F)
pnm
dev.off()

write.csv(melt(as.matrix(n), varnames = c("sp1", "sp2"), value.name = "no"),
          "~/Google Drive/Rudolf_MRE/results/ch2/metabolism/final/pianka.csv")

####HOMOGENEOUS DATA
setwd("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/")
data = read.csv("matrix.csv", header = T)

intra = data %>% 
      filter(mix == "meth85_meth92" |
                   mix == "mr01_meth85" |
                   mix == "mr01_meth92" |
                   mix == "smfr1_spfa2")
inter = data %>% 
      filter(mix != "meth85_meth92" |
                   mix != "mr01_meth85" |
                   mix != "mr01_meth92" |
                   mix != "smfr1_spfa2")

data1ry = data %>% 
      filter(exp == 1 & group == "r_y") %>%
      filter(s1 == "spfa2") %>% 
      pivot_wider(id_cols = resource, names_from=s2, values_from = a1, values_fill = 0) %>% 
      arrange(resource) %>% 
      dplyr::select(resource, meth85, meth92, mr01, smfr1) %>% 
      data.frame()
row.names(data1ry) = data1ry$resource
data1ry = data1ry[,-1]
mat1ry = as.matrix(data1ry)
mds1ry = monoMDS(dist(data1ry))

data2ry = data %>% 
      filter(exp == 2 & group == "r_y") %>% 
      filter(s1 == "spfa2") %>% 
      pivot_wider(id_cols = resource, names_from=s2, values_from = a1, values_fill = 0) %>% 
      arrange(resource) %>% 
      dplyr::select(resource, meth85, meth92, mr01, smfr1) %>% 
      data.frame()
row.names(data2ry) = data2ry$resource
data2ry = data2ry[,-1]
mat2ry = as.matrix(data2ry)
mds2ry = monoMDS(dist(data2ry))

data1yr = data %>% 
      filter(exp == 1 & group == "y_r") %>% 
      filter(s1 == "spfa2") %>% 
      pivot_wider(id_cols = resource, names_from=s2, values_from = a1, values_fill = 0) %>% 
      arrange(resource) %>% 
      dplyr::select(resource, meth85, meth92, mr01, smfr1) %>% 
      data.frame()
row.names(data1yr) = data1yr$resource
data1yr = data1yr[,-1]
mat1yr = as.matrix(data1yr)
mds1yr = monoMDS(dist(data1yr))

data2yr = data %>% 
      filter(exp == 2 & group == "y_r") %>% 
      filter(s1 == "spfa2") %>% 
      pivot_wider(id_cols = resource, names_from=s2, values_from = a1, values_fill = 0) %>% 
      arrange(resource) %>% 
      dplyr::select(resource, meth85, meth92, mr01, smfr1) %>% 
      data.frame()
row.names(data2yr) = data2yr$resource
data2yr = data2yr[,-1]
mat2yr = as.matrix(data2yr)
mds2yr = monoMDS(dist(data2yr))



#     Procrustes test
mdslist = list(exp1_yr = mds1yr, 
               exp1_ry = mds1ry,
               exp2_yr = mds2yr,
               exp2_ry = mds2ry)
names(mdslist)

df = data.frame(mat1 = character(0),
                mat2 = character(0),
                ss = numeric(0),
                t0 = numeric(0),
                p = numeric(0),
                stringsAsFactors = FALSE)
df2 = data.frame(NULL)

for (i in 1:length(mdslist)){
      for (j in 1:length(mdslist)){
            if (i < j){
                  pro = protest(X = mdslist[[i]], Y = mdslist[[j]])
                  df[1,1] = names(mdslist)[i]
                  df[1,2] = names(mdslist)[j]
                  df[1,3] = pro$ss
                  df[1,4] = pro$t0
                  df[1,5] = pro$signif
                  df2 = rbind(df2, df)
            }
      }
}

df2




# Heatmaps
lgd2 = Legend(title = "Fold\nChange",
              title_gap = unit(5, "mm"),
              title_gp = gpar(fontsize=12),
              col=my_palette2,
              labels_gp = gpar(fontsize = 10),
              at = c(-3, 3),
              legend_height = unit(2, "cm"),
              border = TRUE)
h1ry = Heatmap(mat1ry,
               rect_gp = gpar(col= "#7e7e7e"),
               show_heatmap_legend = F,
               col = my_palette2,
               row_title = "Resource",
               cluster_rows = F,
               row_names_side = "left",
               column_title = "Exp1_RY",
               cluster_columns = F,
               column_names_side = "top",
               cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.2f", mat1ry[i, j]), x, y, gp = gpar(fontsize = 15))},
               width = unit(8, "cm"), height = unit(8, "cm"))
h2ry = Heatmap(mat2ry,
               rect_gp = gpar(col= "#7e7e7e"),
               show_heatmap_legend = F,
               col = my_palette2,
               row_title = "Resource",
               cluster_rows = F,
               row_names_side = "left",
               column_title = "Exp2_RY",
               cluster_columns = F,
               column_names_side = "top",
               cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.2f", mat2ry[i, j]), x, y, gp = gpar(fontsize = 15))},
               width = unit(8, "cm"), height = unit(8, "cm"))
h1yr = Heatmap(mat1yr,
               rect_gp = gpar(col= "#7e7e7e"),
               show_heatmap_legend = F,
               col = my_palette2,
               row_title = "Resource",
               cluster_rows = F,
               row_names_side = "left",
               column_title = "Exp1_YR",
               cluster_columns = F,
               column_names_side = "top",
               cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.2f", mat1yr[i, j]), x, y, gp = gpar(fontsize = 15))},
               width = unit(8, "cm"), height = unit(8, "cm"))
h2yr = Heatmap(mat2yr,
               rect_gp = gpar(col= "#7e7e7e"),
               show_heatmap_legend = F,
               col = my_palette2,
               row_title = "Resource",
               cluster_rows = F,
               row_names_side = "left",
               column_title = "Exp2_YR",
               cluster_columns = F,
               column_names_side = "top",
               cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.2f", mat2yr[i, j]), x, y, gp = gpar(fontsize = 15))},
               width = unit(8, "cm"), height = unit(8, "cm"))

hlist = h1ry + h1yr + h2ry + h2yr
hlist

#pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/heatmap_meth85.pdf", 
#pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/heatmap_meth92.pdf", 
#pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/heatmap_mr01.pdf", 
pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/heatmap_spfa2.pdf", 
#pdf("~/Google Drive/Rudolf_MRE/results/ch2/in vitro growth/competition_assay/comp/heatmap_smfr1.pdf",
    width = 15, height = 5, useDingbats=FALSE)
#draw(hlist, column_title = "M85", padding = unit(c(5, 20, 5, 5), "mm"))
#draw(hlist, column_title = "M92", padding = unit(c(5, 20, 5, 5), "mm"))
#draw(hlist, column_title = "Mr", padding = unit(c(5, 20, 5, 5), "mm"))
draw(hlist, column_title = "Sp", padding = unit(c(5, 20, 5, 5), "mm"))
#draw(hlist, column_title = "Sm", padding = unit(c(5, 20, 5, 5), "mm"))
draw(lgd2, x = unit(.5, "cm"), y = unit(.5, "cm"), just = c("left", "bottom"))
dev.off()

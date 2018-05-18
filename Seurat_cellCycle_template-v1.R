# Load libraries and user functions ----------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("Seurat", "dplyr", "colorRamps")

# Load libraries
ipak(packages)


# Define variables --------------------------------------------------------

all_data.files <- c("LSK" = "LSKm2/outs/filtered_gene_bc_matrices/mm10/",
                    "CMP" = "CMPm2/outs/filtered_gene_bc_matrices/mm10/",
                    "MEP" = "MEPm/outs/filtered_gene_bc_matrices/mm10/",
                    "GMP" = "GMPm/outs/filtered_gene_bc_matrices/mm10/")
ProjectName<-"msAggr_ccycle"
genome<- "mm10"
output_prefix <-"20180509_msAggr"
max_pcs <-10
# resolution_list <-c(0.6, 0.8, 1.0, 1.5, 2.0, 2.5)
resolution_list <-0.6
my_palette<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99', "gray40", "black","maroon3", "tan4", "paleturquoise4", "mediumorchid", "moccasin", "lemonchiffon", "darkslategrey", "lightsalmon1", "plum1", "seagreen", "firebrick4", "khaki2", primary.colors(80))

# Define common objects ---------------------------------------------------

color_pallet<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99', "gray40", "black","maroon3", "tan4", "paleturquoise4", "mediumorchid", "moccasin", "lemonchiffon", "darkslategrey", "lightsalmon1", "plum1", "seagreen", "firebrick4", "khaki2", "darkorange3", "darksalmon", "darkslategray3", "cadetblue", "blue2", "aquamarine4", "antiquewhite2", "aliceblue", "cornflowerblue", "chartreuse")

my_transparent_dot<-rgb(0,0,0, alpha = 0)

s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes

# Import data -------------------------------------------------------------

my_object.data<-Read10X(all_data.files)
my_object<- CreateSeuratObject(raw.data = my_object.data, min.cells = 3, min.genes = 200, project = ProjectName)
my_object<-NormalizeData(my_object, normalization.method = "LogNormalize", scale.factor = 10000)
my_object<-FindVariableGenes(my_object,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
my_object<-ScaleData(my_object, vars.to.regress = "nUMI")













# Cell cycle analysis -------------------------------------------------------




# Shows cycle of genes in each cluster
dim14cycle<-CellCycleScoring(dim14, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = F)



# Regress out ALL markers of cell cycle -----------------------------------

# This method means you no longer see difference between cycling cells and non-cycling cells
head(dim14cycle@meta.data)
RidgePlot(dim14cycle, features.plot = c("Pcna", "Top2a", "Mcm6", "Mki67"), nCol = 2)
dim14cycle<-RunPCA(dim14cycle, pc.genes=c(s.genes, g2m.genes), do.print=TRUE)
PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "Phase")
PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "orig.ident", pt.size = 1)
PCAPlot(dim14cycle, cols.use = c(my_transparent_dot, my_transparent_dot, "black", my_transparent_dot), group.by = "orig.ident", pt.size = 1)
PCAPlot(dim14cycle, cols.use = c(my_transparent_dot, "black", my_transparent_dot, my_transparent_dot), group.by = "orig.ident", pt.size = 1)
PCAPlot(dim14cycle, cols.use = c("black", my_transparent_dot, my_transparent_dot, my_transparent_dot), group.by = "orig.ident", pt.size = 1)
PCAPlot(dim14cycle, cols.use = c(my_transparent_dot, my_transparent_dot, my_transparent_dot, "black"), group.by = "orig.ident", pt.size = 1)
PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "res.1.5", pt.size = 0.5)


png(filename = "dim14res15_cellcyclebyID_ridgeplot.png", width = 800, height = 800)
RidgePlot(dim14cycle, features.plot = c("Pcna", "Top2a", "Mcm6", "Mki67"), nCol = 2)
dev.off()

png(filename = "dim14res15_cellcyclebyID_res.png", width = 800, height = 800)
PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "res.1.5", pt.size = 0.5)
dev.off()

png(filename = "dim14res15_cellcyclebyID_phase.png", width = 800, height = 800)
PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "Phase", pt.size = 0.5)
dev.off()

png(filename = "dim14res15_cellcyclebyID_orig.ident.png", width = 800, height = 800)
PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "orig.ident", pt.size = 0.5)
dev.off()

png(filename = "dim14res15_cellcyclebyID_LSKcycle.png", width = 800, height = 800)
PCAPlot(dim14cycle, cols.use = c(my_transparent_dot, my_transparent_dot, "black", my_transparent_dot), group.by = "orig.ident", pt.size =  0.5)
dev.off()

png(filename = "dim14res15_cellcyclebyID_GMPcycle.png", width = 800, height = 800)
PCAPlot(dim14cycle, cols.use = c(my_transparent_dot, "black", my_transparent_dot, my_transparent_dot), group.by = "orig.ident", pt.size =  0.5)
dev.off()

png(filename = "dim14res15_cellcyclebyID_CMPcycle.png", width = 800, height = 800)
PCAPlot(dim14cycle, cols.use = c("black", my_transparent_dot, my_transparent_dot, my_transparent_dot), group.by = "orig.ident", pt.size =  0.5)
dev.off()

png(filename = "dim14res15_cellcyclebyID_MEPcycle.png", width = 800, height = 800)
PCAPlot(dim14cycle, cols.use = c(my_transparent_dot, my_transparent_dot, my_transparent_dot, "black"), group.by = "orig.ident", pt.size =  0.5)
dev.off()


for (i in 0:max(as.integer(dim14cycle@meta.data$res.1.5))){
  vector_length<-length(unique(dim14cycle@meta.data$res.1.5))
  my_color_vector<-rep("black", vector_length)
  my_color_vector<-replace(my_color_vector, i+1, "red")
  my_plot_order<-c(max(as.integer(dim14cycle@meta.data$res.1.5)): min(as.integer(dim14cycle@meta.data$res.1.5)))
  # png(filename = paste("dim14res15_cellcyclebyID_res", i, ".png", sep = ""), width = 800, height = 800)
  myplot<-PCAPlot(dim14cycle, cols.use = my_color_vector, group.by = "res.1.5", pt.size = 0.5, plot.order = my_plot_order, do.return=TRUE)+theme(plot.background = element_rect(fill="transparent"), panel.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))
  # dev.off()
  ggsave(myplot, filename = paste("dim14res15_cellcyclebyID_res", i, ".png", sep = ""), bg="transparent", width = 4, height = 4, units = "in")
}



for (i in 0:max(as.integer(dim14cycle@meta.data$res.1.5))){
  print(i)
  vector_length<-length(unique(dim14cycle@meta.data$res.1.5))
  my_color_vector<-rep(my_transparent_dot, vector_length)
  my_color_vector<-replace(my_color_vector, i+1, "black")
  my_plot_order<-c(max(as.integer(dim14cycle@meta.data$res.1.5)): min(as.integer(dim14cycle@meta.data$res.1.5)))
  myplot<-PCAPlot(dim14cycle, cols.use = my_color_vector, group.by = "res.1.5", pt.size = 0.5, plot.order = my_plot_order, do.return=TRUE)
  png(filename = paste("test/dim14res15_cellcyclebyres_trans", i, ".png", sep = ""), width = 800, height = 800)
  try(plot(myplot))
  dev.off()
  # ggsave(myplot, filename = paste("dim14res15_cellcyclebyID_res_trans", i, ".png", sep = ""), bg="transparent", width = 4, height = 4, units = "in")
}


library(ggplot2)
myplot<-PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "res.1.5")+theme(plot.background = element_rect(fill="transparent"), panel.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent")) + geom_point(size=10)
ggsave(myplot, filename = "dim14res15_cellcyclebyID_res.png", bg="transparent", width = 4, height = 4, units = "in")

myplot<-(PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "Phase"))+theme(plot.background = element_rect(fill="transparent"), panel.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))
ggsave(myplot, filename = "dim14res15_cellcyclebyID_phase.png", bg="transparent", width = 4, height = 4, units = "in")

myplot<-(PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "orig.ident", pt.size = 1))+theme(plot.background = element_rect(fill="transparent"), panel.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))
ggsave(myplot, filename = "dim14res15_cellcyclebyID_orig.ident.png", bg="transparent", width = 4, height = 4, units = "in")

myplot<-(PCAPlot(dim14cycle, cols.use = c(my_transparent_dot, my_transparent_dot, "black", my_transparent_dot), group.by = "orig.ident", pt.size = 1))+theme(plot.background = element_rect(fill="transparent"), panel.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))
ggsave(myplot, filename = "dim14res15_cellcyclebyID_LSKcycle.png", bg="transparent", width = 4, height = 4, units = "in")

myplot<-(PCAPlot(dim14cycle, cols.use = c(my_transparent_dot, "black", my_transparent_dot, my_transparent_dot), group.by = "orig.ident", pt.size = 1))+theme(plot.background = element_rect(fill="transparent"), panel.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))
ggsave(myplot, filename = "dim14res15_cellcyclebyID_GMPcycle.png", bg="transparent", width = 4, height = 4, units = "in")

myplot<-(PCAPlot(dim14cycle, cols.use = c("black", my_transparent_dot, my_transparent_dot, my_transparent_dot), group.by = "orig.ident", pt.size = 1))+theme(plot.background = element_rect(fill="transparent"), panel.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))
ggsave(myplot, filename = "dim14res15_cellcyclebyID_CMPcycle.png", bg="transparent", width = 4, height = 4, units = "in")

myplot<-(PCAPlot(dim14cycle, cols.use = c(my_transparent_dot, my_transparent_dot, my_transparent_dot, "black"), group.by = "orig.ident", pt.size = 1))+theme(plot.background = element_rect(fill="transparent"), panel.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))
ggsave(myplot, filename = "dim14res15_cellcyclebyID_MEPcycle.png", bg="transparent", width = 4, height = 4, units = "in")

myplot<-(PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "res.1.5", pt.size = 0.5))+theme(plot.background = element_rect(fill="transparent"), panel.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))
ggsave(myplot, filename = , width = 800, height = 800)









myplot<-PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "orig.ident", pt.size = 1, do.return=TRUE)+theme(plot.background = element_rect(fill="transparent"), panel.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))

png(filename = "test.png")
# PCAPlot(dim14cycle, cols.use = color_pallet, group.by = "orig.ident", pt.size = 1, do.return=TRUE)+theme(plot.background = element_rect(fill="transparent"))
plot(myplot)
dev.off()

ggsave(myplot, filename = "test2.png", bg="transparent", width = 4, height = 4, units = "in")


# Shows cycle of entire population
# dim14cycle_v2<-CellCycleScoring(dim14, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = T)
# head(dim14cycle_v2@meta.data)
# RidgePlot(dim14cycle_v2, features.plot = c("Pcna", "Top2a", "Mcm6", "Mki67"))
# dim14cycle_v2<-RunPCA(dim14cycle_v2, pc.genes=c(s.genes, g2m.genes), do.print=TRUE)
# PCAPlot(dim14cycle_v2)
# 


dim14_scaledforcellcycle<-ScaleData(dim14cycle, vars.to.regress=c("S.Score", "G2M.Score"), display.progress = TRUE)
dim14_scaledforcellcycle<-RunPCA(dim14_scaledforcellcycle, pc.genes = dim14@var.genes, genes.print = 10)





# Regress out STAGES of cell cycle ----------------------------------------

# This means you see a difference between cycling and non-cycling, but lose clustering based on things like S vs G2M stage

head(dim14cycle@meta.data)


dim14cycle@meta.data$CC.Difference<-dim14cycle@meta.data$S.Score - dim14cycle@meta.data$G2M.Score
dim14cycle<-ScaleData(dim14cycle, vars.to.regress = "CC.Difference", display.progress = TRUE)
dim14cycle<-RunPCA(dim14cycle, pc.genes = dim14cycle@var.genes, genes.print = 10)




dim14cycle.projectpca<-ProjectPCA(dim14cycle)
# saveRDS(dim14cycle.projectpca, file = "mouse_projectpca.rds")
png(filename = "dim14cycle_VizPCA.png", width = 800, height = 800)
VizPCA(dim14cycle.projectpca, use.full = T, font.size = 1, nCol = 3, pcs.use = 1:6)
dev.off()

png(filename = "dim14cycle_ElbowPlot.png", width = 800, height = 800)
PCElbowPlot(dim14cycle.projectpca)
dev.off()


dim14cycle.clustered<-FindClusters(dim14cycle.projectpca, dims.use = 1:14, resolution = 1.5)


####################Run tSNE and find markers####################

dim14cycle.tsne<-RunTSNE(dim14cycle.clustered, dims.use = 1:14)   
TSNEPlot(dim14cycle.tsne, do.label = T, colors.use = color_pallet, pt.size = 0.5, label.size = 5)
saveRDS(dim14cycle.tsne, file="dim14cycle_tsne.rds")


png(filename="dim14cycle_tsne-noLables.png", width=800, heigh=800)
TSNEPlot(dim14cycle.tsne, do.label = F, colors.use = color_pallet, pt.size = 0.5)
dev.off()

png(filename="dim14cycle_tsne-Lables.png", width=800, heigh=800)
TSNEPlot(dim14cycle.tsne, do.label = T, colors.use = color_pallet, pt.size = 0.5, label.size = 5)
dev.off()

png(filename="dim14cycle_tsne-Lables_byID.png", width=800, heigh=800)
TSNEPlot(dim14cycle.tsne, do.label = TRUE, label.size = 5, pt.size=0.5, colors.use = color_pallet, group.by = "orig.ident")
dev.off()

png(filename="dim14cycle_tsne-noLables_byID.png", width=800, heigh=800)
TSNEPlot(dim14cycle.tsne, do.label = F, label.size = 5, pt.size=0.5, colors.use = color_pallet, group.by = "orig.ident")
dev.off()

dim14cycle.tsne<-RunTSNE(dim14cycle.clustered, reduction.use = "pca", dims.use = 1:14, dim.embed=3, reduction.name = "tsne")
TSNEPlot(dim14cycle.tsne, do.label = T, colors.use = color_pallet, pt.size = 0.5, label.size = 5)


dim14cycle.snn<-BuildSNN(dim14cycle.tsne, dims.use = 1:14)
ValidateClusters(dim14cycle.snn, pc.use = 1:14)

dim14cycle.markers<-FindAllMarkers(dim14cycle.tsne, only.pos = F, min.pct = .25, thresh.use = .25)
save(dim14cycle.markers, file="dim14cycle_markers.Robj")

dim14cycle.differential_markers_top100 <- dim14cycle.markers %>% group_by(cluster) %>% top_n(100, avg_diff)
write.table(dim14cycle.differential_markers_top100, file = "~/Desktop/10XGenomicsData/dim14cycle_top100_differentiatlMarkers.txt", sep="\t", quote=F)
dim14cycle.differential_markers_all <- dim14cycle.markers %>% group_by(cluster)
write.table(dim14cycle.differential_markers_all, file = "~/Desktop/10XGenomicsData/dim14cycle_all_differentiatlMarkers.txt", sep="\t", quote=F)

save.image(file = "dim14cycle_tsne.RData")


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

all_data.files <- c("Con3" = "Con3/outs/filtered_gene_bc_matrices/GRCh38/",
                    "DBA1" = "DBA1_2/outs/filtered_gene_bc_matrices/GRCh38/",
                    "DBA3" = "DBA3_2/outs/filtered_gene_bc_matrices/GRCh38/",
                    "DBA4" = "DBA4_2/outs/filtered_gene_bc_matrices/GRCh38/")
ProjectName<-"bmDBA"
genome<- "GRCh38"
Run_mito_filter = FALSE
output_prefix <-"20180510_bmDBA"
max_pcs <-3
resolution_list <-c(0.6, 0.8, 1.0, 1.5, 2.0, 2.5)
my_palette<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99', "gray40", "black","maroon3", "tan4", "paleturquoise4", "mediumorchid", "moccasin", "lemonchiffon", "darkslategrey", "lightsalmon1", "plum1", "seagreen", "firebrick4", "khaki2", primary.colors(20))

# Load data ---------------------------------------------------------------

# Alter output name based on filter parameters

if(Run_mito_filter == TRUE){
  output_prefix<-paste(output_prefix, "_mt", sep = "")
}


con3.data<-Read10X(c("Con3" = "Con3/outs/filtered_gene_bc_matrices/GRCh38/"))
con3<-CreateSeuratObject(raw.data=con3.data, project="bmDBA", min.cells=3, min.genes=200, names.field = 1, names.delim = "_")
# pancreas tutorial has a filter step I don't understand here
con3<-NormalizeData(con3)
con3<-FindVariableGenes(con3, do.plot = FALSE, display.progress = TRUE)
con3<-ScaleData(con3)
con3@meta.data$sample<-"Control3"

dba1.data<-Read10X(c("DBA1" = "DBA1_2/outs/filtered_gene_bc_matrices/GRCh38/"))
dba1<-CreateSeuratObject(raw.data=dba1.data, project="bmDBA", min.cells=3, min.genes=200, names.field = 1, names.delim = "_")
# pancreas tutorial has a filter step I don't understand here
dba1<-NormalizeData(dba1)
dba1<-FindVariableGenes(dba1, do.plot = FALSE, display.progress = TRUE)
dba1<-ScaleData(dba1)
dba1@meta.data$sample<-"DBA1"

dba3.data<-Read10X(c("DBA3" = "DBA3_2/outs/filtered_gene_bc_matrices/GRCh38/"))
dba3<-CreateSeuratObject(raw.data=dba3.data, project="bmDBA", min.cells=3, min.genes=200, names.field = 1, names.delim = "_")
# pancreas tutorial has a filter step I don't understand here
dba3<-NormalizeData(dba3)
dba3<-FindVariableGenes(dba3, do.plot = FALSE, display.progress = TRUE)
dba3<-ScaleData(dba3)
dba3@meta.data$sample<-"DBA3"

dba4.data<-Read10X(c("DBA4" = "DBA4_2/outs/filtered_gene_bc_matrices/GRCh38/"))
dba4<-CreateSeuratObject(raw.data=dba4.data, project="bmDBA", min.cells=3, min.genes=200, names.field = 1, names.delim = "_")
# pancreas tutorial has a filter step I don't understand here
dba4<-NormalizeData(dba4)
dba4<-FindVariableGenes(dba4, do.plot = FALSE, display.progress = TRUE)
dba4<-ScaleData(dba4)
dba4@meta.data$sample<-"DBA4"

object_list<-list(con3, dba1, dba3, dba4)
genes.use<-c()
for(i in 1:length(object_list)){
    genes.use<-c(genes.use, head(rownames(object_list[[i]]@hvg.info), 1000))
}
length(genes.use)

genes.use <-names(which(table(genes.use)>1))
length(genes.use)
for(i in 1:length(object_list)){
  genes.use<-genes.use[genes.use %in% rownames(object_list[[i]]@scale.data)]
}

length(genes.use)
integrated_dba<-RunMultiCCA(object_list, genes.use = genes.use, num.ccs = 15)
MetageneBicorPlot(integrated_dba, grouping.var = "sample", dims.eval = 1:15)
png(filename = "20180516_bmDBA_BicorePlot.png", width = 800, height = 800)
MetageneBicorPlot(integrated_dba, grouping.var = "sample", dims.eval = 1:15)
dev.off()

png(filename = "20180516_bmDBA_DimHeatmap.png", width = 800, height = 800)
DimHeatmap(integrated_dba, reduction.type = "cca", cells.use = 500, dim.use = 1:15, do.balanced = TRUE)
dev.off()



integrated_dba<-CalcVarExpRatio(integrated_dba, reduction.type = "pca", grouping.var = "sample", dims.use = 1:12)
integrated_dba<-SubsetData(integrated_dba, subset.name = "var.ratio.pca", accept.low = 0.5)
integrated_dba<-AlignSubspace(integrated_dba, reduction.type = "cca", grouping.var = "sample", dims.align = 1:12)
integrated_dba<-FindClusters(integrated_dba, reduction.type = "cca.aligned", dims.use = 1:12, save.SNN = TRUE, resolution = 0.4)
integrated_dba<-RunTSNE(integrated_dba, reduction.use = "cca.aligned", dims.use = 1:12)

png(filename = "20180516_bmDBA_tsne.png", width = 800, height = 800)
TSNEPlot(integrated_dba, do.label = TRUE, colors.use = my_palette, label.size = 10)
dev.off()

png(filename = "20180516_bmDBA_tsneID.png", width = 800, height = 800)
TSNEPlot(integrated_dba, do.label = TRUE, group.by = "sample", colors.use = my_palette, label.size = 10)
dev.off()
















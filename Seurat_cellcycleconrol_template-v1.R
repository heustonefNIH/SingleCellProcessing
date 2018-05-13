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

all_data.files <- c("Con3" = "Con3/outs/filtered_gene_bc_matrices/GRCh38/","DBA1b" = "DBA1_2/outs/filtered_gene_bc_matrices/GRCh38/","DBA3b" = "DBA3_2/outs/filtered_gene_bc_matrices/GRCh38/","DBA4b" = "DBA4_2/outs/filtered_gene_bc_matrices/GRCh38/") 
ProjectName<-"OxfordDBA"
genome<- "GRCh38"
Run_mito_filter = FALSE
output_prefix <-"20180509_bmDBA"
max_pcs <-40
resolution_list <-c(0.6, 0.8, 1.0, 1.5, 2.0, 2.5)
my_palette<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99', "gray40", "black","maroon3", "tan4", "paleturquoise4", "mediumorchid", "moccasin", "lemonchiffon", "darkslategrey", "lightsalmon1", "plum1", "seagreen", "firebrick4", "khaki2", primary.colors(20))


# Define cc genes ---------------------------------------------------------

#List of human-formatted cc genes now preloaded into Seurat as cc.genes$s.genes and cc.genes$g2m.genes
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes


# Load data ---------------------------------------------------------------

my_object.data<-Read10X(all_data.files)
my_object<- CreateSeuratObject(raw.data = my_object.data, min.cells = 3, min.genes = 200, project = ProjectName)

VlnPlot(my_object, features.plot = c("nGene", "nUMI"), nCol = 2)
GenePlot(my_object, gene1 = "nUMI", gene2 = "nGene")

png(paste(output_prefix, "_VlnStats.png", sep = ""), width = 800, height = 800)
VlnPlot(my_object, features.plot = c("nGene", "nUMI"), nCol = 2)
dev.off()
png(paste(output_prefix, "_GenePlots.png", sep = ""), width = 800, height = 800)
GenePlot(my_object, gene1 = "nUMI", gene2 = "nGene")
dev.off()



# Normalize and scale -----------------------------------------------------

my_object<-NormalizeData(my_object, normalization.method = "LogNormalize", scale.factor = 10000)
my_object<-FindVariableGenes(my_object,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(my_object@var.genes)
my_object<-ScaleData(my_object, vars.to.regress = "nUMI")


# PCA ---------------------------------------------------------------------

my_object<-RunPCA(my_object, pc.genes = my_object@var.genes, pcs.compute = max_pcs, do.print = TRUE, pcs.print = 1:5, genes.print = 10)
















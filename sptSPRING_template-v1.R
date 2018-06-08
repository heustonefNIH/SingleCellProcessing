# Acknowledge -------------------------------------------------------------

# Code provided by Supat via Beth Psaila. Have to figure out Supat's full name :)

# Load libraries and user functions ----------------------------------------------------------

# source("https://bioconductor.org/biocLite.R")
# biocLite("scran")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("Matrix", "RColorBrewer", "matrixStats", "irlba", "factoextra", "Rtsne", "ggplot2", "gridExtra", "pracma", "cccd", "vegan", "threejs", "scran", "igraph", "statmod", "cluster", "randomcoloR", "dplyr", "cowplot")

# Load libraries
ipak(packages)



use.components<-NULL
output_prefix<-"MEPm_supat_v1"
data.dir<-"~/Desktop/10XGenomicsData/MEPm"
project.folder<-"MEPm_sptSpring_v1"
gene_list <-NULL
max_dims<-12
n.max.cluster<-NULL
k<-NULL #for Kmeans clustering

# Source Supat’s scripts --------------------------------------------------

source("~/Desktop/10XGenomicsData/SupatSpring/SingleCellUtils.R")
source("~/Desktop/10XGenomicsData/SupatSpring/SingleCell10xRNASeq.R")
source("~/Desktop/10XGenomicsData/SupatSpring/SingleCell10xClasses.R")
source("~/Desktop/10XGenomicsData/SupatSpring/SingleCell10xGenerics.R")
source("~/Desktop/10XGenomicsData/SupatSpring/SingleCell10xPlots.R")



# Import data -------------------------------------------------------------

my_obj<-new("scRNAseq10X", dir_path_10x_matrix=data.dir, sample_uniq_id = output_prefix, project.name = output_prefix, project.folder = project.folder)

load_matrices_from_cellranger(my_obj)
my_obj@umi_count@Dimnames[1]<-lapply(my_obj@umi_count@Dimnames[1:length(my_obj@umi_count@Dimnames[1])], toupper)

process_cellsAnnotation(my_obj)

plotCells_Annotation(my_obj)
plotUMIs_Vs_Detected_genes(my_obj)
plotUMI_of_HouseKeepingGenes(object = my_obj, selected_top_cells = 100) #Found this function while surfing the code and threw it in here #20180607: not working right now, thoguh

# UMI and QC --------------------------------------------------------------

png(filename = paste(output_prefix, "_UMIvsgenes.png", sep = ""), height = 800, width = 800)
plotCells_Annotation(my_obj)
plotUMIs_Vs_Detected_genes(my_obj)
# p2<-plotUMI_of_HouseKeepingGenes(my_obj) #Found this function while surfing the code and threw it in here
# plot_grid(p1, p2)
dev.off()

filter_cells_and_genes(my_obj,min_UMIs=1000,max_UMIs=25000,min_detected_genes=1000,max_detected_genes=5000,max_percent_mito=9)
normalizedUMIs(my_obj,use.scaled.factor = T)

# Find variable gene expression -------------------------------------------

getVariableGenesFor10Xdata(my_obj, mean_expr_cutoff = 0.1, disp_zscore_cutoff = 0.1)
plotVariableGenes(my_obj)

png(filename = paste(output_prefix, "_plotVariableGenes.png", sep = ""), height = 800, width = 800)
plotVariableGenes(my_obj)
dev.off()


# Run PCA -----------------------------------------------------------------

runPCA(my_obj, use.components = 50)

plot_screeplot(my_obj)

png(filename = paste(output_prefix, "_plotPCA.png", sep = ""), height = 800, width = 800)
plot_screeplot(my_obj)
dev.off()


# Run TSNE ----------------------------------------------------------------

runFast_TSNE(my_obj, n.dims.use = max_dims)
plot_tsne_label_by_qc(my_obj)
plot_tsne_label_by_genes(my_obj, gene_list = toupper(c("Gata1", "Gata2", "Pf4", "Itga2b", "Cnrip1", "Vwf", "Kit")))

png(filename = paste(output_prefix, "_runTSNE-qcLabel.png", sep = ""), height = 800, width = 800)
plot_tsne_label_by_qc(my_obj)
dev.off()

png(filename = paste(output_prefix, "_runTSNE-geneLabel.png", sep = ""), height = 800, width = 800)
plot_tsne_label_by_genes(my_obj, gene_list = toupper(c("Gata1", "Gata2", "Pf4", "Itga2b", "Cnrip1", "Vwf", "Kit")))
dev.off()


# Density clustering ------------------------------------------------------

DensityClustering(my_obj, n.max.cluster = 11)
plot_tsne_label_by_density_cluster(my_obj, point.size = 2)
png(filename = paste(output_prefix, "_densityClustering.png", sep = ""), height = 800, width = 800)
plot_tsne_label_by_density_cluster(my_obj, point.size = 2)
dev.off()



# KNN-Graph ---------------------------------------------------------------

runKNN_Graph(my_obj, use.components = 12)

saveRDS(my_obj, paste(output_prefix, "_runKNN.rds"))

s.layout<-get_knn_graph.layout(my_obj)
KMeansClusteringOnGraph(my_obj, k=8)
kmeans.info<-get_knn_graph.kmeans.cluster(my_obj)

my.cols<-kmeans.info$kmeans_cluster_color
my.cluster<-kmeans.info$kmeans_cluster


g<-set_vertex_attr(get_knn_graph.graph(my_obj), "color", value = my.cols)
graphjs(g, vertex.size=0.4, edge.color = "gray45", vertex.label = my.cluster, layout = s.layout, fpl = 300)

saveRDS(my_obj, paste(output_prefix, "kclst", k, ".rds", sep = ""))














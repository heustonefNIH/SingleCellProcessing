# Load libraries and user functions ----------------------------------------------------------

# install_github("efheuston/MyPlotTrajectoryPackage")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("Seurat", "plyr", "dplyr", "colorRamps", "monocle", "stringr", "MyPlotTrajectoryPackage")

# Load libraries
ipak(packages)

###Concerned about importing nomalized data incorrectly from Seurat object

# Define variables --------------------------------------------------------

seurat_object_filename<-"20180601_bmDBA_cca_dim21res2.5_tsne.rds"
perform_expression_filtering <- TRUE
color_by_seurat_res = TRUE
order_by_seurat_varGenes = FALSE
UMI_bounded_filtering<- "upper" #Can be "upper", "lower", "both", "none" (or NULL)
cca_variables<- "~nGene + nUMI + orig.ident"
num_dim <- 21
max_components <- 2

# Define script-specific functions ----------------------------------------

# define function to create pngs
png_plotFunction<-function(plot_to_make = x, filename = y, height = 800, width = 800){ # default window is 800x800
  png(filename = filename, height = height, width = width)
  plot_to_make
  dev.off()
}

# define function to cycle through multiple seurat-defined resolutions (if applicable)
res_cycle_plot_clusters<-function(x, y){
  if(y == "cluster"){
  png_plotFunction(plot_cell_clusters(cca_object, 1, 2, color = x), 
                   filename = paste(output_prefix, "_clstr-", x,".png", sep = ""),
                   height = 1600, 
                   width = 1600)
  png_plotFunction(plot_cell_clusters(cca_object, 1, 2, color = x) + 
                     facet_wrap(x), 
                   filename = paste(paste(output_prefix, "_clstr-", x,"FACET.png", sep = "")),
                   height = 1600, 
                   width = 1600)
  }
  if(y == "trajectory")
    png_plotFunction()
                   }

output_prefix<-gsub(x=seurat_object_filename, pattern = ".rds", replacement = "")
output_prefix<-gsub(x = output_prefix, pattern =  "_tsne", replacement =  "")

# specify_root_state<-function(cds){
#   LSK_pops<-table(pData(cds)$State, pData(cds)$Type[,"LSK"])
#   as.numeric(names(LSK_pops)[which(LSK_pops==max(LSK_pops))])
# }

# Load Seurat object as CDS -----------------------------------------------

if(grep(x = seurat_object_filename, pattern = "rds$", ignore.case = T)){
  seurat_object<-readRDS(seurat_object_filename)
} else if (grep(x = seurat_object_filename, pattern = "Robj$", ignore.case = T)){
  seurat_object<-get(load(seurat_object_filename))
} else {
  print("I don't recognize this type of file")
}

seurat_varGenes<-seurat_object@var.genes

if(color_by_seurat_res == TRUE){
  seurat_res<-colnames(seurat_object@meta.data)[grep(x=colnames(seurat_object@meta.data), pattern = "^res")]
}

# Define base colour palette ----------------------------------------------

my_palette<-c("#cb4bbe",
                 "lightskyblue",
                 "grey37",
                 "#53ba48",
                 "moccasin",
                 "#dcca44",
                 "#502d70",
                 "#afd671",
                 "#cb4471",
                 "#69da99",
                 "#d44732",
                 "#6fd3cf",
                 "#5d2539",
                 "#cfd0a0",
                 "blue",
                 "#d2883b",
                 "#6793c0",
                 "#898938",
                 "#c98ac2",
                 "yellow",
                 "#c4c0cc",
                 "#7d3d23",
                 "#00a5ff",
                 "#d68d7a",
                 "#a2c1a3")

# Initialize monocle object -----------------------------------------------

monocle_object<-importCDS(seurat_object, import_all = FALSE)
remove(seurat_object)
gc()

monocle_object<-estimateSizeFactors(monocle_object)
monocle_object<-estimateDispersions(monocle_object)


# Filter low-quality cells ------------------------------------------------

if(perform_expression_filtering == TRUE){
  monocle_object<-detectGenes(monocle_object, min_expr = 0.1)
  expressed_genes<-row.names(subset(fData(monocle_object), num_cells_expressed >=10))
}

# can create a subset of valid cells using columns from pData
# might want to include qc info like Mapped.Fragments, or set valid_cells based on nUMI
# valid_cells<-row.names(subset(pData(monocle_object), nUMI > someNumber))
# monocle_object<-monocle_object[,valid_cells]
 
pData(monocle_object)$Total_nUMI<-Matrix::colSums(exprs(monocle_object))
monocle_object<-monocle_object[,pData(monocle_object)$Total_UMIs < 1e6]

lower_bound<-10^(mean(log10(pData(monocle_object)$Total_nUMI)) - 2*sd(log10(pData(monocle_object)$Total_nUMI)))
upper_bound<-10^(mean(log10(pData(monocle_object)$Total_nUMI)) + 2*sd(log10(pData(monocle_object)$Total_nUMI)))

# Plot UMI distributions and filter --------------------------------------------------

png_plotFunction(qplot(Total_UMIs, data = pData(cca_object), color = orig.ident, geom = "density") + 
                   geom_vline(xintercept = upper_bound) + 
                   geom_vline(xintercept = lower_bound),
                 filename = paste(output_prefix, "_nUMI_2SD_byID.png", sep = ""))

png_plotFunction(qplot(res.0.6, data = pData(cca_object), geom = "density") + 
                   geom_vline(xintercept = upper_bound) + 
                   geom_vline(xintercept = lower_bound),
                 filename = paste(output_prefix, "_nUMI_2SD_byCluster.png", sep = ""))

if(grep("upper", UMI_bounded_filtering, ignore.case = TRUE)){
  monocle_object<-monocle_object[,pData(monocle_object)$Total_nUMI < upper_bound]
  monocle_object<-detectGenes(monocle_object, min_expr = 0.1)
} else {
  if(grep("both", UMI_bounded_filtering, ignore.case = TRUE)){
  monocle_object<-monocle_object[,pData(monocle_object)$Total_nUMI < upper_bound & pData(monocle_object)$Total_nUMI > lower_bound]
  monocle_object<-detectGenes(monocle_object, min_expr = 0.1)
  }
} else {
  if(UMI_bounded_filtering = NULL | UMI_bounded_filtering == "none"){
  print("Not filtering based on UMI counts")
  monocle_object<-detectGenes(monocle_object, min_expr = 0.1)
  }
}



# Cluster cells without marker genes --------------------------------------

monocle_object<-detectGenes(monocle_object)
disp_table<-dispersionTable(monocle_object)

monocle_clustering_genes<-subset(disp_table, mean_expression >= 0.1)
monocle_object<-setOrderingFilter(monocle_object, ordering_genes = monocle_clustering_genes$gene_id)

png_plotFunction(plot_ordering_genes(monocle_object), filename = paste(output_prefix, "_PlotOrderingGenes.png", sep = ""))
png_plotFunction(plot_pc_variance_explained(monocle_object, return_all = FALSE), filename = paste(output_prefix, "_PCVarXplained.png", sep = ""))

savehistory(file=paste(output_prefix, ".Rhistory", sep = ""))
save.image(file=paste(output_prefix, ".RData", sep = ""))


# Perform clustering ------------------------------------------------------

monocle_object<-reduceDimension(monocle_object, 
                                max_components = max_components, 
                                num_dim = num_dim, 
                                reduction_method = 'tSNE', 
                                residualModelFormulaStr = cca_variables, 
                                verbose = TRUE)
monocle_object<-clusterCells(monocle_object)

png_plotFunction(plot_cell_clusters(cca_object, 1, 2, color = "orig.ident") + 
                   facet_wrap(~res.0.6), 
                 filename = paste(output_prefix, "_clstr-orig.identFACET.png", sep = ""),
                 height = 1600, 
                 width = 1600)

png_plotFunction(plot_cell_clusters(cca_object, 1, 2, color = "orig.ident"), 
                 filename = paste(output_prefix, "_clstr-orig.ident.png", sep = ""),
                 height = 1600, 
                 width = 1600)


if(color_by_seurat_res = TRUE){
  lapply(seurat_res[1:length(seurat_res)], res_cycle_plot_clusters )
}
savehistory(file=paste(output_prefix, ".Rhistory", sep = ""))
save.image(file=paste(output_prefix, ".RData", sep = ""))


# Perform trajectory analysis with monocle-defined genes ------------------

clustering_DEG_monocle<- differentialGeneTest(monocle_object[expressed_genes,], fullModelFormulaStr='~Cluster', cores=detectCores()-1)
#current default is to take top 1000 clustering_DEG_monocle genes

monocleBased_orderedgenes<-row.names(clustering_DEG_monocle)[order(clustering_DEG_monocle$qval)][1:1000]
monocle_object<-setOrderingFilter(monocle_object, ordering_genes=monocleBased_orderedgenes)
monocle_object<-reduceDimension(monocle_object, method="DDRTree")
monocle_object<-orderCells(monocle_object)

# HERE we skip defining the root state because we want to leave that as an open question for now


# Plot cell trajectories --------------------------------------------------

try(
)

png(filename="plotCellTrajectory_noGMstate-norm.png", width=800, height = 800)
my_plot_cell_trajectory(cca_monocleOrdered, cell_size = 0.5, point_colors = my_palette)
dev.off()
png(filename="plotCellTrajectory_noGMstate-res.png", width=800, height = 800)
my_plot_cell_trajectory(cca_monocleOrdered, color_by = "res.0.6", cell_size = 0.5, point_colors = my_palette)
dev.off()
png(filename="plotCellTrajectory_noGMstate-orig.ident.png", width=800, height = 800)
my_plot_cell_trajectory(cca_monocleOrdered, color_by = "orig.ident", cell_size = 0.5, point_colors = my_palette)
dev.off()
png(filename="plotCellTrajectory_noGMstate-orig.identFACET.png", width=800, height = 800)
my_plot_cell_trajectory(cca_monocleOrdered, color_by = "orig.ident", cell_size = 0.5, point_colors = my_palette)+facet_wrap(~orig.ident)
dev.off()
png(filename="plotCellTrajectory_noGMstate-resFACET.png", width=800, height = 800)
my_plot_cell_trajectory(cca_monocleOrdered, color_by = "res.0.6", cell_size = 0.5, point_colors = my_palette)+facet_wrap(~res.0.6)
dev.off()
png(filename="plotCellTrajectory_noGMstate-resFACETorig.png", width=800, height = 800)
my_plot_cell_trajectory(cca_monocleOrdered, color_by = "res.0.6", cell_size = 0.5, point_colors = my_palette)+facet_wrap(~orig.ident)
dev.off()
png(filename="plotCellTrajectory_noGMstate-orig.identFACETres.png", width=800, height = 800)
my_plot_cell_trajectory(cca_monocleOrdered, color_by = "orig.ident", cell_size = 0.5, point_colors = my_palette)+facet_wrap(~res.0.6)
dev.off()





  # Create color pallete
  new_length <-0
  if(length(unique(my_object@ident)) > length(my_palette)){
    new_length<-length(unique(my_object@ident)) - length(my_palette)
  }
  new_length
  my_palette<-c(my_palette, primary.colors(new_length))

# 
# 
# WHEN DRAWING IMAGES, set up so can rotate through list of res's
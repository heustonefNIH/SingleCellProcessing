# Load libraries and user functions ----------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("Seurat", "dplyr", "colorRamps", "monocle", "stringr", "plyr")

# Load libraries
ipak(packages)


###Concerned about importing nomalized data incorrectly from Seurat object

# Define variables --------------------------------------------------------
setwd("~/Desktop/10XGenomicsData/")
seurat_object_filename<-"testdim5res0.6_tsne.rds"
perform_filtering<-NULL #currently hard-coded for upper bound only




my_palette<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99', "gray40", "black","maroon3", "tan4", "paleturquoise4", "mediumorchid", "moccasin", "lemonchiffon", "darkslategrey", "lightsalmon1", "plum1", "seagreen", "firebrick4", "khaki2", primary.colors(80))

# Clear environment and define functions ----------------------------------

# rm(list = ls())
get_correct_root_state<-function(cds){
  T0_counts<-table(pData(cds)$State, pData(cds)$Type[,"LSK"])
  as.numeric(names(T0_counts)[which(T0_counts==max(T0_counts))])
}
output_prefix<-gsub(x=seurat_object_filename, pattern = ".rds", replacement = "")
output_prefix<-gsub(x = output_prefix, pattern =  "_tsne", replacement =  "")

# Load Seurat object as CDS -----------------------------------------------

if(grep(x = seurat_object_filename, pattern = "rds$", ignore.case = T)){
  seurat_object<-readRDS(seurat_object_filename)
} else if (grep(x = seurat_object_filename, pattern = "Robj$", ignore.case = T)){
  seurat_object<-get(load(seurat_object_filename))
} else {
  print("I don't recognize this type of object file")
}

monocle_object<-importCDS(seurat_object, import_all = TRUE)
monocle_object<-setOrderingFilter(monocle_object, seurat_object@var.genes)
monocle_object@phenoData@varMetadata$labelDescription[grep(x = row.names(monocle_object@phenoData@varMetadata), pattern = "^res")]<-"seurat_res"
monocle_object<-estimateSizeFactors(monocle_object)
monocle_object<-estimateDispersions(monocle_object)


# Filter low-quality cells ------------------------------------------------

median(seurat_object@data)
monocle_object<-detectGenes(monocle_object, min_expr = 0.1)
print(head(fData(monocle_object)))
expressed_genes<-row.names(subset(fData(monocle_object), num_cells_expressed >=10))



# can create a subset of valid cells using columns from pData
# might want to include qc info like Mapped.Fragments, or set valid_cells based on nUMI
# valid_cells<-row.names(subset(pData(monocle_object), nUMI > someNumber))
# monocle_object<-monocle_object[,valid_cells]

pData(monocle_object)$Total_nUMI<-Matrix::colSums(exprs(monocle_object))
range(monocle_object@phenoData@data$Total_nUMI)

lower_bound<-10^(mean(log10(pData(monocle_object)$Total_nUMI)) - 2*sd(log10(pData(monocle_object)$Total_nUMI)))
upper_bound<-10^(mean(log10(pData(monocle_object)$Total_nUMI)) + 2*sd(log10(pData(monocle_object)$Total_nUMI)))
qplot(Total_nUMI, data = pData(monocle_object), color = res.0.6, geom = "density") + scale_colour_manual(values = my_palette) + geom_vline(xintercept = upper_bound) + geom_vline(xintercept = lower_bound)

png(filename = paste(output_prefix, "_TotalnUMI-upandlow.png"), height = 800, width = 800)
qplot(Total_nUMI, data = pData(monocle_object), color = res.0.6, geom = "density") + scale_colour_manual(values = my_palette) + geom_vline(xintercept = upper_bound) + geom_vline(xintercept = lower_bound)
dev.off()
png(filename = paste(output_prefix, "_TotalnUMI-uponly.png"), height = 800, width = 800)
qplot(Total_nUMI, data = pData(monocle_object), color = res.0.6, geom = "density") + scale_colour_manual(values = my_palette) + geom_vline(xintercept = upper_bound)
dev.off()



# Will filter based on cells < upper_bound ONLY

monocle_object<-monocle_object[,pData(monocle_object)$Total_nUMI < upper_bound]
monocle_object<-detectGenes(monocle_object, min_expr = 0.1)
head(dispersionTable(monocle_object))




# Cluster cells without marker genes --------------------------------------

disp_table<-dispersionTable(monocle_object)
head(seurat_object@var.genes)

filter_list<-subset(disp_table, mean_expression >=0.1)
filter_list<-subset(filter_list, filter_list$gene_id %in% seurat_object@var.genes)
monocle_object<-setOrderingFilter(monocle_object, filter_list$gene_id)
plot_ordering_genes(monocle_object)

plot_pc_variance_explained(monocle_object, return_all = FALSE, max_components = 20)







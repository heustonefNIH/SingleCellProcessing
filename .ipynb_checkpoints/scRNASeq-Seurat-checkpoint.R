# Run analysis on single cell RNA data
# Sample data is from [PandDB](https://hpap.pmacs.upenn.edu/)


# Sample inclusion --------------------------------------------------------




##NB: 

# Global parameters -------------------------------------------------------

rnaProject <- "PancT2D_allSmpls"
regression.vars <- c("sequencerID", "SampleSex", "SampleAge")
cum.var.thresh <- 90
resolution <- 0.5
comp.type <- "biowulf" # one of macbookPro, biowulf, or workPC
do.sctransform <- "each" # one of FALSE, each, pooled

## infrequently modified
do.doubletFinder <- TRUE
run.jackstraw <- FALSE
min.cells <- 3
min.features <- 200
doublet.var.thresh <- 90
predicted.doubletRate <- 0.05

# Directories -------------------------------------------------------------

if(comp.type == "macbookPro"){
	rna.dir <- "/Users/heustonef/Desktop/Obesity/scRNA/"
	path_to_data <- "/Users/heustonef/Desktop/PancDB_data/scRNA_noBams"
	sourceable.functions <- list.files(path = "/Users/heustonef/OneDrive-NIH/SingleCellMetaAnalysis/GitRepository/scMultiomics_MetaAnalysis/RFunctions", pattern = "*.R$", full.names = TRUE)
	metadata.location <- "/Users/heustonef/OneDrive-NIH/SingleCellMetaAnalysis/GitRepository/scMultiomics_MetaAnalysis/"
} else if(comp.type == "biowulf"){
	rna.dir <- "/data/CRGGH/heustonef/hpapdata/cellranger_scRNA/"
	path_to_data <- "/data/CRGGH/heustonef/hpapdata/cellranger_scRNA/scRNA_transfer"
	sourceable.functions <- list.files(path = "/data/CRGGH/heustonef/hpapdata/RFunctions/", pattern = "*.R", full.names = TRUE)
	metadata.location <- "/data/CRGGH/heustonef/hpapdata/"
}

# Load libraries ----------------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)

##load local functions
invisible(sapply(sourceable.functions, source))

# Load data ---------------------------------------------------------------

try(setwd(rna.dir), silent = TRUE)
writeLines(capture.output(sessionInfo()), paste0(rnaProject, "_sessionInfo.txt"))




##load seurat object

sc.data <- sapply(list.dirs(path = path_to_data, recursive = FALSE, full.names = TRUE), 
									basename, 
									USE.NAMES = TRUE)

#remember to drop HPAP-093 b/c it clusters separately
for(i in sc.data){
	if(grepl("HPAP-093", i, ignore.case = TRUE)){
		sc.data <- sc.data[sc.data!=i]
	}
}

metadata <- read.table(file = paste0(metadata.location, "HPAPMetaData.txt"), header = TRUE, sep = "\t")
rownames(metadata) <- metadata$DonorID

# Exclude OW donors
for(i in sc.data){
	x <- strsplit(i, "_")[[1]][1]
	if(!(metadata[x, "SimpDisease"] == "NoDM" & (metadata[x, "BMI"] < 25 | metadata[x, "BMI"] > 30) & metadata[x, "scRNA"] > 0)){
		sc.data <- sc.data[sc.data!=i]}
}


object.list <- c()
for(i in 1:length(sc.data)){
	object.list[[i]] <- Read10X_h5(paste0(names(sc.data)[i], "/outs/filtered_feature_bc_matrix.h5"))
	object.list[[i]] <- CreateSeuratObject(object.list[[i]], 
																				 project = rnaProject, 
																				 min.cells = min.cells, 
																				 min.features = min.features)
	object.list[[i]]$orig.ident <- sc.data[[i]]
	object.list[[i]] <- AssignMetadata(metadata.df = metadata, seurat.object = object.list[[i]])
	object.list[[i]] <- PercentageFeatureSet(object.list[[i]], pattern = "MT-", col.name = "percent.mt")
	
	object.list[[i]] <- subset(object.list[[i]], 
														 subset = nFeature_RNA >= 200 &
														 	nFeature_RNA <= 2500 &
														 	percent.mt <= 5)
	
	print(paste("finished", sc.data[[i]]))
}


	seurat.object <- merge(object.list[[1]], y = object.list[2:length(object.list)], add.cell.ids = names(object.list))

	
	seurat.object$sequencerID <- seurat.object$orig.ident
	seurat.object$sequencerID <- with(seurat.object, stringi::stri_replace_all_fixed(seurat.object$sequencerID, seurat.object$DonorID, ""))
# QC ----------------------------------------------------------------------

##plot qc stats
VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1<- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Normalize and scale data ----------------------------------------------------------


if(do.sctransform == FALSE){ # standard method
	print("Performing standard normalization and scaling")
	if(length(regression.vars) >1){
		print("HEY YOU! You're performing standard scaling on more than 1 regression variable. You should probably be doing SCTransform. Set `do.sctransform` to TRUE")
	}
	
	##normalize
	seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
	
	##find HVG
	seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)
	top10hvg <- head(VariableFeatures(seurat.object), 10)
	plot1 <- VariableFeaturePlot(seurat.object)
	plot2	<- LabelPoints(plot = plot1, points = top10hvg, repel = TRUE)
	plot1 + plot2
	top10hvg
	
	##scale (a linear transformation)
	all.genes <- rownames(seurat.object)
	
	seurat.object <- ScaleData(seurat.object, features = all.genes, vars.to.regress = regression.vars)
	
} else if(do.sctransform == "each"){
	print("Performing SCTransform")
	object.list <- SplitObject(seurat.object, split.by = "orig.ident")
	object.list <- lapply(X = object.list, 
												FUN = SCTransform, assay = "RNA", return.only.var.genes = FALSE, vst.flavor = "v2")

	# RUN DOUBLETFINDER AFTER UMAPhttps://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_doublets.R
	object.list <- lapply(X = object.list, 
												FUN = runDoubletFinder, sctransformed = TRUE, tot.var = doublet.var.thresh, predicted.doubletRate = predicted.doubletRate)
	# object.list <- lapply(X = object.list,
												# FUN = CellCycleScoring, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
	object.list <- lapply(X = object.list,
												FUN = subset, subset = DF.classifications == "Singlet")
	

	integration.features <- SelectIntegrationFeatures(object.list = object.list, verbose = TRUE, nfeatures = 3000)
	object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = integration.features, verbose = TRUE)
	integration.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = integration.features, normalization.method = "SCT", verbose = TRUE)
	seurat.object <- IntegrateData(anchorset = integration.anchors, verbose = TRUE, preserve.order = FALSE, normalization.method = "SCT")
	
} else if(do.sctransform == "pooled") {
	seurat.object <- SCTransform(seurat.object, method = "glsGamPoi", vars.to.regress = regression.vars, verbose = TRUE, return.only.var.genes = FALSE, vst.flavor = "v2")
	seurat.object <- runDoubletFinder(seurat.object = seurat.object, sctransformed = TRUE, tot.var = doublet.var.thresh, predicted.doubletRate = predicted.doubletRate)
	seurat.object <- CellCycleScoring(seurat.object, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
	seurat.object <- subset(seurat.object, subset = DF.classifications == "Singlet")
	
}else {
	print("Must set do.sctransform to one of: FALSE, each, pooled")
}

saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, ".RDS"))
# Linear dimensional reduction --------------------------------------------


seurat.object <- RunPCA(seurat.object, features = VariableFeatures(object = seurat.object))

print(seurat.object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat.object, dims = 1:2, reduction = "pca")
DimPlot(seurat.object, reduction = "pca")
DimHeatmap(seurat.object, dims = 1:2, cells = 500, balanced = TRUE)

saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, ".RDS"))

# Determine dimensionality ------------------------------------------------

if(run.jackstraw == TRUE){
	seurat.object <- JackStraw(seurat.object, num.replicate = 100)
	seurat.object <- ScoreJackStraw(seurat.object, dims = 1:40)
	JackStrawPlot(seurat.object, dims = 1:40)
}

ElbowPlot(seurat.object)

# account for variance
tot.var <- percent.variance(seurat.object@reductions$pca@stdev, plot.var = FALSE, return.val = TRUE)
paste0("Num pcs for 80% variance: ", length(which(cumsum(tot.var) <= 80)))
paste0("Num pcs for 85% variance: ", length(which(cumsum(tot.var) <= 85)))
paste0("Num pcs for 90% variance: ", length(which(cumsum(tot.var) <= 90)))
paste0("Num pcs for 95% variance: ", length(which(cumsum(tot.var) <= 95)))

cluster.dims <- 0
if(cum.var.thresh > 0){
	cluster.dims <- length(which(cumsum(tot.var) <= cum.var.thresh))
}

saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, ".RDS"))

# Louvain cluster ---------------------------------------------------------

##cluster cells
seurat.object <- FindNeighbors(seurat.object, dims = 1:cluster.dims)
seurat.object <- FindClusters(seurat.object, resolution = resolution)
seurat.object <- RunUMAP(seurat.object, dims = 1:cluster.dims)

saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, "-90pctvar.RDS"))


##umap
DimPlot(seurat.object, reduction = "umap", cols = color.palette, label = T, label.size = 7, repel = T)
DimPlot(seurat.object, reduction = "umap", group.by = "Obesity", cols = color.palette, label = T, label.size = 7, repel = T)
DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident", cols = color.palette, label = T, label.size = 7, repel = T)
DimPlot(seurat.object, reduction = "umap", group.by = "DonorID", cols = color.palette, label = T, label.size = 7, repel = T)


FeaturePlot(seurat.object, reduction = "umap", features = "BMI")



# Find cluster biomarkers -------------------------------------------------

seurat.object <- PrepSCTFindMarkers(seurat.object, assay = "SCT")
##find positively expressed markers for all clusters compared to all remaining clusters

markers.seurat.pos <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers.seurat.pos <- markers.seurat.pos %>% 
	group_by(cluster) %>%
	arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.seurat.pos, file = paste0(rnaProject, "-posmarkers-90pctvar.rds"))

markers.seurat.all <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
markers.seurat.all %>%
	group_by(cluster) %>%
	arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.seurat.all, file = paste0(rnaProject, "-allmarkers-90pctvar.rds"))

##create workbook
markers.table <- openxlsx::createWorkbook()


##write positive markers to table
openxlsx::addWorksheet(markers.table, sheetName = "PosMarkers")
openxlsx::writeData(markers.table, sheet = "PosMarkers", x = markers.seurat.pos,startCol = 1, startRow = 1, colNames = TRUE)

##write all markers to table
openxlsx::addWorksheet(markers.table, sheetName = "AllMarkers")
openxlsx::writeData(markers.table, sheet = "AllMarkers", x = markers.seurat.all, startCol = 1, startRow = 1, colNames = TRUE)

##save workbook
openxlsx::saveWorkbook(wb = markers.table, file = paste0(rnaProject, "_seuratMarkers-90pctvar.xlsx"), overwrite = TRUE, returnValue = TRUE)


# Cell Cycle Scoring ------------------------------------------------------

if(!do.sctransform == FALSE){
	DefaultAssay(seurat.object) <- "SCT"
}
seurat.object <- CellCycleScoring(seurat.object, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)

# Local visualization -----------------------------------------------------
# seurat.object <- readRDS("Obesity_scRNA.RDS")

# did doublet removal happen?
colnames(seurat.object@meta.data)

DimPlot(seurat.object, cols = color.palette)
DimPlot(seurat.object, cols = color.palette, group.by = "seurat_clusters", split.by = "Obesity", pt.size = 0.4, ncol = 4)
FeaturePlot(seurat.object, features = "BMI", pt.size = 0.4, cols = c("blue", "red"))
FeaturePlot(seurat.object, features = "BMI", pt.size = 0.4, cols = c("blue", "red"), split.by = "SCT_snn_res.0.5", ncol = 4)



# Visualize after biowulf run ---------------------------------------------

seurat.object <- readRDS("Obesity_scRNA-SCTRegression-NW-OB.RDS")
colnames(seurat.object@meta.data)

cds <- readRDS("Obesity_scRNA-SCTRegression-NW-OB-monocle3CDS.RDS")

# Differential abundance testing ------------------------------------------
library(speckle)
library(limma)
library(ggplot2)


# Get some example data which has two groups, three cell types and two 
# biological replicates in each group
# Run propeller testing for cell type proportion differences between the two 
# groups
propeller(clusters = Idents(seurat.object), sample = seurat.object$DonorID, group = seurat.object$Obesity, transform = "asin")

# Plot cell type proportions
plotCellTypeProps(clusters=seurat.object$Obesity, sample=seurat.object$SCT_snn_res.0.5)




# Trajectory analysis -----------------------------------------------------

library(SeuratWrappers)
library(monocle3)


cds <- as.cell_data_set(seurat.object)
fData(cds)$gene_short_name <- rownames(fData(cds))


recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- seurat.object@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seurat.object@reductions$umap@cell.embeddings
cds <- learn_graph(cds, use_partition = F)










# https://rdrr.io/github/satijalab/seurat-wrappers/f/docs/monocle3.Rmd
# seurat.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)

# From monocle3 tutorial
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds){
	cell_ids <- which(colData(cds)[, "SCT_snn_res.0.5"] %in% c(6, 7))

	closest_vertex <-
		cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
	closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
	root_pr_nodes <-
		igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
																															(which.max(table(closest_vertex[cell_ids,]))))]

	root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
saveRDS(cds, file = paste0(rnaProject, "-monocle3CDS-90pctvar.RDS"))

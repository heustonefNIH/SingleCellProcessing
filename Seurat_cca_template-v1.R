# Load libraries and user functions ----------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("Seurat", "dplyr", "colorRamps", "R.utils")

# Load libraries
suppressMessages(ipak(packages))

# Process commandLine input -----------------------------------------------


##ADD dimensions so can plot in 3d!!


args <- commandArgs(trailingOnly = TRUE, asValues = TRUE, 
                    defaults = c(all_data.files = NULL,
                                 ProjectName = NULL,
                                 genome = NULL,
                                 max_pcs = NULL,
                                 resolution_list = NULL,
                                 Run_mito_filter = FALSE 
                    )
)

if(length(unlist(args)) == 0){
  print("Arguments:")
  print("all_data.files: List of data files to be processed")
  print("ProjectName: Name of project")
  print("genome: can be mm10 or GRCh38")
  print("Run_mito_filter: Logical; To filter, or not to filter (on expression level of mitochondrial genes)")
  print("max_pcs: number of prinicple components to use for clustering")
  print("resolution_list: list of resoultions for cluster analysis iterations")
}else if(length(unlist(args)) < 6){
  print("Must supply Seurat object, number of dimensions, and max components")
}else{
  all_data.files <- args$all_data.files
  ProjectName<-args$ProjectName
  genome<-args$genome
  Run_mito_filter<-args$Run_mito_filter
  max_pcs <- as.integer(args$max_pcs)
  resolution_list <- args$resolution_list
  if(Run_mito_filter = TRUE){
    print("Filtering based on mt-gene level")
  }else{
    print("Not filtering for mt-gene level")
  }
  print(paste("Performing clustering on", all_data.files, sep = " "))
  print(paste("Creating files with the output name", ProjectName, sep = " "))
  print(paste("Running seurat with genome", genome, "on", max_pcs, "principle components at resolutions", resolution_list, sep = " "))
}


# Define basic color palette ----------------------------------------------

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

# Load data ---------------------------------------------------------------

# Alter output name based on filter parameters
ProjectName<-paste(ProjectName, "_cca", sep = "")

if(Run_mito_filter == TRUE){
  ProjectName<-paste(ProjectName, "_mt", sep = "")
  ProjectName<-paste(ProjectName, "_mt", sep = "")
}

multi_object_list<-list()
# Create multi_object_list
create_multi_object_list<-function(x){
  cellranger_files<-paste(x[[1]], "/outs/filtered_gene_bc_matrices/",genome,"/", sep="")
  names(cellranger_files)<-names(x)
  print(cellranger_files)
  my_object.data<-Read10X(cellranger_files)
  my_object<- CreateSeuratObject(raw.data = my_object.data, min.cells = 3, min.genes = 200, project = ProjectName)
  print(head(my_object@cell.names))
  # Apply mito.filter if applicable
  if(Run_mito_filter == TRUE){
    mito.genes<-grep(pattern = "^MT-", x = rownames(x=my_object@data), value = TRUE, ignore.case = TRUE)
    percent.mito <-Matrix::colSums(my_object@raw.data[mito.genes,])/Matrix::colSums(my_object@raw.data)
    my_object<-AddMetaData(my_object, metadata = percent.mito, col.name = "percent.mito")
  } else {
    print("Not running mito filter")
  }
  my_object<-NormalizeData(my_object, normalization.method = "LogNormalize", scale.factor = 10000)
  my_object<-FindVariableGenes(my_object,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  
  # Scale on percent.mito if applicable
  if(Run_mito_filter == TRUE){
    my_object<-ScaleData(my_object, vars.to.regress = c("nUMI", "percent.mito"))
  } else {
    my_object<-ScaleData(my_object, vars.to.regress = "nUMI")
  }
  my_object@meta.data$sample<-names(x)
  multi_object_list[[names(x)]]<-my_object
}
multi_object_list<-lapply(1:length(all_data.files), function(x) create_multi_object_list(all_data.files[x]))
names(multi_object_list)<-names(all_data.files)

genes.use<-c()
for (i in 1:length(multi_object_list)){
  print(i)
  if(length(multi_object_list[[i]]@var.genes)>=1000){
    genes.use<-c(genes.use, head(rownames(multi_object_list[[i]]@hvg.info), 1000))
    print(length(paste("Using", genes.use, "genes to align subspace", sep = " ")))
  } else {
    genes.use<-c(genes.use, head(rownames(multi_object_list[[i]]@hvg.info)))
    print(length(paste("Using", genes.use, "genes to align subspace", sep = " ")))
  }
}

genes.use<-names(which(table(genes.use)>1))
for (i in 1:length(multi_object_list)){
  genes.use<-genes.use[genes.use %in% rownames(multi_object_list[[i]]@scale.data)]
}
head(genes.use)
# Run multi-set CCA


# Align subspace ----------------------------------------------------------

integrated_object<-RunMultiCCA(multi_object_list, genes.use = genes.use, num.ccs = max_pcs)

p1<-DimPlot(integrated_object, reduction.use="cca", group.by="orig.ident", do.return=TRUE)
p2<-VlnPlot(integrated_object, features.plot="CC1", group.by = "orig.ident", do.return = TRUE)

png(filename = paste(ProjectName,"_dim",max_pcs, "_CCFit.png", sep = ""), height = 800, width = 800)
plot_grid(p1, p2)
dev.off()

png(filename = paste(ProjectName,"_dim",max_pcs, "_BicorPlot.png", sep = ""), height = 800, width = 800)
MetageneBicorPlot(integrated_object, grouping.var = "sample", dims.eval = 1:max_pcs)
dev.off()

png(filename = paste(ProjectName,"_dim",max_pcs, "_Heatmap.png", sep = ""), height = 2400, width = 800)
DimHeatmap(integrated_object, reduction.type = "cca", cells.use = 500, dim.use = 1:max_pcs, do.balanced = TRUE)
dev.off()


integrated_object<-CalcVarExpRatio(integrated_object, reduction.type = "pca", grouping.var = "sample", dims.use = 1:max_pcs)
integrated_object<-SubsetData(integrated_object, subset.name = "var.ratio.pca", accept.low = 0.5)
integrated_object<-AlignSubspace(integrated_object, reduction.type = "cca", grouping.var = "sample", dims.align = 1:max_pcs)

# Iterate across multiple resolutions 
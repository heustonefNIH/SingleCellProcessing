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
                    "DBA1_2" = "DBA1_2/outs/filtered_gene_bc_matrices/GRCh38/",
                    "DBA3_2" = "DBA3_2/outs/filtered_gene_bc_matrices/GRCh38/",
                    "DBA4_2" = "DBA4_2/outs/filtered_gene_bc_matrices/GRCh38/")
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

multi_object_list<-list()
list_of_names<-names(all_data.files)
list_of_names


# Create multi_object_list
create_multi_object_list<-function(x){
  my_object.data<-Read10X(x)
  my_object<- CreateSeuratObject(raw.data = my_object.data, min.cells = 3, min.genes = 200, project = ProjectName)
  
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
names(multi_object_list)
names(multi_object_list)<-list_of_names
multi_object_list

genes.use<-c()
for (i in 1:length(multi_object_list)){
  print(i)
  if(length(multi_object_list[[i]]@var.genes)>=1000){
  genes.use<-c(genes.use, head(rownames(multi_object_list[[i]]@hvg.info), 1000))
  } else {
    genes.use<-c(genes.use, head(rownames(multi_object_list[[i]]@hvg.info)))
  }
}

genes.use<-names(which(table(genes.use)>1))
for (i in 1:length(multi_object_list)){
  genes.use<-genes.use[genes.use %in% rownames(multi_object_list[[i]]@scale.data)]
}
head(genes.use)
# Run multi-set CCA

integrated_object<-RunMultiCCA(multi_object_list, genes.use = genes.use, num.ccs = max_pcs)

MetageneBicorPlot(integrated_object, grouping.var = "sample", dims.eval = 1:max_pcs)

integrated_object<-CalcVarExpRatio(integrated_object, reduction.type = "pca", grouping.var = "sample", dims.use = 1:max_pcs)
integrated_object<-SubsetData(integrated_object, subset.name = "var.ratio.pca", accept.low = 0.5)
integrated_object<-AlignSubspace(integrated_object, reduction.type = "cca", grouping.var = "sample", dims.align = 1:max_pcs)
integrated_object<-FindClusters(integrated_object, reduction.type = "cca.aligned", dims.use = 1:max_pcs, save.SNN = T, resolution = 0.6)







































# Load libraries and user functions ----------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("Seurat", "dplyr", "colorRamps", "tools", "ggplot2")

# Load libraries
ipak(packages)


# Define variables --------------------------------------------------------

# all_data.files <- c(LSK = "LSKm2",
#                     CMP = "CMPm2",
#                     MEP = "MEPm",
#                     GMP = "GMPm")
all_data.files <- c(MEP = "MEPm")
ProjectName<-"test"
genome<- "mm10"
max_pcs <-5
Run_mito_filter = FALSE
# resolution_list <-c(0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5)
resolution_list<-0.6
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


# Define script-specific functions ---------------------------------------------------

# define transparnet dot
my_transparent_dot<-rgb(0,0,0, alpha = 0)

# define function to create pngs
png_plotFunction<-function(plot_to_make = x, filename = y, height = 800, width = 800){ # default window is 800x800
  png(filename = filename, height = height, width = width)
  plot_to_make
  dev.off()
}

# define PCA plotting function that colors one variable red and all others black, and cycles through ##can change colors, too
emphasis_plots_PCA<-function(object, file_suffix, iterationSlot, bland_color = "black", emphasis_color = "red", ...){
  iteration_vector<-(unique(FetchData(object = object, vars.all = iterationSlot)[,1]))
  for(x in 1:length(iteration_vector)){
  color_vector <-c(rep(bland_color, length(iteration_vector)))
    color_vector[x]<-emphasis_color
    unique_id<-iteration_vector[x]
    png_plotFunction(PCAPlot(object = object, cols.use = color_vector, group.by = iterationSlot, pt.size = 1), filename = paste(ProjectName, "_", unique_id,"_", file_suffix,".png", sep = ""), ...)
    }
}

# define function to change case of cell cycle genes to match mouse gene nomenclature
if(genome != "GRCh38"){
  proper=function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
  s.genes<-proper(Seurat::cc.genes$s.genes)
  g2m.genes<-proper(Seurat::cc.genes$g2m.genes)
} else{
  s.genes<-Seurat::cc.genes$s.genes
  g2m.genes<-Seurat::cc.genes$g2m.genes
}


# Import data -------------------------------------------------------------

read_10XGenomics_data<-function(x){
  x[[1]]<-paste(x[[1]], "/outs/filtered_gene_bc_matrices/",genome,"/", sep="")
}
data_files.list<-as.character(lapply(1:length(all_data.files), function(x) read_10XGenomics_data(all_data.files[x])))
names(data_files.list)<-names(all_data.files)

my_object.data<-Read10X(data_files.list)
my_object<- CreateSeuratObject(raw.data = my_object.data, min.cells = 3, min.genes = 200, project = ProjectName)
my_object<-NormalizeData(my_object, normalization.method = "LogNormalize", scale.factor = 10000)
my_object<-FindVariableGenes(my_object,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
my_object<-ScaleData(my_object, vars.to.regress = "nUMI")


# Cell cycle analysis -------------------------------------------------------

cc_my_object<-CellCycleScoring(my_object, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)
ridgeplot_genes<- c("Pcna", "Top2a", "Mcm6", "Mki67")


# Regress out ALL markers of cell cycle -----------------------------------

# This method means you no longer see difference between cycling cells and non-cycling cells

if(remove_cellcycle_stages == TRUE){
  cc_my_object<-RunPCA(cc_my_object, pc.genes = c(s.genes, g2m.genes), do.print = TRUE)
  PCAPlot(cc_my_object, group.by = "Phase")
  png_plotFunction(PCAPlot(cc_my_object, group.by = "Phase"), filename = paste(ProjectName, "_PCAPlot_phase.png", sep = ""))
  
  # Shows cycle of genes in each cluster
  cc_my_objectcycle<-CellCycleScoring(cc_my_object, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = F)
  PCAPlot(cc_my_objectcycle, group.by = "Phase")
}



# Regress out STAGES of cell cycle ----------------------------------------

# This means you see a difference between cycling and non-cycling, but lose clustering based on things like S vs G2M stage


cc_my_object@meta.data$CC.Difference<-cc_my_object@meta.data$S.Score - cc_my_object@meta.data$G2M.Score
cc_my_object<-ScaleData(cc_my_object, vars.to.regress = "CC.Difference", display.progress = TRUE)
cc_my_object<-RunPCA(cc_my_object, pc.genes = cc_my_object@var.genes, genes.print = 10)  #Removes differences between G1, S, and G2m

if(separate_cycling_cells == TRUE){
  cc_my_object<-RunPCA(cc_my_object, pc.genes = c(s.genes, g2m.genes), genes.print = 10)  #Separates G1, from c(S, G2m)
}


try(emphasis_plots_PCA(object = cc_my_object, file_suffix = "pcaplot", iterationSlot = "Phase"))
try(emphasis_plots_PCA(object = cc_my_object, file_suffix = "pcaplot", iterationSlot = "orig.ident"))

# cc_my_object.projectpca<-ProjectPCA(cc_my_objectcycle)


png(filename = "cc_my_object.projectpca_VizPCA.png", width = 800, height = 800)
VizPCA(cc_my_object.projectpca, use.full = T, font.size = 1, nCol = 3, pcs.use = 1:6)
dev.off()



RidgePlot(cc_my_object, features.plot = ridgeplot_genes, nCol = 2)
DotPlot(cc_my_object, genes.plot = ridgeplot_genes)
# try(emphasis_plots_PCA(object = cc_my_object, file_suffix = "pcaplot", iterationSlot = grep(pattern = "res", colnames(cc_my_object@meta.data), ignore.case = TRUE)))
dim14_scaledforcellcycle<-ScaleData(cc_my_object, vars.to.regress=c("S.Score", "G2M.Score"), display.progress = TRUE)
dim14_scaledforcellcycle<-RunPCA(cc_my_object, pc.genes = dim14@var.genes, genes.print = 10)




png(filename = "dim14cycle_ElbowPlot.png", width = 800, height = 800)
PCElbowPlot(cc_my_object.projectpca)
dev.off()



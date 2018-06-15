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
remove_cellcycle_stages = FALSE
separate_cycling_cells = FALSE
ridgeplot_genes<- c("Pcna", "Top2a", "Mcm6", "Mki67")
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
if(genome == "mm10"){
  s.genes<-c("Mcm4", "Exo1", "Slbp", "Gmnn", "Cdc45", "Msh2", "Mcm6", "Rrm2", "Pold3", "Blm", "Ubr7", "Mcm5", "Clspn", "Hells", "Nasp", "Rpa2", "Rad51ap1", "Tyms", "Rrm1", "Rfc2", "Prim1", "Brip1", "Usp1", "Ung", "Pola1", "Mcm2", "Fen1", "Tipin", "Pcna", "Cdca7", "Uhrf1", "Casp8ap2", "Cdc6", "Dscc1", "Wdr76", "E2f8", "Dtl", "Ccne2", "Atad2", "Gins2", "Chaf1b", "Pcna-ps2")
  g2m.genes<-c("Nuf2", "Psrc1", "Ncapd2", "Ccnb2", "Smc4", "Lbr", "Tacc3", "Cenpa", "Kif23", "Cdca2", "Anp32e", "G2e3", "Cdca3", "Anln", "Cenpe", "Gas2l3", "Tubb4b", "Cenpf", "Dlgap5", "Hjurp", "Cks1brt", "Gtse1", "Bub1", "Birc5", "Ube2c", "Rangap1", "Hmmr", "Ect2", "Tpx2", "Ckap5", "Cbx5", "Nek2", "Ttk", "Cdca8", "Nusap1", "Ctcf", "Cdc20", "Cks2", "Mki67", "Tmpo", "Ckap2l", "Aurkb", "Kif2c", "Cdk1", "Kif20b", "Top2a", "Aurka", "Ckap2", "Hmgb2", "Cdc25c", "Ndc80", "Kif11")  # list of mouse cell cycle genes provided by leonfodoulian from https://github.com/satijalab/seurat/issues/462; ##is a converted set of the seurat cc.genes list
 } else if (genome == "GRCh38") {
   s.genes<-Seurat::cc.genes$s.genes
   g2m.genes<-Seurat::cc.genes$g2m.genes
 } else {
   print ("No cell cycle genes match this your genome")
 }


# define function to import cellranger file(s)
read_10XGenomics_data<-function(x){
  x[[1]]<-paste(x[[1]], "/outs/filtered_gene_bc_matrices/",genome,"/", sep="")
}

# Import data -------------------------------------------------------------

data_files.list<-as.character(lapply(1:length(all_data.files), function(x) read_10XGenomics_data(all_data.files[x])))
names(data_files.list)<-names(all_data.files)

my_object.data<-Read10X(data_files.list)
my_object<- CreateSeuratObject(raw.data = my_object.data, min.cells = 3, min.genes = 200, project = ProjectName)
my_object<-NormalizeData(my_object, normalization.method = "LogNormalize", scale.factor = 10000)
my_object<-FindVariableGenes(my_object,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)


# Regress out ALL markers of cell cycle -----------------------------------

# This method means you no longer see difference between cycling cells and non-cycling cells

if(remove_cellcycle_stages == TRUE){
  my_object<-ScaleData(my_object, vars.to.regress = "nUMI") # only scale on UMI but PCA on cc.genes
  cc_my_object<-CellCycleScoring(my_object, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)
  cc_my_object<-RunPCA(cc_my_object, pc.genes = c(s.genes, g2m.genes), do.print = TRUE)
  PCAPlot(cc_my_object, group.by = "Phase")
  cc_my_object<-ScaleData(cc_my_object, vars.to.regress = c("nUMI","S.Score", "G2M.Score"), display.progress = TRUE)
  cc_my_object<-RunPCA(cc_my_object, pc.genes = c(s.genes, g2m.genes), do.print = TRUE)
  png_plotFunction(PCAPlot(cc_my_object, group.by = "Phase"), filename = paste(ProjectName, "PCAPlot_phase.png", sep = ""))
  
  # Shows cycle of genes in each cluster
  cc_my_objectcycle<-CellCycleScoring(cc_my_object, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = F)
  PCAPlot(cc_my_objectcycle, group.by = "Phase")
} 


# Regress out STAGES of cell cycle ----------------------------------------

# This means you see a difference between cycling and non-cycling, but lose clustering based on things like S vs G2M stage

if(remove_cellcycle_stages == FALSE){
  cc_my_object<-CellCycleScoring(my_object, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)
  cc_my_object@meta.data$CC.Difference<-cc_my_object@meta.data$S.Score - cc_my_object@meta.data$G2M.Score
  cc_my_object<-ScaleData(cc_my_object, vars.to.regress = c("nUMI","CC.Difference"), display.progress = TRUE)  # scale data works on the @data slot, so if you do this at different steps you will lose the info from the first round of scaling
  cc_my_object<-RunPCA(cc_my_object, pc.genes = cc_my_object@var.genes, genes.print = 10)  #Removes differences between G1, S, and G2m
  
  if(separate_cycling_cells == TRUE){
    cc_my_object<-RunPCA(cc_my_object, pc.genes = c(s.genes, g2m.genes), genes.print = 10)  #Separates G1, from c(S, G2m)
  }
  
  try(emphasis_plots_PCA(object = cc_my_object, file_suffix = "pcaplot", iterationSlot = "Phase"))
  try(emphasis_plots_PCA(object = cc_my_object, file_suffix = "pcaplot", iterationSlot = "orig.ident"))
}
# cc_my_object.projectpca<-ProjectPCA(cc_my_objectcycle)

try(png_plotFunction(VizPCA(cc_my_object, font.size = 1, nCol = 3, pcs.use = 1:6), filename = paste(ProjectName,"VizPlot.png")))
try(png_plotFunction(RidgePlot(cc_my_object, features.plot = ridgeplot_genes, nCol = 2),filename = paste(ProjectName,"Ridgeplot.png")))
try(png_plotFunction(DotPlot(cc_my_object, genes.plot = ridgeplot_genes),filename = paste(ProjectName,"ccDotPlotgenes.png")))
try(png_plotFunction(PCElbowPlot(cc_my_object), filename = paste(ProjectName,"ElbowPlot.png")))




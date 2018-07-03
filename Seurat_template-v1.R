# Load libraries and user functions ----------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("Seurat", "dplyr", "colorRamps", "ggplot2")

# Load libraries
ipak(packages)

# Define variables --------------------------------------------------------

# all_data.files <- c(LSK = "LSKm2",
#                     CMP = "CMPm2",
#                     MEP = "MEPm",
#                     GMP = "GMPm")
all_data.files<-c(MEP="MEPm", LSK= "LSKm2")
ProjectName<-"test"
genome<- "mm10"
Run_mito_filter = FALSE
max_pcs <-5
# resolution_list <-c(0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5)
resolution_list<-0.6

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

if(Run_mito_filter == TRUE){
  ProjectName<-paste(ProjectName, "_mt", sep = "")
}

# Read in 10XGenomics files
read_10XGenomics_data<-function(x){
  x = paste(x[[1]], "/outs/filtered_gene_bc_matrices/",genome,"/", sep="")
}
data_files.list<-sapply(all_data.files, function(x) read_10XGenomics_data(all_data.files))

my_object.data<-Read10X(data_files.list)
my_object<- CreateSeuratObject(raw.data = my_object.data, min.cells = 3, min.genes = 200, project = ProjectName)
# filter_cells_and_genes(my_obj,min_UMIs=1000,max_UMIs=25000,min_detected_genes=1000,max_detected_genes=5000,max_percent_mito=9) <--Supat's filter

# Apply mito.filter if applicable

if(Run_mito_filter == TRUE){
  mito.genes<-grep(pattern = "^MT-", x = rownames(x=my_object@data), value = TRUE, ignore.case = TRUE)
  percent.mito <-Matrix::colSums(my_object@raw.data[mito.genes,])/Matrix::colSums(my_object@raw.data)
  my_object<-AddMetaData(my_object, metadata = percent.mito, col.name = "percent.mito")
  VlnPlot(my_object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  par(mfrow = c(1,2))
  GenePlot(my_object, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(my_object, gene1 = "nUMI", gene2 = "nGene")
  par(mfrow=c(1,1))
  png(paste(ProjectName, "_VlnStats.png", sep = ""), width = 800, height = 800)
  VlnPlot(my_object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  dev.off()
  
  png(paste(ProjectName, "_GenePlots.png", sep = ""), width = 800, height = 800)
  par(mfrow = c(1,2))
  GenePlot(my_object, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(my_object, gene1 = "nUMI", gene2 = "nGene")
  dev.off()
  
} else {
  VlnPlot(my_object, features.plot = c("nGene", "nUMI"), nCol = 2)
  GenePlot(my_object, gene1 = "nUMI", gene2 = "nGene")

  png(paste(ProjectName, "_VlnStats.png", sep = ""), width = 800, height = 800)
  VlnPlot(my_object, features.plot = c("nGene", "nUMI"), nCol = 2)
  dev.off()
  
  png(paste(ProjectName, "_GenePlots.png", sep = ""), width = 800, height = 800)
  GenePlot(my_object, gene1 = "nUMI", gene2 = "nGene")
  dev.off()
}


# Normalize and scale data ----------------------------------------------------------

my_object<-NormalizeData(my_object, normalization.method = "LogNormalize", scale.factor = 10000)


#Variable gene detection helps control for relationship between variability and average expression
#Set to mark visual outliers, potentially ~2000

my_object<-FindVariableGenes(my_object,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
VariableGenePlot(my_object)

png(paste(ProjectName, "_VarGenePlot.png", sep = ""), width = 800, height = 800)
VariableGenePlot(my_object)
dev.off()





if(Run_mito_filter == TRUE){
  my_object<-ScaleData(my_object, vars.to.regress = c("nUMI", "percent.mito"))
  print("I scaled with percent.mito regression")
} else {
  my_object<-ScaleData(my_object, vars.to.regress = "nUMI")
}

my_object<-CellCycleScoring(my_object, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)


# Run PCA -----------------------------------------------------------------
my_object<-RunPCA(my_object, pc.genes = my_object@var.genes, pcs.compute = max_pcs, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

plot_title<-paste(ProjectName, "dim", max_pcs, sep = "")

PrintPCA(my_object, pcs.print = 1:5, genes.print = 5, use.full = FALSE) #Set use.full to TRUE to see projected PCA
VizPCA(my_object, pcs.use = 1:5, use.full = FALSE, font.size = 1)
PCAPlot(my_object, dim.1 = 1, dim.2 = 2, use.full = FALSE)
PCHeatmap(my_object, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, do.return=TRUE)

my_plot<-PCAPlot(my_object, dim.1 = 1, dim.2 = 2, use.full = FALSE, do.return=TRUE)
my_plot<-my_plot+ggtitle(plot_title)
png(paste(ProjectName,"dim",max_pcs, "_PCAPlot.png", sep = ""), width = 800, height = 800)
plot(my_plot)
dev.off()

VizPCA(my_object, pcs.use = 1:5, use.full = FALSE, font.size = 1.5)
png(paste(ProjectName, "dim",max_pcs, "_VizPCA.png", sep = ""), width = 800, height = 1400)
par(mar=c(5,10,5,5))
VizPCA(my_object, pcs.use = 1:5, use.full = FALSE, font.size = 1.5)
dev.off()
PCHeatmap(my_object, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, do.return=FALSE)
png(paste(ProjectName,"dim",max_pcs, "_PCAHeatmap.png", sep = ""), width = 800, height = 1400)
par(mar=c(5,5,5,10))
PCHeatmap(my_object, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, do.return = FALSE)
dev.off()
par(mar=c(0,0,0,0))




# ProjectPCA --------------------------------------------------------------

#Scores each gene in the dataset based on their correlation with calculated components, regarless of whether or not the genes were included in the PCA
#Can be used to ID markers correlated with cellular heterogeneity but didn't pass through variable selection
#Can probably use this to calculate PCA based on 1 sample's ID, then project other data onto it
#Authors say they no longer use the projected PCA in downstream analysis

projected_my_object <-ProjectPCA(my_object, do.print = TRUE)

PrintPCA(projected_my_object, pcs.print = 1:5, genes.print = 5, use.full = TRUE) #Set use.full to TRUE to see projected PCA
VizPCA(projected_my_object, pcs.use = 1:5, use.full = TRUE)
PCAPlot(projected_my_object, dim.1 = 1, dim.2 = 2, use.full = TRUE)
PCHeatmap(my_object, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = TRUE)

png(paste(ProjectName,"dim",max_pcs, "_projectedVizPCA.png", sep = ""), width = 800, height = 800)
VizPCA(projected_my_object, pcs.use = 1:5)
dev.off()

png(paste(ProjectName,"dim",max_pcs, "_projectedHeatmap.png", sep = ""), width = 800, height = 800)
PCHeatmap(projected_my_object, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = TRUE)
dev.off()


# Find statistically signficant principal components -----------------------

my_object<-JackStraw(my_object, num.replicate = 100, num.pc = max_pcs, display.progress = TRUE)

PCElbowPlot(my_object)
JackStrawPlot(my_object, PCs = 1:max_pcs)


png(paste(ProjectName,"dim",max_pcs, "_ElbowPlot.png", sep = ""), width = 800, height = 800)
PCElbowPlot(my_object)
dev.off()
png(paste(ProjectName,"dim",max_pcs, "_JackStraw.png", sep = ""), width = 800, height = 800)
JackStrawPlot(my_object, PCs = 1:max_pcs)
dev.off()


# Clustering --------------------------------------------------------------

for(resolution in resolution_list){
  my_object<-FindClusters(my_object, reduction.type = "pca", dims.use = 1:max_pcs, resolution = resolution, print.output = 0, save.SNN = TRUE)
  PrintFindClustersParams(my_object)
  
  
  
  my_object<-RunTSNE(my_object, dims.use = 1:max_pcs, do.fast = TRUE)
  
  # Create color pallete
  new_length <-0
  if(length(unique(my_object@ident)) > length(my_palette)){
    new_length<-length(unique(my_object@ident)) - length(my_palette)
  }
  new_length
  my_palette<-c(my_palette, primary.colors(new_length))
  

  TSNEPlot(my_object, colors.use = my_palette)
  TSNEPlot(my_object, colors.use = my_palette, do.label = TRUE, label.size = 10, group.by = "orig.ident")
  
  
  plot_title<-paste(ProjectName, "dim",max_pcs, "res",resolution, sep = "")
  my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.return=TRUE)
  my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
  png(paste(ProjectName, "dim",max_pcs, "res",resolution,"_tsne.png", sep = ""), width = 800, height = 800)
  try(plot(my_tsne_plot))
  dev.off()
  my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.label = TRUE, label.size = 10, do.return=TRUE)
  my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
  png(paste(ProjectName, "dim",max_pcs, "res",resolution,"_tsne-Labeled.png", sep = ""), width = 800, height = 800)
  try(plot(my_tsne_plot))
  dev.off()
  my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, group.by = "orig.ident", do.return=TRUE)
  my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
  png(paste(ProjectName, "dim",max_pcs, "res",resolution,"_tsnebyID.png", sep = ""), width = 800, height = 800)
  try(plot(my_tsne_plot))
  dev.off()
  my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.label = TRUE, label.size = 10, group.by = "orig.ident", do.return=TRUE)
  my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
  png(paste(ProjectName, "dim",max_pcs, "res",resolution,"_tsnebyID-Labeled.png", sep = ""), width = 800, height = 800)
  try(plot(my_tsne_plot))
  dev.off()
  my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, group.by = "Phase", do.return=TRUE)
  my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
  png(paste(ProjectName, "dim",max_pcs, "res",resolution,"_tsnebyPhase.png", sep = ""), width = 800, height = 800)
  try(plot(my_tsne_plot))
  dev.off()
  my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.label = TRUE, label.size = 10, group.by = "Phase", do.return=TRUE)
  my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
  png(paste(ProjectName, "dim",max_pcs, "res",resolution,"_tsnebyPhase-Labeled.png", sep = ""), width = 800, height = 800)
  try(plot(my_tsne_plot))
  dev.off()

  
  
  saveRDS(my_object, file = paste(ProjectName, "dim",max_pcs, "res",resolution,"_tsne.rds", sep = ""))
  print(paste("Saved ",ProjectName, "dim",max_pcs, "res",resolution,"_tsne.rds", sep = ""))
  # Plot cell cycle states --------------------------------------------------------------
  
  ridgeplot_genes<- c("Pcna", "Top2a", "Mcm6", "Mki67")

  # Finding markers ---------------------------------------------------------

  my_object_markers<-FindAllMarkers(my_object, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

  my_object_markers_top100 <- my_object_markers %>% group_by(cluster) %>% top_n(100, avg_logFC)
  write.table(my_object_markers_top100, file = paste(ProjectName, "dim",max_pcs, "res",resolution,"_markers_top100.txt", sep = ""), sep="\t", quote=F)
  my_object_markers_all <- my_object_markers %>% group_by(cluster)
  write.table(my_object_markers_all, file = paste(ProjectName, "dim",max_pcs, "res",resolution,"_markers_all.txt", sep = ""), sep="\t", quote=F)

  SaveClusters(my_object, paste(ProjectName, "dim",max_pcs, "res",resolution,"_saveclusters.tsv", sep = ""))

  # Print gene plots --------------------------------------------------------

  drawmyplot<-function(geneList, tsne.obj, name){
    for(gene in geneList){
      png(filename=paste(name, "-Vln_",gene, ".png", sep=""), width=800, heigh=800)
      try(print(VlnPlot(tsne.obj, gene, do.return = T, point.size.use = 0, cols.use = my_palette)))
      dev.off()
      png(filename=paste(name, "-Ftr_",gene, ".png", sep=""), width=800, height=800)
      try(FeaturePlot(tsne.obj, gene,cols.use=c("grey", "blue"), pt.size = 1, do.return = T))
      dev.off()
    }
  }

  geneList = c("Nfe2", "Gata1", "Gata2", "Hbg1", "Hbg2", "Zfpm1", "Hbb", "Hba1", "Hba2", "Hbd", "Hbe1", "Klf1", "Fli1", "Meis1", "Kit", "Vwf", "Pf4", "Mpo", "Runx1", "Csf1", "Tfr2", "Cnrip1", "Myc" ,"Tk1", "Rrm2", "Itga2b", "Lmna", "Nfia", "Gp1bb", "Plek", "Pbx1", "Hes6", "E2f4","Dntt", "Vpreb1", "Id3", "Atf3", "Jchain", "Cd79a", "Satb1", "Sp140", "Tgfbi", "Lgmn", "Irf8", "Irf7", "Tcf4",  "Batf3", "Tcf19", "Sell", "Cd52", "Hoxb5", "Gata3", "Eno1", "Mpo", "Atp8b4", "Spi1", "Mafk", "Hdc", "Prg2", "Lmo4","Ctsg","Elane", "Cebpa", "Lgals1", "Fosb", "Prtn3",  "Tfrc", "Mpl", "Flt3", "Ca1", "Cd177", "Cd180","Cd244","Cd24","Cd27","Cd34","Cd37","Cd47","Cd48","Cd52","Cd53","Cd63","Cd68","Cd69","Cd72","Cd74","Cd81","Cd82","Cd84","Cd9","Cd93", "Mt1a", "Mt1g", "Mt1f", "Mt1e", "Mt1x")
  drawmyplot(geneList, my_object, name = paste(ProjectName, "dim",max_pcs, "res",resolution, sep = ""))

  housekeepingGenes<-c("Actb", "Gapdh", "Rn18s", "Ppia", "Rpl13a", "Rplp0", "B2m", "Hmbs", "Pgk1", "Alas1", "Gusb", "Sdha", "Tbp", "Tubb", "Ywhaz")
  drawmyplot(geneList, my_object, name = paste(ProjectName, "dim",max_pcs, "res",resolution, sep = ""))

  # Tabulate data -----------------------------------------------------------
  
  #Number of cells in each cluster
  table(my_object@ident)
  write.table(table(my_object@ident), file = paste(ProjectName, "dim",max_pcs, "res",resolution,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
  
  
  #Number of cells in each ID
  table(my_object@meta.data$orig.ident)
  write.table(table(my_object@meta.data$orig.ident), file = paste(ProjectName, "dim",max_pcs, "res",resolution,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  
  #Proportion of cells in each cluster
  prop.table(x = table(my_object@ident))
  write.table(prop.table(x = table(my_object@ident)), file = paste(ProjectName, "dim",max_pcs, "res",resolution,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  
  #Number of cells from each ID in each cluster
  table(my_object@ident, my_object@meta.data$orig.ident)
  write.table(table(my_object@ident, my_object@meta.data$orig.ident), file = paste(ProjectName, "dim",max_pcs, "res",resolution,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  
  #Proportion of cells of each ID in each cluster
  prop.table(x = table(my_object@ident, my_object@meta.data$orig.ident))
  write.table(prop.table(x = table(my_object@ident, my_object@meta.data$orig.ident)), file = paste(ProjectName, "dim",max_pcs, "res",resolution,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  
  #Average expression of each gene in each cluster
  avgexp<-AverageExpression(my_object)
  write.table(avgexp, file = paste(ProjectName, "dim",max_pcs, "res",resolution,"_AvgXprsn.txt", sep = ""), row.names = TRUE, quote = FALSE, sep = "\t")

}








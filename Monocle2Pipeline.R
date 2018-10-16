# Process commandLine input -----------------------------------------------
args<-commandArgs(trailingOnly = TRUE)
print(args)

main<-function(){
  args<-paste(unlist(args), collapse = " ")
  args<-unlist(strsplit(args, "--"))[-1]
  option_arguments<-sapply(args, function(x){
    unlist(strsplit(x, " "))[-1]
  })
  option_names<-sapply(args, function(x){
    unlist(strsplit(x, " "))[1]
  })
  names(option_arguments) <- unlist(option_names)
  required_keys<-c("seurat_object_filename", "num_dim")
  
  if(!(all(required_keys %in% unlist(args)))){
    print(paste("Missing:", paste(required_list[!(required_list%in%names(option_arguments))], collapse = " "), sep = " "))
    stop("Required arguments:
           --seurat_object_filename: name of seurat object to import (note: file must have .Robj or .rds extension)
           --num_dim: dimensionality of the reduced space
           Additional arguments:
           --max_components: Dimensions to plot (default = 2) note: currently I don't think this does anything
           --perform_expression_filtering: should analyzed RNAs be limited to expression of at least 0.1 and in at least 10 cells (default = TRUE) 
           --color_by_seurat_res: (default = TRUE)
           --order_by_seurat_varGenes: (default = FALSE)
           --UMI_bounded_filtering: (default = \"upper\", can be \"upper\", \"lower\", \"both\", or \"none\")
           --cca_variables: (default = \"~nUMI + nGene\")")
  } else{
<<<<<<< HEAD
    print("all good!")
    args<-paste(unlist(args), collapse = " ")
    args<-unlist(strsplit(args, "--"))[-1]
    option_arguments<-sapply(args, function(x){
      unlist(strsplit(x, " "))[-1]
    })
    option_names<-sapply(args, function(x){
      unlist(strsplit(x, " "))[1]
    })
    names(option_arguments) <- unlist(option_names)


=======
>>>>>>> cd142d8... Fixes to clarify error handeling
    if(!("max_components" %in% names(option_arguments))){
      print("Not filtering expression values")
      option_arguments$max_components<-2
    }
    if(!("perform_expression_filtering" %in% names(option_arguments))){
      print("Not filtering expression values")
      option_arguments$perform_expression_filtering<-TRUE
    }
    if(!("UMI_bounded_filtering" %in% names(option_arguments))){
      print("Using upper UMI_ounded filtering")
      option_arguments$UMI_bounded_filtering<-"upper"
    }
    if(!("order_by_seurat_varGenes" %in% names(option_arguments))){
      print("Not ordering by seurat_varGenes")
      option_arguments$order_by_seurat_varGenes<-FALSE
    }
    if(!("color_by_seurat_res" %in% names(option_arguments))){
      print("Coloring by seurat resolutions")
      option_arguments$color_by_seurat_res<-TRUE
    }
    if(!("cca_variables" %in% names(option_arguments))){
      print("Not ordering by seurat_varGenes")
      option_arguments$cca_variables<- "~nUMI + nGene"
    }

    return(option_arguments)

  }
}

option_arguments<-main()
option_arguments$num_dim<-as.numeric(option_arguments$num_dim)
option_arguments$max_components<-as.numeric(option_arguments$max_components)


print(paste("Setting UMI bounded filtering to", option_arguments$UMI_bounded_filtering, sep = " "))
print(paste("Setting correction variables as", list(option_arguments$cca_variables), sep = " "))
print(paste("Will process", option_arguments$seurat_object_filename, "using", option_arguments$num_dim, "dimensions and", option_arguments$max_components, "components", sep = " "))


# Load libraries ----------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}
packages<-c("Seurat", "monocle", "plyr", "dplyr", "colorRamps", "stringr", "future")
suppressMessages(ipak(packages))

if(!("MyPlotTrajectoryPackage" %in% installed.packages())){
  print("Installing MyPlotTrajectoryPackage from github")
  library(devtools)
  withr::with_libpaths(new = "./", install_github("efheuston/MyPlotTrajectoryPackage", "https://github.com/efheuston/MyPlotTrajectoryPackage.git"))
  library("MyPlotTrajectoryPackage", lib.loc = "./")
}else{
  library(MyPlotTrajectoryPackage)
}

###Concerned about importing nomalized data incorrectly from Seurat object



# Define script-specific functions ----------------------------------------

output_prefix<-gsub(x=option_arguments$seurat_object_filename, pattern = ".rds", replacement = "")
output_prefix<-gsub(x = output_prefix, pattern =  "_tsne", replacement =  "")

# define function to create pngs
png_plotFunction<-function(plot_to_make = plot_to_make, filename = filename, height = 800, width = 800){ # default window is 800x800
  png(filename = filename, height = height, width = width)
  plot(plot_to_make)
  dev.off()
}

cycle_plot_param<-function(plotting_function, cycle_parameter, the_object){
  # Create color pallete
  basic_color_palette<-basic_color_palette[1:25]
  new_length <-0
  if(length(cycle_parameter) > length(basic_color_palette)){
    new_length<-length(unique(my_object@ident)) - length(basic_color_palette)
  }
  new_length
  basic_color_palette<-c(basic_color_palette, primary.colors(new_length))

  if(plotting_function == "cluster"){
    png_plotFunction(my_plot_cell_clusters(the_object, 1, 2, color = cycle_parameter, point_colors = basic_color_palette) +
                       ggtitle(paste(output_prefix, "-", cycle_parameter, sep = "")),
                     filename = paste(output_prefix, "_clstr-", cycle_parameter,".png", sep = ""),
                     height = 1600,
                     width = 1600)
    png_plotFunction(my_plot_cell_clusters(the_object, 1, 2, color = cycle_parameter, point_colors = basic_color_palette) +
                       facet_wrap(as.vector(cycle_parameter)) +
                       ggtitle(paste(paste(output_prefix, "-", cycle_parameter, sep = ""))),
                     filename = paste(paste(output_prefix, "_clstr-", cycle_parameter,"FACET.png", sep = "")),
                     height = 1600,
                     width = 1600)
    png_plotFunction(plot_cell_clusters(the_object, 1, 2, color = "orig.ident") +
                       facet_wrap(as.vector(cycle_parameter)) +
                       ggtitle(paste(output_prefix, "-orig.identFACET", sep = "")),
                     filename = paste(output_prefix, "_clstr-orig.identFACET.png", sep = ""),
                     height = 1600,
                     width = 1600)

  }
  if(plotting_function == "trajectory"){
    png_plotFunction(my_plot_cell_trajectory(the_object,
                                             color_by = cycle_parameter,
                                             cell_size = 1,
                                             point_colors = basic_color_palette) + ggtitle(paste(output_prefix, "-", cycle_parameter, sep = "")),
                     filename = paste(output_prefix, "_trajectory-", cycle_parameter,".png", sep = ""),
                     height = 1600,
                     width = 1600)
    png_plotFunction(my_plot_cell_trajectory(the_object,
                                             color_by = cycle_parameter,
                                             cell_size = 1,
                                             point_colors = basic_color_palette) +
                       facet_wrap(as.vector(cycle_parameter)) + ggtitle(paste(output_prefix, "-", cycle_parameter, sep = "")),
                     filename = paste(output_prefix, "_trajectory-", cycle_parameter,"FACET.png", sep = ""),
                     height = 1600,
                     width = 1600)
    png_plotFunction(my_plot_cell_trajectory(the_object,
                                             color_by = cycle_parameter,
                                             cell_size = 1,
                                             point_colors = basic_color_palette) +
                       facet_wrap(~orig.ident) + ggtitle(paste(output_prefix, "-", cycle_parameter,"FACETorig.ident", sep = "")),
                     filename = paste(output_prefix, "_trajectory-", cycle_parameter,"FACETorig.ident.png", sep = ""),
                     height = 1600,
                     width = 1600)
    png_plotFunction(my_plot_cell_trajectory(the_object,
                                             color_by = "orig.ident",
                                             cell_size = 1,
                                             point_colors = basic_color_palette) +
                       facet_wrap(as.vector(cycle_parameter)) + ggtitle(paste(output_prefix, "-orig.identFACET", cycle_parameter, sep = "")),
                     filename = paste(output_prefix, "_trajectory-orig.identFACET", cycle_parameter,".png", sep = ""),
                     height = 1600,
                     width = 1600)
  }
  if(plotting_function == "qplot"){
    png_plotFunction(qplot(Total_nUMI, data = pData(the_object), color = as.factor(cycle_parameter), geom = "density") +
                       geom_vline(xintercept = upper_bound) +
                       geom_vline(xintercept = lower_bound) +
                       ggtitle(paste(output_prefix, "-", cycle_parameter, sep = "")),
                     filename = paste(output_prefix, "_nUMI_2SD-by", cycle_parameter,".png", sep = ""))

  }
}



# Load Seurat object as CDS -----------------------------------------------

if(grepl(x = option_arguments$seurat_object_filename, pattern = "rds$", ignore.case = T)){
  seurat_object<-readRDS(option_arguments$seurat_object_filename)
} else if (grepl(x = option_arguments$seurat_object_filename, pattern = "Robj$", ignore.case = T)){
  seurat_object<-get(load(option_arguments$seurat_object_filename))
} else {
  print("I don't recognize this type of file")
}

seurat_varGenes<-seurat_object@var.genes

if(option_arguments$color_by_seurat_res == TRUE){
  seurat_res<-colnames(seurat_object@meta.data)[grep(x=colnames(seurat_object@meta.data), pattern = "^res")]
}




# Define base colour palette ----------------------------------------------

# Create basic colour palette

basic_color_palette<-c("#cb4bbe",
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

# Adjust colour palette
adjust_palette_size <- function(object_length, basic_color_palette = basic_color_palette){
  if(length(unique(object_length)) > length(basic_color_palette)){
    new_length <- length(unique(object_length)) - length(basic_color_palette)
    basic_color_palette <- c(basic_color_palette, primary.colors(new_length))
  } else {
    basic_color_palette <- basic_color_palette
  }
}

# Initialize monocle object -----------------------------------------------

monocle_object<-importCDS(seurat_object, import_all = FALSE)


if(option_arguments$color_by_seurat_res == TRUE){
  adjust_palette_size(length(unique(seurat_object@ident)), basic_color_palette = basic_color_palette)
} else{
  basic_color_palette = basic_color_palette
}

remove(seurat_object)
gc()

monocle_object<-estimateSizeFactors(monocle_object)
monocle_object<-estimateDispersions(monocle_object)


# Filter low-quality cells ------------------------------------------------

if(option_arguments$perform_expression_filtering == TRUE){
  monocle_object<-detectGenes(monocle_object, min_expr = 0.1)
  expressed_genes<-row.names(subset(fData(monocle_object), num_cells_expressed >=10))
}

# can create a subset of valid cells using columns from pData
# might want to include qc info like Mapped.Fragments, or set valid_cells based on nUMI
# valid_cells<-row.names(subset(pData(monocle_object), nUMI > someNumber))
# monocle_object<-monocle_object[,valid_cells]

pData(monocle_object)$Total_nUMI<-Matrix::colSums(exprs(monocle_object))
monocle_object<-monocle_object[,pData(monocle_object)$Total_nUMI < 1e6]

lower_bound<-10^(mean(log10(pData(monocle_object)$Total_nUMI)) - 2*sd(log10(pData(monocle_object)$Total_nUMI)))
upper_bound<-10^(mean(log10(pData(monocle_object)$Total_nUMI)) + 2*sd(log10(pData(monocle_object)$Total_nUMI)))

# Plot UMI distributions and filter --------------------------------------------------

png_plotFunction(qplot(Total_nUMI, data = pData(monocle_object), color = orig.ident, geom = "density") +
                   geom_vline(xintercept = upper_bound) +
                   geom_vline(xintercept = lower_bound),
                 filename = paste(output_prefix, "_nUMI_2SD-byID.png", sep = ""))

try(
  if(option_arguments$color_by_seurat_res == TRUE){
    lapply(seurat_res[1:length(seurat_res)], plotting_function = "qplot", the_object = monocle_object, cycle_plot_param)
  }
)

if(grep("upper", option_arguments$UMI_bounded_filtering, ignore.case = TRUE)){
  monocle_object<-monocle_object[,pData(monocle_object)$Total_nUMI < upper_bound]
} else if(grep("both", option_arguments$UMI_bounded_filtering, ignore.case = TRUE)){
  monocle_object<-monocle_object[,pData(monocle_object)$Total_nUMI < upper_bound & pData(monocle_object)$Total_nUMI > lower_bound]
} else if(is.null(option_arguments$UMI_bounded_filtering) | grep("none", option_arguments$UMI_bounded_filtering, ignore.case = TRUE)){
  print("Not filtering based on UMI counts")
} else {
  print("No boundries requested")
}

monocle_object<-detectGenes(monocle_object, min_expr = 0.1)



# Cluster cells without marker genes --------------------------------------

disp_table<-dispersionTable(monocle_object)

monocle_clustering_genes<-subset(disp_table, mean_expression >= 0.1)
monocle_object<-setOrderingFilter(monocle_object, ordering_genes = monocle_clustering_genes$gene_id)

png_plotFunction(plot_ordering_genes(monocle_object), filename = paste(output_prefix, "_PlotOrderingGenes.png", sep = ""))
png_plotFunction(plot_pc_variance_explained(monocle_object, return_all = FALSE), filename = paste(output_prefix, "_PCVarXplained.png", sep = ""))

# Perform clustering ------------------------------------------------------

monocle_object<-reduceDimension(monocle_object,
                                max_components = option_arguments$max_components,
                                num_dim = option_arguments$num_dim,
                                reduction_method = 'tSNE',
                                residualModelFormulaStr = paste("~", paste(unlist(option_arguments$cca_variables), collapse = " + "), sep = ""),
                                verbose = TRUE)
monocle_object<-clusterCells(monocle_object)

png_plotFunction(plot_cell_clusters(monocle_object, 1, 2, color = "orig.ident"),
                 filename = paste(output_prefix, "_clstr-orig.ident.png", sep = ""),
                 height = 1600,
                 width = 1600)


try(
  if(option_arguments$color_by_seurat_res == TRUE){
    lapply(seurat_res[cycle_parameter = 1:length(seurat_res)], plotting_function = "cluster", the_object = monocle_object, cycle_plot_param)
  }
)

save.image(file=paste(output_prefix, ".RData", sep = ""))


# Perform trajectory analysis with monocle-defined genes ------------------

clustering_DEG_monocle<- differentialGeneTest(monocle_object[expressed_genes,], fullModelFormulaStr='~Cluster')
#current default is to take top 1000 clustering_DEG_monocle genes

monocleBased_orderedgenes<-row.names(clustering_DEG_monocle)[order(clustering_DEG_monocle$qval)][1:1000]
monocle_object<-setOrderingFilter(monocle_object, ordering_genes=monocleBased_orderedgenes)
monocle_object<-reduceDimension(monocle_object, method="DDRTree")
monocle_object<-orderCells(monocle_object)

# HERE we skip defining the root state because we want to leave that as an open question for now


# specify_root_state<-function(cds){
#   LSK_pops<-table(pData(cds)$State, pData(cds)$Type[,"LSK"])
#   as.numeric(names(LSK_pops)[which(LSK_pops==max(LSK_pops))])
# }



# Plot cell trajectories --------------------------------------------------



try(
  png_plotFunction(my_plot_cell_trajectory(monocle_object,
                                           cell_size = 1,
                                           point_colors = basic_color_palette),
                   filename = paste(output_prefix, "_trajectory.png", sep = ""))
)
try(
  png_plotFunction(my_plot_cell_trajectory(monocle_object,
                                           color_by = "orig.ident",
                                           cell_size = 1,
                                           point_colors = basic_color_palette),
                   filename = paste(output_prefix, "_trajectory-orig.ident.png", sep = ""))
)
try(
  png_plotFunction(my_plot_cell_trajectory(monocle_object,
                                           color_by = "orig.ident",
                                           cell_size = 1,
                                           point_colors = basic_color_palette) +
                     facet_wrap(~orig.ident),
                   filename = paste(output_prefix, "_trajectory-orig.identFACET.png", sep = ""))
)

try(
  if(option_arguments$cca_variables == TRUE){
    lapply(seurat_res[cycle_parameter = 1:length(seurat_res)], plotting_function = "trajectory", the_object = monocle_object, cycle_plot_param)
  }
)


saveRDS(monocle_object, file = paste(output_prefix, "_UnsupClustMonocle.rds", sep = ""))
save.image(file=paste(output_prefix, ".RData", sep = ""))

<<<<<<< HEAD
=======
#

# Basic differential expression analysis ----------------------------------

if(perform_de ==TRUE){
  # Hard-coding to use Seurat vargenes as markers... ?
  marker_genes <-row.names(subset(fData(monocle_object), gene_short_name %in% seurat_varGenes))
  diff_test_res <-differentialGeneTest(monocle_object[marker_genes,],
                                       cores = future::availableCores(), 
                                       fullModelFormulaStr = "~Cluster") # ASSUMING CLUSTER IS TEH CORRECT FULL MODEL FORMULA STRING
  sig_genes<-subset(diff_test_res, qval < 0.1) # Set for FDR < 10%
  sig_genes[,c("gene_short_name", "pval", "qval")]
  write.table(sig_genes, file = paste(output_prefix, "_differentiallyExpressed_seuratVar.txt", sep = ""), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}



>>>>>>> cd142d8... Fixes to clarify error handeling



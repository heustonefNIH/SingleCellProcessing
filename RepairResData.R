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



# Define repair_output function -------------------------------------------

repair_output<-function(x){
  my_object<-readRDS(x)
  print(paste("Working on ", x, sep = ""))
  output_prefix<-gsub(pattern = "_tsne.rds", replacement = "", x = x)
  
  # # Create color pallete
  # new_length <-0
  # if(length(unique(integrated_dba@ident)) > length(my_palette)){
  #   new_length<-length(unique(integrated_dba@ident)) - length(my_palette)
  # }
  # new_length
  # my_palette<-c(my_palette, primary.colors(new_length))
  # 
  # plot_title<-paste(output_prefix)
  # 
  # my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.return=TRUE)
  # my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
  # png(paste(output_prefix,"_tsne.png", sep = ""), width = 800, height = 800)
  # try(plot(my_tsne_plot))
  # dev.off()
  # 
  # my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.label = TRUE, label.size = 10, do.return=TRUE)
  # my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
  # png(paste(output_prefix,"_tsne-Labeled.png", sep = ""), width = 800, height = 800)
  # try(plot(my_tsne_plot))
  # dev.off()
  # 
  # my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, group.by = "orig.ident", do.return=TRUE)
  # my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
  # png(paste(output_prefix,"_tsnebyID.png", sep = ""), width = 800, height = 800)
  # try(plot(my_tsne_plot))
  # dev.off()
  # 
  # my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.label = TRUE, label.size = 10, group.by = "orig.ident", do.return=TRUE)
  # my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
  # png(paste(output_prefix,"_tsnebyID-Labeled.png", sep = ""), width = 800, height = 800)
  # try(plot(my_tsne_plot))
  # dev.off()

  my_object_markers<-FindAllMarkers(my_object, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  my_object_markers_top100 <- my_object_markers %>% group_by(cluster) %>% top_n(100, avg_logFC)
  write.table(my_object_markers_top100, file = paste(output_prefix, "_markers_top100.txt", sep = ""), sep="\t", quote=F)
  my_object_markers_all <- my_object_markers %>% group_by(cluster)
  write.table(my_object_markers_all, file = paste(output_prefix,"_markers_all.txt", sep = ""), sep="\t", quote=F)

  # Tabulate data -----------------------------------------------------------

  #Number of cells in each cluster
  table(my_object@ident)
  write.table(table(my_object@ident), file = paste(output_prefix,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")


  #Number of cells in each ID
  table(my_object@meta.data$orig.ident)
  write.table(table(my_object@meta.data$orig.ident), file = paste(output_prefix,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)

  #Proportion of cells in each cluster
  prop.table(x = table(my_object@ident))
  write.table(prop.table(x = table(my_object@ident)), file = paste(output_prefix,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)

  #Number of cells from each ID in each cluster
  table(my_object@ident, my_object@meta.data$orig.ident)
  write.table(table(my_object@ident, my_object@meta.data$orig.ident), file = paste(output_prefix,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)

  #Proportion of cells of each ID in each cluster
  prop.table(x = table(my_object@ident, my_object@meta.data$orig.ident))
  write.table(prop.table(x = table(my_object@ident, my_object@meta.data$orig.ident)), file = paste(output_prefix,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  
  #Average expression of each gene in each cluster
  avgexp<-AverageExpression(my_object)
  write.table(table(my_object@ident), file = paste(output_prefix, "_AvgXprsn.txt", sep = ""), row.names = TRUE, quote = FALSE, sep = "\t")

}



my_rdsdir<-"."

rds_fileList<-file.path(my_rdsdir, list.files(my_rdsdir, pattern=".rds$"))
lapply(rds_fileList[1:length(rds_fileList)], repair_output)



q()









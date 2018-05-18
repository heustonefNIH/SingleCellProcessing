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
  my_object_markers<-FindAllMarkers(my_object, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  my_object_markers_top100 <- my_object_markers %>% group_by(cluster) %>% top_n(100, avg_logFC)
  write.table(my_object_markers_top100, file = paste(output_prefix, "dim",max_pcs, "res",resolution,"_markers_top100.txt", sep = ""), sep="\t", quote=F)
  my_object_markers_all <- my_object_markers %>% group_by(cluster)
  write.table(my_object_markers_all, file = paste(output_prefix, "dim",max_pcs, "res",resolution,"_markers_all.txt", sep = ""), sep="\t", quote=F)
  
  # Tabulate data -----------------------------------------------------------
  
  #Number of cells in each cluster
  table(my_object@ident)
  write.table(table(my_object@ident), file = "test_table.txt", row.names = FALSE, quote = FALSE, sep = "\t")
  
  
  #Number of cells in each ID
  table(my_object@meta.data$orig.ident)
  write.table(table(my_object@meta.data$orig.ident), file = "test_table.txt", row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  
  #Proportion of cells in each cluster
  prop.table(x = table(my_object@ident))
  write.table(prop.table(x = table(my_object@ident)), file = "test_table.txt", row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  
  #Number of cells from each ID in each cluster
  table(my_object@ident, my_object@meta.data$orig.ident)
  write.table(table(my_object@ident, my_object@meta.data$orig.ident), file = "test_table.txt", row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  
  #Proportion of cells of each ID in each cluster
  prop.table(x = table(my_object@ident, my_object@meta.data$orig.ident))
  write.table(prop.table(x = table(my_object@ident, my_object@meta.data$orig.ident)), file = "test_table.txt", row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  
}



my_rdsdir<-"."

rds_fileList<-file.path(my_rdsdir, list.files(my_rdsdir, pattern=".rds$"))
lapply(rds_fileList[1:length(rds_fileList)], repair_output)



q()









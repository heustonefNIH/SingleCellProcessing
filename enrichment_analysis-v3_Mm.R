source("https://bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("org.Mm.eg.db")
# biocLite("org.Hs.eg.db")
# biocLite("GO.db")
library(limma)
library("org.Mm.eg.db")
# library("org.Hs.eg.db")
library(stringr)
library(GO.db)
##############################################Input file format matches Seurat output format##############################################

# #####NOTE: This program requires internet access

enrichment_function<-function(my_folder, number_of_terms, output_name_augment){
  create_GO_list<-function(Seurat_table, number_of_terms, base_name){
    cluster.range<-c(0:max(Seurat_table["cluster"]))
    xprsn_subset<-subset.data.frame(Seurat_table, avg_logFC>0, select = c(6:7))
    KEGG_and_GO_up_dir<-paste("./",base_name,"-KEGGandGO_up", sep="")
    dir.create(KEGG_and_GO_up_dir)
    for(i in cluster.range){
      cluster_name<-paste("cluster_", i, sep="")
      print(cluster_name)
      genelist<-subset.data.frame(xprsn_subset, cluster==i, select=2)
      table_list_entrezIDs <- unlist(mget(x=as.character(genelist$gene), envir=org.Mm.egALIAS2EG, ifnotfound=NA))
      table_list_goana <- goana(table_list_entrezIDs, species="Mm")
      table_list_GO <- topGO(table_list_goana, ontology=c("BP", "MF"), number=number_of_terms)
      try(write.table(table_list_GO, file=paste(KEGG_and_GO_up_dir,"/",base_name,"_", cluster_name, "-GO_up.txt", sep=""), quote = F, row.names = F, sep = "\t"))
      table_list_gokeg <- kegga(table_list_entrezIDs, species="Mm")
      table_list_keg <- topKEGG(table_list_gokeg)
      try(write.table(table_list_keg, file=paste(KEGG_and_GO_up_dir,"/",base_name,"_", cluster_name, "-KEGG_up.txt", sep=""), quote = F, row.names = F, sep = "\t"))
      remove(genelist, cluster_name, table_list_entrezIDs, table_list_goana, table_list_GO, table_list_gokeg, table_list_keg)
    }
  }
  create_GOdn_list<-function(Seurat_table, number_of_terms, base_name){
    cluster.range<-c(0:max(Seurat_table["cluster"]))
    xprsn_subset<-subset.data.frame(Seurat_table, avg_logFC<0, select = c(6:7))
    KEGG_and_GO_dn_dir<-paste("./",base_name,"-KEGGandGO_dn", sep="")
    dir.create(KEGG_and_GO_dn_dir)
    for(i in cluster.range){
      cluster_name<-paste("cluster_", i, sep="")
      print(cluster_name)
      genelist<-subset.data.frame(xprsn_subset, cluster==i, select=2)
      table_list_entrezIDs <- unlist(mget(x=as.character(genelist$gene), envir=org.Mm.egALIAS2EG, ifnotfound=NA))
      if(length(table_list_entrezIDs)<1){
        print(paste("No entries in ", cluster_name, sep = ""))
        next
      }
      table_list_goana <- goana(table_list_entrezIDs, species="Mm")
      table_list_GO <- topGO(table_list_goana, ontology=c("BP", "MF"), number=number_of_terms)
      try(write.table(table_list_GO, file=paste(KEGG_and_GO_dn_dir,"/",base_name,"_", cluster_name, "-GO_dn.txt", sep=""), quote = F, row.names = F, sep = "\t"))
      table_list_gokeg <- kegga(table_list_entrezIDs, species="Mm")
      table_list_keg <- topKEGG(table_list_gokeg)
      try(write.table(table_list_keg, file=paste(KEGG_and_GO_dn_dir,"/",base_name,"_", cluster_name, "-KEGG_dn.txt", sep=""), quote = F, row.names = F, sep = "\t"))
      remove(genelist, cluster_name, table_list_entrezIDs, table_list_goana, table_list_GO, table_list_gokeg, table_list_keg)
    }
  }
  for(my_file in my_folder){
    print(basename(my_file))
    file_name<-paste(gsub(pattern = ".txt", replacement = "", x = basename(my_file)))
    print(paste("Creating files with names ", file_name, sep=""))
    Seurat_gene_list<-read.delim(file = my_file, sep = "\t", header=T, na.strings = c("", "NA"))
    Seurat_gene_list<-na.omit(Seurat_gene_list)
    create_GO_list(Seurat_gene_list, number_of_terms, file_name)
    create_GOdn_list(Seurat_gene_list, number_of_terms, file_name)
    remove(Seurat_gene_list, file_name, my_file)
  }
}


my_folder<-list.files(path = "/Users/heustonef/Desktop/10XGenomicsData/msAggr/msAggr_mtdim34/",pattern="_top100.txt$", full.names = T)
enrichment_function(my_folder, output_name_augment="", number_of_terms=100L)

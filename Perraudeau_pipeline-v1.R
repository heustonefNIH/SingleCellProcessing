# 20180831 - Running Perraudeau analysis on msAggr
# Raw data are imported from 
# Load libraries ----------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}
packages<-c("BiocParallel", "clusterExperiment", "scone", "zinbwave", "slingshot", "doParallel", "gam", "RColorBrewer", "Seurat", "dplyr", "colorRamps")
suppressMessages(ipak(packages))

msaggr<-readRDS("/data/heustonef/10XGenomics/msAggr/20180517_msAggrdim36res2.5_tsne.rds")

slotNames(msaggr)
colnames(msaggr@meta.data)


msaggr_counts<-as.matrix(msaggr@raw.data)
dim(msaggr_counts)
msaggr_counts<-na.omit(msaggr_counts)
msaggr_counts<-msaggr_counts[rowSums(msaggr_counts)>0,]
dim(msaggr_counts)
summary(rowSums(msaggr_counts))

colnames(msaggr@meta.data)
range(msaggr@meta.data$nGene)
range(msaggr@meta.data$nUMI)

# qc currently available is basically nGene and nUMI from meta.data
# In future could try using nUMI boundries like in Monocle, mt data

# Instead of "Batch", use "orig.ident"
msaggr_origIdent<-msaggr@meta.data$orig.ident
msaggr_ClusterLabels<-data.frame(cbind(names(msaggr@ident), msaggr@ident), row.names = 1:length(msaggr@ident))

msaggr_matched<-match(colnames(msaggr_counts), msaggr_ClusterLabels)

msaggr_MetaData<-data.frame("OrigIdent" = msaggr_origIdent,
                        "nGenes" = msaggr@meta.data$nGene,
                        "nUMI" = msaggr@meta.data$nUMI,
                        "res0.6" = msaggr@meta.data$res.0.6,
                        "res0.8" = msaggr@meta.data$res.0.8,
                        "res1.0" = msaggr@meta.data$res.1,
                        "res1.5" = msaggr@meta.data$res.1.5,
                        "res2.0" = msaggr@meta.data$res.2,
                        "res2.5" = msaggr@meta.data$res.2.5)

msaggr_se<-SummarizedExperiment(assays = list(counts = msaggr_counts), 
                                colData = msaggr_MetaData)
msaggr_se

data("housekeeping")

msaggr_hk<-rownames(msaggr_se)[toupper(rownames(msaggr_se)) %in% housekeeping$V1]
msaggr_mfilt <-metric_sample_filter(assay(msaggr_se), 
                                    pos_controls = rownames(msaggr_se) %in% msaggr_hk,
                                    zcut = 3, mixture = FALSE, 
                                    plot = FALSE)
png(filename = "20180831_msaggr_mfilt.png", height = 1600, width = 1600)
msaggr_mfilt <-metric_sample_filter(assay(msaggr_se), 
                                    pos_controls = rownames(msaggr_se) %in% msaggr_hk,
                                    zcut = 3, mixture = FALSE, 
                                    plot = TRUE)
dev.off()

# Will need to check results to see if applying this filter makes sense
# Looks like it's having the same problem as before, with a threshhold being drawn @ nREADS = 4.4
# Until I figure this out, not going to apply this filter.
# Does it make sesne to fake nreads = nGene?


msaggr_vars <-rowVars(log1p(assay(msaggr_se)))
names(msaggr_vars) <-rownames(msaggr_se)
msaggr_vars<-sort(msaggr_vars, decreasing = TRUE)
summary(msaggr_vars)

min(msaggr_vars[1:1000])

# looks like 1000 is fine

msaggr_core <- msaggr_se[names(msaggr_vars)[1:1000],]
msaggr_core

msaggr_core_origIdent<-colData(msaggr_core)$OrigIdent
msaggr_core_sClusID<-colData(msaggr_core)$res2.5

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

msaggr_col_sClusID<-c(basic_color_palette, primary.colors(13))
names(msaggr_col_sClusID) = unique(msaggr_core_sClusID)
table(msaggr_core_sClusID)

msaggr_col_origIdent<-basic_color_palette[1:4]
names(msaggr_col_origIdent) = unique(msaggr_core_origIdent)

remove(msaggr)
register(MulticoreParam(detectCores()))


print(system.time(msaggr_se <- zinbwave(msaggr_core, K = 36,
                                 residuals = TRUE,
                                 normalizedValues = TRUE)))

save(msaggr_se, file = "20180831_msaggr_se_afterZinbwave.rds")
save.image(file = "20180831_msaggr_se_afterZinbwave.RData")

msaggr_norm<-assays(msaggr_se)$normalizedValues

msaggr_norm_order<-msaggr_norm[,order(as.numeric(msaggr_core_sClusID))]
msaggr_colsClust_order<-msaggr_col_sClusID[msaggr_core_sClusID(as.numeric(msaggr_core_sClusID))]
png(filename = "20180831_msaggr_normBoxplot_sClustID.png", height = 1600, width = 1600)
boxplot(msaggr_norm_order, col = msaggr_colsClust_order, staplewex = 0, outline = 0, border = msaggr_colsClust_order, xaxt = "n", ylab = "Expression measure")
abline(h=0)
dev.off()

msaggr_norm_orderOI<-msaggr_norm_orderOI[,order(as.numeric(msaggr_core_origIdent))]
msaggr_colOI_order<-msaggr_col_origIdent[msaggr_core_origIdent(as.numeric(msaggr_core_origIdent))]
png(filename = "20180831_msaggr_normBoxplot_origIdent.png", height = 1600, width = 1600)
boxplot(msaggr_norm_orderOI, col = msaggr_colOI_order, staplewex = 0, outline = 0, border = msaggr_colOI_order, xaxt = "n", ylab = "Expression measure")
abline(h=0)
dev.off()

msaggr_pca<-prcomp(t(msaggr_norm))
par(mfrow = c(1,2))
png(filename = "20180831_msaggr_pca.png", height = 1600, width = 1600)
plot(msaggr_pca$x, col = msaggr_col_origIdent[msaggr_core_origIdent], pch = 20, main = "")
plot(msaggr_pca$x, col = msaggr_col_sClusID[msaggr_core_sClusID], pch = 20, main = "")
dev.off()
par(mfrow = c(1,1))

msaggr_redDim<-reducedDim(msaggr_se, "zinbwave")
msaggr_redDim<-as.matrix(msaggr_redDim)
msaggr_dist<-dist(msaggr_redDim)
msaggr_fit<-cmdscale(msaggr_dist, eig = TRUE, k = 2) # where k = maximum dimensions of space that data will be represented in

png(filename = "20180831_msaggr_DimRedctn_36_sClust.png", height = 1600, width = 1600)
plot(msaggr_fit$points, col = msaggr_col_sClusID[as.character(msaggr_core_sClusID)], main = "", 
     pch = 20, xlab = "Component 1", ylab = "Component 2")
legend(x = "topleft", levend = unique(names(msaggr_col_sClusID)), cex = .5, 
       fill = unique(msaggr_col_sClusID), title = "Colored by seurat msaggr_dim36res2.5_ClusterID")
dev.off()

png(filename = "20180831_msaggr_DimRedctn_36_OrigIdent.png", height = 1600, width = 1600)
plot(msaggr_fit$points, col = msaggr_col_origIdent[as.character(msaggr_core_origIdent)], main = "", 
     pch = 20, xlab = "Component 1", ylab = "Component 2")
legend(x = "topleft", levend = unique(names(msaggr_col_origIdent)), cex = .5, 
       fill = unique(msaggr_col_origIdent), title = "Colored by seurat msaggr_dim36res2.5_OrigIdent")
dev.off()

save.history("20180831_Perraudeau.Rhistory")
save.image("20180831_Perraudeau.RData")















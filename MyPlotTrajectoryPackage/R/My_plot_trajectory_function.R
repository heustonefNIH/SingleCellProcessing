# Modified from Monocle2 in torder to change geom_point colours
library(ggplot2)
library(reshape2)
library(igraph)
library(plyr)
library(dplyr)
library(viridis)
# Monocle_theme_opts ------------------------------------------------------

#' @title Loads theme
#'
#' @description Accepts a subset of a CellDataSet and an attribute to group cells by,
#' and produces one or more ggplot2 objects that plots the level of expression for
#' each group of cells.
#' @export


monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}




# my_plot_cell_trajectory ------------------------------------------------------
#' @title Plots cell trajectory
#'
#' @description Accepts a subset of a CellDataSet and an attribute to group cells by,
#' and produces one or more ggplot2 objects that plots the level of expression for
#' each group of cells.
#' @export

my_plot_cell_trajectory <- function(cds,
            x=1,
            y=2,
            color_by="State",
            show_tree=TRUE,
            show_backbone=TRUE,
            backbone_color="black",
            markers=NULL,
            markers_linear = FALSE,
            show_cell_names=FALSE,
            show_state_number = FALSE,
            cell_size=1.5,
            cell_link_size=0.75,
            cell_name_size=2,
            state_number_size = 2.9,
            show_branch_points=TRUE,
            point_colors = c("#cb4bbe", "lightskyblue","grey37", "#53ba48","moccasin","#dcca44","#502d70","#afd671","#cb4471","#69da99","#d44732","#6fd3cf","#5d2539","#cfd0a0","blue","#d2883b","#6793c0","#898938","#c98ac2","yellow","#c4c0cc","#7d3d23","#00a5ff","#d68d7a","#a2c1a3"),
            theta = 0){
  gene_short_name <- NA
  sample_name <- NA
  sample_state <- pData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA

  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)

  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }

  if (cds@dim_reduce_type == "ICA"){
    reduced_dim_coords <- reducedDimS(cds)
  }else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree") ){
    reduced_dim_coords <- reducedDimK(cds)
  }else {
    stop("Error: unrecognized dimensionality reduction method.")
  }

  ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(x,y),]))
  colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")

  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_df$sample_state <- row.names(ica_space_df)
  #ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  #print(ica_space_with_state_df)
  dp_mst <- minSpanningTree(cds)

  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }

  edge_list <- as.data.frame(get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")

  edge_df <- merge(ica_space_df, edge_list, by.x="sample_name", by.y="source", all=TRUE)
  #edge_df <- ica_space_df
  edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="source_prin_graph_dim_1", "prin_graph_dim_2"="source_prin_graph_dim_2"))
  edge_df <- merge(edge_df, ica_space_df[,c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2")], by.x="target", by.y="sample_name", all=TRUE)
  edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="target_prin_graph_dim_1", "prin_graph_dim_2"="target_prin_graph_dim_2"))

  S_matrix <- reducedDimS(cds)
  data_df <- data.frame(t(S_matrix[c(x,y),]))
  data_df <- cbind(data_df, sample_state)
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")

  return_rotation_mat <- function(theta) {
    theta <- theta / 180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  }

  tmp <- return_rotation_mat(theta) %*% t(as.matrix(data_df[, c(2, 3)]))
  data_df$data_dim_1 <- tmp[1, ]
  data_df$data_dim_2 <- tmp[2, ]

  tmp <- return_rotation_mat(theta = theta) %*% t(as.matrix(edge_df[, c('source_prin_graph_dim_1', 'source_prin_graph_dim_2')]))
  edge_df$source_prin_graph_dim_1 <- tmp[1, ]
  edge_df$source_prin_graph_dim_2 <- tmp[2, ]

  tmp <- return_rotation_mat(theta) %*% t(as.matrix(edge_df[, c('target_prin_graph_dim_1', 'target_prin_graph_dim_2')]))
  edge_df$target_prin_graph_dim_1 <- tmp[1, ]
  edge_df$target_prin_graph_dim_2 <- tmp[2, ]

  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
      colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
    if(markers_linear){
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size= (value * 0.1))) + facet_wrap(~feature_label)
    } else {
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size=log10(value + 0.1))) + facet_wrap(~feature_label)
    }
  }else{
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
  }
  if (show_tree){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=cell_link_size, linetype="solid", na.rm=TRUE, data=edge_df)
  }

  # FIXME: setting size here overrides the marker expression funtionality.
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0) {
    g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
  } else {
    g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE)
  }


  if (show_branch_points && cds@dim_reduce_type == 'DDRTree'){
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[,c("sample_name", "source_prin_graph_dim_1", "source_prin_graph_dim_2")]
    branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, mst_branch_nodes)
    branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx), ]

    g <- g + geom_point(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2"),
                        size=5, na.rm=TRUE, data=branch_point_df) +
      scale_color_manual(values = point_colors) +
      geom_text(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", label="branch_point_idx"),
                size=4, color="white", na.rm=TRUE, data=branch_point_df)
  }
  if (show_cell_names){
    g <- g +geom_text(aes(label=sample_name), size=cell_name_size)
  }
  if (show_state_number){
    g <- g + geom_text(aes(label = sample_state), size = state_number_size)
  }

  g <- g +
    # scale_color_brewer(palette="Set1") +
    monocle_theme_opts() +
    xlab(paste("Component", x)) +
    ylab(paste("Component", y)) +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
  g
}



# my_plot_genes_violin -------------------------------------------------------

#' @title Plots expression for one or more genes as a violin plot
#'
#' @description Accepts a subset of a CellDataSet and an attribute to group cells by,
#' and produces one or more ggplot2 objects that plots the level of expression for
#' each group of cells.
#'
#' @param cds_subset CellDataSet for the experiment
#' @param grouping the cell attribute (e.g. the column of pData(cds)) to group cells by on the horizontal axis
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to be used to color each cell
#' @param plot_trend whether to plot a trendline tracking the average expression across the horizontal axis.
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr Whether to transform expression into relative values
#' @param log_scale a boolean that determines whether or not to scale data logarithmically
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics sizeFactors
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' my_genes <- HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("ACTA1", "ID1", "CCNB2"))),]
#' plot_genes_violin(my_genes, grouping="Hours", ncol=2, min_expr=0.1)
#' }
my_plot_genes_violin <- function (cds_subset, grouping = "State", min_expr = NULL, cell_size = 0.75,
                              nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL,
                              plot_trend = FALSE, label_by_short_name = TRUE, relative_expr = TRUE,
                              log_scale = TRUE, point_colors = c("#cb4bbe", "grey37", "#53ba48","moccasin","#dcca44","#502d70","#afd671","#cb4471","#69da99","#d44732","#6fd3cf","#5d2539","#cfd0a0","blue","#d2883b","#6793c0","#898938","#c98ac2","yellow","#c4c0cc","#7d3d23","#00a5ff","#d68d7a","#a2c1a3"))
{
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial",
                                                 "negbinomial.size")) {
    integer_expression = TRUE
  }
  else {
    integer_expression = FALSE
    relative_expr = TRUE
  }
  if (integer_expression) {
    cds_exprs = exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs = Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    #cds_exprs = reshape2::melt(round(as.matrix(cds_exprs)))
    cds_exprs = reshape2::melt(as.matrix(cds_exprs))
  }
  else {
    cds_exprs = exprs(cds_subset)
    cds_exprs = reshape2::melt(as.matrix(cds_exprs))
  }
  if (is.null(min_expr)) {
    min_expr = cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) = c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] = min_expr
  cds_pData = pData(cds_subset)

  #
  # # Custom bit for adding in a group for
  # if(! is.null(show_combined)) {
  #   for(combine_gene in show_combined) {
  #     cds_pData_all = subset(cds_pData, gene == combine_gene)
  #     cds_pData_all[, grouping] = paste("All", combine_gene)
  #     cds_pData = rbind(cds_pData, cds_pData_all)
  #   }
  # }

  cds_fData = fData(cds_subset)
  cds_exprs = merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs = merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  cds_exprs$adjusted_expression = log10(cds_exprs$expression)




  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label = cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] = cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label = cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label = cds_exprs$f_id
  }
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label = factor(cds_exprs$feature_label,
                                     levels = panel_order)
  }
  q = ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs)
  if (is.null(color_by) == FALSE) {
    q = q + geom_violin(aes_string(fill = color_by)) + scale_fill_manual(values = point_colors)
  }
  else {
    q = q + geom_violin() + scale_fill_manual(values = point_colors)
  }
  if (plot_trend == TRUE) {
    q = q + stat_summary(fun.data = "mean_cl_boot",
                         size = 0.2)
    q = q + stat_summary(aes_string(x = grouping, y = "expression",
                                    group = color_by), fun.data = "mean_cl_boot",
                         size = 0.2, geom = "line")
  }
  q = q + facet_wrap(~feature_label, nrow = nrow,
                     ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
     q = q + expand_limits(y = c(min_expr, 1))
  }


  q = q + ylab("Expression") + xlab(grouping)

  if (log_scale == TRUE){

    q = q + scale_y_log10()
  }
  q
}

# my_plot_genes_jitter -------------------------------------------------------

#' Plots expression for one or more genes as a jittered, grouped points
#'
#' @description Accepts a subset of a CellDataSet and an attribute to group cells by,
#' and produces one or more ggplot2 objects that plots the level of expression for
#' each group of cells.
#'
#' @param cds_subset CellDataSet for the experiment
#' @param grouping the cell attribute (e.g. the column of pData(cds)) to group cells by on the horizontal axis
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to be used to color each cell
#' @param plot_trend whether to plot a trendline tracking the average expression across the horizontal axis.
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr Whether to transform expression into relative values
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics sizeFactors
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' my_genes <- HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MYOG", "ID1", "CCNB2"))),]
#' plot_genes_jitter(my_genes, grouping="Media", ncol=2)
#' }
my_plot_genes_jitter <- function(cds_subset,
                              grouping = "State",
                              min_expr=NULL,
                              cell_size=0.75,
                              nrow=NULL,
                              ncol=1,
                              panel_order=NULL,
                              color_by=NULL,
                              plot_trend=FALSE,
                              label_by_short_name=TRUE,
                              point_colors = c("#cb4bbe", "lightskyblue","grey37", "#53ba48","moccasin","#dcca44","#502d70","#afd671","#cb4471","#69da99","#d44732","#6fd3cf","#5d2539","#cfd0a0","blue","#d2883b","#6793c0","#898938","#c98ac2","yellow","#c4c0cc","#7d3d23","#00a5ff","#d68d7a","#a2c1a3"),
                              relative_expr=TRUE){

  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){



    integer_expression <- TRUE
  }else{
    integer_expression <- FALSE
    relative_expr <- TRUE
  }

  if (integer_expression)
  {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset)))
      {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs) / sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  }else{
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  }
  if (is.null(min_expr)){
    min_expr <- cds_subset@lowerDetectionLimit
  }

  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)

  cds_exprs <- merge(cds_exprs, cds_fData, by.x="f_id", by.y="row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x="Cell", by.y="row.names")

  cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  #cds_exprs$adjusted_expression <- log10(cds_exprs$adjusted_expression + abs(rnorm(nrow(cds_exprs), min_expr, sqrt(min_expr))))

  if (label_by_short_name == TRUE){
    if (is.null(cds_exprs$gene_short_name) == FALSE){
      cds_exprs$feature_label <- cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)]  <- cds_exprs$f_id
    }else{
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }else{
    cds_exprs$feature_label <- cds_exprs$f_id
  }

  #print (head(cds_exprs))

  if (is.null(panel_order) == FALSE)
  {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label, levels=panel_order)
  }

  q <- ggplot(aes_string(x=grouping, y="expression"), data=cds_exprs) +
    scale_color_manual(values = point_colors)


  if (is.null(color_by) == FALSE){
    q <- q + geom_jitter(aes_string(color=color_by), size=I(cell_size))
  }else{
    q <- q + geom_jitter(size=I(cell_size))
  }
  if (plot_trend == TRUE){
    q <- q + stat_summary(aes_string(color=color_by), fun.data = "mean_cl_boot", size=0.35)
    q <- q + stat_summary(aes_string(x=grouping, y="expression", color=color_by, group=color_by), fun.data = "mean_cl_boot", size=0.35, geom="line")
  }

  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow=nrow, ncol=ncol, scales="free_y")

  # Need this to guard against plotting failures caused by non-expressed genes
  if (min_expr < 1)
  {
    q <- q + expand_limits(y=c(min_expr, 1))
  }

  q <- q + ylab("Expression") + xlab(grouping)
  q <- q + monocle_theme_opts()
  q
}


# my_plot_genes_positive_cells -----------------------------------------------

#' Plots the number of cells expressing one or more genes as a barplot
#'
#'  @description Accetps a CellDataSet and a parameter,"grouping", used for dividing cells into groups.
#'  Returns one or more bar graphs (one graph for each gene in the CellDataSet).
#'  Each graph shows the percentage of cells that express a gene in the in the CellDataSet for
#'  each sub-group of cells created by "grouping".
#'
#'  Let's say the CellDataSet passed in included genes A, B, and C and the "grouping parameter divided
#'  all of the cells into three groups called X, Y, and Z. Then three graphs would be produced called A,
#'  B, and C. In the A graph there would be three bars one for X, one for Y, and one for Z. So X bar in the
#'  A graph would show the percentage of cells in the X group that express gene A.
#'
#' @param cds_subset CellDataSet for the experiment
#' @param grouping the cell attribute (e.g. the column of pData(cds)) to group cells by on the horizontal axis
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param plot_as_fraction whether to show the percent instead of the number of cells expressing each gene
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr Whether to transform expression into relative values
#' @param plot_limits A pair of number specifying the limits of the y axis. If NULL, scale to the range of the data.
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics sizeFactors
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' MYOG_ID1 <- HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MYOG", "ID1"))),]
#' plot_genes_positive_cells(MYOG_ID1, grouping="Media", ncol=2)
#' }
my_plot_genes_positive_cells <- function(cds_subset,
                                      grouping = "State",
                                      min_expr=0.1,
                                      nrow=NULL,
                                      ncol=1,
                                      panel_order=NULL,
                                      plot_as_fraction=TRUE,
                                      label_by_short_name=TRUE,
                                      relative_expr=TRUE,
                                      point_colors = c("#cb4bbe", "lightskyblue","grey37", "#53ba48","moccasin","#dcca44","#502d70","#afd671","#cb4471","#69da99","#d44732","#6fd3cf","#5d2539","#cfd0a0","blue","#d2883b","#6793c0","#898938","#c98ac2","yellow","#c4c0cc","#7d3d23","#00a5ff","#d68d7a","#a2c1a3"),
                                      plot_limits=c(0,100)){

  percent <- NULL

  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
    integer_expression <- TRUE
  }else{
    integer_expression <- FALSE
    relative_expr <- TRUE
  }

  if (integer_expression)
  {
    marker_exprs <- exprs(cds_subset)
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset)))
      {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      marker_exprs <- Matrix::t(Matrix::t(marker_exprs) / sizeFactors(cds_subset))
    }
    marker_exprs_melted <- reshape2::melt(round(as.matrix(marker_exprs)))
  }else{
    marker_exprs_melted <- reshape2::melt(exprs(marker_exprs))
  }

  colnames(marker_exprs_melted) <- c("f_id", "Cell", "expression")

  marker_exprs_melted <- merge(marker_exprs_melted, pData(cds_subset), by.x="Cell", by.y="row.names")
  marker_exprs_melted <- merge(marker_exprs_melted, fData(cds_subset), by.x="f_id", by.y="row.names")

  if (label_by_short_name == TRUE){
    if (is.null(marker_exprs_melted$gene_short_name) == FALSE){
      marker_exprs_melted$feature_label <- marker_exprs_melted$gene_short_name
      marker_exprs_melted$feature_label[is.na(marker_exprs_melted$feature_label)]  <- marker_exprs_melted$f_id
    }else{
      marker_exprs_melted$feature_label <- marker_exprs_melted$f_id
    }
  }else{
    marker_exprs_melted$feature_label <- marker_exprs_melted$f_id
  }

  if (is.null(panel_order) == FALSE)
  {
    marker_exprs_melted$feature_label <- factor(marker_exprs_melted$feature_label, levels=panel_order)
  }

  marker_counts <- plyr::ddply(marker_exprs_melted, c("feature_label", grouping), function(x) {
    data.frame(target=sum(x$expression > min_expr),
               target_fraction=sum(x$expression > min_expr)/nrow(x)) } )

  #print (head(marker_counts))
  if (plot_as_fraction){
    marker_counts$target_fraction <- marker_counts$target_fraction * 100
    qp <- ggplot(aes_string(x=grouping, y="target_fraction", fill=grouping), data=marker_counts) +
      scale_color_manual(values = point_colors) +
      ylab("Cells (percent)")
    if (is.null(plot_limits) == FALSE)
      qp <- qp + scale_y_continuous(limits=plot_limits)
  }else{
    qp <- ggplot(aes_string(x=grouping, y="target", fill=grouping), data=marker_counts) +
      scale_color_manual(values = point_colors) +
      ylab("Cells")
  }

  qp <- qp + facet_wrap(~feature_label, nrow=nrow, ncol=ncol, scales="free_y") + scale_color_manual(values = point_colors)
  qp <-  qp + geom_bar(stat="identity") + monocle_theme_opts()

  return(qp)
}




# my_plot_genes_in_pseudotime ------------------------------------------------

#' Plots expression for one or more genes as a function of pseudotime
#'
#' @description Plots expression for one or more genes as a function of pseudotime.
#' Plotting allows you determine if the ordering produced by orderCells() is correct
#' and it does not need to be flipped using the "reverse" flag in orderCells
#'
#' @param cds_subset CellDataSet for the experiment
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to be used to color each cell
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr Whether to transform expression into relative values
#' @param vertical_jitter A value passed to ggplot to jitter the points in the vertical dimension. Prevents overplotting, and is particularly helpful for rounded transcript count data.
#' @param horizontal_jitter A value passed to ggplot to jitter the points in the horizontal dimension. Prevents overplotting, and is particularly helpful for rounded transcript count data.
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply .
#' @importFrom reshape2 melt
#' @importFrom ggplot2 Position
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' my_genes <- row.names(subset(fData(HSMM), gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
#' cds_subset <- HSMM[my_genes,]
#' plot_genes_in_pseudotime(cds_subset, color_by="Time")
#' }
my_plot_genes_in_pseudotime <-function(cds_subset,
                                    min_expr=NULL,
                                    cell_size=0.75,
                                    nrow=NULL,
                                    ncol=1,
                                    panel_order=NULL,
                                    color_by="State",
                                    trend_formula="~ sm.ns(Pseudotime, df=3)",
                                    label_by_short_name=TRUE,
                                    relative_expr=TRUE,
                                    point_colors = c("#cb4bbe", "lightskyblue","grey37", "#53ba48","moccasin","#dcca44","#502d70","#afd671","#cb4471","#69da99","#d44732","#6fd3cf","#5d2539","#cfd0a0","blue","#d2883b","#6793c0","#898938","#c98ac2","yellow","#c4c0cc","#7d3d23","#00a5ff","#d68d7a","#a2c1a3"),
                                    vertical_jitter=NULL,
                                    horizontal_jitter=NULL){

  f_id <- NA
  Cell <- NA
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
        relative_expr <- TRUE
    }
    if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
        }
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    }
    else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    #cds_exprs$f_id <- as.character(cds_exprs$f_id)
    #cds_exprs$Cell <- as.character(cds_exprs$Cell)

    if (integer_expression) {
        cds_exprs$adjusted_expression <- cds_exprs$expression
    }
    else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
    # trend_formula <- paste("adjusted_expression", trend_formula,
    #     sep = "")
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)

    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula,
                        relative_expr = T, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    expectation <- ddply(cds_exprs, .(f_id, Cell), function(x) data.frame("expectation"=model_expectation[x$f_id, x$Cell]))
    cds_exprs <- merge(cds_exprs, expectation)
    #cds_exprs$expectation <- expectation#apply(cds_exprs,1, function(x) model_expectation[x$f_id, x$Cell])

    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    if (is.null(panel_order) == FALSE) {
      cds_exprs$feature_label <- factor(cds_exprs$feature_label,
            levels = panel_order)
    }
    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    if (is.null(color_by) == FALSE) {
        q <- q + geom_point(aes_string(color = color_by), size = I(cell_size), position=position_jitter(horizontal_jitter, vertical_jitter)) +
          scale_color_manual(values = point_colors)
    }
    else {
        q <- q + geom_point(size = I(cell_size), position=position_jitter(horizontal_jitter, vertical_jitter)) +
          scale_color_manual(values = point_colors)
    }

    q <- q + geom_line(aes(x = Pseudotime, y = expectation), data = cds_exprs) + scale_color_manual(values = point_colors)

    q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow,
        ncol = ncol, scales = "free_y")
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    if (relative_expr) {
        q <- q + ylab("Relative Expression")
    }
    else {
        q <- q + ylab("Absolute Expression")
    }
    q <- q + xlab("Pseudo-time")
    q <- q + monocle_theme_opts()
    q
}


# my_plot_genes_branched_pseudotime ------------------------------------------

#' Plot the branch genes in pseduotime with separate branch curves.
#'
#' @description Works similarly to plot_genes_in_psuedotime esceptit shows
#' one kinetic trend for each lineage.
#'
#' @details This plotting function is used to make the branching plots for a branch dependent gene goes through the progenitor state
#' and bifurcating into two distinct branchs (Similar to the pitch-fork bifurcation in dynamic systems). In order to make the
#' bifurcation plot, we first duplicated the progenitor states and by default stretch each branch into maturation level 0-100.
#' Then we fit two nature spline curves for each branchs using VGAM package.
#'
#' @param cds CellDataSet for the experiment
#' @param branch_states The states for two branching branchs
#' @param branch_point The ID of the branch point to analyze. Can only be used when reduceDimension is called with method = "DDRTree".
#' @param branch_labels The names for each branching branch
#' @param method The method to draw the curve for the gene expression branching pattern, either loess ('loess') or VGLM fitting ('fitting')
#' @param min_expr The minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size The size (in points) of each cell used in the plot
#' @param nrow Number of columns used to layout the faceted cluster panels
#' @param ncol Number of columns used to layout the faceted cluster panels
#' @param panel_order The a character vector of gene short names (or IDs, if that's what you're using), specifying order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param color_by The cell attribute (e.g. the column of pData(cds)) to be used to color each cell
#' @param expression_curve_linetype_by The cell attribute (e.g. the column of pData(cds)) to be used for the linetype of each branch curve
#' @param trend_formula The model formula to be used for fitting the expression trend over pseudotime
#' @param reducedModelFormulaStr A formula specifying a null model. If used, the plot shows a p value from the likelihood ratio test that uses trend_formula as the full model
#' @param label_by_short_name Whether to label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr Whether or not the plot should use relative expression values (only relevant for CellDataSets using transcript counts)
#' @param ... Additional arguments passed on to branchTest. Only used when reducedModelFormulaStr is not NULL.
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics sizeFactors
#' @export
my_plot_genes_branched_pseudotime <- function (cds,
                                            branch_states = NULL,
                                            branch_point=1,
                                            branch_labels = NULL,
                                            method = "fitting",
                                            min_expr = NULL,
                                            cell_size = 0.75,
                                            nrow = NULL,
                                            ncol = 1,
                                            panel_order = NULL,
                                            color_by = "State",
                                            expression_curve_linetype_by = "Branch",
                                            trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch",
                                            reducedModelFormulaStr = NULL,
                                            label_by_short_name = TRUE,
                                            relative_expr = TRUE,
                                            point_colors = c("#cb4bbe", "lightskyblue","grey37", "#53ba48","moccasin","#dcca44","#502d70","#afd671","#cb4471","#69da99","#d44732","#6fd3cf","#5d2539","#cfd0a0","blue","#d2883b","#6793c0","#898938","#c98ac2","yellow","#c4c0cc","#7d3d23","#00a5ff","#d68d7a","#a2c1a3"),
                                            #gene_pairs = NULL,
                                            ...)
{
  Branch <- NA
  if (is.null(reducedModelFormulaStr) == FALSE) {
    pval_df <- branchTest(cds,
                          branch_states=branch_states,
                          branch_point=branch_point,
                          fullModelFormulaStr = trend_formula,
                          reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)",
                          ...)
    fData(cds)[, "pval"] <- pval_df[row.names(cds), 'pval']
  }
  if("Branch" %in% all.vars(terms(as.formula(trend_formula)))) { #only when Branch is in the model formula we will duplicate the "progenitor" cells
    cds_subset <- buildBranchCellDataSet(cds = cds,
                                         branch_states = branch_states,
                                         branch_point=branch_point,
                                         branch_labels = branch_labels,
                                         progenitor_method = 'duplicate',
                                         ...)
  }
  else {
    cds_subset <- cds
    pData(cds_subset)$Branch <- pData(cds_subset)$State
  }
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
  }
  if (integer_expression) {
    CM <- exprs(cds_subset)
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      CM <- Matrix::t(Matrix::t(CM)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(CM)))
  }
  else {
    cds_exprs <- reshape2::melt(exprs(cds_subset))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_pData <- pData(cds_subset)

  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  if (integer_expression) {
    cds_exprs$adjusted_expression <- round(cds_exprs$expression)
  }
  else {
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  }
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$feature_label <- as.factor(cds_exprs$feature_label)
  # trend_formula <- paste("adjusted_expression", trend_formula,
  #     sep = "")
  cds_exprs$Branch <- as.factor(cds_exprs$Branch)

  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Branch = pData(cds_subset)$Branch)

  full_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula,
                                            relative_expr = T, new_data = new_data)
  colnames(full_model_expectation) <- colnames(cds_subset)

  cds_exprs$full_model_expectation <- apply(cds_exprs,1, function(x) full_model_expectation[x[2], x[1]])
  if(!is.null(reducedModelFormulaStr)){
    reduced_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = reducedModelFormulaStr,
                                                 relative_expr = T, new_data = new_data)
    colnames(reduced_model_expectation) <- colnames(cds_subset)
    cds_exprs$reduced_model_expectation <- apply(cds_exprs,1, function(x) reduced_model_expectation[x[2], x[1]])
  }

  # FIXME: If you want to show the bifurcation time for each gene, this function
  # should just compute it. Passing it in as a dataframe is just too complicated
  # and will be hard on the user.
  # if(!is.null(bifurcation_time)){
  #     cds_exprs$bifurcation_time <- bifurcation_time[as.vector(cds_exprs$gene_short_name)]
  # }
  if (method == "loess")
    cds_exprs$expression <- cds_exprs$expression + cds@lowerDetectionLimit
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  cds_exprs$expression[is.na(cds_exprs$expression)] <- min_expr
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$full_model_expectation[is.na(cds_exprs$full_model_expectation)] <- min_expr
  cds_exprs$full_model_expectation[cds_exprs$full_model_expectation < min_expr] <- min_expr

  if(!is.null(reducedModelFormulaStr)){
    cds_exprs$reduced_model_expectation[is.na(cds_exprs$reduced_model_expectation)] <- min_expr
    cds_exprs$reduced_model_expectation[cds_exprs$reduced_model_expectation < min_expr] <- min_expr
  }

  cds_exprs$State <- as.factor(cds_exprs$State)
  cds_exprs$Branch <- as.factor(cds_exprs$Branch)

  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs) + scale_color_manual(values = point_colors)
  # if (!is.null(bifurcation_time)) {
  #   q <- q + geom_vline(aes(xintercept = bifurcation_time),
  #                       color = "black", linetype = "longdash")
  # }
  if (is.null(color_by) == FALSE) {
    q <- q + geom_point(aes_string(color = color_by), size = I(cell_size))
  }
  if (is.null(reducedModelFormulaStr) == FALSE)
    q <- q + scale_y_log10() + facet_wrap(~feature_label +
                                            pval, nrow = nrow, ncol = ncol, scales = "free_y")
  else q <- q + scale_y_log10() + facet_wrap(~feature_label,
                                             nrow = nrow, ncol = ncol, scales = "free_y")
  if (method == "loess")
    q <- q + stat_smooth(aes(fill = Branch, color = Branch),
                         method = "loess")
  else if (method == "fitting") {
    q <- q + geom_line(aes_string(x = "Pseudotime", y = "full_model_expectation",
                                  linetype = "Branch"), data = cds_exprs) #+ scale_color_manual(name = "Type", values = c(colour_cell, colour), labels = c("Pre-branch", "AT1", "AT2", "AT1", "AT2")
  }

  if(!is.null(reducedModelFormulaStr)) {
    q <- q + geom_line(aes_string(x = "Pseudotime", y = "reduced_model_expectation"),
                       color = 'black', linetype = 2, data =  cds_exprs)
  }

  q <- q + ylab("Expression") + xlab("Pseudotime (stretched)")

  q <- q + monocle_theme_opts()
  q + expand_limits(y = min_expr)
}


# my_plot_cell_clusters ------------------------------------------------------

#' Plots clusters of cells .
#'
#' @param cds CellDataSet for the experiment
#' @param x the column of reducedDimS(cds) to plot on the horizontal axis
#' @param y the column of reducedDimS(cds) to plot on the vertical axis
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to map to each cell's color
#' @param markers a gene name or gene id to use for setting the size of each cell in the plot
#' @param show_cell_names draw the name of each cell in the plot
#' @param cell_size The size of the point for each cell
#' @param cell_name_size the size of cell name labels
#' @param ... additional arguments passed into the scale_color_viridis function
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom viridis scale_color_viridis
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' HSMM <- reduceD
#' plot_cell_clusters(HSMM)
#' plot_cell_clusters(HSMM, color_by="Pseudotime")
#' plot_cell_clusters(HSMM, markers="MYH3")
#' }
my_plot_cell_clusters <- function(cds,
                               x=1,
                               y=2,
                               color_by="Cluster",
                               markers=NULL,
                               show_cell_names=FALSE,
                               cell_size=1.5,
                               cell_name_size=2,
                               point_colors = c("#cb4bbe", "lightskyblue","grey37", "#53ba48","moccasin","#dcca44","#502d70","#afd671","#cb4471","#69da99","#d44732","#6fd3cf","#5d2539","#cfd0a0","blue","#d2883b","#6793c0","#898938","#c98ac2","yellow","#c4c0cc","#7d3d23","#00a5ff","#d68d7a","#a2c1a3"),
                               ...){
  if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) == 0){
    stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
  }

  gene_short_name <- NULL
  sample_name <- NULL
  data_dim_1 <- NULL
  data_dim_2 <- NULL

  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info <- pData(cds)

  tSNE_dim_coords <- reducedDimA(cds)
  data_df <- data.frame(t(tSNE_dim_coords[c(x,y),]))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- colnames(cds)
  data_df <- merge(data_df, lib_info, by.x="sample_name", by.y="row.names")

  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
      colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")

    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + facet_wrap(~feature_label)
  }else{
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
  }

  # FIXME: setting size here overrides the marker expression funtionality.
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    g <- g + geom_point(aes(color=log10(value + 0.1)), size=I(cell_size), na.rm = TRUE) +
      scale_color_manual(values = point_colors) +
      scale_color_viridis(name = paste0("log10(value + 0.1)"), ...)
  }else {
    g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE)+
      scale_color_manual(values = point_colors)
  }

  g <- g +
    monocle_theme_opts() +
    xlab(paste("Component", x)) +
    ylab(paste("Component", y)) +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(text = element_text(size = 15))
  g
}

# my_plot_rho_delta ----------------------------------------------------------

#' Plots the decision map of density clusters .
#'
#' @param cds CellDataSet for the experiment after running clusterCells_Density_Peak
#' @param rho_threshold The threshold of local density (rho) used to select the density peaks for plotting
#' @param delta_threshold The threshold of local distance (delta) used to select the density peaks for plotting
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' plot_rho_delta(HSMM)
#' }

my_plot_rho_delta <- function(cds, rho_threshold = NULL, delta_threshold = NULL){
    if(!is.null(cds@auxClusteringData[["tSNE"]]$densityPeak)
    & !is.null(pData(cds)$Cluster)
    & !is.null(pData(cds)$peaks)
    & !is.null(pData(cds)$halo)
    & !is.null(pData(cds)$delta)
    & !is.null(pData(cds)$rho)) {
    rho <- NULL
    delta <- NULL

    # df <- data.frame(rho = as.numeric(levels(pData(cds)$rho))[pData(cds)$rho],
    #   delta = as.numeric(levels(pData(cds)$delta))[pData(cds)$delta])
    if(!is.null(rho_threshold) & !is.null(delta_threshold)){
      peaks <- pData(cds)$rho >= rho_threshold & pData(cds)$delta >= delta_threshold
    }
    else
      peaks <- pData(cds)$peaks

    df <- data.frame(rho = pData(cds)$rho, delta = pData(cds)$delta, peaks = peaks)

    g <- qplot(rho, delta, data = df, alpha = I(0.5), color = peaks) +  monocle_theme_opts() +
      theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
      scale_color_manual(values=c("grey","black")) +
      theme(legend.key = element_blank()) +
      theme(panel.background = element_rect(fill='white'))
  }
  else {
    stop('Please run clusterCells_Density_Peak before using this plotting function')
  }
  g
}


# my_plot_complex_cell_trajectory --------------------------------------------

#' Plots the minimum spanning tree on cells.
#' @description Plots the minimum spanning tree on cells.
#' @param cds CellDataSet for the experiment
#' @param x the column of reducedDimS(cds) to plot on the horizontal axis
#' @param y the column of reducedDimS(cds) to plot on the vertical axis
#' @param root_states the state used to set as the root of the graph
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to map to each cell's color
#' @param show_tree whether to show the links between cells connected in the minimum spanning tree
#' @param show_backbone whether to show the diameter path of the MST used to order the cells
#' @param backbone_color the color used to render the backbone.
#' @param markers a gene name or gene id to use for setting the size of each cell in the plot
#' @param show_cell_names draw the name of each cell in the plot
#' @param cell_size The size of the point for each cell
#' @param cell_link_size The size of the line segments connecting cells (when used with ICA) or the principal graph (when used with DDRTree)
#' @param cell_name_size the size of cell name labels
#' @param show_branch_points Whether to show icons for each branch point (only available when reduceDimension was called with DDRTree)
#' @param ... Additional arguments passed to the scale_color_viridis function
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom igraph V get.edgelist layout_as_tree
#' @importFrom reshape2 melt
#' @importFrom viridis scale_color_viridis
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' plot_complex_cell_trajectory(HSMM)
#' plot_complex_cell_trajectory(HSMM, color_by="Pseudotime", show_backbone=FALSE)
#' plot_complex_cell_trajectory(HSMM, markers="MYH3")
#' }
my_plot_complex_cell_trajectory <- function(cds,
                                         x=1,
                                         y=2,
                                         root_states = NULL,
                                         color_by="State",
                                         show_tree=TRUE,
                                         show_backbone=TRUE,
                                         backbone_color="black",
                                         markers=NULL,
                                         show_cell_names=FALSE,
                                         cell_size=1.5,
                                         cell_link_size=0.75,
                                         cell_name_size=2,
                                         show_branch_points=TRUE,
                                         point_colors = c("#cb4bbe", "lightskyblue","grey37", "#53ba48","moccasin","#dcca44","#502d70","#afd671","#cb4471","#69da99","#d44732","#6fd3cf","#5d2539","#cfd0a0","blue","#d2883b","#6793c0","#898938","#c98ac2","yellow","#c4c0cc","#7d3d23","#00a5ff","#d68d7a","#a2c1a3"),
                                         ...){
  gene_short_name <- NA
  sample_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA

  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)

  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }

  if (cds@dim_reduce_type == "ICA"){
    reduced_dim_coords <- reducedDimS(cds)
  }else if (cds@dim_reduce_type %in% c("SimplePPT", "DDRTree", "SGL-tree") ){
    reduced_dim_coords <- reducedDimK(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }else {
    stop("Error: unrecognized dimensionality reduction method.")
  }

  if (is.null(reduced_dim_coords)){
    stop("You must first call reduceDimension() before using this function")
  }

  dp_mst <- minSpanningTree(cds)


  if(is.null(root_states)) {
    if(is.null(lib_info_with_pseudo$Pseudotime)){
      root_cell <- row.names(lib_info_with_pseudo)[degree(dp_mst) == 1][1]
    }
    else
      root_cell <- row.names(subset(lib_info_with_pseudo, Pseudotime == 0))

    if(cds@dim_reduce_type != "ICA")
      root_cell <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_cell, ]]

  }
  else {
    candidate_root_cells <- row.names(subset(pData(cds), State %in% root_states))
    if(cds@dim_reduce_type == "ICA") {
      root_cell <- candidate_root_cells[which(degree(dp_mst, candidate_root_cells) == 1)]
    } else {
      Y_candidate_root_cells <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[candidate_root_cells, ]]
      root_cell <- Y_candidate_root_cells[which(degree(dp_mst, Y_candidate_root_cells) == 1)]
    }

  }

  # #root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  # root_state <- pData(cds)[root_cell,]$State
  # #root_state <- V(pr_graph_cell_proj_mst)[root_cell,]$State

  # pr_graph_root <- subset(pData(cds), State == root_state)

  # closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  # root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root),]
  tree_coords <- layout_as_tree(dp_mst, root=root_cell)

  #ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(x,y),]))
  ica_space_df <- data.frame(tree_coords)
  row.names(ica_space_df) <- colnames(reduced_dim_coords)
  colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")

  ica_space_df$sample_name <- row.names(ica_space_df)
  #ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  #print(ica_space_with_state_df)


  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }

  edge_list <- as.data.frame(get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")

  edge_df <- merge(ica_space_df, edge_list, by.x="sample_name", by.y="source", all=TRUE)
  #edge_df <- ica_space_df
  edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="source_prin_graph_dim_1", "prin_graph_dim_2"="source_prin_graph_dim_2"))
  edge_df <- merge(edge_df, ica_space_df[,c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2")], by.x="target", by.y="sample_name", all=TRUE)
  edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="target_prin_graph_dim_1", "prin_graph_dim_2"="target_prin_graph_dim_2"))

  #S_matrix <- reducedDimS(cds)
  #data_df <- data.frame(t(S_matrix[c(x,y),]))

  if(cds@dim_reduce_type == "ICA"){
    S_matrix <- tree_coords[,] #colnames(cds)

  } else if(cds@dim_reduce_type %in% c("DDRTree", "SimplePPT", "SGL-tree")){
    S_matrix <- tree_coords[closest_vertex,]
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }

  data_df <- data.frame(S_matrix)
  row.names(data_df) <- colnames(reducedDimS(cds))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")

  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
      colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
    #print (head(edge_df))
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, I(cell_size))) + facet_wrap(~feature_label)
  }else{
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
  }
  if (show_tree){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=cell_link_size, linetype="solid", na.rm=TRUE, data=edge_df)
  }

  # FIXME: setting size here overrides the marker expression funtionality.
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    if(class(data_df[, color_by]) == 'numeric') {
      g <- g + geom_jitter(aes_string(color = paste0("log10(", color_by, " + 0.1)")), size=I(cell_size), na.rm = TRUE, height=5) +
                             scale_color_viridis(name = paste0("log10(", color_by, ")"), ...)
    } else {
      g <- g + geom_jitter(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE, height=5)
    }
  }else {
    if(class(data_df[, color_by]) == 'numeric') {
      g <- g + geom_jitter(aes_string(color = paste0("log10(", color_by, " + 0.1)")), size=I(cell_size), na.rm = TRUE, height=5) +
        scale_color_viridis(name = paste0("log10(", color_by, " + 0.1)"), ...)
    } else {
      g <- g + geom_jitter(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE, height=5)
    }
  }

  if (show_branch_points && cds@dim_reduce_type == 'DDRTree'){
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[,c("sample_name", "source_prin_graph_dim_1", "source_prin_graph_dim_2")]
    branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, mst_branch_nodes)
    branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx), ]

    g <- g + geom_point(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2"),
                        size=2 * cell_size, na.rm=TRUE, data=branch_point_df) +
      scale_color_manual(values = point_colors) +
      geom_text(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", label="branch_point_idx"),
                size=1.5 * cell_size, color="white", na.rm=TRUE, data=branch_point_df)
  }
  if (show_cell_names){
    g <- g +geom_text(aes(label=sample_name), size=cell_name_size)
  }
  g <- g +
    #scale_color_brewer(palette="Set1") +
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    # theme(axis.line.x = element_line(size=0.25, color="black")) +
    # theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank()) +
    xlab('') +
    ylab('') +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
  g
}














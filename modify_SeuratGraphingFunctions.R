# Load libraries ----------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("Seurat", "dplyr", "colorRamps", "ggplot2", "tidyr")

ipak(packages)

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
                 "#a2c1a3",
                  primary.colors(8))


# Required background functions from Seuratv2.3.1 -------------------------

PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

# Modify DotPlot (Seuratv2.3.1) -------------------------------------------

#Modify dotplot to increase minimum size of pct.exp points

my_dotplot<-function (object, genes.plot, cols.use = c("lightgrey", "blue"), 
  col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
  scale.by = "radius", scale.min = NA, scale.max = NA, group.by, 
  plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE) 
{
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
    radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (!missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  colnames(x = data.to.plot) <- genes.plot
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot <- data.to.plot %>% gather(key = genes.plot, 
    value = expression, -c(cell, id))
  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
    summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression, 
      threshold = 0))
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
    mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale, 
    max = col.max, min = col.min))
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
    levels = rev(x = genes.plot))
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
    y = id)) + geom_point(mapping = aes(size = pct.exp, 
    color = avg.exp.scale)) + scale.func(range = c(2, dot.scale+2), 
    limits = c(0, scale.max)) + theme(axis.title.x = element_blank(), 
    axis.title.y = element_blank())
  if (length(x = cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  }
  else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, 
      vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}

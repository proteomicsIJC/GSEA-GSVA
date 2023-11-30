################ 
## treemaping ##
###############

## soup = minestrone result
## graph_title = graph title

treemaping <- function(soup, graph_title){
  treemap::treemap(soup, index = c("parent", "child"), 
                   vSize = "size", type = "index", title = graph_title,
                   ## palette = gg_color_hue(length(unique(soup1$parent))),
                   palette = hue_pal()(length(unique(soup$parent))),
                   fontcolor.labels = c("#FFFFFFDD", "#00000080"), bg.labels = 0, 
                   border.col = "#00000080", fontsize.labels = 10)}

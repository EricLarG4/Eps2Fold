
#PCA plotter----
pca.plotR <- function(pca.data, dim.1 = input$dim.cd[1], dim.2 = input$dim.cd[1], color = 'GBA', shape = 'GBA'){
  
  dim.1.cap <- gsub('Dim.', 'Dimension ', dim.1)
  dim.2.cap <- gsub('Dim.', 'Dimension ', dim.2)
  
  pca.data %>% 
    mutate(
      km = factor(km),
      tetrad = factor(tetrad),
      tetrad.loop = if_else(
        loop != 'Z DNA',
        paste0(tetrad, loop),
        loop
      )
    ) %>% 
    rename(
      'k-means' = km,
      'Topology' = topo, 
      'GBA' = gba, 
      'Tetrads' = tetrad, 
      'Tetrad combination' = tetrad.id, 
      'Loop progression' = loop, 
      'Tetrads x Loops' = tetrad.loop,
      'Tetrad handedness' = plus.minus, 
      'Grooves' = groove, 
      'Cation' = salt
    ) %>% 
    ggplot(
      aes(
        x = .data[[dim.1]], y = .data[[dim.2]], 
        color = .data[[color]], fill = .data[[color]],
        shape = .data[[shape]]
      )
    ) +
    scale_shape_manual(values = c(0:15)) +
    geom_polygon(
      data = . %>%
        slice(chull(.data[[dim.1]], .data[[dim.2]]), .by = color),
      alpha = 0.25,
      show.legend = FALSE
    ) +
    geom_point(
      size = 6,
      stroke = 1.25,
      show.legend = TRUE
    ) +
    geom_text_repel(
      aes(label = oligo),
      force = 20,
      size = 8,
      point.padding = 1,
      min.segment.length = 0.75,
      show.legend = FALSE
    ) +
    theme(
      panel.background = element_blank(),
      strip.background = element_blank(),
      legend.position = 'bottom',
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 22),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 22),
      strip.text = element_text(size = 20),
      axis.line = element_line(size = 0.75),
      axis.ticks = element_line(size = 0.75)
    )  +
    labs(
      y = dim.2.cap,
      x = dim.1.cap
    ) 
}

#correlation circle----
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(cir.x = xx, cir.y = yy))
}
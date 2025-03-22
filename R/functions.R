## Dependencies----
library(tidyverse)
library(ggtext)

## Theme----
custom.theme.markdown <- ggplot2::theme(
  panel.background = element_blank(),
  strip.background = element_blank(),
  legend.position = 'bottom',
  legend.background = element_blank(),
  legend.box.background = element_blank(),
  legend.key = element_blank(),
  legend.text = element_markdown(
    size = 18,
  ),
  legend.title = element_markdown(
    size = 22
  ),
  axis.text = element_markdown(
    size = 18
  ),
  axis.title.x = element_markdown(
    size = 22
  ),
  axis.title.y = element_markdown(size = 22, angle = 90),
  strip.text = element_markdown(
    size = 20
  ),
  axis.line = element_line(
    size = 0.75
  ),
  axis.ticks = element_line(
    size = 0.75
  )
)

#PCA plotter----
pca.plotR <- function(
  pca.data,
  dim.1 = input$dim.cd[1],
  dim.2 = input$dim.cd[1],
  color = 'GBA',
  shape = 'GBA'
) {
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
      'Conformer' = conformer,
      'GBA' = gba,
      'GBA stacks' = gba.stacks,
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
        x = .data[[dim.1]],
        y = .data[[dim.2]],
        color = .data[[color]],
        fill = .data[[color]],
        shape = .data[[shape]]
      )
    ) +
    scale_shape_manual(values = 0:15) +
    geom_polygon(
      data = . %>%
        filter(.data[[color]] != "User") %>%
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
      legend.text = element_markdown(size = 18),
      legend.title = element_markdown(size = 22),
      axis.text = element_markdown(size = 18),
      axis.title = element_markdown(size = 22),
      strip.text = element_markdown(size = 20),
      axis.line = element_line(size = 0.75),
      axis.ticks = element_line(size = 0.75)
    ) +
    labs(
      y = dim.2.cap,
      x = dim.1.cap
    )
}

#correlation circle----
circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100) {
  r = diameter / 2
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  data.frame(cir.x = xx, cir.y = yy)
}

# CD and IDS plots----

## Plotting function----
render_ref_plot <- function(df, input, norm_col, axis_label_text) {
  # Initialize color and sd columns
  df <- df %>%
    add_column(color = NA, sd.y = NA)

  # Assign y values based on normalization method
  df$y <- df[[norm_col]]

  # Define color mapping
  color_map <- list(
    "Topology" = df$topo,
    "Conformer" = df$conformer,
    "GBA" = df$gba,
    "GBA stacks" = df$gba.stacks,
    "Tetrads" = factor(df$tetrad),
    "Tetrad combination" = df$tetrad.id,
    "Loop progression" = df$loop,
    "Tetrad handedness" = df$plus.minus,
    "Grooves" = df$groove,
    "Salt" = df$salt
  )
  df$color <- color_map[[input$ref.color]]

  #Initialize plot
  p <- df %>%
    ggplot(aes(x = wl, y = y, color = color)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = c(295, 275, 262, 245), linetype = 'dashed') +
    labs(
      x = "&lambda; (nm)",
      y = axis_label_text,
      color = input$ref.color,
      fill = input$ref.color
    ) +
    custom.theme.markdown +
    scale_y_continuous(n.breaks = 3) +
    scale_x_continuous(expand = c(0, 0))

  # If panel mode is "Mean", add geom_ribbon and geom_line of mean
  if (input$ref.panel == "Mean") {
    # Calculate mean
    df_mean <- df %>%
      group_by(color, wl) %>%
      summarise(
        sd.y = sd(y, na.rm = TRUE),
        y = mean(y, na.rm = TRUE),
        .groups = "drop"
      )

    # Add geom_ribbon and geom_line
    p <- p +
      geom_ribbon(
        data = df_mean,
        inherit.aes = FALSE,
        aes(x = wl, y = y, ymin = y - sd.y, ymax = y + sd.y, fill = color),
        alpha = 0.5,
        show.legend = FALSE
      ) +
      geom_line(
        data = df_mean,
        size = 1,
        show.legend = TRUE
      )
  } else {
    #else add lines of individuals
    p <- p +
      geom_line(aes(group = oligo), size = 1, show.legend = TRUE)
  }

  # If panel mode is "Panels", add facet_wrap
  if (input$ref.panel == "Panels") {
    p <- p + facet_wrap(~oligo)
  }

  p
}

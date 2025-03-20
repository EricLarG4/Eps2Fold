shiny::bootstrapLib(bslib::bs_theme())

# operators----
`%notin%` <- Negate(`%in%`)

#server----
server <- shinyServer(function(input, output, session) {
  # Theme handling
  observe({
    # Update the theme when the selector changes
    session$setCurrentTheme(
      bslib::bs_theme(
        preset = input$theme_select,
        font_scale = input$font_size
      )
    )
  })

  # Extraction of primary theme color for discrete ggplot theming
  current_primary <- reactive({
    theme <- session$getCurrentTheme()
    bslib::bs_get_variables(theme, "primary")
  })

  # Enable thematic
  thematic::thematic_shiny(font = "auto")
  thematic::thematic_rmd(font = "auto")

  # theme----
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

  # Change ggplot2's default "gray" theme
  theme_set(custom.theme.markdown)

  # Reference input----

  ## Data files----
  # input.file <- reactive({
  #   input$ref.data
  # })

  input.file <- reactive({
    if (isFALSE(input$data_source) && !is.null(input$ref.data)) {
      return(input$ref.data$datapath)
    } else if (isTRUE(input$data_source)) {
      "data/Reference_data.xlsx"
    } else {
      return(NULL)
    }
  })

  ### File import toggle----
  file.toggle <- reactive({
    if (!is.null(input.file())) {
      return('yes')
    } else {
      return('no')
    }
  })

  # outputs the value of the file toggle
  # output$file_toggle_value <- renderPrint({
  #   file.toggle()
  # })

  ### Reference set----
  ref.set <- reactive({
    if (isFALSE(input$data_source)) {
      ref.set <- read_excel(input.file())
    } else {
      ref.set <- read_excel(input.file())
    }

    if (input$ids.ref.select == TRUE) {
      ref.set <- ref.set %>%
        mutate(
          delta.eps.ids = delta.eps.ids.th,
          uv.eps = if_else(
            salt == 'none',
            th.eps,
            uv.eps
          )
        )
    }

    ref.set %>%
      select(-c(delta.eps.ids.th, th.eps)) %>%
      filter(
        wl <= max(input$wl),
        wl >= min(input$wl),
        wl %% 1 == 0
      ) %>%
      group_by(oligo, salt) %>%
      mutate(
        plus.minus = gsub("\"", "", plus.minus),
        salt = gsub('<sup>', '', salt),
        salt = gsub('</sup>', '', salt),
        norm.ids = case_when(
          input$ids.norm == "-1/+1" ~
            2 *
              (delta.eps.ids - min(delta.eps.ids)) /
              (max(delta.eps.ids) - min(delta.eps.ids)) -
              1,
          TRUE ~ delta.eps.ids
        ),
        norm.cd = case_when(
          input$cd.norm == "Δε/ε" ~ delta.eps.cd / eps,
          input$cd.norm == "-1/+1" ~
            2 *
              (delta.eps.cd - min(delta.eps.cd, na.rm = TRUE)) /
              (max(delta.eps.cd, na.rm = TRUE) -
                min(delta.eps.cd, na.rm = TRUE)) -
              1,
          TRUE ~ delta.eps.cd
        )
      )
  })

  #### Sequences----
  ref.seq <- reactive({
    ref.set() %>%
      select(
        c(
          "oligo.number",
          "oligo",
          "description",
          "sequence",
          "nt",
          "eps",
          "salt",
          "topo",
          "conformer",
          "feature",
          "gba",
          "gba.stacks",
          "tetrad",
          "tetrad.id",
          "loop",
          "plus.minus",
          "groove",
          "conc"
        )
      ) %>%
      filter(salt != 'none') %>%
      unique()
  })

  #### Picker inputs: filtering source----
  output$ref.seq.topo.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.topo.0",
        label = "Topology",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$topo))

      selectizeInput(
        "ref.seq.topo.0",
        label = "Topology",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.confomer.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.conformer.0",
        label = "Conformer",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$conformer))

      selectizeInput(
        "ref.seq.conformer.0",
        label = "Conformer",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.gba.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.gba.0",
        label = "GBA",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$gba))

      selectizeInput(
        "ref.seq.gba.0",
        label = "GBA",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.gba.stacks.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.gba.stacks.0",
        label = "GBA stacks",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$gba.stacks))

      selectizeInput(
        "ref.seq.gba.stacks.0",
        label = "GBA stacks",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.tetrad.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.tetrad.0",
        label = "Tetrads",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$tetrad))

      selectizeInput(
        "ref.seq.tetrad.0",
        label = "Tetrads",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.tetrad.id.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.tetrad.id.0",
        label = "Tetrad combination",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$tetrad.id))

      selectizeInput(
        "ref.seq.tetrad.id.0",
        label = "Tetrad combination",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.loop.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.loop.0",
        label = "Loop progression",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$loop))

      selectizeInput(
        "ref.seq.loop.0",
        label = "Loop progression",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.plus.minus.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.plus.minus.0",
        label = "Tetrad handedness",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$plus.minus))

      selectizeInput(
        "ref.seq.plus.minus.0",
        label = "Tetrad handedness",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.groove.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.groove.0",
        label = "Grooves",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$groove))

      selectizeInput(
        "ref.seq.groove.0",
        label = "Grooves",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.cation.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.cation.0",
        label = "Cation",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$salt))

      selectizeInput(
        "ref.seq.cation.0",
        label = "Cation",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.oligo.0_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.oligo.0",
        label = "Oligonucleotide",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(ref.seq()$oligo))

      selectizeInput(
        "ref.seq.oligo.0",
        label = "Oligonucleotide",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  #### Picker inputs: filtering plots----
  output$ref.seq.topo_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.topo",
        label = "Topology",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- unique(input$ref.seq.topo.0)

      selectizeInput(
        "ref.seq.topo",
        label = "Topology",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.conformer_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.conformer",
        label = "Conformer",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- unique(input$ref.seq.conformer.0)

      selectizeInput(
        "ref.seq.conformer",
        label = "Conformer",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.gba_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.gba",
        label = "GBA",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- unique(input$ref.seq.gba.0)

      selectizeInput(
        "ref.seq.gba",
        label = "GBA",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.gba.stacks_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.gba.stacks",
        label = "GBA stacks",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- unique(input$ref.seq.gba.stacks.0)

      selectizeInput(
        "ref.seq.gba.stacks",
        label = "GBA stacks",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.tetrad_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.tetrad",
        label = "Tetrads",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- unique(input$ref.seq.tetrad.0)

      selectizeInput(
        "ref.seq.tetrad",
        label = "Tetrads",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.tetrad.id_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.tetrad.id",
        label = "Tetrad combination",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- unique(input$ref.seq.tetrad.id.0)

      selectizeInput(
        "ref.seq.tetrad.id",
        label = "Tetrad combination",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.loop_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.loop",
        label = "Loop progression",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- unique(input$ref.seq.loop.0)

      selectizeInput(
        "ref.seq.loop",
        label = "Loop progression",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.plus.minus_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.plus.minus",
        label = "Tetrad handedness",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- unique(input$ref.seq.plus.minus.0)

      selectizeInput(
        "ref.seq.plus.minus",
        label = "Tetrad handedness",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.groove_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.groove",
        label = "Grooves",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- unique(input$ref.seq.groove.0)

      selectizeInput(
        "ref.seq.groove",
        label = "Grooves",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  output$ref.seq.cation_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.cation",
        label = "Cation",
        choices = "upload data first",
        multiple = T
      )
    } else {
      selectizeInput(
        "ref.seq.cation",
        label = "Cation",
        choices = unique(input$ref.seq.cation.0),
        selected = unique(input$ref.seq.cation.0),
        multiple = T
      )
    }
  })

  output$ref.seq.oligo_ui <- renderUI({
    if (file.toggle() == 'no') {
      selectizeInput(
        "ref.seq.oligo",
        label = "Oligonucleotide",
        choices = "upload data first",
        multiple = T
      )
    } else {
      selectizeInput(
        "ref.seq.oligo",
        label = "Oligonucleotide",
        choices = unique(input$ref.seq.oligo.0),
        selected = unique(input$ref.seq.oligo.0),
        multiple = T
      )
    }
  })

  ## Outputs of inputs----

  ### Reference sequences----
  output$seq <- renderDT(server = FALSE, {
    if (file.toggle() == 'no') {
      return(NULL)
    } else {
      ref.seq() %>%
        filter(
          salt %in% input$ref.seq.cation.0,
          topo %in% input$ref.seq.topo.0,
          conformer %in% input$ref.seq.conformer.0,
          gba %in% input$ref.seq.gba.0,
          gba.stacks %in% input$ref.seq.gba.stacks.0,
          oligo %in% input$ref.seq.oligo.0,
          tetrad %in% input$ref.seq.tetrad.0,
          tetrad.id %in% input$ref.seq.tetrad.id.0,
          loop %in% input$ref.seq.loop.0,
          plus.minus %in% input$ref.seq.plus.minus.0,
          groove %in% input$ref.seq.groove.0
        ) %>%
        mutate(conc = round(conc, 2)) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          colnames = c(
            "#" = "oligo.number",
            "Name" = "oligo",
            "&epsilon;<sub>260nm</sub> (M<sup>-1</sup>cm<sup>-1</sup>)" = "eps",
            "Cation" = "salt",
            "C (µM)" = "conc",
            "Sequence" = "sequence",
            "Number of nucleotides" = "nt",
            "Topology" = "topo",
            "Conformer" = "conformer",
            "Other feature" = "feature",
            "Description" = "description",
            "GBA" = "gba",
            "GBA stacks" = "gba.stacks",
            "Tetrads" = "tetrad",
            "Tetrad combination" = "tetrad.id",
            "Loop progression" = "loop",
            "Tetrad handedness" = "plus.minus",
            "Grooves" = "groove"
          ),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 750,
            scroller = F,
            pageLength = 50,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "reference oligos"),
              list(
                extend = 'excel',
                title = NULL,
                filename = "reference oligos"
              ),
              list(extend = 'colvis')
            ),
            title = NULL,
            columnDefs = list(list(visible = FALSE, targets = c(0, 5, 6, 17)))
          )
        )
    }
  })

  ### Reference UV----

  output$p.ref.uv <- renderPlot({
    if (file.toggle() == 'no') {
      ggplot()
    } else {
      ref.set() %>%
        filter(
          topo %in% input$ref.seq.topo,
          conformer %in% input$ref.seq.conformer,
          gba %in% input$ref.seq.gba,
          gba.stacks %in% input$ref.seq.gba.stacks,
          salt %in% c(input$ref.seq.cation, "none"),
          oligo %in% input$ref.seq.oligo,
          tetrad %in% input$ref.seq.tetrad,
          tetrad.id %in% input$ref.seq.tetrad.id,
          loop %in% input$ref.seq.loop,
          plus.minus %in% input$ref.seq.plus.minus,
          groove %in% input$ref.seq.groove
        ) %>%
        ggplot(
          aes(
            x = wl,
            y = uv.eps,
            color = salt,
            group = salt
          )
        ) +
        geom_line(
          size = 1,
          show.legend = T
        ) +
        facet_wrap(
          ~oligo
          # ncol = 4
        ) +
        labs(
          x = "&lambda; (nm)",
          y = "&epsilon; (M<sup>-1</sup>cm<sup>-1</sup>)",
          color = "Cation"
        ) +
        custom.theme.markdown +
        scale_y_continuous(n.breaks = 3) +
        scale_x_continuous(expand = c(0, 0))
    }
  })

  ### Reference UV/CD/IDS tables----
  output$ref.uv <- renderDT(server = FALSE, {
    if (file.toggle() == 'no') {
      return(NULL)
    } else {
      ref.set() %>%
        mutate(
          conc = round(conc, 2),
          abs = round(abs, 3)
        ) %>%
        select(c(
          "oligo.number",
          "oligo",
          "description",
          "sequence",
          "nt",
          "eps",
          "salt",
          "topo",
          "conformer",
          "feature",
          "gba",
          "gba.stacks",
          "tetrad",
          "tetrad.id",
          "loop",
          "plus.minus",
          "groove",
          "conc",
          "wl",
          "abs",
          "uv.eps"
        )) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          colnames = c(
            "Wavelength (nm)" = "wl",
            "Name" = "oligo",
            "A" = "abs",
            "#" = "oligo.number",
            "&epsilon;<sub>260nm</sub> (M<sup>-1</sup>cm<sup>-1</sup>)" = "eps",
            "Experimental &epsilon; (M<sup>-1</sup>cm<sup>-1</sup>)" = "uv.eps",
            "Cation" = "salt",
            "C (µM)" = "conc",
            "Sequence" = "sequence",
            "Number of nucleotides" = "nt",
            "Topology" = "topo",
            "Conformer" = "conformer",
            "Other feature" = "feature",
            "Description" = "description",
            "GBA" = "gba",
            "GBA stacks" = "gba.stacks",
            "Tetrads" = "tetrad",
            "Loop progression" = "loop",
            "Tetrad combination" = "tetrad.id",
            "Loop progression" = "loop",
            "Tetrad handedness" = "plus.minus",
            "Grooves" = "groove"
          ),
          rownames = F,
          escape = FALSE,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "reference uv"),
              list(extend = 'excel', title = NULL, filename = "reference uv"),
              list(extend = 'colvis')
            ),
            title = NULL,
            columnDefs = list(list(visible = FALSE, targets = c(2:6, 9, 17)))
          )
        )
    }
  })

  output$ref.ids <- renderDT(server = FALSE, {
    if (file.toggle() == 'no') {
      return(NULL)
    } else {
      ref.set() %>%
        ungroup() %>%
        filter(
          salt != 'none',
          !is.na(delta.eps.ids)
          # !is.na(delta.eps.ids.th)
        ) %>%
        mutate(
          conc = round(conc, 2),
          delta.eps.ids = round(delta.eps.ids, 0)
        ) %>%
        select(c(
          "oligo.number",
          "oligo",
          "description",
          "sequence",
          "nt",
          "eps",
          "salt",
          "topo",
          "conformer",
          "feature",
          "gba",
          "gba.stacks",
          "tetrad",
          "loop",
          "conc",
          "wl",
          "delta.eps.ids",
          "norm.ids"
        )) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          colnames = c(
            "#" = "oligo.number",
            "Wavelength (nm)" = "wl",
            "Name" = "oligo",
            "&epsilon;<sub>260nm</sub> (M<sup>-1</sup>cm<sup>-1</sup>)" = "eps",
            "&Delta;&epsilon;<sub>260nm</sub> (M<sup>-1</sup>cm<sup>-1</sup>)" = "delta.eps.ids",
            "Cation" = "salt",
            "C (µM)" = "conc",
            "Normalized IDS" = "norm.ids",
            "Sequence" = "sequence",
            "Number of nucleotides" = "nt",
            "Topology" = "topo",
            "Conformer" = "conformer",
            "Other feature" = "feature",
            "Description" = "description",
            "GBA" = "gba",
            "GBA stacks" = "gba.stacks",
            "Tetrads" = "tetrad",
            "Loop progression" = "loop"
          ),
          rownames = F,
          escape = FALSE,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "reference ids"),
              list(extend = 'excel', title = NULL, filename = "reference ids"),
              list(extend = 'colvis')
            ),
            title = NULL,
            columnDefs = list(list(
              visible = FALSE,
              targets = c(2:6, 9, 14, 17)
            ))
          )
        )
    }
  })

  output$ref.cd <- renderDT(server = FALSE, {
    if (file.toggle() == 'no') {
      return(NULL)
    } else {
      ref.set() %>%
        mutate(
          conc = round(conc, 2),
          delta.eps.cd = round(delta.eps.cd, 2),
          norm.cd = round(norm.cd, 2)
        ) %>%
        select(c(
          "oligo.number",
          "oligo",
          "description",
          "sequence",
          "nt",
          "eps",
          "salt",
          "topo",
          "conformer",
          "feature",
          "gba",
          "gba.stacks",
          "tetrad",
          "loop",
          "conc",
          "wl",
          "delta.eps.cd",
          "norm.cd"
        )) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          colnames = c(
            "Wavelength (nm)" = "wl",
            "Name" = "oligo",
            "Molar ellipticity (M<sup>-1</sup>cm<sup>-1</sup>)" = "delta.eps.cd",
            "#" = "oligo.number",
            "&epsilon;<sub>260nm</sub> (M<sup>-1</sup>cm<sup>-1</sup>)" = "eps",
            "C (µM)" = "conc",
            "Cation" = "salt",
            "Sequence" = "sequence",
            "Number of nucleotides" = "nt",
            "Topology" = "topo",
            "Conformer" = "conformer",
            "Other feature" = "feature",
            "Description" = "description",
            "GBA" = "gba",
            "GBA stacks" = "gba.stacks",
            "Tetrads" = "tetrad",
            "Loop progression" = "loop",
            "Normalized CD" = "norm.cd"
          ),
          rownames = F,
          escape = FALSE,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "reference cd"),
              list(extend = 'excel', title = NULL, filename = "reference cd"),
              list(extend = 'colvis')
            ),
            title = NULL,
            columnDefs = list(list(
              visible = FALSE,
              targets = c(2:6, 9, 14, 17)
            ))
          )
        )
    }
  })

  ### Reference CD/IDS plots----

  output$p.ref.ids <- renderPlot({
    ref.ids <- ref.set() %>%
      ungroup() %>%
      filter(salt != 'none') %>%
      filter(
        topo %in% input$ref.seq.topo,
        conformer %in% input$ref.seq.conformer,
        gba %in% input$ref.seq.gba,
        gba.stacks %in% input$ref.seq.gba.stacks,
        salt %in% c(input$ref.seq.cation, "none"),
        oligo %in% input$ref.seq.oligo,
        tetrad %in% input$ref.seq.tetrad,
        tetrad.id %in% input$ref.seq.tetrad.id,
        loop %in% input$ref.seq.loop,
        plus.minus %in% input$ref.seq.plus.minus,
        groove %in% input$ref.seq.groove
      )

    norm_col <- if (input$ids.norm == "Δε") "delta.eps.ids" else "norm.ids"
    axis_label_text <- if (input$ids.norm == "Δε")
      "&Delta;&epsilon; (M<sup>-1</sup>cm<sup>-1</sup>)" else "Normalized IDS"

    render_ref_plot(ref.ids, input, norm_col, axis_label_text)
  })

  output$p.ref.cd <- renderPlot({
    ref.cd <- ref.set() %>%
      filter(
        topo %in% input$ref.seq.topo,
        conformer %in% input$ref.seq.conformer,
        gba %in% input$ref.seq.gba,
        gba.stacks %in% input$ref.seq.gba.stacks,
        salt %in% input$ref.seq.cation,
        oligo %in% input$ref.seq.oligo,
        tetrad %in% input$ref.seq.tetrad,
        tetrad.id %in% input$ref.seq.tetrad.id,
        loop %in% input$ref.seq.loop,
        plus.minus %in% input$ref.seq.plus.minus,
        groove %in% input$ref.seq.groove
      )

    norm_col <- if (input$cd.norm == "Δε") "delta.eps.cd" else "norm.cd"
    axis_label_text <- if (input$cd.norm == "Δε")
      "&Delta;&epsilon; (M<sup>-1</sup>cm<sup>-1</sup>)" else "Normalized CD"

    render_ref_plot(ref.cd, input, norm_col, axis_label_text)
  })

  # PCA----

  ## training set----

  ### CD----

  training.cd <- reactive({
    ref.set() %>%
      ungroup() %>%
      filter(
        salt %in% input$ref.seq.cation.0,
        topo %in% input$ref.seq.topo.0,
        conformer %in% input$ref.seq.conformer.0,
        gba %in% input$ref.seq.gba.0,
        gba.stacks %in% input$ref.seq.gba.stacks.0,
        oligo %in% input$ref.seq.oligo.0,
        tetrad %in% input$ref.seq.tetrad.0,
        tetrad.id %in% input$ref.seq.tetrad.id.0,
        loop %in% input$ref.seq.loop.0,
        plus.minus %in% input$ref.seq.plus.minus.0,
        groove %in% input$ref.seq.groove.0
      ) %>%
      mutate(
        delta.eps = case_when(
          input$cd.norm == "Δε" ~ delta.eps.cd,
          input$cd.norm == "Δε/ε" ~ norm.cd,
          input$cd.norm == "-1/+1" ~ norm.cd
        )
      ) %>%
      select(oligo, wl, delta.eps) %>%
      mutate(wl = paste0("cd-", wl)) %>%
      pivot_wider(
        names_from = wl,
        values_from = delta.eps
      )
  })

  output$training.cd <- renderDT(server = FALSE, {
    if (file.toggle() == 'no') {
      return(NULL)
    } else {
      training.cd() %>%
        mutate_at(., vars(2:last_col()), round, 2) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "training cd"),
              list(extend = 'excel', title = NULL, filename = "training cd"),
              list(extend = 'colvis')
            ),
            title = NULL
            # columnDefs = list(list(visible=FALSE, targets=c(0,2:6,9,12)))
          )
        )
      # formatStyle(
      #   columns = 0:17,
      #   target = "cell",
      #   backgroundColor = "#272c30"
      # )
    }
  })

  ### IDS----

  training.ids <- reactive({
    ref.set() %>%
      ungroup() %>%
      filter(
        !is.na(delta.eps.ids),
        salt %in% input$ref.seq.cation.0,
        topo %in% input$ref.seq.topo.0,
        conformer %in% input$ref.seq.conformer.0,
        gba %in% input$ref.seq.gba.0,
        gba.stacks %in% input$ref.seq.gba.stacks.0,
        oligo %in% input$ref.seq.oligo.0,
        tetrad %in% input$ref.seq.tetrad.0,
        tetrad.id %in% input$ref.seq.tetrad.id.0,
        loop %in% input$ref.seq.loop.0,
        plus.minus %in% input$ref.seq.plus.minus.0,
        groove %in% input$ref.seq.groove.0
      ) %>%
      mutate(
        delta.eps = case_when(
          input$ids.norm == "Δε" ~ delta.eps.ids,
          input$ids.norm == "-1/+1" ~ norm.ids
        )
      ) %>%
      select(oligo, wl, delta.eps) %>%
      mutate(wl = paste0("ids-", wl)) %>%
      pivot_wider(
        names_from = wl,
        values_from = delta.eps
      )
  })

  output$training.ids <- renderDT(server = FALSE, {
    if (file.toggle() == 'no') {
      return(NULL)
    } else {
      training.ids() %>%
        # mutate_at(., vars(2:last_col()), round, 2) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "training ids"),
              list(extend = 'excel', title = NULL, filename = "training ids"),
              list(extend = 'colvis')
            ),
            title = NULL
            # columnDefs = list(list(visible=FALSE, targets=c(0,2:6,9,12)))
          )
        )
      # formatStyle(
      #   columns = 0:17,
      #   target = "cell",
      #   backgroundColor = "#272c30"
      # )
    }
  })

  training.cd.ids <- reactive({
    training.cd() %>%
      left_join(
        training.ids(),
        by = "oligo"
      )
  })

  ## Calculation----

  ### CD----

  #### PCA----
  pca.cd <- reactive({
    pca.cd <- FactoMineR::PCA(
      X = training.cd()[, -1],
      ncp = input$ncp,
      scale.unit = input$scale.unit,
      graph = FALSE
    )

    # save(pca.cd, file = "pca.cd.RData")

    return(pca.cd)
  })

  #### Clusters----

  #Optimal number of clusters

  pca.cd.twss <- reactive({
    data.frame(
      clusters = 1:10
    ) %>%
      group_by(clusters) %>%
      mutate(
        totwithin = kmeans(
          pca.cd()$ind$coord,
          algorithm = input$k.mean.algo,
          centers = clusters,
          nstart = 100
        )$tot.withinss
      ) %>%
      ggplot(
        aes(x = clusters, y = totwithin)
      ) +
      geom_point(size = 4) +
      geom_line(size = 1) +
      scale_x_continuous(breaks = 1:10) +
      labs(
        x = "clusters",
        y = "total within sum of square"
      ) +
      custom.theme.markdown
  })

  output$pca.cd.twss <- renderPlot({
    pca.cd.twss()
  })

  output$pca.cd.twss.2 <- renderPlot({
    pca.cd.twss()
  })

  pca.cd.gap <- reactive({
    clusGap(
      pca.cd()$ind$coord,
      FUN = kmeans,
      nstart = 25,
      K.max = 10,
      B = 50
    )$Tab %>%
      as_tibble() %>%
      rownames_to_column(var = "cluster") %>%
      mutate(cluster = as.numeric(cluster)) %>%
      ggplot(
        aes(x = cluster, y = gap)
      ) +
      geom_vline(
        data = . %>% filter(gap == max(gap)),
        aes(xintercept = cluster),
        linetype = 'dashed',
        size = 1,
        color = "tomato"
      ) +
      geom_errorbar(
        aes(x = cluster, ymin = gap - SE.sim, ymax = gap + SE.sim),
        width = 0
      ) +
      geom_point(size = 5) +
      geom_line(size = 1) +
      scale_x_continuous(breaks = 1:10) +
      labs(
        x = "clusters",
        y = "gap statistic"
      ) +
      custom.theme.markdown
  })

  output$pca.cd.gap <- renderPlot({
    pca.cd.gap()
  })

  output$pca.cd.gap.2 <- renderPlot({
    pca.cd.gap()
  })

  #k-means clustering
  pca.cd.km <- reactive({
    set.seed(42)

    kmeans(
      pca.cd()$ind$coord,
      algorithm = input$k.mean.algo,
      centers = input$cluster.center,
      nstart = 100
    )
  })

  #### Analytics ----

  ##### Scree----

  pca.cd.scree <- reactive({
    pca.cd.scree <- data.frame(pca.cd()$eig) %>%
      mutate(dim = 1:nrow(data.frame(pca.cd()$eig))) %>%
      # mutate(dim.cum.sum = cumsum(dim)) %>%
      filter(dim <= 10) %>%
      ggplot(aes(x = dim, y = percentage.of.variance)) +
      geom_bar(
        stat = 'identity',
        aes(fill = -dim),
        show.legend = FALSE
      ) +
      geom_point(
        aes(
          x = dim,
          y = cumsum(percentage.of.variance)
        ),
        size = 5
      ) +
      geom_line(
        aes(
          x = dim,
          y = cumsum(percentage.of.variance)
        ),
        size = 1
      ) +
      custom.theme.markdown +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
      scale_y_continuous(expand = c(0, NA), limits = c(0, 105)) +
      labs(
        y = 'Variance contribution (%)',
        x = 'Dimension'
      )

    #diagnostics
    # save(pca.cd.scree, file = "scree.rdata")

    return(pca.cd.scree)
  })

  output$pca.cd.scree <- renderPlot({
    pca.cd.scree()
  })

  output$pca.cd.scree.2 <- renderPlot({
    pca.cd.scree()
  })

  ##### Factor map----
  pca.cd.fac.map <- reactive({
    pca.cd.fac.map <- data.frame(pca.cd()$var$cos2) %>%
      select(input$dim.cd) %>%
      mutate(
        var = as.numeric(substr(rownames(data.frame(pca.cd()$var$cos2)), 4, 6))
      ) %>%
      pivot_longer(
        cols = 1:length(input$dim.cd),
        names_to = 'dim',
        values_to = 'cos2'
      ) %>%
      mutate(
        dim = substr(dim, 5, 5)
      ) %>%
      ggplot(
        aes(
          x = dim,
          y = var,
          size = cos2
        )
      ) +
      geom_point(
        color = "#2AA198",
        alpha = 0.5
      ) +
      labs(
        y = 'Variable',
        x = 'Dimension'
      ) +
      custom.theme.markdown +
      scale_size_continuous(
        limits = c(0, 1),
        range = c(0, 25),
        name = "cos<sup>2</sup>"
      )

    #diagnostics
    # save(pca.cd.fac.map, file = "factormap.rdata")

    return(pca.cd.fac.map)
  })

  output$pca.cd.fac.map <- renderPlot({
    pca.cd.fac.map()
  })
  output$pca.cd.fac.map.2 <- renderPlot({
    pca.cd.fac.map()
  })

  ##### Correlation circle----

  circle <- circleFun(c(0, 0), 2, npoints = 100)

  pca.cd.var.cor <- reactive({
    data.frame(pca.cd()$var$cor) %>%
      select(input$dim.cd[1], input$dim.cd[2]) %>%
      mutate(var = rownames(data.frame(pca.cd()$var$cor))) %>%
      magrittr::set_colnames(c("Dim.1", "Dim.2", "var")) %>%
      left_join(
        data.frame(pca.cd()$var$cos2) %>%
          select(input$dim.cd[1], input$dim.cd[2]) %>%
          magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
          mutate(sum.cos2 = Dim.1 + Dim.2) %>%
          select(sum.cos2) %>%
          mutate(var = rownames(data.frame(pca.cd()$var$cos2))),
        by = 'var'
      ) %>%
      ggplot(
        aes(
          x = Dim.1,
          y = Dim.2,
          color = sum.cos2
        )
      ) +
      geom_path(
        data = circle,
        aes(x = cir.x, y = cir.y),
        inherit.aes = FALSE
      ) +
      geom_segment(
        aes(
          x = 0,
          y = 0,
          xend = Dim.1,
          yend = Dim.2
        ),
        arrow = arrow(length = unit(0.1, "inches")),
        show.legend = TRUE
      ) +
      geom_text_repel(
        aes(label = stringr::str_extract(var, "\\d+") %>% as.numeric()),
        show.legend = FALSE,
        size = 5,
        fontface = 'bold'
      ) +
      scale_x_continuous(limits = c(-1, 1)) +
      scale_y_continuous(limits = c(-1, 1)) +
      labs(
        y = input$dim.cd[2],
        x = input$dim.cd[1]
      ) +
      custom.theme.markdown +
      scale_color_continuous(
        type = 'viridis',
        limits = c(0, 1),
        name = "sum cos<sup>2</sup>",
        guide = guide_colorbar(barwidth = 15)
      )
  })

  output$pca.cd.var.cor <- renderPlot({
    pca.cd.var.cor()
  })
  output$pca.cd.var.cor.2 <- renderPlot({
    pca.cd.var.cor()
  })

  #### Variable coordinates-------
  pca.cd.var.coord <- reactive({
    data.frame(pca.cd()$var$coord) %>%
      select(input$dim.cd[1], input$dim.cd[2]) %>%
      magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
      mutate(var = rownames(data.frame(pca.cd()$var$coord))) %>%
      left_join(
        data.frame(pca.cd()$var$cos2) %>%
          select(input$dim.cd[1], input$dim.cd[2]) %>%
          magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
          mutate(sum.cos2 = Dim.1 + Dim.2) %>%
          select(sum.cos2) %>%
          mutate(var = rownames(data.frame(pca.cd()$var$cos2))),
        by = 'var'
      ) %>%
      ggplot(
        aes(
          x = Dim.1,
          y = Dim.2,
          color = sum.cos2
        ),
        show.legend = FALSE
      ) +
      geom_segment(
        aes(
          x = 0,
          y = 0,
          xend = Dim.1,
          yend = Dim.2
        ),
        arrow = arrow(length = unit(0.1, "inches")),
        show.legend = TRUE
      ) +
      geom_text_repel(
        aes(label = stringr::str_extract(var, "\\d+") %>% as.numeric()),
        show.legend = FALSE,
        size = 5,
        fontface = 'bold'
      ) +
      custom.theme.markdown +
      scale_color_continuous(
        type = 'viridis',
        limits = c(0, 1),
        name = "sum cos<sup>2</sup>",
        guide = guide_colorbar(barwidth = 15)
      ) +
      labs(
        y = input$dim.cd[2],
        x = input$dim.cd[1]
      )
  })

  output$pca.cd.var.coord <- renderPlot({
    pca.cd.var.coord()
  })
  output$pca.cd.var.coord.2 <- renderPlot({
    pca.cd.var.coord()
  })

  #### Individuals coordinates--------

  output$pca.cd.table <- renderDT(server = FALSE, {
    if (file.toggle() == 'no') {
      return(NULL)
    } else {
      data.frame(pca.cd()$ind$coord) %>%
        cbind(training.cd() %>% select(oligo)) %>%
        add_column(km = pca.cd.km()$cluster) %>%
        left_join(
          ref.seq() %>%
            select(
              oligo,
              topo,
              conformer,
              feature,
              gba,
              gba.stacks,
              tetrad,
              tetrad.id,
              loop,
              plus.minus,
              groove,
              salt
            ),
          by = 'oligo'
        ) %>%
        setcolorder(
          c(
            'oligo',
            'topo',
            'conformer',
            'feature',
            'gba',
            'gba.stacks',
            'tetrad',
            'tetrad.id',
            'loop',
            'plus.minus',
            'groove',
            'salt',
            'km'
          )
        ) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          colnames = c(
            "Name" = "oligo",
            "Cation" = "salt",
            "Topology" = "topo",
            "Conformer" = "conformer",
            "Other feature" = "feature",
            "GBA" = "gba",
            "GBA stacks" = "gba.stacks",
            "Tetrads" = "tetrad",
            "Tetrad combination" = "tetrad.id",
            "Loop progression" = "loop",
            "Tetrad handedness" = "plus.minus",
            "Grooves" = "groove",
            "Cluster" = "km"
          ),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "pca cd"),
              list(extend = 'excel', title = NULL, filename = "pca cd"),
              list(extend = 'colvis')
            ),
            title = NULL,
            columnDefs = list(list(visible = FALSE, targets = c(2:11)))
          )
        ) %>%
        formatSignif(columns = 13:(13 + input$ncp), digits = 5)
    }
  })

  pca.cd.coord <- reactive({
    data.frame(pca.cd()$ind$coord) %>%
      select(input$dim.cd[1], input$dim.cd[2]) %>%
      cbind(training.cd()) %>%
      add_column(km = pca.cd.km()$cluster) %>%
      left_join(
        ref.seq() %>%
          select(
            oligo,
            topo,
            conformer,
            gba,
            gba.stacks,
            tetrad,
            tetrad.id,
            loop,
            plus.minus,
            groove,
            salt
          ),
        by = 'oligo'
      ) %>%
      pca.plotR(
        .,
        dim.1 = input$dim.cd[1],
        dim.2 = input$dim.cd[2],
        color = input$pca.color.cd,
        shape = input$pca.shape.cd
      )
  })

  output$pca.cd.coord <- renderPlot({
    pca.cd.coord()
  })

  pca.cd.coord.2 <- reactive({
    data.frame(pca.cd()$ind$coord) %>%
      select(input$dim.cd.2[1], input$dim.cd.2[2]) %>%
      cbind(training.cd()) %>%
      add_column(km = pca.cd.km()$cluster) %>%
      left_join(
        ref.seq() %>%
          select(
            oligo,
            topo,
            conformer,
            gba,
            gba.stacks,
            tetrad,
            tetrad.id,
            loop,
            plus.minus,
            groove,
            salt
          ),
        by = 'oligo'
      ) %>%
      pca.plotR(
        .,
        dim.1 = input$dim.cd.2[1],
        dim.2 = input$dim.cd.2[2],
        color = input$pca.color.cd.2,
        shape = input$pca.shape.cd.2
      )
  })

  output$pca.cd.coord.2 <- renderPlot({
    pca.cd.coord.2()
  })

  output$dwn.pca.cd <- downloadHandler(
    filename = function() {
      paste("PCA - CD -- view 1", input$device, sep = ".")
    },
    content = function(file) {
      if (input$device == "png") {
        ggsave(
          file,
          plot = pca.cd.coord(),
          device = input$device,
          dpi = input$resolution,
          width = input$width,
          height = input$height,
          scaling = input$scaling
        )
      } else {
        ggsave(
          file,
          plot = pca.cd.coord(),
          device = input$device,
          dpi = input$resolution,
          width = input$width,
          height = input$height
        )
      }
    }
  )

  output$dwn.pca.cd.2 <- downloadHandler(
    filename = function() {
      paste("PCA - CD -- view 2", input$device, sep = ".")
    },
    content = function(file) {
      if (input$device == "png") {
        ggsave(
          file,
          plot = pca.cd.coord.2(),
          device = input$device,
          dpi = input$resolution,
          width = input$width,
          height = input$height,
          scaling = input$scaling
        )
      } else {
        ggsave(
          file,
          plot = pca.cd.coord.2(),
          device = input$device,
          dpi = input$resolution,
          width = input$width,
          height = input$height
        )
      }
    }
  )

  #### Parameter table----

  output$param.cd.table <- renderDT(server = FALSE, {
    combined_list <- list(
      'Cations' = input$ref.seq.cation.0,
      'Topologies' = input$ref.seq.topo.0,
      'Conformers' = input$ref.seq.conformer.0,
      'GBA' = input$ref.seq.gba.0,
      'GBA stacks' = input$ref.seq.gba.stacks.0,
      'Oligonucleotides' = input$ref.seq.oligo.0,
      'Tetrad #' = input$ref.seq.tetrad.0,
      'Tetrad combination' = input$ref.seq.tetrad.id.0,
      'Loops' = input$ref.seq.loop.0,
      '+/-' = input$ref.seq.plus.minus.0,
      'Grooves' = input$ref.seq.groove.0,
      'Max wavelength' = max(input$wl),
      'Min wavelength' = min(input$wl),
      'Theoretical UV reference' = input$ids.ref.select,
      'IDS normalization' = input$ids.norm,
      'CD normalization' = input$cd.norm,
      'PCA dimensions' = input$ncp,
      'Variance scaling' = input$scale.unit,
      'k-means algorithm' = input$k.mean.algo,
      'k-means centers' = input$cluster.center
    )

    max.length <- max(sapply(combined_list, length))

    for (i in 1:length(combined_list)) {
      length(combined_list[[i]]) <- max.length
    }

    do.call(cbind, combined_list) %>%
      datatable(
        style = "bootstrap",
        extensions = c('Buttons', 'Responsive', 'Scroller'),
        rownames = F,
        escape = T,
        filter = 'top',
        autoHideNavigation = T,
        options = list(
          deferRender = TRUE,
          scrollY = 200,
          scroller = F,
          pageLength = 25,
          autoWidth = F,
          dom = 'Bfrtip',
          buttons = list(
            list(extend = 'copy'),
            list(extend = 'csv', title = NULL, filename = "pca parameters"),
            list(extend = 'excel', title = NULL, filename = "pca parameters"),
            list(extend = 'colvis')
          ),
          title = NULL,
          columnDefs = list(list(visible = FALSE, targets = c(0:8)))
        )
      )
  })

  ### IDS----

  ####PCA----

  pca.ids <- reactive({
    FactoMineR::PCA(
      X = training.ids()[, -1],
      # select(-salt),
      ncp = input$ncp,
      scale.unit = input$scale.unit,
      graph = FALSE
    )
  })

  #### Clusters----

  #Optimal number of clusters
  pca.ids.twss <- reactive({
    data.frame(
      clusters = 1:10
    ) %>%
      group_by(clusters) %>%
      mutate(
        totwithin = kmeans(
          pca.ids()$ind$coord,
          algorithm = input$k.mean.algo,
          centers = clusters,
          nstart = 100
        )$tot.withinss
      ) %>%
      ggplot(
        aes(x = clusters, y = totwithin)
      ) +
      geom_point(size = 4) +
      geom_line(size = 1) +
      scale_x_continuous(breaks = 1:10) +
      labs(
        x = "clusters",
        y = "total within sum of square"
      ) +
      custom.theme.markdown
  })

  output$pca.ids.twss <- renderPlot({
    pca.ids.twss()
  })

  output$pca.ids.twss.2 <- renderPlot({
    pca.ids.twss()
  })

  pca.ids.gap <- reactive({
    clusGap(
      pca.ids()$ind$coord,
      FUN = kmeans,
      nstart = 25,
      K.max = 10,
      B = 50
    )$Tab %>%
      as_tibble() %>%
      rownames_to_column(var = "cluster") %>%
      mutate(cluster = as.numeric(cluster)) %>%
      ggplot(
        aes(x = cluster, y = gap)
      ) +
      geom_vline(
        data = . %>% filter(gap == max(gap)),
        aes(xintercept = cluster),
        linetype = 'dashed',
        size = 1,
        color = "tomato"
      ) +
      geom_errorbar(
        aes(x = cluster, ymin = gap - SE.sim, ymax = gap + SE.sim),
        width = 0
      ) +
      geom_point(size = 5) +
      geom_line(size = 1) +
      scale_x_continuous(breaks = 1:10) +
      labs(
        x = "clusters",
        y = "gap statistic"
      ) +
      custom.theme.markdown
  })

  output$pca.ids.gap <- renderPlot({
    pca.ids.gap()
  })

  output$pca.ids.gap.2 <- renderPlot({
    pca.ids.gap()
  })

  #k-means clustering
  #done on all 5 kept dimensions
  pca.ids.km <- reactive({
    set.seed(42)

    kmeans(
      pca.ids()$ind$coord,
      algorithm = input$k.mean.algo,
      centers = input$cluster.center,
      nstart = 100
    )
  })

  #### Analytics----

  ##### Scree----
  pca.ids.scree <- reactive({
    pca.ids.scree <- data.frame(pca.ids()$eig) %>%
      mutate(dim = 1:nrow(data.frame(pca.ids()$eig))) %>%
      filter(dim <= 10) %>%
      ggplot(aes(x = dim, y = percentage.of.variance)) +
      geom_bar(
        stat = 'identity',
        aes(fill = -dim),
        show.legend = FALSE
      ) +
      geom_point(
        aes(
          x = dim,
          y = cumsum(percentage.of.variance)
        ),
        size = 5
      ) +
      geom_line(
        aes(
          x = dim,
          y = cumsum(percentage.of.variance)
        ),
        size = 1
      ) +
      custom.theme.markdown +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
      scale_y_continuous(expand = c(0, NA), limits = c(0, 105)) +
      labs(
        y = 'Variance contribution (%)',
        x = 'Dimension'
      )

    #diagnostics
    # save(pca.ids.scree, file = "scree.rdata")

    return(pca.ids.scree)
  })

  output$pca.ids.scree <- renderPlot({
    pca.ids.scree()
  })

  output$pca.ids.scree.2 <- renderPlot({
    pca.ids.scree()
  })

  ##### Factor map----
  pca.ids.fac.map <- reactive({
    pca.cd.fac.map <- data.frame(pca.ids()$var$cos2) %>%
      select(input$dim.ids) %>%
      mutate(
        var = as.numeric(substr(rownames(data.frame(pca.ids()$var$cos2)), 5, 7))
      ) %>%
      pivot_longer(
        cols = 1:length(input$dim.ids),
        names_to = 'dim',
        values_to = 'cos2'
      ) %>%
      mutate(
        dim = substr(dim, 5, 5)
      ) %>%
      ggplot(
        aes(
          x = dim,
          y = var,
          size = cos2
        )
      ) +
      geom_point(
        color = "#2AA198",
        alpha = 0.5
      ) +
      labs(
        y = 'Variable',
        x = 'Dimension'
      ) +
      custom.theme.markdown +
      scale_size_continuous(
        limits = c(0, 1),
        range = c(0, 25),
        name = "cos<sup>2</sup>",
      )

    #diagnostics
    # save(pca.cd.fac.map, file = "factormap.rdata")

    return(pca.cd.fac.map)
  })

  output$pca.ids.fac.map <- renderPlot({
    pca.ids.fac.map()
  })
  output$pca.ids.fac.map.2 <- renderPlot({
    pca.ids.fac.map()
  })

  ##### Correlation circle----
  pca.ids.var.cor <- reactive({
    data.frame(pca.ids()$var$cor) %>%
      select(input$dim.ids[1], input$dim.ids[2]) %>%
      mutate(var = rownames(data.frame(pca.ids()$var$cor))) %>%
      magrittr::set_colnames(c("Dim.1", "Dim.2", "var")) %>%
      left_join(
        data.frame(pca.ids()$var$cos2) %>%
          select(input$dim.ids[1], input$dim.ids[2]) %>%
          magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
          mutate(sum.cos2 = Dim.1 + Dim.2) %>%
          select(sum.cos2) %>%
          mutate(var = rownames(data.frame(pca.ids()$var$cos2))),
        by = 'var'
      ) %>%
      ggplot(
        aes(
          x = Dim.1,
          y = Dim.2,
          color = sum.cos2
        )
      ) +
      geom_path(
        data = circle,
        aes(x = cir.x, y = cir.y),
        inherit.aes = FALSE
      ) +
      geom_segment(
        aes(
          x = 0,
          y = 0,
          xend = Dim.1,
          yend = Dim.2
        ),
        arrow = arrow(length = unit(0.1, "inches")),
        show.legend = TRUE
      ) +
      geom_text_repel(
        aes(label = stringr::str_extract(var, "\\d+") %>% as.numeric()),
        show.legend = FALSE,
        size = 5,
        fontface = 'bold'
      ) +
      scale_x_continuous(limits = c(-1, 1)) +
      scale_y_continuous(limits = c(-1, 1)) +
      labs(
        y = input$dim.ids[2],
        x = input$dim.ids[1]
      ) +
      custom.theme.markdown +
      scale_color_continuous(
        type = 'viridis',
        limits = c(0, 1),
        name = "sum cos<sup>2</sup>",
        guide = guide_colorbar(barwidth = 15)
      )
  })

  output$pca.ids.var.cor <- renderPlot({
    pca.ids.var.cor()
  })

  output$pca.ids.var.cor.2 <- renderPlot({
    pca.ids.var.cor()
  })

  #### Variable coordinates-----

  pca.ids.var.coord <- reactive({
    data.frame(pca.ids()$var$coord) %>%
      select(input$dim.ids[1], input$dim.ids[2]) %>%
      magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
      mutate(var = rownames(data.frame(pca.ids()$var$coord))) %>%
      left_join(
        data.frame(pca.ids()$var$cos2) %>%
          select(input$dim.ids[1], input$dim.ids[2]) %>%
          magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
          mutate(sum.cos2 = Dim.1 + Dim.2) %>%
          select(sum.cos2) %>%
          mutate(var = rownames(data.frame(pca.ids()$var$cos2))),
        by = 'var'
      ) %>%
      ggplot(
        aes(
          x = Dim.1,
          y = Dim.2,
          color = sum.cos2
        ),
        show.legend = FALSE
      ) +
      geom_segment(
        aes(
          x = 0,
          y = 0,
          xend = Dim.1,
          yend = Dim.2
        ),
        arrow = arrow(length = unit(0.1, "inches")),
        show.legend = TRUE
      ) +
      geom_text_repel(
        aes(label = stringr::str_extract(var, "\\d+") %>% as.numeric()),
        show.legend = FALSE,
        size = 5,
        fontface = 'bold'
      ) +
      custom.theme.markdown +
      labs(
        y = input$dim.ids[2],
        x = input$dim.ids[1]
      ) +
      scale_color_continuous(
        type = 'viridis',
        limits = c(0, 1),
        name = "sum cos<sup>2</sup>",
        guide = guide_colorbar(barwidth = 15)
      )
  })

  output$pca.ids.var.coord <- renderPlot({
    pca.ids.var.coord()
  })
  output$pca.ids.var.coord.2 <- renderPlot({
    pca.ids.var.coord()
  })

  #### Individuals coordinates----

  output$pca.ids.table <- renderDT(server = FALSE, {
    if (file.toggle() == 'no') {
      return(NULL)
    } else {
      data.frame(pca.ids()$ind$coord) %>%
        cbind(training.ids() %>% select(oligo)) %>%
        add_column(km = pca.ids.km()$cluster) %>%
        left_join(
          ref.seq() %>%
            select(
              oligo,
              topo,
              conformer,
              feature,
              gba,
              gba.stacks,
              tetrad,
              tetrad.id,
              loop,
              plus.minus,
              groove,
              salt
            ),
          by = 'oligo'
        ) %>%
        setcolorder(c(
          'oligo',
          'topo',
          'conformer',
          'feature',
          'gba',
          'gba.stacks',
          'tetrad',
          'tetrad.id',
          'loop',
          'plus.minus',
          'groove',
          'salt',
          'km'
        )) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          colnames = c(
            "Name" = "oligo",
            "Cation" = "salt",
            "Topology" = "topo",
            "Conformer" = "conformer",
            "Other feature" = "feature",
            "GBA" = "gba",
            "GBA stacks" = "gba.stacks",
            "Tetrads" = "tetrad",
            "Tetrad combination" = "tetrad.id",
            "Loop progression" = "loop",
            "Tetrad handedness" = "plus.minus",
            "Grooves" = "groove",
            "Cluster" = "km"
          ),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "pca ids"),
              list(extend = 'excel', title = NULL, filename = "pca ids"),
              list(extend = 'colvis')
            ),
            title = NULL,
            columnDefs = list(list(visible = FALSE, targets = c(2:11)))
          )
        ) %>%
        formatRound(columns = 13:(13 + input$ncp), digits = 2)
    }
  })

  pca.ids.coord <- reactive({
    data.frame(pca.ids()$ind$coord) %>%
      select(input$dim.ids[1], input$dim.ids[2]) %>%
      cbind(training.ids()) %>%
      add_column(km = pca.ids.km()$cluster) %>%
      left_join(
        ref.seq() %>%
          select(
            oligo,
            topo,
            conformer,
            gba,
            gba.stacks,
            tetrad,
            tetrad.id,
            loop,
            plus.minus,
            groove,
            salt
          ),
        by = 'oligo'
      ) %>%
      pca.plotR(
        .,
        dim.1 = input$dim.ids[1],
        dim.2 = input$dim.ids[2],
        color = input$pca.color.ids,
        shape = input$pca.shape.ids
      )
  })

  output$pca.ids.coord <- renderPlot({
    pca.ids.coord()
  })

  pca.ids.coord.2 <- reactive({
    data.frame(pca.ids()$ind$coord) %>%
      select(input$dim.ids.2[1], input$dim.ids.2[2]) %>%
      cbind(training.ids()) %>%
      add_column(km = pca.ids.km()$cluster) %>%
      left_join(
        ref.seq() %>%
          select(
            oligo,
            topo,
            conformer,
            gba,
            gba.stacks,
            tetrad,
            tetrad.id,
            loop,
            plus.minus,
            groove,
            salt
          ),
        by = 'oligo'
      ) %>%
      pca.plotR(
        .,
        dim.1 = input$dim.ids.2[1],
        dim.2 = input$dim.ids.2[2],
        color = input$pca.color.ids.2,
        shape = input$pca.shape.ids.2
      )
  })

  output$pca.ids.coord.2 <- renderPlot({
    pca.ids.coord.2()
  })

  output$dwn.pca.ids <- downloadHandler(
    filename = function() {
      paste("PCA - IDS -- view 1", input$device, sep = ".")
    },
    content = function(file) {
      if (input$device == "png") {
        ggsave(
          file,
          plot = pca.ids.coord(),
          device = input$device,
          dpi = input$resolution,
          width = input$width,
          height = input$height,
          scaling = input$scaling
        )
      } else {
        ggsave(
          file,
          plot = pca.ids.coord(),
          device = input$device,
          dpi = input$resolution,
          width = input$width,
          height = input$height
        )
      }
    }
  )

  output$dwn.pca.ids.2 <- downloadHandler(
    filename = function() {
      paste("PCA - IDS -- view 2", input$device, sep = ".")
    },
    content = function(file) {
      if (input$device == "png") {
        ggsave(
          file,
          plot = pca.ids.coord.2(),
          device = input$device,
          dpi = input$resolution,
          width = input$width,
          height = input$height,
          scaling = input$scaling
        )
      } else {
        ggsave(
          file,
          plot = pca.ids.coord.2(),
          device = input$device,
          dpi = input$resolution,
          width = input$width,
          height = input$height
        )
      }
    }
  )

  #### Parameter table----

  output$param.ids.table <- renderDT(server = FALSE, {
    combined_list <- list(
      'Cations' = input$ref.seq.cation.0,
      'Topologies' = input$ref.seq.topo.0,
      'GBA' = input$ref.seq.gba.0,
      'Oligonucleotides' = input$ref.seq.oligo.0,
      'Tetrad #' = input$ref.seq.tetrad.0,
      'Tetrad combination' = input$ref.seq.tetrad.id.0,
      'Loops' = input$ref.seq.loop.0,
      '+/-' = input$ref.seq.plus.minus.0,
      'Grooves' = input$ref.seq.groove.0,
      'Max wavelength' = max(input$wl),
      'Min wavelength' = min(input$wl),
      'Theoretical UV reference' = input$ids.ref.select,
      'IDS normalization' = input$ids.norm,
      'CD normalization' = input$cd.norm,
      'PCA dimensions' = input$ncp,
      'Variance scaling' = input$scale.unit,
      'k-means algorithm' = input$k.mean.algo,
      'k-means centers' = input$cluster.center
    )

    max.length <- max(sapply(combined_list, length))

    for (i in 1:length(combined_list)) {
      length(combined_list[[i]]) <- max.length
    }

    do.call(cbind, combined_list) %>%
      datatable(
        style = "bootstrap",
        extensions = c('Buttons', 'Responsive', 'Scroller'),
        rownames = F,
        escape = T,
        filter = 'top',
        autoHideNavigation = T,
        options = list(
          deferRender = TRUE,
          scrollY = 200,
          scroller = F,
          pageLength = 25,
          autoWidth = F,
          dom = 'Bfrtip',
          buttons = list(
            list(extend = 'copy'),
            list(extend = 'csv', title = NULL, filename = "pca parameters"),
            list(extend = 'excel', title = NULL, filename = "pca parameters"),
            list(extend = 'colvis')
          ),
          title = NULL,
          columnDefs = list(list(visible = FALSE, targets = c(0:8)))
        )
      )
  })

  ### Investigations----

  #### CD----

  pca.cd.invest <- observeEvent(input$button.cd.invest, {
    # thematic::thematic_off()

    withProgress(
      message = 'Investigating CD data',
      detail = 'Please wait',
      value = 0,
      {
        incProgress(amount = 1 / 2)

        FactoInvestigate::Investigate(pca.cd())

        incProgress(amount = 2 / 2)
      }
    )

    # thematic::thematic_on()
  })

  #### IDS----
  pca.ids.invest <- observeEvent(input$button.ids.invest, {
    # thematic::thematic_off()

    withProgress(
      message = 'Investigating IDS data',
      detail = 'Please wait',
      value = 0,
      {
        incProgress(amount = 1 / 2)

        FactoInvestigate::Investigate(pca.ids())

        incProgress(amount = 2 / 2)
      }
    )

    # thematic::thematic_on()
  })

  pca.cd.ids.invest <- observeEvent(input$button.cd.ids.invest, {
    # thematic::thematic_off()

    withProgress(
      message = 'Investigating CD+IDS data',
      detail = 'Please wait',
      value = 0,
      {
        incProgress(amount = 1 / 2)

        FactoInvestigate::Investigate(pca.cd.ids())

        incProgress(amount = 2 / 2)
      }
    )

    # thematic::thematic_on()
  })

  # User prediction----

  ## input----

  ### data files----
  # input.file <- reactive({
  #   input$ref.data
  # })

  input.file.user <- reactive({
    input$user.data
  })

  #### file import toggle----
  # file.toggle <- reactive({
  #   if(is.null(input.file())){
  #     return('no')
  #   } else {
  #     return('yes')
  #   }
  # })

  file.toggle.user <- reactive({
    if (is.null(input.file.user())) {
      return('no')
    } else {
      return('yes')
    }
  })

  ##### sequences----

  user.seq <- reactive({
    read_excel(
      input.file.user()$datapath,
      sheet = 'Sample information'
    )
  })

  #####picker inputs----
  ###### source----
  output$user.seq.oligo.0_ui <- renderUI({
    if (file.toggle.user() == 'no') {
      selectizeInput(
        "user.seq.oligo.0",
        label = "User oligonucleotide",
        choices = "upload data first",
        multiple = T
      )
    } else {
      choices <- sort(unique(user.seq()$oligo))

      selectizeInput(
        "user.seq.oligo.0",
        label = "User oligonucleotide",
        choices = choices,
        selected = choices,
        multiple = T
      )
    }
  })

  #######plot-----
  output$user.seq.oligo_ui <- renderUI({
    if (file.toggle.user() == 'no') {
      selectizeInput(
        "user.seq.oligo",
        label = "Oligonucleotide",
        choices = "upload data first",
        multiple = T
      )
    } else {
      selectizeInput(
        "user.seq.oligo",
        label = "Oligonucleotide",
        choices = unique(input$user.seq.oligo.0),
        selected = unique(input$user.seq.oligo.0),
        multiple = T
      )
    }
  })

  ##### CD----

  user.cd.input <- reactive({
    user_cd_input <- read_excel(
      input.file.user()$datapath,
      sheet = "CD"
    ) %>%
      pivot_longer(
        cols = 2:(ncol(.)),
        names_to = 'oligo',
        values_to = 'cd'
      ) %>%
      left_join(
        user.seq(),
        by = "oligo"
      ) %>%
      filter(wl %% 1 == 0) %>% #removing non integer wl to have the same number of data points than in IDS
      group_by(oligo)

    # writexl::write_xlsx(
    #   user_cd_input,
    #   path = "user_cd.xlsx"
    # )

    user_cd_input %>%
      left_join(
        user.uv.input() %>%
          filter(cation == "none") %>%
          select(wl, oligo, eps) %>%
          unique() %>%
          set_colnames(c("wl", "oligo", "eps.uv")),
        by = c("oligo", "wl")
      ) %>%
      mutate(
        cd = cd - mean(cd[wl > 320]),
        delta.eps = cd / (32980 * l * Concentration / 1E6)
      ) %>%
      mutate(
        norm.cd = case_when(
          input$cd.norm == "Δε/ε" ~ delta.eps / eps.uv,
          input$cd.norm == "-1/+1" ~
            2 *
              (delta.eps - min(delta.eps)) /
              (max(delta.eps) - min(delta.eps)) -
              1,
          TRUE ~ delta.eps
        )
      ) %>%
      ungroup() %>%
      filter(
        wl <= max(input$wl),
        wl >= min(input$wl),
        wl %% 1 == 0,
        oligo %in% input$user.seq.oligo.0
      ) %>%
      arrange(wl)
  })

  ##### UV----

  #data with cation
  user.uv.input.cation <- reactive({
    read_excel(
      input.file.user()$datapath,
      sheet = "UV"
    ) %>%
      add_column(cation = 'cation') %>%
      arrange(wl)
  })

  #data without cation
  user.uv.input.no.cation <- reactive({
    read_excel(
      input.file.user()$datapath,
      sheet = "UV - no cation"
    ) %>%
      add_column(cation = 'no cation')
  })

  #merge data
  user.uv.input <- reactive({
    #data rbinding
    user.uv.input.cation() %>%
      pivot_longer(
        cols = 2:(ncol(.) - 1),
        values_to = 'abs',
        names_to = 'oligo'
      ) %>%
      rbind(
        user.uv.input.no.cation() %>%
          pivot_longer(
            cols = 2:(ncol(.) - 1),
            values_to = 'abs',
            names_to = 'oligo'
          )
      ) %>%
      left_join(
        user.seq(),
        by = "oligo"
      ) %>%
      group_by(oligo, cation) %>%
      mutate(
        abs = abs - mean(abs[wl > 320]),
        eps = abs / (1 * Concentration / 1E6)
      ) %>%
      ungroup() %>%
      filter(
        wl <= max(input$wl),
        wl >= min(input$wl),
        wl %% 1 == 0,
        oligo %in% input$user.seq.oligo.0
      )
  })

  # user.uv.th <- reactive({
  #   read_excel(
  #     input.file.user()$datapath,
  #     sheet = "th.spectra"
  #   )
  # })

  # user.uv.th <- reactive({
  #   user.uv.input() %>%
  #     filter(cation == 'no cation') %>%
  #     group_by(oligo) %>%
  #     select(oligo, wl) %>%
  #     unique() %>%
  #     left_join(
  #       user.seq() %>%
  #         select(oligo, seq)
  #     )
  # })
  #

  ##### IDS----

  user.ids.input <- reactive({
    if (isFALSE(input$ids.ref.select)) {
      user.ids.input <- user.uv.input()
    } else {
      withProgress(
        message = 'Calculating theoretical spectra',
        detail = 'Please wait',
        value = 0,
        {
          incProgress(amount = 1 / 3)

          user.ids.input <- user.uv.input() %>%
            filter(cation == 'no cation') %>%
            group_by(oligo, wl) %>%
            mutate(
              eps = dt.spec.calcR(
                input.contrib = dt.contributR(input.seq = seq),
                input.wl = wl
              )
            ) %>%
            rbind(
              user.uv.input() %>%
                filter(cation == 'cation')
            )

          incProgress(amount = 3 / 3)
        }
      )
    }

    user.ids.input %>%
      group_by(oligo, wl) %>%
      mutate(
        delta.eps = eps[cation == 'no cation'] - eps[cation == 'cation']
      ) %>%
      filter(cation == 'cation') %>%
      group_by(oligo) %>%
      mutate(
        norm.ids = case_when(
          input$ids.norm == "-1/+1" ~
            2 *
              (delta.eps - min(delta.eps)) /
              (max(delta.eps) - min(delta.eps)) -
              1,
          TRUE ~ delta.eps
        )
      )
  })

  #### outputs of inputs----

  ##### CD-----

  output$user.cd.input <- renderDT(server = FALSE, {
    if (file.toggle.user() == 'no') {
      return(NULL)
    } else {
      user.cd.input() %>%
        mutate(
          Concentration = round(Concentration, 2),
          cd = round(cd, 2),
          delta.eps = round(delta.eps, 2),
          norm.cd = round(norm.cd, 2)
        ) %>%
        setcolorder(c(
          "oligo",
          "Concentration",
          "wl",
          "cd",
          "delta.eps",
          "norm.cd"
        )) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          colnames = c(
            "Wavelength (nm)" = "wl",
            "Name" = "oligo",
            "Molar ellipticity" = "delta.eps",
            "Ellipticity (mdeg)" = "cd",
            "Concentration (µM)" = "Concentration",
            "Normalized CD" = "norm.cd"
          ),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "user cd"),
              list(extend = 'excel', title = NULL, filename = "user cd"),
              list(extend = 'colvis')
            ),
            title = NULL
            # columnDefs = list(list(visible=FALSE, targets=c(0,2:6,9,12)))
          )
        )
      # formatStyle(columns = 0:6, target = "cell", backgroundColor = "#272c30")
    }
  })

  output$p.user.cd <- renderPlot({
    user.cd <- user.cd.input() %>%
      filter(oligo %in% input$user.seq.oligo)

    if (input$cd.norm == "Δε") {
      user.cd$y <- user.cd$delta.eps
      axis.label <- "&Delta;&eps; (M<sup>-1</sup> cm<sup>-1</sup>)"
    } else {
      user.cd$y <- user.cd$norm.cd
      axis.label <- "Normalized CD"
    }

    if (file.toggle.user() == 'no') {
      ggplot()
    } else {
      p.user.cd <- user.cd %>%
        ggplot(
          aes(
            x = wl,
            y = y,
            color = oligo
          )
        ) +
        geom_hline(
          yintercept = 0,
          color = "white"
        ) +
        geom_line(
          aes(group = oligo),
          size = 1,
          show.legend = T
        ) +
        labs(
          x = "&lambda; (nm)",
          y = axis.label
          # color = input$user.color,
          # fill = input$user.color
        ) +
        custom.theme.markdown +
        scale_y_continuous(n.breaks = 4) +
        scale_x_continuous(expand = c(0, 0))
    }

    if (input$user.panel == "Panels") {
      p.user.cd <- p.user.cd +
        facet_wrap(~oligo, ncol = 4)

      return(p.user.cd)
    } else {
      return(p.user.cd)
    }
  })

  ##### UV----

  output$user.uv.input <- renderDT(server = FALSE, {
    if (file.toggle.user() == 'no') {
      return(NULL)
    } else {
      user.uv.input() %>%
        mutate(
          Concentration = round(Concentration, 2),
          abs = round(abs, 2)
        ) %>%
        setcolorder(c(
          "oligo",
          "cation",
          "Concentration",
          "wl",
          "abs",
          "eps"
        )) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          colnames = c(
            "Wavelength (nm)" = "wl",
            "Name" = "oligo",
            "A" = "abs",
            "Cation?" = "cation",
            "Concentration (µM)" = "Concentration",
            "Molar extinction coefficient" = "eps"
          ),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "user uv"),
              list(extend = 'excel', title = NULL, filename = "user uv"),
              list(extend = 'colvis')
            ),
            title = NULL
            # columnDefs = list(list(visible=FALSE, targets=c(0,2:6,9,12)))
          )
        )
      # formatStyle(columns = 0:5, target = "cell", backgroundColor = "#272c30")
    }
  })

  output$p.user.uv <- renderPlot({
    if (file.toggle.user() == 'no') {
      ggplot()
    } else {
      user.uv.input() %>%
        filter(oligo %in% input$user.seq.oligo) %>%
        ggplot(
          aes(
            x = wl,
            y = eps,
            color = cation,
            group = cation
          )
        ) +
        geom_line(
          size = 1,
          show.legend = T
        ) +
        facet_wrap(
          ~oligo,
          ncol = 4
        ) +
        labs(
          x = "&lambda; (nm)",
          y = "&epsilon; (M<sup>-1</sup> cm<sup>-1</sup>)",
          color = "Cation"
        ) +
        custom.theme.markdown +
        scale_y_continuous(n.breaks = 3) +
        scale_x_continuous(expand = c(0, 0))
    }
  })

  ##### IDS----

  output$user.ids.input <- renderDT(server = FALSE, {
    if (file.toggle.user() == 'no') {
      return(NULL)
    } else {
      user.ids.input() %>%
        select(-c("cation", "abs")) %>%
        mutate(
          Concentration = round(Concentration, 2),
          delta.eps = round(delta.eps, 0),
          norm.ids = round(norm.ids, 2)
        ) %>%
        setcolorder(c(
          "oligo",
          "Concentration",
          "wl",
          "delta.eps",
          "norm.ids"
        )) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          colnames = c(
            # "#" = "oligo.number",
            "Wavelength (nm)" = "wl",
            "Name" = "oligo",
            "IDS" = "delta.eps",
            "C (µM)" = "Concentration",
            "Normalized IDS" = "norm.ids",
            "Molar extinction coefficient" = "eps"
          ),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "user ids"),
              list(extend = 'excel', title = NULL, filename = "user ids"),
              list(extend = 'colvis')
            ),
            title = NULL
            # columnDefs = list(list(visible=FALSE, targets=c(0,2:6,9,12)))
          )
        )
      # formatStyle(
      #   columns = 0:15,
      #   target = "cell",
      #   backgroundColor = "#272c30"
      # )
    }
  })

  output$p.user.ids <- renderPlot({
    user.ids <- user.ids.input() %>%
      filter(oligo %in% input$user.seq.oligo) %>%
      add_column(
        y = NA
      )

    if (input$ids.norm == "Δε") {
      user.ids$y <- user.ids$delta.eps
      axis.label <- "&Delta;&epsilon; (M<sup>-1</sup> cm<sup>-1</sup>)"
    } else {
      user.ids$y <- user.ids$norm.ids
      axis.label <- "Normalized IDS"
    }

    if (file.toggle.user() == 'no') {
      ggplot()
    } else {
      p.user.ids <- user.ids %>%
        ggplot(
          aes(
            x = wl,
            y = y,
            color = oligo
          )
        ) +
        geom_vline(
          xintercept = 295,
          linetype = 'dashed',
          color = "white"
        ) +
        geom_vline(
          xintercept = 275,
          linetype = 'dashed',
          color = "white"
        ) +
        geom_vline(
          xintercept = 262,
          linetype = 'dashed',
          color = "white"
        ) +
        geom_vline(
          xintercept = 245,
          linetype = 'dashed',
          color = "white"
        ) +
        geom_hline(
          yintercept = 0,
          color = "white"
        ) +
        geom_line(
          aes(group = oligo),
          size = 1,
          show.legend = TRUE
        ) +
        labs(
          x = "&lambda; (nm)",
          y = axis.label
          # color = input$user.color,
          # fill = input$user.color
        ) +
        custom.theme.markdown +
        scale_y_continuous(n.breaks = 3) +
        scale_x_continuous(expand = c(0, 0))
    }

    if (input$user.panel == "Panels") {
      p.user.ids <- p.user.ids +
        facet_wrap(~oligo, ncol = 4)

      return(p.user.ids)
    } else {
      return(p.user.ids)
    }
  })

  ## PCA prediction----

  ### distance function----
  closest.cluster <- function(user.data, km.centers) {
    cluster.dist <- apply(
      km.centers,
      1,
      function(y) sqrt(sum((user.data - y)^2))
    )
    return(which.min(cluster.dist)[1])
  }

  ### CD----

  #### Prediction set prep----

  user.cd <- reactive({
    user.cd.input() %>%
      mutate(
        delta.eps = case_when(
          input$cd.norm == "Δε" ~ delta.eps,
          input$cd.norm == "Δε/ε" ~ norm.cd,
          input$cd.norm == "-1/+1" ~ norm.cd
        )
      ) %>%
      select(oligo, wl, delta.eps) %>%
      filter(wl %% 1 == 0) %>%
      mutate(wl = paste0("cd-", wl)) %>%
      pivot_wider(
        names_from = wl,
        values_from = delta.eps
      )
  })

  output$user.cd <- renderDT(server = FALSE, {
    if (file.toggle.user() == 'no') {
      return(NULL)
    } else {
      user.cd() %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "user cd"),
              list(extend = 'excel', title = NULL, filename = "user cd"),
              list(extend = 'colvis')
            ),
            title = NULL
          )
        )
      # formatStyle(
      #   columns = 0:50,
      #   target = "cell",
      #   backgroundColor = "#272c30"
      # )
    }
  })

  #### PCA predict----

  cd.predict <- reactive({
    #prediction
    cd.predict.0 <- FactoMineR::predict.PCA(
      pca.cd(),
      user.cd(),
      scale.unit = input$scale.unit
    )$coord

    #cluster assignment
    km.predict <- apply(
      cd.predict.0,
      1,
      closest.cluster,
      km.centers = pca.cd.km()$centers
    )

    #merging
    cd.predict <- cd.predict.0 %>%
      as.data.frame() %>%
      select(input$dim.cd[1], input$dim.cd[2]) %>%
      magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
      cbind(user.cd()) %>%
      cbind(km = km.predict) %>%
      mutate(topo = "User", gba = "User", salt = "User")

    return(cd.predict)
  })

  cd.predict.2 <- reactive({
    #prediction
    cd.predict.0 <- FactoMineR::predict.PCA(
      pca.cd(),
      user.cd(),
      scale.unit = input$scale.unit
    )$coord

    #cluster assignment
    km.predict <- apply(
      cd.predict.0,
      1,
      closest.cluster,
      km.centers = pca.cd.km()$centers
    )

    #merging
    cd.predict <- cd.predict.0 %>%
      as.data.frame() %>%
      select(input$dim.cd.2[1], input$dim.cd.2[2]) %>%
      magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
      cbind(user.cd()) %>%
      cbind(km = km.predict) %>%
      mutate(topo = "User", gba = "User", salt = "User")

    return(cd.predict)
  })

  ##### plots----

  output$predict.cd.table <- renderDT(server = FALSE, {
    if (file.toggle.user() == 'no') {
      return(NULL)
    } else {
      cd.predict() %>%
        rbind(
          data.frame(pca.cd()$ind$coord) %>%
            select(input$dim.cd[1], input$dim.cd[2]) %>%
            magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
            cbind(training.cd()) %>%
            add_column(km = pca.cd.km()$cluster) %>%
            left_join(
              ref.seq() %>%
                select(oligo, topo, gba, salt),
              by = 'oligo'
            )
        ) %>%
        add_column(color = NA, shape = NA) %>%
        left_join(
          ref.seq() %>%
            select(oligo, tetrad, tetrad.id, loop, plus.minus, groove),
          by = "oligo"
        ) %>%
        mutate(
          topo = if_else(is.na(topo), "User", topo),
          tetrad = if_else(is.na(tetrad), "User", as.character(tetrad)),
          tetrad.id = if_else(is.na(tetrad.id), "User", tetrad.id),
          loop = if_else(is.na(loop), "User", loop),
          plus.minus = if_else(is.na(plus.minus), "User", plus.minus),
          groove = if_else(is.na(groove), "User", groove)
        ) %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          colnames = c(
            "Name" = "oligo",
            "Cation" = "salt",
            "Topology" = "topo",
            "GBA" = "gba",
            "Tetrads" = "tetrad",
            "Tetrad combination" = "tetrad.id",
            "Loop progression" = "loop",
            "Tetrad handedness" = "plus.minus",
            "Grooves" = "groove",
            "Cluster" = "km"
          ),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "pca cd"),
              list(extend = 'excel', title = NULL, filename = "pca cd"),
              list(extend = 'colvis')
            ),
            title = NULL
            # columnDefs = list(list(visible=FALSE, targets=c(2:9)))
          )
        )
      # formatStyle(columns = 0:(11+input$ncp),  target = "cell",
      #             backgroundColor = "#272c30") %>%
      # formatSignif(columns = 11:(11+input$ncp), digits = 5)
    }
  })

  output$conditional_predict_cd_table <- renderUI({
    if (file.toggle.user() == 'no') {
      #text output if no file is uploaded
      h4("No user file uploaded. Please upload file in the Data tab")
    } else {
      DTOutput("predict.cd.table")
    }
  })

  predict.cd.coord <- reactive({
    pca.cd.coord <- cd.predict() %>%
      rbind(
        data.frame(pca.cd()$ind$coord) %>%
          ungroup() %>%
          select(input$dim.cd[1], input$dim.cd[2]) %>%
          magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
          cbind(training.cd()) %>%
          add_column(km = pca.cd.km()$cluster) %>%
          left_join(
            ref.seq() %>%
              select(oligo, topo, gba, salt),
            by = 'oligo'
          )
      ) %>%
      left_join(
        ref.seq() %>%
          ungroup() %>%
          select(
            oligo,
            tetrad,
            tetrad.id,
            loop,
            plus.minus,
            groove,
            conformer,
            gba.stacks
          ),
        by = "oligo"
      ) %>%
      mutate(
        topo = if_else(is.na(topo), "User", topo),
        conformer = if_else(is.na(conformer), "User", conformer),
        gba.stacks = if_else(is.na(gba.stacks), "User", gba.stacks),
        tetrad = if_else(is.na(tetrad), "User", as.character(tetrad)),
        tetrad.id = if_else(is.na(tetrad.id), "User", tetrad.id),
        loop = if_else(is.na(loop), "User", loop),
        plus.minus = if_else(is.na(plus.minus), "User", plus.minus),
        groove = if_else(is.na(groove), "User", groove),
        salt = if_else(is.na(salt), "User", groove) # added
      )

    pca.cd.coord %>%
      pca.plotR(
        .,
        dim.1 = input$dim.cd[1],
        dim.2 = input$dim.cd[2],
        color = input$pca.color.cd,
        shape = input$pca.shape.cd
      )
  })

  output$predict.cd.coord <- renderPlot({
    predict.cd.coord()
  })

  output$conditional_pca_plot_cd <- renderUI({
    plotOutput(
      if (is.null(input$user.data)) {
        "pca.cd.coord"
      } else {
        "predict.cd.coord"
      },
      height = "800px"
    )
  })

  predict.cd.coord.2 <- reactive({
    pca.cd.coord <- cd.predict.2() %>%
      rbind(
        data.frame(pca.cd()$ind$coord) %>%
          ungroup() %>%
          select(input$dim.cd.2[1], input$dim.cd.2[2]) %>%
          magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
          cbind(training.cd()) %>%
          add_column(km = pca.cd.km()$cluster) %>%
          left_join(
            ref.seq() %>%
              select(oligo, topo, gba, salt),
            by = 'oligo'
          )
      ) %>%
      left_join(
        ref.seq() %>%
          ungroup() %>%
          select(
            oligo,
            tetrad,
            tetrad.id,
            loop,
            plus.minus,
            groove,
            conformer,
            gba.stacks
          ),
        by = "oligo"
      ) %>%
      mutate(
        topo = if_else(is.na(topo), "User", topo),
        conformer = if_else(is.na(conformer), "User", conformer),
        gba.stacks = if_else(is.na(gba.stacks), "User", gba.stacks),
        tetrad = if_else(is.na(tetrad), "User", as.character(tetrad)),
        tetrad.id = if_else(is.na(tetrad.id), "User", tetrad.id),
        loop = if_else(is.na(loop), "User", loop),
        plus.minus = if_else(is.na(plus.minus), "User", plus.minus),
        groove = if_else(is.na(groove), "User", groove),
        salt = if_else(is.na(salt), "User", groove) # added
      )

    pca.cd.coord %>%
      pca.plotR(
        .,
        dim.1 = input$dim.cd.2[1],
        dim.2 = input$dim.cd.2[2],
        color = input$pca.color.cd.2,
        shape = input$pca.shape.cd.2
      )
  })

  output$predict.cd.coord.2 <- renderPlot({
    predict.cd.coord.2()
  })

  output$conditional_pca_plot_cd_2 <- renderUI({
    plotOutput(
      if (is.null(input$user.data)) {
        "pca.cd.coord.2"
      } else {
        "predict.cd.coord.2"
      },
      height = "800px"
    )
  })

  #### IDS----

  ##### prep----

  user.ids <- reactive({
    user.ids.input() %>%
      mutate(
        delta.eps = case_when(
          input$ids.norm == "Δε" ~ delta.eps,
          input$ids.norm == "295 nm" ~ norm.ids,
          input$ids.norm == "-1/+1" ~ norm.ids
        )
      ) %>%
      select(oligo, wl, delta.eps) %>%
      mutate(wl = paste0("ids-", wl)) %>%
      pivot_wider(
        names_from = wl,
        values_from = delta.eps
      )
  })

  output$user.ids <- renderDT(server = FALSE, {
    if (file.toggle.user() == 'no') {
      return(NULL)
    } else {
      user.ids() %>%
        datatable(
          style = "bootstrap",
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          rownames = F,
          escape = T,
          filter = 'top',
          autoHideNavigation = T,
          options = list(
            deferRender = TRUE,
            scrollY = 200,
            scroller = F,
            pageLength = 25,
            autoWidth = F,
            dom = 'Bfrtip',
            buttons = list(
              list(extend = 'copy'),
              list(extend = 'csv', title = NULL, filename = "user ids"),
              list(extend = 'excel', title = NULL, filename = "user ids"),
              list(extend = 'colvis')
            ),
            title = NULL
            # columnDefs = list(list(visible=FALSE, targets=c(0,2:6,9,12)))
          )
        )
      # formatStyle(
      #   columns = 0:17,
      #   target = "cell",
      #   backgroundColor = "#272c30"
      # )
    }
  })

  #####predict----

  ids.predict <- reactive({
    #prediction
    ids.predict.0 <- FactoMineR::predict.PCA(
      pca.ids(),
      user.ids(),
      scale.unit = input$scale.unit
    )$coord

    #cluster assignment
    km.predict <- apply(
      ids.predict.0,
      1,
      closest.cluster,
      km.centers = pca.ids.km()$centers
    )

    #merging
    ids.predict <- ids.predict.0 %>%
      as.data.frame() %>%
      select(input$dim.ids[1], input$dim.ids[2]) %>%
      magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
      cbind(user.ids()) %>%
      cbind(km = km.predict) %>%
      mutate(topo = "User", gba = "User", salt = "User")

    return(ids.predict)
  })

  ids.predict.2 <- reactive({
    #prediction
    ids.predict.0 <- FactoMineR::predict.PCA(
      pca.ids(),
      user.ids(),
      scale.unit = input$scale.unit
    )$coord

    #cluster assignment
    km.predict <- apply(
      ids.predict.0,
      1,
      closest.cluster,
      km.centers = pca.ids.km()$centers
    )

    #merging
    ids.predict <- ids.predict.0 %>%
      as.data.frame() %>%
      select(input$dim.ids.2[1], input$dim.ids.2[2]) %>%
      magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
      cbind(user.ids()) %>%
      cbind(km = km.predict) %>%
      mutate(topo = "User", gba = "User", salt = "User")

    return(ids.predict)
  })

  #####plots----

  predict.ids.coord <- reactive({
    pca.ids.coord <- ids.predict() %>%
      rbind(
        data.frame(pca.ids()$ind$coord) %>%
          ungroup() %>%
          select(input$dim.ids[1], input$dim.ids[2]) %>%
          magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
          cbind(training.ids()) %>%
          add_column(km = pca.ids.km()$cluster) %>%
          left_join(
            ref.seq() %>%
              select(oligo, topo, gba, salt),
            by = 'oligo'
          )
      ) %>%
      left_join(
        ref.seq() %>%
          ungroup() %>%
          select(
            oligo,
            tetrad,
            tetrad.id,
            loop,
            plus.minus,
            groove,
            conformer,
            gba.stacks
          ),
        by = "oligo"
      ) %>%
      mutate(
        topo = if_else(is.na(topo), "User", topo),
        conformer = if_else(is.na(conformer), "User", conformer),
        gba_stacks = if_else(is.na(gba.stacks), "User", gba.stacks),
        tetrad = if_else(is.na(tetrad), "User", as.character(tetrad)),
        tetrad.id = if_else(is.na(tetrad.id), "User", tetrad.id),
        loop = if_else(is.na(loop), "User", loop),
        plus.minus = if_else(is.na(plus.minus), "User", plus.minus),
        groove = if_else(is.na(groove), "User", groove),
        salt = if_else(is.na(salt), "User", groove) # added
      )

    pca.ids.coord %>%
      pca.plotR(
        .,
        dim.1 = input$dim.ids[1],
        dim.2 = input$dim.ids[2],
        color = input$pca.color.ids,
        shape = input$pca.shape.ids
      )
  })

  output$predict.ids.coord <- renderPlot({
    predict.ids.coord()
  })

  output$conditional_pca_plot_ids <- renderUI({
    plotOutput(
      if (is.null(input$user.data)) {
        "pca.ids.coord"
      } else {
        "predict.ids.coord"
      },
      height = "800px"
    )
  })

  predict.ids.coord.2 <- reactive({
    pca.ids.coord <- ids.predict.2() %>%
      rbind(
        data.frame(pca.ids()$ind$coord) %>%
          ungroup() %>%
          select(input$dim.ids.2[1], input$dim.ids.2[2]) %>%
          magrittr::set_colnames(c("Dim.1", "Dim.2")) %>%
          cbind(training.ids()) %>%
          add_column(km = pca.ids.km()$cluster) %>%
          left_join(
            ref.seq() %>%
              select(oligo, topo, gba, salt),
            by = 'oligo'
          )
      ) %>%
      left_join(
        ref.seq() %>%
          ungroup() %>%
          select(
            oligo,
            tetrad,
            tetrad.id,
            loop,
            plus.minus,
            groove,
            conformer,
            gba.stacks
          ),
        by = "oligo"
      ) %>%
      mutate(
        topo = if_else(is.na(topo), "User", topo),
        conformer = if_else(is.na(conformer), "User", conformer),
        gba_stacks = if_else(is.na(gba.stacks), "User", gba.stacks),
        tetrad = if_else(is.na(tetrad), "User", as.character(tetrad)),
        tetrad.id = if_else(is.na(tetrad.id), "User", tetrad.id),
        loop = if_else(is.na(loop), "User", loop),
        plus.minus = if_else(is.na(plus.minus), "User", plus.minus),
        groove = if_else(is.na(groove), "User", groove),
        salt = if_else(is.na(salt), "User", groove) # added
      )

    pca.ids.coord %>%
      pca.plotR(
        .,
        dim.1 = input$dim.ids.2[1],
        dim.2 = input$dim.ids.2[2],
        color = input$pca.color.ids.2,
        shape = input$pca.shape.ids.2
      )
  })

  output$predict.ids.coord.2 <- renderPlot({
    predict.ids.coord.2()
  })

  output$conditional_pca_plot_ids_2 <- renderUI({
    plotOutput(
      if (is.null(input$user.data)) {
        "pca.ids.coord.2"
      } else {
        "predict.ids.coord.2"
      },
      height = "800px"
    )
  })

  # UV calculator------

  ## User oligos----
  user_oligos <- reactiveVal(
    data.frame(
      oligo = "Oligo 1",
      seq = sample(
        c("TTAGGG", "TTAGGC", "TTTAGGG", "TTGGG", "TTACAGG", "TGGG", "TATGGG"),
        1
      ),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        eps_260 = dt.spec.calcR(
          input.contrib = dt.contributR(input.seq = seq),
          input.wl = 260
        )
      )
  )

  ## Render editable DataTable----
  output$spectra_user_input <- renderDT({
    datatable(
      user_oligos() %>%
        group_by(oligo, seq) %>%
        mutate(
          eps_260 = dt.spec.calcR(
            input.contrib = dt.contributR(input.seq = seq),
            input.wl = 260
          )
        ),
      colnames = c(
        "Oligonucleotide",
        "Sequence",
        "&epsilon;<sub>260nm</sub> (M<sup>-1</sup>cm<sup>-1</sup>)"
      ),
      editable = TRUE,
      escape = FALSE #HTML rendering
    )
  })

  ### Track user edits and update reactiveVal----
  observeEvent(input$spectra_user_input_cell_edit, {
    info <- input$spectra_user_input_cell_edit
    new_data <- user_oligos()
    new_data[info$row, info$col] <- info$value # Use correct column index
    user_oligos(new_data) # Update reactiveVal
  })

  ### Add a new row when button is clicked
  observeEvent(input$add_row, {
    new_data <- user_oligos()
    row_number <- nrow(new_data) + 1
    new_row <- data.frame(
      oligo = paste0("Oligo ", row_number),
      seq = sample(
        c("TTAGGG", "TTAGGC", "TTTAGGG", "TTGGG", "TTACAGG", "TGGG", "TATGGG"),
        1
      ),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        eps_260 = dt.spec.calcR(
          input.contrib = dt.contributR(input.seq = seq),
          input.wl = 260
        )
      )
    user_oligos(rbind(new_data, new_row)) # Append new row
  })

  ## Perform calculations based on updated data----
  user_oligos_calc <- reactive({
    as.data.table(user_oligos()) %>%
      group_by(oligo, seq) %>%
      expand(wl = 220:310) %>%
      mutate(
        contrib = map(seq, dt.contributR),
        eps = map2_dbl(wl, contrib, dt.spec.calcR)
      )
  })

  ## Render updated plot----
  output$p_user_oligos_calc <- renderPlot({
    user_oligos_calc() %>%
      ggplot(aes(x = wl, y = eps, color = oligo)) +
      geom_vline(
        xintercept = input$calc_wavelength,
        color = current_primary(),
        linewidth = 0.75,
        linetype = 'dashed',
        inherit.aes = FALSE
      ) +
      geom_label_repel(
        data = . %>% filter(wl == input$calc_wavelength),
        aes(label = round(eps, 0), color = oligo),
        size = 6,
        fontface = "bold",
        show.legend = FALSE
      ) +
      geom_line(linewidth = 0.8) +
      labs(
        x = "&lambda; (nm)",
        y = "&epsilon; (M<sup>-1</sup>cm<sup>-1</sup>)",
        color = "Oligonucleotide"
      ) +
      custom.theme.markdown +
      scale_x_continuous(
        limits = c(220, 310),
        breaks = c(220, 240, 260, 280, 300, 310, input$calc_wavelength),
        expand = c(0, 0)
      )
  })
})

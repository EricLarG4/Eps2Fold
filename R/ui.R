library(shiny)
library(bslib)
library(bsicons)
library(DT)

shiny::bootstrapLib(bslib::bs_theme())

link_git <- tags$a(
  shiny::icon("github"),
  "GitHub",
  href = "https://github.com/EricLarG4/Eps2Fold",
  target = "_blank"
)

link_wiki <- tags$a(
  shiny::icon("book"),
  "Wiki",
  href = "https://github.com/EricLarG4/Eps2Fold/wiki",
  target = "_blank"
)

link_template_ref <- tags$a(
  shiny::icon("download"),
  href = "https://github.com/EricLarG4/G4_IDS/raw/refs/heads/main/data/Reference_data.xlsx",
  "Reference Data",
  target = "_blank",
  # style = "font-size: 14px; color: coral; text-decoration: underline;",
  download = NA # Optional, allows browsers to download instead of opening
)

link_template_user <- tags$a(
  shiny::icon("download"),
  href = "https://github.com/EricLarG4/G4_IDS/raw/refs/heads/main/data/User%20data.xlsx",
  "User Data Template",
  target = "_blank",
  # style = "font-size: 14px; color: coral; text-decoration: underline;",
  download = NA # Optional, allows browsers to download instead of opening
)

shiny::bootstrapLib() # Ensure Bootstrap resources are included

ui <- bslib::page_navbar(
  title = 'Eps2Fold',
  id = "nav",
  window_title = "Eps2Fold - a Shiny app for Eps2Fold analysis",
  theme = bslib::bs_theme(
    preset = 'zephyr',
    font_scale = 0.9
  ),
  #Ensure that all accordions are open on page load, on all pages
  tags$head(
    tags$script(HTML(
      '
    $(document).ready(function() {
      // Open accordions by default on page load
      $(".accordion-collapse").collapse("show");
    });
  '
    ))
  ),
  tags$head(tags$style(".shiny-output-error{visibility: hidden}")),
  tags$head(tags$style(
    ".shiny-output-error:after{content: 'Please wait...';visibility: visible}"
  )),
  ##sidebar------
  sidebar = sidebar(
    accordion(
      accordion_panel(
        title = "Theming",
        selectInput(
          "theme_select",
          label = "Theme:",
          choices = c(
            "default",
            "cerulean",
            "cosmo",
            "cyborg",
            "darkly",
            "flatly",
            "journal",
            "litera",
            "lumen",
            "lux",
            "materia",
            "minty",
            "morph",
            "pulse",
            "quartz",
            "sandstone",
            "simplex",
            "sketchy",
            "slate",
            "solar",
            "spacelab",
            "superhero",
            "united",
            "vapor",
            "yeti",
            "zephyr"
          ),
          selected = "zephyr",
          multiple = FALSE,
          width = "100%"
        ),
        #slider for font size
        sliderInput(
          inputId = "font_size",
          label = "Font size:",
          min = 0.1,
          max = 2,
          value = 0.9,
          step = 0.1,
          width = "100%"
        )
      ),
      conditionalPanel(
        condition = "input.nav === 'Data'",
        accordion_panel(
          id = 'tabinput',
          title = "Data input",
          icon = bsicons::bs_icon("menu-app"),
          open = TRUE,
          class = "collapse in",
          multiple = TRUE,
          input_switch(
            id = "data_source",
            label = "Online data source",
            value = TRUE
          ),
          fileInput(
            "ref.data",
            "Select reference data file",
            multiple = FALSE,
            accept = c(".xlsx", ".xls")
          ),
          fileInput(
            "user.data",
            "Select user data file",
            multiple = FALSE,
            accept = c(".xlsx", ".xls")
          ),
          link_template_ref,
          link_template_user
        ),
        accordion_panel(
          title = 'Data filters',
          icon = icon('filter'),
          open = TRUE,
          multiple = TRUE,
          uiOutput("ref.seq.topo.0_ui"),
          uiOutput("ref.seq.confomer.0_ui"),
          uiOutput("ref.seq.gba.0_ui"),
          uiOutput("ref.seq.gba.stacks.0_ui"),
          uiOutput("ref.seq.tetrad.0_ui"),
          uiOutput("ref.seq.tetrad.id.0_ui"),
          uiOutput("ref.seq.loop.0_ui"),
          uiOutput("ref.seq.plus.minus.0_ui"),
          uiOutput("ref.seq.groove.0_ui"),
          uiOutput("ref.seq.cation.0_ui"),
          uiOutput("ref.seq.oligo.0_ui"),
          uiOutput("user.seq.oligo.0_ui")
        )
      ),
      conditionalPanel(
        condition = "input.nav === 'Data plots'",
        accordion_panel(
          id = 'tabplot',
          title = 'Plot options',
          icon = icon("compass-drafting"),
          open = TRUE,
          multiple = TRUE,
          input_switch(
            id = "ids.ref.select",
            label = "Theoretical UV reference",
            value = TRUE
          ),
          sliderInput(
            inputId = "wl",
            label = "Wavelength (nm)",
            min = 220,
            max = 350,
            value = c(220, 310),
            step = 5
          ),
          selectInput(
            inputId = "ids.norm",
            label = "Eps2Fold normalization",
            choices = c("Δε", "-1/+1"),
            selected = "Δε"
          ),
          selectInput(
            inputId = "cd.norm",
            label = "CD normalization",
            choices = c("Δε", "Δε/ε", "-1/+1"),
            selected = "Δε"
          ),
          selectInput(
            inputId = "ref.panel",
            label = "Layout",
            choices = c("Panels", "Superimposed", "Mean"),
            selected = "Panels"
          ),
          selectInput(
            inputId = "ref.color",
            label = "Coloring",
            choices = c(
              "Topology",
              "Conformer",
              "GBA",
              "GBA stacks",
              "Tetrads",
              "Tetrad combination",
              "Loop progression",
              "Tetrad handedness",
              "Grooves",
              "Cation"
            ),
            selected = "Topology"
          )
        ),
        accordion_panel(
          title = "Data filters",
          icon = icon('filter'),
          open = TRUE,
          multiple = TRUE,
          uiOutput("ref.seq.topo_ui"),
          uiOutput("ref.seq.conformer_ui"),
          uiOutput("ref.seq.gba_ui"),
          uiOutput("ref.seq.gba.stacks_ui"),
          uiOutput("ref.seq.tetrad_ui"),
          uiOutput("ref.seq.tetrad.id_ui"),
          uiOutput("ref.seq.loop_ui"),
          uiOutput("ref.seq.plus.minus_ui"),
          uiOutput("ref.seq.groove_ui"),
          uiOutput("ref.seq.cation_ui"),
          uiOutput("ref.seq.oligo_ui")
        ),
        accordion_panel(
          title = 'User data',
          icon = icon('user-friends'),
          open = TRUE,
          multiple = TRUE,
          uiOutput("user.seq.oligo_ui"),
          selectInput(
            inputId = "user.panel",
            label = "Layout",
            choices = c("Panels", "Superimposed"),
            selected = "Panels"
          )
        )
      ),
      conditionalPanel(
        condition = "input.nav === 'PCA'",
        accordion_panel(
          id = "tabpca",
          class = "collapse in",
          title = "PCA",
          icon = icon("object-ungroup"),
          open = TRUE,
          multiple = TRUE,
          sliderInput(
            inputId = "ncp",
            label = "PCA dimensions",
            min = 2,
            max = 10,
            value = 3,
            step = 1
          ),
          input_switch(
            id = "scale.unit",
            label = "Scale to unit variance",
            value = FALSE
          ),
          selectInput(
            inputId = "k.mean.algo",
            label = "k-means algorithm",
            choices = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
            selected = "Hartigan-Wong"
          ),
          sliderInput(
            inputId = "cluster.center",
            label = "k-means centers",
            min = 2,
            max = 10,
            value = 4,
            step = 1
          ),
          actionButton(
            inputId = "button.cd.invest",
            label = "Investigate CD",
            icon = icon("magnifying-glass-plus"),
          ),
          actionButton(
            inputId = "button.ids.invest",
            label = "Investigate Eps2Fold",
            icon = icon("magnifying-glass-plus"),
          )
        )
      ),
      conditionalPanel(
        condition = "input.nav === 'UV spectrum calculator'",
        accordion_panel(
          title = "UV spectrum calculation",
          icon = icon('calculator'),
          open = TRUE,
          multiple = TRUE,
          actionButton("add_row", "Add oligonucleotide", icon = icon("plus")),
        ),
        accordion_panel(
          title = "Plot parameters",
          icon = icon('compass-drafting'),
          open = TRUE,
          multiple = TRUE,
          sliderInput(
            inputId = "calc_wavelength",
            label = "Wavelength (nm)",
            min = 220,
            max = 310,
            value = 260,
            step = 1
          )
        )
      )
    )
  ),
  nav_panel(
    id = "Data",
    title = "Data",
    icon = icon("table"),
    navset_card_tab(
      nav_panel(
        "Reference",
        icon = icon("table-cells-row-lock"),
        navset_card_tab(
          nav_panel(
            title = 'Oligonucleotides',
            DTOutput("seq", height = "1200px")
          ),
          nav_panel(
            title = "UV",
            DTOutput("ref.uv", height = "1200px")
          ),
          nav_panel(
            title = "Eps2Fold",
            DTOutput("ref.ids", height = "1200px")
          ),
          nav_panel(
            title = "CD",
            DTOutput("ref.cd", height = "1200px")
          ),
          nav_panel(
            title = "training Eps2Fold",
            DTOutput("training.ids", height = "1200px")
          ),
          nav_panel(
            title = "training CD",
            DTOutput("training.cd", height = "1200px")
          )
        )
      ),
      nav_panel(
        "User",
        icon = icon('user-plus'),
        navset_card_tab(
          nav_panel(
            title = "UV",
            DTOutput("user.uv.input")
          ),
          nav_panel(
            title = "Eps2Fold",
            DTOutput("user.ids.input")
          ),
          nav_panel(
            title = "CD",
            DTOutput("user.cd.input")
          ),
          nav_panel(
            title = "Eps2Fold set for PCA",
            DTOutput("user.ids")
          ),
          nav_panel(
            title = "CD set for PCA",
            DTOutput("user.cd")
          )
        )
      )
    )
  ),
  nav_panel(
    title = "Data plots",
    icon = icon("chart-line"),
    navset_card_tab(
      nav_panel(
        title = 'Reference',
        icon = icon("table-cells-row-lock"),
        navset_card_tab(
          nav_panel(
            title = "UV",
            plotOutput(
              "p.ref.uv",
              height = "1200px"
            )
          ),
          nav_panel(
            title = "Eps2Fold",
            plotOutput(
              "p.ref.ids",
              height = "1200px"
            )
          ),
          nav_panel(
            title = "CD",
            plotOutput(
              "p.ref.cd",
              height = "1200px"
            )
          )
        )
      ),
      nav_panel(
        title = "User",
        icon = icon('user-plus'),
        navset_card_tab(
          nav_panel(
            title = "UV",
            plotOutput(
              "p.user.uv",
              height = "1200px"
            )
          ),
          nav_panel(
            title = "Eps2Fold",
            plotOutput(
              "p.user.ids",
              height = "1200px"
            )
          ),
          nav_panel(
            title = "CD",
            plotOutput(
              "p.user.cd",
              height = "1200px"
            )
          )
        )
      )
    )
  ),
  nav_panel(
    title = "PCA",
    icon = icon("circle-nodes"),
    navset_card_tab(
      nav_panel(
        title = 'CD',
        icon = icon("circle-half-stroke"),
        layout_column_wrap(
          width = 1 / 2,
          # width = "600px",
          # height = "600px",
          card(
            height = "800px",
            full_screen = TRUE,
            card_header(
              "View 1",
              tooltip(
                bs_icon("info-circle"),
                "Expand the view by clicking on the bottom-right button."
              )
            ),
            layout_sidebar(
              fillable = TRUE,
              sidebar = sidebar(
                position = 'right',
                selectInput(
                  inputId = "dim.cd",
                  label = "Dimensions",
                  choices = paste0("Dim.", 1:10),
                  multiple = TRUE,
                  selected = c("Dim.1", "Dim.2")
                ),
                selectInput(
                  inputId = "pca.color.cd",
                  label = "PCA colors",
                  choices = c(
                    "k-means",
                    "Topology",
                    "Conformer",
                    "GBA",
                    "GBA stacks",
                    "Tetrads",
                    "Tetrad combination",
                    "Loop progression",
                    'Tetrads x Loops',
                    "Tetrad handedness",
                    "Grooves",
                    "Cation"
                  ),
                  selected = "Topology"
                ),
                selectInput(
                  inputId = "pca.shape.cd",
                  label = "PCA shapes",
                  choices = c(
                    "k-means",
                    "Topology",
                    "Conformer",
                    "GBA",
                    "GBA stacks",
                    "Tetrads",
                    "Tetrad combination",
                    "Loop progression",
                    'Tetrads x Loops',
                    "Tetrad handedness",
                    "Grooves",
                    "Cation"
                  ),
                  selected = "Topology"
                ),
                downloadButton("dwn.pca.cd", "Download")
              ),
              #if user input supplied, show prediction, otherwise ref. panel PCA only
              uiOutput("conditional_pca_plot_cd")
            )
          ),
          card(
            height = "800px",
            full_screen = TRUE,
            card_header(
              "View 2",
              tooltip(
                bs_icon("info-circle"),
                "Expand the view by clicking on the bottom-right button."
              )
            ),
            layout_sidebar(
              fillable = TRUE,
              sidebar = sidebar(
                position = 'right',
                selectInput(
                  inputId = "dim.cd.2",
                  label = "Dimensions",
                  choices = paste0("Dim.", 1:10),
                  multiple = TRUE,
                  selected = c("Dim.1", "Dim.2")
                ),
                selectInput(
                  inputId = "pca.color.cd.2",
                  label = "PCA colors",
                  choices = c(
                    "k-means",
                    "Topology",
                    "Conformer",
                    "GBA",
                    "GBA stacks",
                    "Tetrads",
                    "Tetrad combination",
                    "Loop progression",
                    'Tetrads x Loops',
                    "Tetrad handedness",
                    "Grooves",
                    "Cation"
                  ),
                  selected = "GBA"
                ),
                selectInput(
                  inputId = "pca.shape.cd.2",
                  label = "PCA shapes",
                  choices = c(
                    "k-means",
                    "Topology",
                    "Conformer",
                    "GBA",
                    "GBA stacks",
                    "Tetrads",
                    "Tetrad combination",
                    "Loop progression",
                    'Tetrads x Loops',
                    "Tetrad handedness",
                    "Grooves",
                    "Cation"
                  ),
                  selected = "GBA"
                ),
                downloadButton(
                  "dwn.pca.cd.2",
                  "Download"
                )
              ),
              #if user input supplied, show prediction, otherwise ref. panel PCA only
              uiOutput("conditional_pca_plot_cd_2")
            )
          ),
          card(
            height = "800px",
            full_screen = TRUE,
            card_header(
              "Diagnostics",
              tooltip(
                bs_icon("info-circle"),
                "Expand the view by clicking on the bottom-right button."
              )
            ),
            navset_card_tab(
              selected = 'Variable contributions',
              nav_panel(
                title = "Scree",
                plotOutput(
                  "pca.cd.scree",
                  height = "800px"
                )
              ),
              nav_panel(
                title = "Factor map",
                plotOutput(
                  "pca.cd.fac.map",
                  height = "800px"
                )
              ),
              nav_panel(
                title = "Correlation circle",
                plotOutput(
                  "pca.cd.var.cor",
                  height = "800px"
                )
              ),
              nav_panel(
                title = "Variable contributions",
                plotOutput(
                  "pca.cd.var.coord",
                  height = "800px"
                )
              ),
              nav_panel(
                title = "Total within sum of squares",
                plotOutput(
                  "pca.cd.twss",
                  height = "800px"
                )
              ),
              nav_panel(
                title = "Gap statistic",
                plotOutput(
                  "pca.cd.gap",
                  height = "800px"
                )
              )
            )
          ),
          card(
            height = "800px",
            full_screen = TRUE,
            card_header(
              "PCA data",
              tooltip(
                bs_icon("info-circle"),
                "Expand the view by clicking on the bottom-right button."
              )
            ),
            navset_card_tab(
              nav_panel(
                title = "Data table",
                DTOutput("pca.cd.table")
              ),
              nav_panel(
                title = "Parameters table",
                DTOutput("param.cd.table")
              ),
              nav_panel(
                title = "Predict data table",
                uiOutput("conditional_predict_cd_table")
              )
            )
          )
        )
      ),
      nav_panel(
        title = 'Eps2Fold',
        icon = icon("paper-plane"),
        layout_column_wrap(
          width = 1 / 2,
          # width = "600px",
          # height = "600px",
          card(
            height = "800px",
            full_screen = TRUE,
            card_header(
              "View 1",
              tooltip(
                bs_icon("info-circle"),
                "Expand the view by clicking on the bottom-right button."
              )
            ),
            layout_sidebar(
              fillable = TRUE,
              sidebar = sidebar(
                position = 'right',
                selectInput(
                  inputId = "dim.ids",
                  label = "Dimensions",
                  choices = paste0("Dim.", 1:10),
                  multiple = TRUE,
                  selected = c("Dim.1", "Dim.2")
                ),
                selectInput(
                  inputId = "pca.color.ids",
                  label = "PCA colors",
                  choices = c(
                    "k-means",
                    "Topology",
                    "Conformer",
                    "GBA",
                    "GBA stacks",
                    "Tetrads",
                    "Tetrad combination",
                    "Loop progression",
                    'Tetrads x Loops',
                    "Tetrad handedness",
                    "Grooves",
                    "Cation"
                  ),
                  selected = "Topology"
                ),
                selectInput(
                  inputId = "pca.shape.ids",
                  label = "PCA shapes",
                  choices = c(
                    "k-means",
                    "Topology",
                    "Conformer",
                    "GBA",
                    "GBA stacks",
                    "Tetrads",
                    "Tetrad combination",
                    "Loop progression",
                    'Tetrads x Loops',
                    "Tetrad handedness",
                    "Grooves",
                    "Cation"
                  ),
                  selected = "Topology"
                ),
                downloadButton("dwn.pca.ids", "Download")
              ),
              #if user input supplied, show prediction, otherwise ref. panel PCA only
              uiOutput("conditional_pca_plot_ids")
            )
          ),
          card(
            height = "800px",
            full_screen = TRUE,
            card_header(
              "View 2",
              tooltip(
                bs_icon("info-circle"),
                "Expand the view by clicking on the bottom-right button."
              )
            ),
            layout_sidebar(
              fillable = TRUE,
              sidebar = sidebar(
                position = 'right',
                selectInput(
                  inputId = "dim.ids.2",
                  label = "Dimensions",
                  choices = paste0("Dim.", 1:10),
                  multiple = TRUE,
                  selected = c("Dim.1", "Dim.2")
                ),
                selectInput(
                  inputId = "pca.color.ids.2",
                  label = "PCA colors",
                  choices = c(
                    "k-means",
                    "Topology",
                    "Conformer",
                    "GBA",
                    "GBA stacks",
                    "Tetrads",
                    "Tetrad combination",
                    "Loop progression",
                    'Tetrads x Loops',
                    "Tetrad handedness",
                    "Grooves",
                    "Cation"
                  ),
                  selected = "GBA"
                ),
                selectInput(
                  inputId = "pca.shape.ids.2",
                  label = "PCA shapes",
                  choices = c(
                    "k-means",
                    "Topology",
                    "Conformer",
                    "GBA",
                    "GBA stacks",
                    "Tetrads",
                    "Tetrad combination",
                    "Loop progression",
                    'Tetrads x Loops',
                    "Tetrad handedness",
                    "Grooves",
                    "Cation"
                  ),
                  selected = "GBA"
                ),
                downloadButton(
                  "dwn.pca.ids.2",
                  "Download"
                )
              ),
              #if user input supplied, show prediction, otherwise ref. panel PCA only
              uiOutput("conditional_pca_plot_ids_2")
            )
          ),
          card(
            height = "800px",
            full_screen = TRUE,
            card_header(
              "Diagnostics",
              tooltip(
                bs_icon("info-circle"),
                "Expand the view by clicking on the bottom-right button."
              )
            ),
            navset_card_tab(
              selected = 'Variable contributions',
              nav_panel(
                title = "Scree",
                plotOutput(
                  "pca.ids.scree",
                  height = "800px"
                )
              ),
              nav_panel(
                title = "Factor map",
                plotOutput(
                  "pca.ids.fac.map",
                  height = "800px"
                )
              ),
              nav_panel(
                title = "Correlation circle",
                plotOutput(
                  "pca.ids.var.cor",
                  height = "800px"
                )
              ),
              nav_panel(
                title = "Variable contributions",
                plotOutput(
                  "pca.ids.var.coord",
                  height = "800px"
                )
              ),
              nav_panel(
                title = "Total within sum of squares",
                plotOutput(
                  "pca.ids.twss",
                  height = "800px"
                )
              ),
              nav_panel(
                title = "Gap statistic",
                plotOutput(
                  "pca.ids.gap",
                  height = "800px"
                )
              )
            )
          ),
          card(
            height = "800px",
            full_screen = TRUE,
            card_header(
              "PCA data",
              tooltip(
                bs_icon("info-circle"),
                "Expand the view by clicking on the bottom-right button."
              )
            ),
            navset_card_tab(
              nav_panel(
                title = "Data table",
                DTOutput("pca.ids.table")
              ),
              nav_panel(
                title = "Parameters table",
                DTOutput("param.ids.table")
              ),
              nav_panel(
                title = "Predict data table",
                DTOutput("predict.ids.table")
              )
            )
          )
        )
      ),
    )
  ),
  nav_panel(
    "UV spectrum calculator",
    icon = icon("calculator"),
    layout_column_wrap(
      width = 1,
      card(
        card_header(
          "Oligonucleotide table",
          tooltip(
            bs_icon("info-circle"),
            "Input  oligonucleotide sequences.\nUse unique names.\nClick on the 'Add oligonucleotide' button to add more oligonucleotides."
          )
        ),
        DTOutput("spectra_user_input")
      ),
      card(
        card_header(
          "Theoretical UV spectra",
          tooltip(
            bs_icon("info-circle"),
            "Spectra expressed in molar extinction coefficients.\nThe label indicates the values at 260 nm by default.\nUse the slider to change the wavelength."
          )
        ),
        full_screen = TRUE,
        plotOutput("p_user_oligos_calc"),
      )
    )
  ),
  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(link_git),
    nav_item(link_wiki),
    nav_item(link_template_ref),
    nav_item(link_template_user)
  )
)

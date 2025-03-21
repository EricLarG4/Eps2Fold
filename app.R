# libraries----

# Define required packages
required_packages <- c(
    "tidyverse", "readxl", "data.table", "magrittr",
    
    # Shiny
    "shiny", "bslib", "bsicons",
    
    # Theming
    "thematic",
    
    # Data display
    "ggrepel", "ggthemes", "ggtext", "patchwork", "plotly", "DT",
    
    # PCA
    "FactoMineR", "cluster", "FactoInvestigate"
  )
  
  # Function to install missing packages
  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  
  # Install missing packages
  invisible(lapply(required_packages, install_if_missing))
  
  # Load all packages
  lapply(required_packages, library, character.only = TRUE)
  

## data manipulation----
library('tidyverse')
library('readxl')
library('data.table')
library('magrittr')

## shiny----
library('shiny')
library('bslib')
library('bsicons')

## Theming----
library('thematic')

## data display----
library('ggrepel')
library('ggthemes')
library('ggtext')
library('patchwork')
library('plotly')
library('DT')

## PCA----
library('FactoMineR')
library('cluster')
library('FactoInvestigate')

# sources----
source('R/functions.R')
source('R/theoretical spectra.R')

source('R/ui.R')
source('R/server.R')


# options----
# Set options to suppress warnings in the app
options(warn = -1)
# also see 'error message management' in UI

# App----

shinyApp(ui = ui, server = server)

# libraries----

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
#to solve text not rendering with shinylive. Not necesarry for connect.posit.cloud deployement
webr::install("ggtext")
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

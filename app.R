# libraries----

## data manipulation----
library('tidyverse')
library('readxl')
library('data.table')
library('magrittr')

## shiny----
library('shiny')
library('shinydashboard')
library('shinydashboardPlus')
library('shinyWidgets')

## Theming----
library('thematic')
library('fresh')
library('bsicons')
library('bslib')

## data display----
library('ggrepel')
library('ggdark')
library('ggthemes')
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
source('ui.R')
source('server.R')


# options----
# Set options to suppress warnings in the app
options(warn = -1)
# also see 'error message management' in UI

# App----

shinyApp(ui = ui, server = server)

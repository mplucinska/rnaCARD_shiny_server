library(readr)
library(ggplot2)
library(ggthemes)  
library(shiny)
library(plotly)
library(shinyWidgets)
library(shinyjs)
library(future)
library(promises)
library(shinycssloaders)
library(DT)


if (future::supportsMulticore()) {
  future::plan(future::multicore)
} else {
  future::plan(future::multisession)
}

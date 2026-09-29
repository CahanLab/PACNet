# Loaded once before ui.R and server.R (locally and under Shinylive).
suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinycssloaders)
  library(shinyWidgets)
  library(shinyjs)
  library(shinyalert)
  library(shinybusy)
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
  library(pheatmap)
  library(plotly)
  library(dplyr)
  library(bslib)
})
# randomForest's namespace must be loaded so predict() dispatches on the
# stored classifiers, but it is not attached: its margin() would mask
# ggplot2::margin() used in plotting.R. cnlib.R calls randomForest::randomForest().
loadNamespace("randomForest")
source("cnlib.R")
source("data_access.R")
source("pacnet_utils.R")
source("pipeline.R")
source("plotting.R")
source("ui_helpers.R")

options(shiny.maxRequestSize = 50*1024^2)

# plotly's JavaScript, loaded with the page (see ui.R). Under Shinylive,
# plotly widgets arriving in the same update otherwise rendered before the
# dynamically injected plotly.js had loaded ("Plotly is not defined").
plotly_page_deps <- htmltools::renderTags(
  plot_ly(x = 1, y = 1, type = "scatter", mode = "markers"))$dependencies

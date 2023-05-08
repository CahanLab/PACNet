# CellNet WebApp ui for top-pairs, rank-based agnostic analysis of expression profiles

# On AWS image:
# sudo su - -c "R -e \"install.packages('BiocManager')\""
# sudo su - -c "R -e \"BiocManager::install('AnnotationDbi')\""
# sudo su - -c "R -e \"BiocManager::install('GO.db')\""
# sudo su - -c "R -e \"BiocManager::install('org.Hs.eg.db')\""
# sudo su - -c "R -e \"devtools::install_github('pcahan1/CellNet', ref='master')\""
# sudo su - -c "R -e \"devtools::install_github('pcahan1/cancerCellNet')\""
# sudo su - -c "R -e \"install.packages('shinythemes')\""
# sudo su - -c "R -e \"install.packages('shinycssloaders')\""
# sudo su - -c "R -e \"install.packages('shinyWidgets')\""
# sudo su - -c "R -e \"install.packages('shinyjs')\""
# sudo su - -c "R -e \"install.packages('shinyalert')\""
# sudo su - -c "R -e \"install.packages('plotly')\""

suppressPackageStartupMessages({
  library(CellNet)
  library(cancerCellNet)
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
  library(igraph)
  library(plotly)
  source("plotting.R")
  source("pacnet_utils.R")
})

# ui <- fluidPage(
#    #theme = shinytheme("cosmo"),
#    #theme = shinytheme("flatly"),
#    
#    useShinyjs(),
#    
#    titlePanel(windowTitle="CellNet Web",
#               title = div("CellNet Web App", img(src="images/logo.png", height=60))),
#    hr(),
#    h2("Introduction"),
#    p("CellNet is a network-biology-based, computational platform that assesses the fidelity of
#      cellular engineering and generates hypotheses for improving cell derivations.
#      CellNet is based on the reconstruction of cell type-specific gene regulatory networks (GRNs),
#      which we performed using publicly available RNA-Seq data of 16 human cell and tissue types.
#      Below, we describe how to apply CellNet to your RNA-Seq data."),
#    p("See more at ", tags$a("https://github.com/pcahan1/CellNet")),
#    hr(),
#    
#    h2("How to use this application"),
#    h4("Inputs: "),
#    p("1) ", strong("Study name"), ": desired study name"),
#    p("2) ", strong("Sample metadata table"), "with the following columns at minimum: sample_name, description1"),
#    p("3) ", strong("non-normalized"), "counts, TPM, or FPKM expression matrix with gene", strong("symbols"), "as row names and sample names as column names.
#      Column names in counts matrix must match sample_name column in sample metadata table."),
#    br(),
#    p("Alternatively, choose 'Browse reference data' to browse available engineered reference datasets rather analyzing an input dataset."),
#    hr(),
#    
#    h4("Outputs: "),
#    p("1) ", strong("Classification heatmap"),": Columns represent query samples, and rows represent tissue types of the training data.
#      Each square is colored by the query sample's classification score for the given tissue type.",
#      strong("Scores range from 0 (distinct from the tissue type of the training data) to 1 (indistinguishable from the tissue type of the training data).")),
#    p("2) ", strong("GRN Status"),": GRN status indicates the extent to which a tissue GRN is established in the training and query samples.
#      The raw GRN status is computed as the mean z-score of all genes in a tissue GRN, weighted by their importance to the associated tissue classifier.
#      The raw GRN status is then normalized to the mean raw GRN status of the training data samples of the given tissue.
#      Error bars represent mean ± 1 s.d."),
#    p("3) ", strong("Transcription Factor Scores"),": The transcriptional regulators of the tissue GRN are shown on the y axis,
#      with the Network Influence Score (NIS) of each regulator on the x axis.
#      The NIS prioritizes transcription factors (TFs) such that their experimental perturbation is predicted to improve the target tissue classification,
#      with", strong("negative scores suggesting that upregulation of a given TF's expression is needed to achieve the target cell type, and positive scores suggesting that downregulation is needed."),
#      "The NIS of a TF is computed based on three components:
#      The first component is the extent to which the TF is dysregulated as compared with its expected value in the target tissue type.
#      The second component is the number of predicted targets of the TF. The third component is the extent to which the target genes are dysregulated."),
#    
#    #h3("Protocol"),
#    #p("Step 1: "),
#    hr(),
#    
#    sidebarLayout(
#       sidebarPanel(h4("Uploading your inputs"),
#                    selectInput(inputId="actionChoice", label="What do you want to do?",
#                                choices = c("CellNet analysis", "Browse reference data"),
#                                selected = "CellNet"),
#                    uiOutput(outputId = "appOptions"),
#                    actionButton(inputId = "submit", label= "Submit"),
#                    width=3),
#       
#       mainPanel(
#          h3("CellNet Results"),
#          hr(),
#          #shinyjs for hidden, shinyWidgets for progressBar
#          hidden(div(id="progressDiv", progressBar(id = "progress", value = 0, total = 100, status = "info",
#                                                   display_pct = TRUE, striped = TRUE, title = "Progress..."))),
#          tabsetPanel(type = "tabs",
#                      tabPanel("Classification",
#                               br(),
#                               hidden(div(id="classDiv",
#                                          shinycssloaders::withSpinner(
#                                             plotlyOutput(outputId="classHm", width="900px", height="500px"),
#                                             type=7, color = "#9A50BD"),
#                                          br(),
#                                          textOutput("engineeredDescrip"),
#                                          br()
#                               )
#                               ),
#                               hidden(div(id="engineeredRefClassDiv",
#                                          shinycssloaders::withSpinner(
#                                             plotlyOutput("engineeredRef", width="900px", height="1200px"),
#                                             type=7, color = "#9A50BD"),
#                                          br()
#                               )
#                               )
#                               
#                      ),
#                      tabPanel("GRN Status",
#                               br(),
#                               # hidden(div(id="grnStatDiv",
#                               #            shinycssloaders::withSpinner(plotOutput("GRNstatus", height="4000px"), type=7, color = "#9A50BD")))
#                               # ),
#                               hidden(div(id="grnStatDiv",
#                                          shinycssloaders::withSpinner(
#                                             plotlyOutput(outputId = "GRNstatus", height="3000px", width="900px"),
#                                             type=7, color = "#9A50BD")
#                               ))
#                      ),
#                      tabPanel("Network Scores",
#                               br(),
#                               hidden(div(id="nisDiv",
#                                          shinycssloaders::withSpinner(
#                                             plotlyOutput(outputId="NIS", height="4000px", width="900px"),
#                                             type=7, color = "#9A50BD")
#                               )
#                               ),
#                               hidden(div(id="nisDiv2",
#                                          textOutput("noNIS")
#                               )
#                               )
#                      ),
#                      hidden(div(id="validat_plots",
#                                 tabPanel("New Training Validation",
#                                          h4("Newly trained classifier performance"),
#                                          textOutput("iGenes"),
#                                          plotOutput("newTrainHm"),
#                                          br(),
#                                          plotOutput("newTrainPR"))))
#          ),
#          width=9
#       )
#    )
#    )

htmlTemplate("www/index.html",
             cellNetApp = sidebarLayout(
                            sidebarPanel(h4("Uploading your inputs"),
                                         selectInput(inputId="actionChoice", label="What do you want to do?",
                                                     choices = c("CellNet analysis", "Browse reference data"),
                                                     selected = "CellNet"),
                                         uiOutput(outputId = "appOptions"),
                                         actionButton(inputId = "submit", label= "Submit"),
                                         width=3),
                            
                            mainPanel(
                               h3("CellNet Results"),
                               hr(),
                               use_busy_spinner(spin_id = "overall_spinner", position = "top-right", margins = c(20, 20),
                                                spin = "bounce", color="#9A50BD"),
                               #shinyjs for hidden, shinyWidgets for progressBar
                               hidden(div(id="progressDiv", progressBar(id = "progress", value = 0, total = 100, status = "info",
                                                                        display_pct = TRUE, striped = FALSE, title = "Progress..."))),
                               tabsetPanel(id="analysis_tabs", type = "tabs",
                                           tabPanel(title="Classification",
                                                    value="classification",
                                                    textOutput="classification",
                                                    br(),
                                                    hidden(downloadButton(outputId="classDownload", label="Download classification results")),
                                                    hidden(div(id="classDiv",
                                                                  shinycssloaders::withSpinner(
                                                                     plotlyOutput(outputId="classHm", width="900px", height="1200px"),
                                                                     type=7, color = "#9A50BD"),
                                                                  br(),
                                                                  textOutput("engineeredDescrip"),
                                                                  br()
                                                    )),
                                                    hidden(div(id="engineeredRefClassDiv",
                                                                  shinycssloaders::withSpinner(
                                                                     plotlyOutput("engineeredRef", width="900px", height="1500px"),
                                                                     type=7, color = "#9A50BD"),
                                                                  br()
                                                    ))
                                                    
                                           ),
                                           # tabPanel(title="GRN Status",
                                           #          value="grnStatus",
                                           #          textOutput="grnStatus",
                                           #          br(),
                                           #          hidden(downloadButton(outputId="gsDownload", label="Download GRN Status Results")),
                                           #          br(),
                                           #          # hidden(div(id="grnStatDiv",
                                           #          #            shinycssloaders::withSpinner(plotOutput("GRNstatus", height="4000px"), type=7, color = "#9A50BD")))
                                           #          # ),
                                           #          hidden(div(id="grnStatDiv",
                                           #                     hidden(div(id="grnBrowseDiv", 
                                           #                                uiOutput(outputId = "browseGRN"),
                                           #                                actionButton(inputId = "submitGRN", label= "Go"),
                                           #                                shinycssloaders::withSpinner(
                                           #                                   plotlyOutput(outputId = "GRNstatusBrowse", height="4000px", width="900px"),
                                           #                                   type=7, color = "#9A50BD"
                                           #                                )
                                           #                                )),
                                           #                     hidden(div(id="GRNstatQueryDiv", 
                                           #                                shinycssloaders::withSpinner(
                                           #                                 plotlyOutput(outputId = "GRNstatus", height="3000px", width="900px"),
                                           #                                 type=7, color = "#9A50BD"
                                           #                                 )
                                           #                                ))
                                           #          ))
                                           # ),
                                           tabPanel(title="NIS",
                                                    value="networkScores",
                                                    textOutput="networkScores",
                                                    br(),
                                                    hidden(downloadButton(outputId="nisDownload", label="Download Network Influence Scores for all TFs")),
                                                    br(),
                                                    hidden(div(id="nisDiv",
                                                               shinycssloaders::withSpinner(
                                                                  plotlyOutput(outputId="NIS", height="4000px", width="900px"),
                                                                  type=7, color = "#9A50BD")
                                                    )),
                                                    hidden(div(id="nisDiv2",
                                                               textOutput("noNIS")
                                                    ))
                                           ),
                                           tabPanel(title="New Training Validation",
                                                    value="newTrainingValidation",
                                                    textOutput="newTrainingValidation",
                                                    h4("Newly trained classifier performance"),
                                                    textOutput("iGenes"),
                                                    shinycssloaders::withSpinner(
                                                       plotOutput("newTrainHm", width="auto", height="auto"),
                                                       type=7, color = "#9A50BD"),
                                                    br(),
                                                    shinycssloaders::withSpinner(
                                                       plotOutput("newTrainPR", width="auto", height="auto"),
                                                       type=7, color = "#9A50BD"))
                               ), width=9
                            )
                         )
             )
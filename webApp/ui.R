# PACNet web app UI. Packages and helpers are loaded in global.R.

page_head <- tags$head(
  # Google Analytics
  tags$script(async = NA, src = "https://www.googletagmanager.com/gtag/js?id=G-836SPJH7DH"),
  tags$script(HTML(
    "window.dataLayer = window.dataLayer || [];
     function gtag(){dataLayer.push(arguments);}
     gtag('js', new Date());
     gtag('config', 'G-836SPJH7DH');")),
  tags$link(rel = "icon", type = "image/jpeg", href = "images/PACNet_icon.jpg"),
  tags$style(HTML(
    ".doc-page table { margin: 0.5rem 0 1.25rem; border-collapse: collapse; }
     .doc-page th, .doc-page td { padding: 0.4rem 0.9rem; border-bottom: 1px solid #dee2e6;
                                  vertical-align: top; }
     .doc-page th { border-bottom-width: 2px; }"))
)

inputs_panel <- tagList(
  selectInput(inputId = "actionChoice",
              label = label_tip("What do you want to do?", TIPS$actionChoice),
              choices = c("CellNet analysis", "Browse reference data"),
              selected = "CellNet"),
  uiOutput(outputId = "appOptions"),
  bslib::tooltip(
    actionButton(inputId = "submit", label = "Submit", class = "btn-primary w-100"),
    TIPS$submit, placement = "right")
)

results_panel <- tagList(
  use_busy_spinner(spin_id = "overall_spinner", position = "top-right", margins = c(20, 20),
                   spin = "bounce", color = "#9A50BD"),
  #shinyjs for hidden, shinyWidgets for progressBar
  hidden(div(id = "progressDiv", progressBar(id = "progress", value = 0, total = 100, status = "info",
                                             display_pct = TRUE, striped = FALSE, title = "Progress..."))),
  tabsetPanel(id = "analysis_tabs", type = "tabs",
              tabPanel(title = "Classification",
                       value = "classification",
                       div(id = "introDiv", class = "text-muted py-4",
                           p("Upload your sample metadata table and expression matrix, choose the target cell/tissue type, and click Submit. Results will appear here."),
                           p("New to PACNet? See ", tags$b("Help"), " for file formats and example data you can try.")),
                       br(),
                       hidden(downloadButton(outputId = "classDownload", label = "Download classification results")),
                       hidden(div(id = "classDiv",
                                  shinycssloaders::withSpinner(
                                    plotlyOutput(outputId = "classHm", width = "900px", height = "1200px"),
                                    type = 7, color = "#9A50BD"),
                                  br(),
                                  textOutput("engineeredDescrip"),
                                  br()
                       )),
                       hidden(div(id = "engineeredRefClassDiv",
                                  shinycssloaders::withSpinner(
                                    plotlyOutput("engineeredRef", width = "900px", height = "1500px"),
                                    type = 7, color = "#9A50BD"),
                                  br()
                       ))
              ),
              tabPanel(title = "NIS",
                       value = "networkScores",
                       br(),
                       hidden(downloadButton(outputId = "nisDownload", label = "Download Network Influence Scores for all TFs")),
                       br(),
                       hidden(div(id = "nisDiv",
                                  shinycssloaders::withSpinner(
                                    plotlyOutput(outputId = "NIS", height = "4000px", width = "900px"),
                                    type = 7, color = "#9A50BD")
                       )),
                       hidden(div(id = "nisDiv2",
                                  textOutput("noNIS")
                       ))
              ),
              tabPanel(title = "New Training Validation",
                       value = "newTrainingValidation",
                       h4("Newly trained classifier performance"),
                       textOutput("iGenes"),
                       shinycssloaders::withSpinner(
                         plotOutput("newTrainHm", width = "auto", height = "auto"),
                         type = 7, color = "#9A50BD"),
                       br(),
                       shinycssloaders::withSpinner(
                         plotOutput("newTrainPR", width = "auto", height = "auto"),
                         type = 7, color = "#9A50BD"))
  )
)

doc_page <- function(md_file) {
  div(class = "container py-4 doc-page", style = "max-width: 900px;",
      shiny::markdown(readLines(md_file, encoding = "UTF-8", warn = FALSE)))
}

page_footer <- div(
  class = "text-center text-muted small py-3 border-top",
  "© 2019–2026 ",
  tags$a(href = "https://cahanlab.org/", target = "_blank", "Cahan Lab"),
  " · ",
  tags$a(href = "https://github.com/CahanLab/PACNet", target = "_blank", "GitHub"),
  " · ",
  tags$a(href = "mailto:patrick.cahan@jhmi.edu", "Contact")
)

page_navbar(
  title = tags$span(tags$img(src = "images/PACNet_icon.jpg", height = "32", class = "me-2",
                             alt = "PACNet logo"),
                    "PACNet"),
  window_title = "PACNet",
  # Pages flow normally so long Help/About text scrolls and the footer follows it.
  fillable = FALSE,
  theme = bs_theme(version = 5, bootswatch = "cosmo", primary = "#9A50BD"),
  header = tagList(page_head, useShinyjs(), plotly_page_deps),
  footer = page_footer,
  nav_panel("Analyze",
            layout_sidebar(
              sidebar = sidebar(title = "Inputs", width = 330, inputs_panel),
              results_panel)),
  nav_panel("Help", doc_page("help.md")),
  nav_panel("About", doc_page("about.md")),
  nav_spacer(),
  nav_item(tags$a(shiny::icon("flask"), "Cahan Lab", href = "https://cahanlab.org/", target = "_blank")),
  nav_item(tags$a(shiny::icon("github"), "GitHub", href = "https://github.com/CahanLab/PACNet",
                  target = "_blank"))
)

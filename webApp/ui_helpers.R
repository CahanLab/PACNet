# Small UI helpers shared by ui.R and server.R (renderUI inputs).

# Input label followed by an info icon that shows `tip` on hover or focus.
label_tip <- function(label, tip) {
  tags$span(
    label,
    bslib::tooltip(
      tags$span(shiny::icon("circle-info"), class = "text-muted ms-1", tabindex = "0"),
      tip, placement = "right"))
}

TIPS <- list(
  actionChoice = "CellNet analysis scores your own samples. Browse reference data shows the published engineered-cell panels without uploading anything.",
  sampTab = "CSV with one row per sample. Required columns: sample_name (must match the expression matrix column names) and description1 (group label shown in the plots). Other columns are ignored.",
  expMat = "CSV of counts, TPM or FPKM. First column: gene symbols (e.g. ALB); then one column per sample, named as in sample_name. Genes with missing values are dropped.",
  species = "The web app supports human data only.",
  tissueType = "The cell or tissue type your protocol aims to make. Heart, HSPC, intestine/colon, liver, lung, neuron and skeletal muscle also show published engineered reference samples. Choose none_of_the_above if yours is not listed (network scores are then skipped).",
  tissueTypeRef = "The engineered reference panel to display.",
  submit = "Runs the analysis in your browser; your data never leaves your computer. Usually under a minute. If your genes do not cover the classifier, it is retrained first, which takes several minutes."
)

# CellNet WebApp server for top-pairs, rank-based agnostic analysis of expression profiles
# Packages and helper sources are loaded in global.R.

server <- function(input, output, session) {
  
  ########################################################################################
  ##### FUNCTION DEFINITIONS #####
  
  show_validation <- function(res) {
    updateProgressBar(session = session, id = "progress",
                      title = "Creating validation plots...",
                      value = 35, total = 100)

    output$iGenes <- renderText({
      paste("Number of intersecting genes: ", length(res$iGenes))
    })
    outputOptions(output, "iGenes", suspendWhenHidden = FALSE)

    # Classification Validation Heatmap
    Sys.setlocale("LC_COLLATE", "C") # So that sort() is case sensitive
    output$newTrainHm <- renderPlot({
      acn_queryClassHm(res$classMatrix, grps=res$grps, fontsize_row=10)
    }, height=400, width=900)
    outputOptions(output, "newTrainHm", suspendWhenHidden = FALSE)

    # Classifier Assessment
    output$newTrainPR <- renderPlot({
      plot_class_PRs(res$assessment)
    }, height=720, width=720)
    outputOptions(output, "newTrainPR", suspendWhenHidden = FALSE)
  }
  
  queryClassifier <- function(classMatrixQuery, querySampTab, groupBy="description1") {
    # Adjust labels
    grp_names <- c(as.character(querySampTab[,groupBy]), rep("random", 3))
    names(grp_names) <- c(as.character(querySampTab[,"sample_name"]), "rand_1", "rand_2", "rand_3")
    #names(grp_names) <- c(as.character(querySampTab$sample_name), "rand_1", "rand_2", "rand_3")
    # Re-order classMatrixQuery to match order of rows in querySampTab
    classMatrixQuery <- classMatrixQuery[,names(grp_names)]
    
    # acn_queryClassHm(classMatrixQuery, main = paste0("Classification Heatmap, ", studyName),
    #                  grps = grp_names, is.Big = TRUE,
    #                  fontsize_row=9, fontsize_col = 10)
    heatmapPlotly(classMatrixQuery, querySampTab)
  }
  
  queryClassifierRef <- function(broadReturn, querySampTab, queryExpDat, groupBy="citation") {
    cnProc_broad <- broadReturn$cnProc
    
    # Query the classifier:
    classMatrixQuery <- broadClass_predict(cnProc = cnProc_broad, expDat = queryExpDat, nrand = 5)
    
    # Adjust labels
    grp_names <- c(as.character(querySampTab[,groupBy]), rep("random", 5))
    names(grp_names) <- c(as.character(querySampTab[,"sample_name"]), 
                          "rand_1", "rand_2", "rand_3", "rand_4", "rand_5")
    
    # Re-order classMatrixQuery to match order of rows in querySampTab
    classMatrixQuery <- classMatrixQuery[,names(grp_names)]
    
    # acn_queryClassHm(classMatrixQuery, main = paste0("Classification Heatmap, ", studyName),
    #                  grps = grp_names, is.Big = TRUE,
    #                  fontsize_row=9, fontsize_col = 10)
    heatmapPlotlyRef(classMatrixQuery, querySampTab)
  }
  
  # queryGRN <- function(GRN_statusQuery, querySampTab, ncols=2) {
  #    cell_types <- rownames(GRN_statusQuery)
  #    GRN_statusQuery <- GRN_statusQuery[,querySampTab[,"sample_name"]]
  #    
  #    plot_list <- list()
  #    i <- 1
  #    for(type in cell_types) {
  #       plot_df <-  data.frame(
  #          #"SampleNames" = paste(colnames(GRN_statusQuery), querySampTab[,"description1"]),
  #          "SampleNames" = colnames(GRN_statusQuery),
  #          "GRN_Status" = as.vector(GRN_statusQuery[type, ]),
  #          "Description" = querySampTab[,"description1"]
  #       )
  #       plot_df$SampleNames <- factor(plot_df$SampleNames, levels=plot_df$SampleNames)
  #       type_plot <- ggplot(plot_df, aes(text=Description)) + geom_bar(stat="identity", data = plot_df, 
  #                                                                      aes(x=SampleNames, y=GRN_Status), width = 0.7) +
  #          #ggtitle(paste0(type, " Network GRN Status")) + 
  #          xlab("Samples") + ylab("GRN Status") + theme_bw() +
  #          theme(text = element_text(size=10), 
  #                legend.position="none",
  #                axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  #          geom_hline(yintercept=1, linetype="dashed", color = "steelblue")
  #       
  #       plot_list[[i]] <- ggplotly(type_plot, tooltip="all") %>% layout(annotations = list(
  #          text = paste0(type, " Network GRN Status"),
  #          xref = "paper",
  #          yref = "paper",
  #          yanchor = "bottom",
  #          xanchor = "left",
  #          x = 0,
  #          y = 1,
  #          showarrow = FALSE
  #       ),
  #       xaxis = list(automargin=TRUE, tickangle=75),
  #       margin=list(t=50, b=30)
  #       )
  #       i <- i+1
  #    }
  #    
  #    output$GRNstatus <- renderPlotly({
  #       subplot(plot_list, nrows=ceiling(length(plot_list)/ncols), 
  #               titleX=TRUE, titleY=TRUE,
  #               margin=c(0.05,0.05,0.02,
  #                        floor(mean(nchar(colnames(GRN_statusQuery))))/(20*length(plot_list))
  #               )
  #       )
  #    })
  #    outputOptions(output, "GRNstatus", suspendWhenHidden = FALSE) # Essential for making tabs work in htmlTemplate!
  # }
  
  queryNIS <- function(TF_scores, querySampTab, tissueType, tissueTypeTopTFs) {
    updateProgressBar(session = session, id = "progress", title="Finishing up...",
                      value = 97, total = 100)
    
    # sample_names <- querySampTab$sra_id
    sample_names <- querySampTab[,"sample_name"]
    
    iTFs <- intersect(tissueTypeTopTFs, rownames(TF_scores))
    print(which(iTFs %in% rownames(TF_scores)) == 1:length(iTFs))
    
    plot_list <- list()
    title_list <- list()
    i <- 1
    for(sample in sample_names) {
      print(sample)
      descript <- querySampTab$description1[which(querySampTab$sample_name == sample)]
      #descript <- querySampTab$description3[which(querySampTab$sample_name == sample)]
      plot_df <- data.frame("TFs" = iTFs,
                            "Scores" = as.vector(TF_scores[iTFs,sample]))
      sample_TFplot <- ggplot(plot_df, aes(x = reorder(TFs,Scores,mean) , y = Scores)) + geom_bar(stat="identity") + #aes(fill = medVal)) +
        theme_bw() + #scale_fill_gradient2(low = "purple", 
        # mid = "white", 
        # high = "orange") + 
        #ggtitle(paste0(sample, ", ", descript, ", ", tissueType, " transcription factor scores")) +
        ylab("Network influence score") + xlab("Transcriptional regulator") + 
        theme(legend.position = "none", axis.text = element_text(size = 8)) +
        theme(text = element_text(size=10), 
              legend.position="none",
              axis.text.x = element_text(angle = 45, vjust=0.5),
              axis.title.y=element_blank())
      #coord_flip()
      plot_list[[i]] <- ggplotly(sample_TFplot) %>% 
        layout(annotations = list(
          text = paste(strwrap(
            paste0(sample, ", ", descript, ", ", tissueType, " TF scores"), width=45), 
            collapse="\n"),
          xref = "paper",
          yref = "paper",
          yanchor = "bottom",
          xanchor = "left",
          x = 0,
          y = 1,
          showarrow = FALSE
        ),
        margin=list(t=100, b=50))
      i <- i+1
    }
    
    output$NIS <- renderPlotly({
      subplot(plot_list, nrows=ceiling(length(plot_list)/2), 
              titleX=TRUE, titleY=TRUE,
              margin=c(0.05,0.05,0.02,0.02)
      )
      
    })
    outputOptions(output, "NIS", suspendWhenHidden = FALSE)
  }
  
  # queryGRNBrowse <- function(GRN_statusQuery, querySampTab, targetTissueType, ncols=2) {
  #    study_names <- unique(querySampTab$citation)
  #    plot_list <- list()
  #    i <- 1
  #    for(study in study_names) {
  #       temp_df <- querySampTab[which(querySampTab$citation == study),]
  #       plot_df <-  data.frame(#"SampleNames" = paste(colnames(GRN_statusQuery), querySampTab[,"description1"]),
  #          "SampleNames" = temp_df$sample_name,
  #          "GRN_Status" = as.vector(GRN_statusQuery[targetTissueType, temp_df$sample_name]),
  #          "StudyName" = rep(study, nrow(temp_df)),
  #          "Description" = temp_df[,"description1"]
  #       )
  #       plot_df$SampleNames <- factor(plot_df$SampleNames, levels=unique(plot_df$SampleNames))
  #       study_plot <- ggplot(plot_df, aes(text=Description)) + 
  #          geom_bar(stat="identity", data = plot_df, 
  #                   aes(x=SampleNames, y=GRN_Status), width=0.7) +
  #          xlab("Samples") + ylab("GRN Status") + theme_bw() +
  #          theme(text = element_text(size=10), 
  #                legend.position="none",
  #                axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  #          geom_hline(yintercept=1, linetype="dashed", color = "steelblue")
  #       
  #       plot_list[[i]] <- ggplotly(study_plot, tooltip="all") %>% layout(annotations = list(
  #          text = study,
  #          xref = "paper",
  #          yref = "paper",
  #          yanchor = "bottom",
  #          xanchor = "left",
  #          x = 0,
  #          y = 1,
  #          showarrow = FALSE
  #       ),
  #       xaxis = list(automargin=TRUE, tickangle=75),
  #       margin=list(t=50, b=30)
  #       )
  #       i <- i+1
  #    }
  #    
  #    output$GRNstatusBrowse <- renderPlotly({
  #       subplot(plot_list, nrows=ceiling(length(plot_list)/ncols), titleX=TRUE, titleY=TRUE,
  #               #margin=c(0.05,0.05,0,floor(mean(nchar(colnames(GRN_statusQuery))))/180)
  #               margin=c(0.05,0.05,0.01, 1/(length(plot_list) + 10) )
  #       )
  #    })
  #    outputOptions(output, "GRNstatusBrowse", suspendWhenHidden = FALSE) # Essential for making tabs work in htmlTemplate!
  # }
  
  
  
  ##### END OF FUNCTION DEFINITIONS #####
  ########################################################################################
  ########################################################################################
  
  ########################################################################################
  ########################################################################################
  
  
  ########################################################################################
  ########################################################################################
  
  
  ########################################################################################
  ########################################################################################
  
  ########################################################################################
  ########################################################################################
  ##### PROCESS INPUTS #####
  
  actionChoice <- eventReactive(input$actionChoice, {
    return(input$actionChoice)
  })
  
  
  output$appOptions <- renderUI({
    if (actionChoice() == "Browse reference data") {
      tagList(
        selectInput(inputId="species", label="Species", choices = c("Human")),#, "Mouse")),
        uiOutput("tissueTypeRef")
      )
    } else {
      tagList(
        fileInput(inputId = "sampTabUploadCSV", label= "Upload .csv sample metadata table", 
                  accept = ".csv", multiple=FALSE),
        fileInput(inputId = "expMatUpload", label = "Upload .csv expression matrix", 
                  accept = ".csv", multiple=FALSE),
        selectInput(inputId="species", label="Species", choices = c("Human")), #, "Mouse")),
        uiOutput("tissueType")
      )
    }
  })
  
  # Read in file input for sample metadata table
  querySampTab <- eventReactive(input$sampTabUploadCSV, { 
    st <- tryCatch(
      {
        read.csv(input$sampTabUploadCSV$datapath, stringsAsFactors = FALSE, header=TRUE)
      },
      error=function(condition) {
        shinyalert("Oops!", paste0("Something went wrong with your sample metadata table upload: ",
                                   condition), type="error")
        return(NULL)
      }
    )
    rownames(st) <- st$sample_name
    return(st)
  })
  
  # Read in file input for expression matrix
  queryExpDat <- eventReactive(input$expMatUpload, {
    tryCatch(
      read_expression_csv(input$expMatUpload$datapath),
      error=function(condition) {
        shinyalert("Oops!", paste0("Something went wrong with your expression matrix upload: ",
                                   condition), type="error")
        return(NULL)
      })
  })

  species <- eventReactive(input$species, {
    input$species
  })
  
  output$tissueType <- renderUI({
    if (species() == "Human") {
      selectInput(inputId = "tissueType", label = "Target Cell/Tissue Type",
                  choices = c("heart","hspc","intestine_colon",
                              "liver","lung","neuron","skeletal_muscle",
                              "esc","endothelial_cell","fibroblast","kidney",
                              "b_cell","t_cell","monocyte_macrophage",
                              "none_of_the_above"
                  ))
      
    } 
    # else if (species() == "Mouse") {
    #    selectInput(inputId = "tissueType", label = "Target Cell/Tissue Type",
    #                choices = c("t_cell","lung","nk_cell","hspc","heart",
    #                            "macrophage","skeletal_muscle","neuron","wat",
    #                            "kidney","intestine_colon","dendritic_cell",
    #                            "b_cell","fibroblast","esc","liver"))
    # }   
  })
  
  tissueType <- eventReactive(input$tissueType, { 
    return(input$tissueType)
  })
  
  output$tissueTypeRef <- renderUI({
    if (species() == "Human") {
      selectInput(inputId = "tissueTypeRef", label = "Target Cell/Tissue Type",
                  choices = c("heart","hspc","intestine_colon",
                              "liver","lung","neuron","skeletal_muscle"
                  ))
      
    } 
    # else if (species() == "Mouse") {
    #    selectInput(inputId = "tissueType", label = "Target Cell/Tissue Type",
    #                choices = c("t_cell","lung","nk_cell","hspc","heart",
    #                            "macrophage","skeletal_muscle","neuron","wat",
    #                            "kidney","intestine_colon","dendritic_cell",
    #                            "b_cell","fibroblast","esc","liver"))
    # }   
  })
  
  tissueTypeRef <- eventReactive(input$tissueTypeRef, { 
    return(input$tissueTypeRef)
  })
  
  # output$browseGRN <- renderUI({
  #    selectInput(inputId="grnTissueType", label="Tissue Type", 
  #                choices = c("esc","hspc","endothelial_cell","fibroblast",
  #                            "lung","intestine_colon","kidney",
  #                            "heart","liver","neuron","skeletal_muscle",
  #                            "b_cell","t_cell","monocyte_macrophage"),
  #                selected = "heart"
  #    )
  # })
  # outputOptions(output, "browseGRN", suspendWhenHidden = FALSE)
  # 
  # grnTissueType <- eventReactive(input$grnTissueType, {
  #    return(input$grnTissueType)
  # })
  
  ########################################################################################
  ########################################################################################
  
  ########################################################################################
  ########################################################################################
  
  ########################################################################################
  ########################################################################################
  
  ########################################################################################
  ########################################################################################
  
  ########################################################################################
  ########################################################################################
  
  ########################################################################################
  ########################################################################################
  ####### BEGIN ANALYSIS #######
  
  # Upon clicking submit:
  observeEvent(input$submit, {
    hideTab(inputId="analysis_tabs", target="newTrainingValidation")
    shinyjs::show("progressDiv") #shinyjs
    shinyjs::hide("classDownload")
    # shinyjs::hide("gsDownload")
    shinyjs::hide("nisDownload")
    
    show_spinner(spin_id = "overall_spinner", session=session)
    
    ###################################
    # For CellNet analysis of user data
    
    if (actionChoice() == "CellNet analysis") {
      if (is.null(querySampTab()) || is.null(queryExpDat())) {
        shinyalert("Oops!", "Please upload proper sample metadata table and expression matrix.", type = "error")
        return(NULL)
      }
      
      if (!("sample_name" %in% colnames(querySampTab()))) {
        shinyalert("Oops!", "Metadata table missing sample_name column. 
                       Please make sure columns are named correctly.", type = "error")
        return(NULL)
      }
      
      if (!("description1" %in% colnames(querySampTab()))) {
        shinyalert("Oops!", "Metadata table missing description1 column. 
                       Please make sure columns are named correctly.", type = "error")
        return(NULL)
      }
      
      if (!all(colnames(queryExpDat()) %in% querySampTab()[,"sample_name"])) {
        shinyalert("Oops!", "Column names of expression matrix DO NOT match sample_name column of metadata table. 
                       Please re-submit correctly formatted tables and try again.", type = "error")
        return(NULL)
      }
      
      # shinyjs::hide("grnBrowseDiv")
      shinyjs::hide("nisDiv2")
      #shinyjs::hide("validat_plots")
      hideTab(inputId="analysis_tabs", target="newTrainingValidation")
      
      queryExpDat_new <- normalize_query(queryExpDat(), querySampTab())

      updateProgressBar(session = session, id = "progress", title = "Loaded in query data...",
                        value = 10, total = 100)

      # Load TISSUE-SPECIFIC classifier and iGenes (other tissues use none_of_the_above)
      clf <- load_classifier(tissueType())
      broadReturn <- clf$broadReturn
      iGenes <- clf$iGenes

      # Check if retraining is necessary
      if (!needs_retrain(broadReturn, queryExpDat_new)) {
        updateProgressBar(session = session, id = "progress",
                          title = "No need to retrain! Starting classification...",
                          value = 35, total = 100)

        output$iGenes <- renderText({
          "Retraining not necessary for this study!"
        })
        outputOptions(output, "iGenes", suspendWhenHidden = FALSE)
      } else {
        updateProgressBar(session = session, id = "progress",
                          title = "Retraining classifier based on intersecting genes...This may take several minutes",
                          value = 20, total = 100)
        iGenes <- intersect(iGenes, rownames(queryExpDat_new))

        if (length(iGenes) == 0) {
          shinyalert("Error!", "Your expression matrix does not share any genes with the training data.
                          Please make sure the row names of your expression matrix are GENE SYMBOLS.", type="error")
        } else if (length(iGenes) < 1000) {
          shinyalert("Warning!", "Your expression matrix shares fewer than 1000 genes with the training data.
                          The most likely cause is that the row names of your expression matrix are not GENE SYMBOLS.
                          If you have made sure that your row names are properly gene symbols and would like to proceed,
                          please note that a classifier trained on fewer than 1000 genes may not be reliable.
                          Check out the New Training Validation tab for information on the classifier's performance.",
                     type="warning",
                     confirmButtonText="I understand")
        }
        showTab(inputId="analysis_tabs", target="newTrainingValidation")
        retrained <- tryCatch(
          retrain_classifier(queryExpDat_new, iGenes),
          error=function(condition) {
            shinyalert("Error", paste0("Something went wrong with retraining: ",
                                       condition), type="error")
            return(NULL)
          }
        )
        if (is.null(retrained)) {
          hide_spinner(spin_id = "overall_spinner", session=session)
          return(NULL)
        }
        broadReturn <- retrained$broadReturn
        show_validation(retrained)
        updateProgressBar(session = session, id = "progress",
                          title = "Done retraining! Starting classification...",
                          value = 40, total = 100)
      }

      # Load tissueType-specific engineered reference data
      if (tissueType() %in% PANEL_TISSUES) {
        ref <- load_engineered_ref(tissueType())
        tissueTypeRefDat <- ref$expDat
        tissueTypeRefSt <- ref$sampTab
        tissueTypeTopTFs <- ref$topTFs
      }

      classMatrixQuery <- classify(broadReturn, queryExpDat_new, nrand = 3)
      
      shinyjs::show("classDiv")
      
      ### Plot classification heatmap for user-inputted study
      output$classHm <- renderPlotly({
        tryCatch(
          {
            queryClassifier(classMatrixQuery, querySampTab(), groupBy="description1")
          },
          error=function(condition) {
            shinyalert("Error", paste0("Something went wrong with classification of your samples: ",
                                       condition), type="error")
            return(NULL)
          }
        )
        
      })
      outputOptions(output, "classHm", suspendWhenHidden = FALSE)
      
      
      ### Plot classification heatmap for engineered reference studies
      if(tissueType() %in% PANEL_TISSUES) {
        shinyjs::show("engineeredRefClassDiv")
        
        output$engineeredDescrip <- renderText({
          "How does your protocol compare to publicly available data from related studies? \n \n"   
        })
        
        output$engineeredRef <- renderPlotly({
          queryClassifierRef(broadReturn, tissueTypeRefSt, tissueTypeRefDat, groupBy="citation")
        })
        outputOptions(output, "engineeredRef", suspendWhenHidden = FALSE)
        
      } else {
        shinyjs::hide("engineeredRefClassDiv")
        # shinyjs::hide("nisDownload")
        
        output$engineeredDescrip <- renderText({
          "Comparison to engineered reference samples unavailable for selected cell/tissue type."   
        })
        outputOptions(output, "engineeredDescrip", suspendWhenHidden = FALSE)
      }
      
      output$classDownload <- downloadHandler(
        filename = function() {
          paste('classification_results_', Sys.Date(), '.csv', sep='')
        },
        content = function(con) {
          write.csv(classMatrixQuery, con, quote=FALSE)
        }
      )
      
      shinyjs::show("classDownload")
      outputOptions(output, "classDownload", suspendWhenHidden = FALSE)
      
      # updateProgressBar(session = session, id = "progress", title="Now starting GRN analysis...",
      #                   value = 50, total = 100)
      
      
      #####
      # NIS

      updateProgressBar(session = session, id = "progress", title="Now analyzing transcriptional regulators...may take several minutes",
                        value = 85, total = 100)

      if(tissueType() != "none_of_the_above") {

        TF_scores <- tryCatch(
          run_nis(queryExpDat_new, querySampTab(), iGenes, tissueType()),
          error=function(condition) {
            shinyalert("Error", paste0("Something went wrong with NIS analysis of your samples: ",
                                       condition), type="error")
            return(NULL)
          })

        output$nisDownload <- downloadHandler(
          filename = function() {
            paste('network_results_', Sys.Date(), '.csv', sep='')
          },
          content = function(con) {
            write.csv(TF_scores, con, quote=FALSE)
          }
        )
        
        if (tissueType() %in% OTHER_TISSUES) {
          keep_rows <- order(rowSums(abs(TF_scores)), decreasing = TRUE)[1:25]
          tissueTypeTopTFs <- rownames(TF_scores)[keep_rows]
        }
        
        shinyjs::hide("nisDiv2")
        shinyjs::show("nisDiv")
        tryCatch(
          {
            queryNIS(TF_scores, querySampTab(), tissueType(), tissueTypeTopTFs)   
          },
          error=function(condition) {
            shinyalert("Error", paste0("Something went wrong with NIS plotting of your samples: ",
                                       condition), type="error")
            return(NULL)
          }
        )
        
        # shinyjs::show("gsDownload")
        # outputOptions(output, "gsDownload", suspendWhenHidden = FALSE)
        
        shinyjs::show("nisDownload")
        outputOptions(output, "nisDownload", suspendWhenHidden = FALSE)
        
      } else {
        shinyjs::hide("nisDiv")
        shinyjs::hide("nisDownload")
        shinyjs::show("nisDiv2")
        output$noNIS <- renderText({
          "Network scores not available for cell/tissue types when \"none of the above\" is chosen."
        })
        outputOptions(output, "noNIS", suspendWhenHidden = FALSE)
      }
      
      
      graphics.off() # Essential for pheatmap to render for some reason...
      
      hide_spinner(spin_id = "overall_spinner", session=session)
      
      delay(ms=5000,
            updateProgressBar(session = session, id = "progress", title="Analysis Complete!",
                              value = 100, total = 100))
      
    } else if (actionChoice() == "Browse reference data") { 
      
      ###################################
      # For just browsing reference data
      shinyjs::hide("classDownload")
      # shinyjs::hide("gsDownload")
      shinyjs::hide("nisDownload")
      shinyjs::hide("classDiv")
      shinyjs::hide("nisDiv")
      #shinyjs::hide("validat_plots")
      hideTab(inputId="analysis_tabs", target="newTrainingValidation")
      
      # Load in TISSUE-SPECIFIC broadReturn and iGenes
      ##### Fix naming scheme for classifiers and iGenes
      broadReturn <- load_classifier(tissueTypeRef())$broadReturn
      ref <- load_engineered_ref(tissueTypeRef())
      tissueTypeRefDat <- ref$expDat
      tissueTypeRefSt <- ref$sampTab
      
      
      updateProgressBar(session = session, id = "progress", title = "Loaded in reference data...",
                        value = 15, total = 100)
      
      shinyjs::show("engineeredRefClassDiv")
      
      output$engineeredDescrip <- renderText({
        paste0("Reference datasets for engineered ", tissueTypeRef(), " studies")   
      })
      outputOptions(output, "engineeredDescrip", suspendWhenHidden = FALSE)
      
      output$engineeredRef <- renderPlotly({
        queryClassifierRef(broadReturn, tissueTypeRefSt, tissueTypeRefDat, groupBy="citation")
      })
      outputOptions(output, "engineeredRef", suspendWhenHidden = FALSE)
      
      updateProgressBar(session = session, id = "progress", title="Finished classification!", #, now starting GRN analysis...",
                        value = 30, total = 100)
      
      ###### GRN Analysis ######
      
      
      ### GRN Status
      # shinyjs::show("grnStatDiv")
      # shinyjs::hide("GRNstatQueryDiv")
      # shinyjs::show("grnBrowseDiv")
      # 
      # updateProgressBar(session = session, id = "progress", title="Now analyzing GRN status...",
      #                   value = 70, total = 100)
      # 
      # 
      # observeEvent(input$submitGRN, {
      #    queryGRNBrowse(tissueTypeRefGRNStat, tissueTypeRefSt, grnTissueType(), ncols=2)
      # })
      
      ### NIS
      shinyjs::show("nisDiv2")
      output$noNIS <- renderText({
        "Network scores unavailable in browse mode."
      })
      outputOptions(output, "noNIS", suspendWhenHidden = FALSE)
      
      hide_spinner(spin_id = "overall_spinner", session=session)
      
      delay(ms=5000,
            updateProgressBar(session = session, id = "progress", title="Analysis Complete!",
                              value = 100, total = 100))
      
    } #else {
    #         
    #      }
    
    
  })
  
}
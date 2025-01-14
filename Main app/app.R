# Load required packages -----------------------------------------------------
library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(plotly)
library(tidyverse)
library(DT)
library(colourpicker) 
library(shinyWidgets)
library(fgsea)
library(gplots)
library(RColorBrewer)
library(rio)
library(ggbeeswarm)
library(GSEABase)
library(org.Hs.eg.db)
library(biomaRt)
library(stringr)

#-------------------------------------
# UI Section ---------------------------------------------------------------

ui <- dashboardPage(
  
  # Header -----------------------------------------------------------------
  dashboardHeader(title = "BF591finalproject"),
  
  # Sidebar ----------------------------------------------------------------
  dashboardSidebar(
    sidebarMenu(
      menuItem("Samples", icon = icon("vial"), tabName = "samples"),
      menuItem("Counts", icon = icon("chart-bar"), startExpanded = TRUE,
               menuSubItem("Summary", tabName = "counts_summary_tab", icon = icon("info")),
               menuSubItem("Diagnostic Plots", tabName = "counts_diag_tab", icon = icon("stethoscope")),
               menuSubItem("Clustered Heatmap", tabName = "counts_heatmap_tab", icon = icon("fire")),
               menuSubItem("PCA", tabName = "counts_pca_tab", icon = icon("project-diagram"))
      ),
      menuItem("DE", icon = icon("fire"), startExpanded = TRUE,
               menuSubItem("Tables", tabName = "de_tables_tab", icon = icon("table")),
               menuSubItem("Volcano Plot", tabName = "de_volcano_tab", icon = icon("volcano")),
               menuSubItem("Plot Table", tabName = "de_plot_table_tab", icon = icon("list"))
      ),
      menuItem("GSEA", icon = icon("dna"), startExpanded = TRUE,
               menuSubItem("Top Results", tabName = "gsea_top_tab", icon = icon("bar-chart")),
               menuSubItem("Table", tabName = "gsea_table_tab", icon = icon("table")),
               menuSubItem("Plots", tabName = "gsea_plots_tab", icon = icon("chart-line"))
      )
    )
  ),
  
  # Body -------------------------------------------------------------------
  dashboardBody(
    tabItems(
      
      # ========== Samples Tab =========================================================
      tabItem(tabName = "samples",
              fluidPage(
                titlePanel("BF591: RShiny App - Samples"),
                h2("--Xinyu Li"),
                sidebarLayout(
                  sidebarPanel(
                    fileInput(inputId ="file_input",  label= "Upload metadata file (.csv)"), 
                    actionButton("submit_button", "Submit")
                  ),
                  mainPanel(
                    tabsetPanel(
                      
                      # Summary of sample data
                      tabPanel("Summary",
                               hr(),
                               # Display basic statistical information of each column
                               DT::dataTableOutput("summary_output")
                      ),
                      
                      # Raw sample data table
                      tabPanel("Table",
                               hr(),
                               # Show the original sample data
                               dataTableOutput("table_output")
                      ),
                      
                      # Diagnostic plot for samples (e.g., violin plots)
                      tabPanel("Diagnostic Plot",
                               hr(),
                               # Use plotly for interactive plots
                               plotlyOutput("plots_output")
                      ),
                      
                      # Distribution plot (histogram or density)
                      tabPanel("Distribution Plot",
                               hr(),
                               # User selects column and plot type
                               selectInput("sample_plot_var", "Select Variable to Plot:", choices = NULL),
                               radioButtons("sample_plot_type", "Plot Type:", c("Histogram" = "hist", "Density" = "density")),
                               plotOutput("sample_dist_plot")
                      )
                    )
                  )
                )
              )
      ),
      
      # ========== Counts Tabs =========================================================
      tabItem(tabName = "counts_summary_tab",
              fluidRow(
                box(width=3, title="Counts Input & Filters", status="primary", solidHeader=TRUE,
                    
                    # Upload normalized counts file
                    fileInput(inputId = "counts",
                              label = "Upload Normalized Counts File (.csv):", 
                              multiple = FALSE, 
                              accept = ".csv",
                              buttonLabel = "Browse",
                              placeholder = "norm_counts.csv"
                    ),
                    
                    # Set variance threshold percentile
                    numericInput(inputId = "count_percentile_slider", 
                                 label = "Genes With at Least `X` Percentile of Variance:", 
                                 min = 1, 
                                 max = 100, 
                                 value = 80
                    ),
                    
                    # Set non-zero sample count threshold
                    numericInput(inputId = "count_nonzero_slider", 
                                 label = "Genes With at Least `X` Samples That Are Non-Zero:", 
                                 min = 1, 
                                 max = 30, 
                                 value = 5
                    ),
                    
                    # Reset filter values
                    actionButton(inputId = "reset_counts_main", 
                                 label = "Reset Values",
                                 icon = icon("undo"))
                ),
                box(width=9, status="primary", solidHeader=TRUE,
                    h3("Get overview of data"),
                    h5("wait a few seconds"),
                    # Display a summary table after filtering counts
                    DTOutput("counts_summary_table")
                )
              )
      ),
      
      tabItem(tabName = "counts_diag_tab",
              fluidRow(
                box(width=12, status="primary", solidHeader=TRUE,
                    h5("Please wait a few seconds for the plots to load."),
                    h5("Plot of median count vs variance"),
                    # Plot median counts vs variance
                    plotOutput("plot_count_vs_variance"),
                    h5("Plot of median count vs number of zeros"),
                    # Plot median counts vs number of zeros
                    plotOutput("plot_count_vs_zeros")
                )
              )
      ),
      
      tabItem(tabName = "counts_heatmap_tab",
              fluidRow(
                box(width=2, status="primary", solidHeader=TRUE,
                    # Whether to log-transform the data for the heatmap
                    radioButtons(inputId = "count_heatmap_transformation", 
                                 label = "Log Transform Data?", 
                                 choices = c("TRUE", "FALSE"), 
                                 selected = "TRUE"
                    ),
                    # Button to display heatmap
                    actionButton(inputId = "count_heatmap_button", 
                                 label = "Display Heatmap",
                                 icon = icon("fire"))
                ),
                box(width=10, status="primary", solidHeader=TRUE,
                    tags$h5("Please wait a few seconds for the plots to load."),
                    # Display hierarchical clustering heatmap of gene expression
                    plotOutput("count_heatmap")
                )
              )
      ),
      
      tabItem(tabName = "counts_pca_tab",
              fluidRow(
                box(width=3, status="primary", solidHeader=TRUE,
                    # Choose PCs for PCA plot axes
                    selectInput("count_pca_x", "PC for X-axis:", choices = paste0("PC",1:10), selected = "PC1"),
                    selectInput("count_pca_y", "PC for Y-axis:", choices = paste0("PC",1:10), selected = "PC2"),
                    # Button to reset PCA
                    actionButton(inputId = "reset_count_pca_button", 
                                 label = "Default PCA",
                                 icon = icon("refresh"))
                ),
                box(width=9, status="primary", solidHeader=TRUE,
                    # Display PCA plot
                    plotOutput("count_pca_plot")
                )
              )
      ),
      
      # ========== DE Tabs =============================================================
      tabItem(tabName = "de_tables_tab",
              fluidRow(
                box(width=3, status="primary", solidHeader=TRUE,
                    
                    # Upload differential expression results
                    fileInput(inputId = "de", 
                              label = "Upload Differential Expression File (.csv):", 
                              multiple = FALSE, 
                              accept = ".csv",
                              buttonLabel = "Browse",
                              placeholder = "deg.csv"
                    ),
                    
                    # Choose X-axis column for Volcano plot
                    radioButtons("x_name", "X-axis Column:",
                                 choices = c("baseMean" = "baseMean",
                                             "log2FoldChange" = "log2FoldChange",
                                             "lfcSE" = "lfcSE",
                                             "stat" = "stat",
                                             "pvalue" = "pvalue",
                                             "padj" = "padj"),
                                 selected = "log2FoldChange"),
                    
                    # Choose Y-axis column for Volcano plot
                    radioButtons("y_name", "Y-axis Column:",
                                 choices = c("baseMean" = "baseMean",
                                             "log2FoldChange" = "log2FoldChange",
                                             "lfcSE" = "lfcSE",
                                             "stat" = "stat",
                                             "pvalue" = "pvalue",
                                             "padj" = "padj"),
                                 selected = "padj"),
                    
                    # Set base and highlight colors for points
                    colourInput("base_color", "Base point color",
                                value = "#000000"),
                    colourInput("highlight_color", "Highlight point color",
                                value = "#FF0000"),
                    
                    # Slider for p-value cutoff in Volcano plot
                    sliderInput("volcano_slider", "P adjusted coloring magnitude:",
                                min = -300, max = 0, value = -10, step = 1),
                    
                    # Reset slider button
                    actionButton(inputId = "reset_de_volcano", 
                                 label = "Reset Slider",
                                 icon = icon("undo")),
                    
                    tags$hr(),
                    
                    # Button to draw Volcano plot
                    actionButton("plot_button", "Plot", icon = icon("paint-brush"))
                ),
                box(width=9, status="primary", solidHeader=TRUE,
                    # Display DE results table
                    DTOutput("deg_table")
                )
              )
      ),
      
      tabItem(tabName = "de_volcano_tab",
              fluidRow(
                box(width=10, status="primary", solidHeader=TRUE,
                    # Display Volcano plot
                    plotOutput("volcano_plot")
                )
              )
      ),
      
      tabItem(tabName = "de_plot_table_tab",
              fluidRow(
                box(width=10, status="primary", solidHeader=TRUE,
                    # Display table of genes highlighted in the Volcano plot
                    DTOutput("volcano_plot_table")
                )
              )
      ),
      
      # ========== GSEA Tabs ===========================================================
      tabItem(tabName = "gsea_top_tab",
              fluidRow(
                box(width=2, status="primary", solidHeader=TRUE,
                    p("Adjust the p-value threshold and click 'Update'"),
                    
                    # Upload fgsea results file
                    fileInput('fgsea_file', label='Load fgsea results.', buttonLabel = 'Browse...', accept=c('.csv', '.tsv'), placeholder = 'fgsea_res.csv'),
                    
                    # Set adjusted p-value threshold for GSEA results
                    sliderInput('slider4', min = 0.00000001, max = 0.05, 'Select adjusted p-value threshold', value=0.01),
                    
                    # Update GSEA bar plot
                    actionButton('gsea_update_button', 'Update', width='100%')
                ),
                box(width=10, status="primary", solidHeader=TRUE,
                    # Display GSEA results as a bar plot
                    plotOutput('gsea_bar', height = 2500, width = 900)
                )
              )
      ),
      
      tabItem(tabName = "gsea_table_tab",
              fluidRow(
                box(width=2, status="primary", solidHeader=TRUE,
                    
                    # GSEA table filtering options
                    radioButtons('radio_gsea', 'Choose an option', choices = c('All genes', "Negative NES", "Positive NES"), selected='All genes'),
                    
                    # Button to filter GSEA table
                    actionButton('gsea_filter_button', 'Filter', width='100%'),
                    
                    # Button to download filtered GSEA results
                    downloadButton('download_gsea', 'Download Filtered GSEA Results')
                ),
                box(width=9, status="primary", solidHeader=TRUE,
                    # Display GSEA results table. Clicking a row shows gene details.
                    DT::dataTableOutput('gsea')
                )
              )
      ),
      
      tabItem(tabName = "gsea_plots_tab",
              fluidRow(
                box(width=3, status="primary", solidHeader=TRUE,
                    # Button to draw GSEA scatter plot
                    actionButton('gsea_plot_button', 'Plot Scatter', width='100%')
                ),
                box(width=9, status="primary", solidHeader=TRUE,
                    # Display GSEA scatter plot
                    plotOutput('gsea_scatter')
                )
              )
      )
    )
  )
)

#-------------------------------------
# Server Section -----------------------------------------------------------

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 50*1024^2)
  
  #---------------------------------------------------------------------------
  # Samples Section -----------------------------------------------------------
  #---------------------------------------------------------------------------
  
  # Load user-uploaded sample data
  load_data <- reactive({
    req(input$file_input) 
    if (!grepl("\\.(csv)$", input$file_input$name, ignore.case = TRUE)) {
      stop("Please upload a CSV file.")
    }
    tryCatch(
      {
        read.csv(input$file_input$datapath, header = TRUE)
      },
      error = function(e) {
        stop("Failed to read the file. Please check if it's well-formatted.")
      }
    )
  })
  
  # Summarize sample data columns (mean, sd or distinct values)
  samples_summary <- function(samples_data) {
    names <- colnames(samples_data)
    types <- sapply(samples_data, class)
    means_or_values <- sapply(samples_data, function(x) {
      if (class(x) %in% c("integer", "numeric", "double")) {
        mean_value <- round(mean(na.omit(x)))
        sd_value <- round(sd(na.omit(x)))
        paste0(mean_value, " (+/-", sd_value, ")")
      } else {
        paste0(paste(unique(x), collapse = ", "))
      }
    })
    summary_df <- data.frame("Column name" = names,
                             "Type" = types,
                             "Mean (sd) or Distinct Values" = means_or_values, check.names = FALSE)
    return(summary_df)
  }
  
  # Diagnostic plot for samples (e.g., violin plot by diagnosis)
  samples_diagnostic_plot <- function(samples_data) {
    ggplot(samples_data, aes(x = Diagnosis , y = age_of_death, fill= Diagnosis)) +
      geom_violin() +
      labs(title = "Age of death based on diagnosis", x = "Diagnosis", y = "Age of Death") 
  }
  
  # After clicking Submit, display data and plots
  observeEvent(input$submit_button, {
    # Show summary stats table
    output$summary_output <- renderDataTable({
      DT::datatable(samples_summary(load_data()), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    # Show original sample data
    output$table_output <- renderDataTable({
      DT::datatable(load_data(), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    # Show diagnostic plot
    output$plots_output <- renderPlotly({
      samples_diagnostic_plot(load_data())
    })
    
    # Update available numeric columns for distribution plot
    observe({
      req(load_data())
      numeric_cols <- names(load_data())[sapply(load_data(), is.numeric)]
      updateSelectInput(session, "sample_plot_var", choices = numeric_cols)
    })
  })
  
  # Show distribution plot according to selected column and type
  output$sample_dist_plot <- renderPlot({
    req(load_data(), input$sample_plot_var, input$sample_plot_type)
    df <- load_data()
    validate(
      need(input$sample_plot_var %in% names(df), "Selected variable not found in data"),
      need(is.numeric(df[[input$sample_plot_var]]), "Selected variable is not numeric")
    )
    if(input$sample_plot_type == "hist") {
      ggplot(df, aes_string(x = input$sample_plot_var)) +
        geom_histogram(fill = "blue", color = "white", bins = 30) +
        theme_classic()
    } else {
      ggplot(df, aes_string(x = input$sample_plot_var)) +
        geom_density(fill = "red", alpha = 0.5) +
        theme_classic()
    }
  })
  
  #---------------------------------------------------------------------------
  # Counts Section ------------------------------------------------------------
  #---------------------------------------------------------------------------
  
  # Reset counts filter thresholds
  observeEvent(input$reset_counts_main, {
    updateNumericInput(session = session, inputId = "count_percentile_slider", value = 80)
    updateNumericInput(session = session, inputId = "count_nonzero_slider", value = 5)
  })
  
  # Filter counts data (based on variance and non-zero counts threshold)
  filtered_counts <- reactive({
    req(input$counts)
    counts_data <- rio::import(input$counts$datapath)
    variances <- apply(counts_data[, -1], 1, var)
    variance_threshold <- quantile(variances, input$count_percentile_slider / 100)
    filtered_by_variance <- counts_data[variances >= variance_threshold, ]
    nonzero_counts <- rowSums(filtered_by_variance[, -1] != 0)
    filtered_by_nonzero <- filtered_by_variance[nonzero_counts >= input$count_nonzero_slider, ]
    filtered_data <- filtered_by_nonzero
    
    num_samples <- ncol(counts_data) - 1
    num_genes <- nrow(counts_data)
    num_passed_filter <- nrow(filtered_data)
    percent_passed_filter <- (num_passed_filter / num_genes) * 100
    num_not_passed_filter <- num_genes - num_passed_filter
    percent_not_passed_filter <- (num_not_passed_filter / num_genes) * 100
    
    summary_table <- tibble::tribble(
      ~Measure, ~Summary,
      "Number of Samples", num_samples,
      "Number of Genes", num_genes,
      "Number of Genes Passed Filter", num_passed_filter,
      "% Passed Filter", percent_passed_filter,
      "Number of Genes Not Passed Filter", num_not_passed_filter,
      "% Not Passed Filter", percent_not_passed_filter
    )
    return(list(counts_data = counts_data, filtered_by_nonzero = filtered_data, 
                summary_table = summary_table, variances = variances, variance_threshold = variance_threshold))
  })
  
  # Display counts summary table
  output$counts_summary_table <- renderDT({
    DT::datatable(filtered_counts()$summary_table, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Median counts vs variance plot
  output$plot_count_vs_variance <- renderPlot({
    req(input$counts)
    counts_data <- filtered_counts()$counts_data
    variances <- filtered_counts()$variances
    variance_threshold <- filtered_counts()$variance_threshold
    
    medians <- apply(counts_data[, -1], 1, median)
    df <- data.frame(medians, variances) %>%
      dplyr::mutate(Threshold = ifelse(variances < variance_threshold, "Failed Threshold", "Passed Threshold"))
    
    ggplot(df, aes(x = medians, y = variances, color = Threshold)) +
      geom_point() +
      scale_color_manual(values = c("Failed Threshold" = "black", "Passed Threshold" = "red")) +
      scale_x_log10() +
      scale_y_log10() +
      labs(x = "Log 10 Median Count", y = "Log 10 Variance") +
      theme_classic() +
      theme(legend.position = "bottom")
  })
  
  # Median counts vs number of zeros plot
  output$plot_count_vs_zeros <- renderPlot({
    req(input$counts)
    counts_data <- filtered_counts()$counts_data
    
    medians <- apply(counts_data[, -1], 1, median)
    num_zeros <- apply(counts_data[, -1], 1, function(x) sum(x == 0))
    df <- data.frame(medians, num_zeros) %>%
      dplyr::mutate(Threshold = ifelse(num_zeros > input$count_nonzero_slider, "Failed Filter", "Passed Filter"))
    
    ggplot(df, aes(x = medians, y = num_zeros, color = Threshold)) +
      geom_point() +
      scale_color_manual(values = c("Failed Filter" = "black", "Passed Filter" = "red")) +
      scale_x_log10() +
      scale_y_log10() +
      ylim(0, 70) +
      labs(x = "Log 10 Median Count", y = "Log 10 Number of Samples with Zero Counts") +
      theme_classic() +
      theme(legend.position = "bottom")
  })
  
  # Reactive to choose if log transform data for heatmap
  count_heatmap_radio <- reactive({
    if (is.null(input$count_heatmap_transformation)) {
      return("TRUE")
    }
    input$count_heatmap_transformation
  })
  
  # Display heatmap
  output$count_heatmap <- renderPlot({
    req(input$count_heatmap_button, input$counts)
    data <- filtered_counts()$filtered_by_nonzero %>%
      data.frame() %>%
      tibble::remove_rownames() %>%
      dplyr::select(-gene) %>%
      as.matrix()
    if(count_heatmap_radio() == "TRUE") {
      data <- log(data + 1)
    } 
    isolate({
      heatmap.2(data, trace = "none", col = bluered(75), scale = "row")
    })
  })
  
  # Reset PCA selection
  observeEvent(input$reset_count_pca_button, {
    updateSelectInput(session = session, inputId = "count_pca_x", selected = "PC1")
    updateSelectInput(session = session, inputId = "count_pca_y", selected = "PC2")
  })
  
  # Display PCA plot
  output$count_pca_plot <- renderPlot({
    req(input$counts)
    data <- filtered_counts()$filtered_by_nonzero %>%
      data.frame() %>%
      tibble::remove_rownames() %>%
      dplyr::select(-gene) %>%
      as.matrix()
    pca <- prcomp(data)
    pca_df <- as.data.frame(pca$x)
    ggplot(pca_df, aes_string(x=input$count_pca_x, y=input$count_pca_y)) +
      geom_point() +
      xlab(input$count_pca_x) +
      ylab(input$count_pca_y) +
      ggtitle("PCA Plot") +
      theme_classic()
  })
  
  #---------------------------------------------------------------------------
  # DE Section ---------------------------------------------------------------
  #---------------------------------------------------------------------------
  
  # Reset volcano plot p-value threshold
  observeEvent(input$reset_de_volcano, {
    updateSliderInput(session, "volcano_slider", value = -10)
  })
  
  # Load DE data
  de_data <- reactive({
    req(input$de)
    deg <- rio::import(input$de$datapath)
    deg
  })
  
  # Display DE results table
  output$deg_table <- renderDT({
    req(input$de)
    data <- de_data()
    DT::datatable(data, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  volcano_data <- reactiveVal(NULL)
  observeEvent(input$plot_button, {
    req(input$de)
    df <- de_data()
    volcano_data(df)
  })
  
  # Display Volcano plot
  output$volcano_plot <- renderPlot({
    req(volcano_data())
    df <- volcano_data()
    req(input$x_name, input$y_name)
    validate(
      need(input$x_name %in% names(df), paste("Column", input$x_name, "not found in data")),
      need(input$y_name %in% names(df), paste("Column", input$y_name, "not found in data"))
    )
    
    df$highlight <- df[[input$y_name]] <= 10^(input$volcano_slider)
    
    ggplot(df, aes_string(x = input$x_name, y = paste0("-log10(", input$y_name, ")"), color = "highlight")) +
      geom_point(size = 1) +
      scale_color_manual(values = c("FALSE" = input$base_color, "TRUE" = input$highlight_color)) +
      labs(x = input$x_name, y = input$y_name) + 
      theme_classic() + 
      theme(legend.position = "bottom")
  })
  
  # Display table of significant genes from Volcano plot
  output$volcano_plot_table <- renderDT({
    req(volcano_data())
    df <- volcano_data()
    req(input$y_name %in% names(df))
    filtered_df <- df %>% dplyr::filter(!!sym(input$y_name) <= 10^(input$volcano_slider))
    DT::datatable(filtered_df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  #---------------------------------------------------------------------------
  # GSEA Section -------------------------------------------------------------
  #---------------------------------------------------------------------------
  
  # Load GSEA data
  load_fgsea <- reactive({
    req(input$fgsea_file)
    data <- read_csv(input$fgsea_file$datapath, col_names = TRUE)
    data
  })
  
  # Function to plot GSEA bar graph
  plot_gsea_bar <- function(fgseaData, slider4){
    fgseaData <- fgseaData %>%
      dplyr::filter(padj<slider4) %>%
      dplyr::mutate(pathway=str_replace_all(pathway, '_', ' ')) %>%
      dplyr::mutate(pathway = str_wrap(pathway, width = 100)) %>%
      dplyr::mutate(pathway=factor(pathway, levels=pathway))
    ggplot(fgseaData, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=NES<0)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title="C2 Canonical Pathways NES from GSEA") + 
      scale_fill_discrete('NES', labels=c('Positive', 'Negative')) +
      theme_minimal()
  }
  
  # Filter GSEA results table based on user selections
  fgsea_res_table <- function(fgseaData, slider4, radio_gsea){
    fgseaData <- dplyr::filter(fgseaData, padj < slider4) %>%
      dplyr::arrange(padj) 
    if (radio_gsea == 'Positive NES') {
      fgseaData <- dplyr::filter(fgseaData, NES>0)
    } else if (radio_gsea == 'Negative NES') {
      fgseaData <- dplyr::filter(fgseaData, NES<0)
    }
    fgseaData
  }
  
  # Function to plot GSEA scatter plot
  plot_gsea_scatter <- function(fgseaData, slider4){
    fgseaData <- fgseaData %>%
      dplyr::mutate(neglog10padj = -log(padj, base=10))
    ggplot(fgseaData, aes(x=NES, y=neglog10padj, color=padj<slider4)) +
      geom_point() +
      labs(title = 'Scatter Plot of GSEA Results', 
           x = 'Normalized Enrichment Score', 
           y='-log10(padj)') +
      scale_color_manual('Adjusted P-value', 
                         labels = c(paste('>', slider4), paste('<', slider4)),
                         values = c('#22577A', '#FFCF56'))
  }
  
  gsea_data <- reactiveVal(NULL)
  
  # Update GSEA data
  observeEvent(input$gsea_update_button, {
    req(input$fgsea_file)
    data <- load_fgsea()
    gsea_data(data)
  })
  
  # Display GSEA bar plot
  output$gsea_bar <- renderPlot({
    req(gsea_data())
    plot_gsea_bar(gsea_data(), input$slider4)
  })
  
  # Filter GSEA table on button click
  observeEvent(input$gsea_filter_button, {
    req(gsea_data())
    output$gsea <- DT::renderDataTable({
      DT::datatable(fgsea_res_table(gsea_data(), input$slider4, input$radio_gsea), 
                    options = list(pageLength=10, scrollX=TRUE), selection = 'single')
    }, striped = TRUE)
  }, ignoreInit = TRUE)
  
  # Download filtered GSEA results
  output$download_gsea <- downloadHandler(
    filename = 'filtered_gsea.csv',
    content = function(file) {
      write_csv(fgsea_res_table(gsea_data(), input$slider4, input$radio_gsea), file)
    }
  )
  
  # Plot GSEA scatter plot
  observeEvent(input$gsea_plot_button, {
    req(gsea_data())
    output$gsea_scatter <- renderPlot({
      plot_gsea_scatter(gsea_data(), input$slider4)
    })
  }, ignoreInit = TRUE)
  
  # Show leadingEdge genes in a modal when selecting a row in GSEA table
  observeEvent(input$gsea_rows_selected, {
    req(gsea_data())
    data_filtered <- fgsea_res_table(gsea_data(), input$slider4, input$radio_gsea)
    sel <- input$gsea_rows_selected
    if(length(sel) == 1) {
      if("leadingEdge" %in% colnames(data_filtered)) {
        genes <- data_filtered$leadingEdge[[sel]]
        genes_text <- paste(genes, collapse=", ")
        showModal(modalDialog(
          title = paste("Genes in pathway:", data_filtered$pathway[sel]),
          p(genes_text),
          easyClose = TRUE,
          footer = NULL
        ))
      }
    }
  })
  
}

#-------------------------------------
# Run the App --------------------------------------------------------------
shinyApp(ui, server)

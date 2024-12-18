library(shiny)
library(tidyverse)
library(lme4)
library(emmeans)
library(DT)
library(patchwork)
library(ggpubr)

# Source helper functions
source("R/data_preparation.R")
source("R/analysis.R")
source("R/subject_analysis.R")

# UI
ui <- fluidPage(
  titlePanel("Brain Region Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      # First level selection
      selectInput("analysisType",
                  "Select Analysis Type",
                  choices = c(
                    "Software Comparison" = "software",
                    "Coefficient of Variation" = "cv",
                    "Distance Metrics" = "distances",
                    "Discrimination Ratios" = "ratios"
                  )),
      
      # Software comparison options
      conditionalPanel(
        condition = "input.analysisType == 'software'",
        selectInput("region", 
                    "Select Brain Region",
                    choices = NULL)
      ),
      
      # Distance metrics options
      conditionalPanel(
        condition = "input.analysisType == 'distances'",
        selectInput("distanceMetric",
                    "Select Distance Metric",
                    choices = c(
                      "Jensen-Shannon Divergence" = "mean_js_div",
                      "Kendall's Tau" = "mean_kendall_tau",
                      "Wilcoxon Effect Size" = "mean_wilcox_paired"
                    ))
      ),
      
      # Ratio options
      conditionalPanel(
        condition = "input.analysisType == 'ratios'",
        selectInput("ratioMetric",
                    "Select Ratio Metric",
                    choices = c(
                      "Jensen-Shannon Ratio" = "js_ratio",
                      "Kendall's Tau Ratio" = "kendall_ratio",
                      "Wilcoxon Ratio" = "wilcox_paired_ratio"
                    ))
      ),
      
      # Region type selector for everything except software comparison
      conditionalPanel(
        condition = "input.analysisType != 'software'",
        selectInput("regionType",
                    "Select Region Type",
                    choices = c(
                      "Left Hemisphere" = "lh",
                      "Right Hemisphere" = "rh",
                      "Subcortical & Other" = "other"
                    ))
      ),
      
      downloadButton("downloadPlot", "Download Plot")
    ),
    
    mainPanel(
      tabsetPanel(id = "mainTabs",
                  tabPanel("Plots", 
                           plotOutput("mainPlot", height = "800px")
                  ),
                  
                  # Statistical Results tab
                  tabPanel("Statistical Results",
                           conditionalPanel(
                             condition = "input.analysisType == 'software'",
                             h4("Pairwise Comparisons"),
                             DTOutput("pairwiseTable"),
                             br(),
                             h4("Correlations between Software"),
                             DTOutput("correlationTable"),
                             br(),
                             h4("Mixed Model Results"),
                             verbatimTextOutput("mixedModelSummary")
                           )
                  ),
                  
              
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Read and process data
  data <- read.csv("data/synthseg_FS_concat_all_versions.csv")
  
  processed_data <- reactive({
    prepare_data(data, normalize_by_icv = FALSE, scale_volumes = FALSE)
  })
  
  # Calculate subject differences once
  subject_analysis <- reactive({
    analyze_subject_differences(processed_data())
  })
  
  # Update region choices for software comparison
  observe({
    req(input$analysisType == "software")
    regions <- unique(processed_data()$region)
    updateSelectInput(session, "region",
                      choices = regions)
  })
  
  # Main plot output
  output$mainPlot <- renderPlot({
    if(input$analysisType == "software") {
      req(input$region)
      # Original software comparison plot
      results <- analyze_region(processed_data(),
                                input$region,
                                p_adjust_method = "bonferroni")
      results[[4]]
      
    } else if(input$analysisType == "cv") {
      # CV analysis
      create_single_metric_plot(subject_analysis()$results, 
                                "mean_within_cv",
                                "Within-Subject Coefficient of Variation",
                                input$regionType)
      
    } else if(input$analysisType == "distances") {
      # Distance metrics
      create_single_metric_plot(subject_analysis()$results,
                                input$distanceMetric,
                                get_metric_title(input$distanceMetric),
                                input$regionType)
      
    } else {
      # Ratio metrics
      create_single_metric_plot(subject_analysis()$results,
                                input$ratioMetric,
                                get_metric_title(input$ratioMetric),
                                input$regionType)
    }
  })
  
  # Model and statistical outputs
  output$pairwiseTable <- renderDT({
    req(input$analysisType == "software", input$region)
    results <- analyze_region(processed_data(),
                              input$region,
                              p_adjust_method = "bonferroni")
    datatable(results$pairwise_tests,
              options = list(pageLength = 5),
              rownames = FALSE) %>%
      formatRound(columns = c("statistic", "p_value", "p_adjusted"), 
                  digits = 4)
  })
  
  output$correlationTable <- renderDT({
    req(input$analysisType == "software", input$region)
    results <- analyze_region(processed_data(),
                              input$region,
                              p_adjust_method = "bonferroni")
    datatable(results$correlations,
              options = list(pageLength = 5),
              rownames = FALSE) %>%
      formatRound(columns = c("correlation", "cor_pvalue"), 
                  digits = 4)
  })
  
  output$mixedModelSummary <- renderPrint({
    req(input$analysisType == "software", input$region)
    results <- analyze_region(processed_data(),
                              input$region,
                              p_adjust_method = "bonferroni")
    results$mixed_model
  })
  
  output$modelFormula <- renderPrint({
    req(input$analysisType == "software", input$region)
    results <- analyze_region(processed_data(),
                              input$region,
                              p_adjust_method = "bonferroni")
    formula(results$mixed_model$object)
  })
  
  output$randomEffects <- renderPrint({
    req(input$analysisType == "software", input$region)
    results <- analyze_region(processed_data(),
                              input$region,
                              p_adjust_method = "bonferroni")
    print(VarCorr(results$mixed_model$object), comp = "Variance")
  })
  
  output$fixedEffects <- renderPrint({
    req(input$analysisType == "software", input$region)
    results <- analyze_region(processed_data(),
                              input$region,
                              p_adjust_method = "bonferroni")
    coef(summary(results$mixed_model$object))
  })
  
  output$modelFit <- renderPrint({
    req(input$analysisType == "software", input$region)
    results <- analyze_region(processed_data(),
                              input$region,
                              p_adjust_method = "bonferroni")
    model <- results$mixed_model$object
    c(AIC = AIC(model),
      BIC = BIC(model),
      logLik = logLik(model))
  })
  
  # Download handler
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("brain_analysis_", input$analysisType, ".pdf")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), width = 12, height = 16)
    }
  )
}

# Helper functions for plotting
create_single_metric_plot <- function(results, metric, title, region_type) {
  create_horizontal_barplot(
    results, 
    metric,
    paste0(title, if(metric == "mean_within_cv") " ↓" else " ↑"),
    region_type,
    sort_by = "metric"
  ) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.text.y = element_text(size = 10)
    ) +
    labs(
      subtitle = case_when(
        region_type == "lh" ~ "Left Hemisphere Regions",
        region_type == "rh" ~ "Right Hemisphere Regions",
        region_type == "other" ~ "Subcortical & Other Regions"
      )
    )
}

get_metric_title <- function(metric) {
  case_when(
    metric == "mean_js_div" ~ "Jensen-Shannon Divergence",
    metric == "mean_kendall_tau" ~ "Kendall's Tau Distance",
    metric == "mean_wilcox_paired" ~ "Wilcoxon Effect Size (Paired)",
    metric == "js_ratio" ~ "Jensen-Shannon Discrimination Ratio",
    metric == "kendall_ratio" ~ "Kendall's Tau Discrimination Ratio",
    metric == "wilcox_paired_ratio" ~ "Wilcoxon Effect Size Discrimination Ratio",
    TRUE ~ gsub("_", " ", metric)
  )
}

# Run the application
shinyApp(ui = ui, server = server)
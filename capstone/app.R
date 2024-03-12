library(shiny)
library(ggplot2)
library(dplyr)
library(plotly)
library(readr)
library(ggrepel)

cancer_type_mapping <- c(
  "Melanoma" = "Melanoma",
  "Breast Carcinoma" = "breast_carcinoma",
  "Ovarian Cancer" = "ovarian_cancer",
  "Prostate Cancer" = "prostate_cancer"
)


ui <- fluidPage(
  titlePanel("Transcriptome-Wide Association Studies Application on Cancers"),
  
  
  h3("Manhattan Plot for TWAS Results"),
  sidebarLayout(
    sidebarPanel(
      selectInput("cancerType", "Select Cancer Type:", choices = cancer_type_mapping)
    ),
    mainPanel(
      plotlyOutput("manhattanPlot")
    )
  ),
  

  hr(),
  
  h3("Explore Data by Gene ID"),
  sidebarLayout(
    sidebarPanel(
      textInput("geneID", "Enter Gene ID:", value = "ExampleGeneID"),
      helpText("Example genes to try: SPEN, TRIT1, PHAX")
    ),
    mainPanel(
      tableOutput("resultsTable")
    )
  )
)


server <- function(input, output) {
  combined_twas_results <- read_csv("full_twas_susie.csv", show_col_types = FALSE)
  combined_twas_results$CHR <- as.factor(combined_twas_results$CHR)
  
  output$resultsTable <- renderTable({
    if (input$geneID != "") {
      filteredData <- dplyr::filter(combined_twas_results, ID == input$geneID) %>%
        dplyr::select(CHR, HSQ, MODEL, TWAS.Z, TWAS.P, cancer_type) 
      if (nrow(filteredData) == 0) {
        return(data.frame(Message = "No results found for the input gene ID. Please try another ID."))
      }
      filteredData
    } else {
      return(data.frame(Message = "Please enter a gene ID to see the results. Example genes to try: SPEN, TRIT1, PHAX"))
    }
  })
  
  output$manhattanPlot <- renderPlotly({
    filtered_data <- combined_twas_results %>%
      filter(cancer_type == input$cancerType)
    
    significance_threshold <- -log10(0.05 / nrow(filtered_data))
    
    manhattan_plot <- ggplot(filtered_data, aes(x = POS, y = -log10(TWAS.P), colour = CHR, 
                                                text = paste("Gene ID:", ID, "<br>P-value:", TWAS.P))) +
      geom_point(alpha = 0.6, size = 1.2) +
      scale_colour_manual(values = rainbow(length(unique(filtered_data$CHR)))) +
      geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(title = paste0("Manhattan Plot for ", input$cancerType), 
           x = "Position", y = "-log10(p-value)", colour = "Chromosome") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    significant_snps <- subset(filtered_data, -log10(TWAS.P) > significance_threshold)
    if(nrow(significant_snps) > 0) {
      manhattan_plot <- manhattan_plot +
        geom_text_repel(data = significant_snps,
                        aes(label = ID),  
                        nudge_y = 0.1,  
                        size = 3)  
    }
    
    ggplotly(manhattan_plot, tooltip = "text")
  })
  
}

shinyApp(ui = ui, server = server)

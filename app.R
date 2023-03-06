# Load libraries
library(shiny)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(Seurat)
library(viridis) #color palettes
library(cowplot) #plot_grid
library(Hmisc)

# Load helper functions
lapply(list.files(pattern = "helpers"), source)
# c("CD14","KRT18","KRT19","EGFR","EPCAM","VIM","FOXA1","LTF","GATA3","TP63","ITGA6","FOLR1")
# Import proteome data
allData <- list(
  "Human Proteome" = readRDS("data/human_proteome.rds"),
  "Mouse Proteome" = readRDS("data/mouse_proteome.rds")
)

# Heatmap color palettes
hm_colors <- rownames(RColorBrewer::brewer.pal.info)

# Color palettes 
ui <- navbarPage("Khokha Lab MEC Explorer",
                 
                 # Tab 1
                 tabPanel("Omics data",  
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            # Sidebar panel for inputs ----
                            sidebarPanel( 
                              # CSS
                              id="omicsPanel", 
                              style="overflow-y:scroll; max-height:800px; position=absolute",
                              
                              # Scaling
                              selectInput("species", label = "Select dataset", choices = c("Human Proteome","Mouse Proteome"), selected = "Human Proteome"),
                              
                              # Title for plot
                              textInput("file_label", h4("Title for Plot:"),  value = ""),
                              
                              # Make drop down menu with all gene symbols in proteome
                              selectizeInput("geneIDs", multiple = T, label =  h4('Select proteins of interest'), choices = NULL, options = list(create = FALSE)),
                              
                              h4("Data Transformations"),
                              # z-score
                              checkboxGroupInput("transf_checkbox", label = "", choices = list("z-score (row)" = "z_row", "z-score (column)" = "z_column"), selected = 0),
                              
                              h4("Heatmap Preferences"),
                              # Scaling
                              selectInput("scale", label = "Scale by: ", choices = c("Row", "Column", "None"), selected = "Row"),
                              # Clustering/showing name
                              checkboxGroupInput("row", label = "Rows", choices = list("cluster" = "cluster", "show_name" = "show_name"), selected = c("show_name", "cluster")),
                              checkboxGroupInput("column", label = "Columns", choices = list("cluster" = "cluster", "show_name" = "show_name"), selected = 0),
                              # Cell size
                              numericInput("cell_h", "Cell height", value = "", width = 100),
                              numericInput("cell_w", "Cell width", value = "", width = 100),
                              
                              # Gradient color picker
                              selectInput("color_pal", label = "Pick heatmap color", choices = hm_colors, selected = hm_colors[3]),
                              radioButtons("color_rev", label = "Reverse palette?",  choices = c("No", "Yes"), selected = "No"),
                              
                              # Download plot and data
                              selectInput("file_format", label = h4("Save plot as"), choices = c("png", "jpeg", "pdf", "svg", "tiff"), selected = "png"),
                              downloadButton('downloadPlot', 'Download Plot'),
                              downloadButton('downloadData', 'Download Data')
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              plotOutput(outputId = "heatmap")
                            ))),
                 
                 # Tab 1
                 tabPanel("Single-cell data",
                          
                          uiOutput("sc_Shiny")
                          )
)

# Define server logic ----
server <- function(input, output, session) {
  # Load data
  dta <- reactive({
    allData[[as.character(input$species)]]
  })
  
  # For drop down menu with gene names, update and filter server-side to speed up program
  observeEvent(input$species, {
    updateSelectizeInput(session, 'geneIDs', 
                         choices = rownames(dta()$vals), 
                         selected = dta()$selected_markers,
                         server = TRUE)
    }
  )
  
  # Process data
  df_mat_s <- reactive({
    df_s <- dta()$vals
    if(length(input$geneIDs) != 0)
      # Subset
      df_s <-
        get_to_keep(rownames(dta()$vals), list_to_match = input$geneIDs) %>%
        subset_table(dta()$vals, rows=.)
    
    # Transform data
    if(is.null(input$transf_checkbox))
      return(df_s)
    transform_data(df_s, input$transf_checkbox)
  })
  # 
  # Color of heatmap
  color <- reactive({
    get_col_pal(brew_pal = input$color_pal, rev =  ifelse(input$color_rev == "No", F, T)) %>%
      get_col_grad(n=50) %>%
      return
  })
  
  # Make heatmap
  heatMap <- reactive({
    pheatmap(mat = df_mat_s(),
             scale = tolower(input$scale),
             main = input$file_label,
             # Colors
             color = color(),
             # Cell size
             cellwidth = ifelse(input$cell_w == "", NA, input$cell_w),
             cellheight = ifelse(input$cell_h == "", NA, input$cell_h),
             # Annotations
             annotation_row = NA, 
             annotation_col = dta()$colAnn,
             annotation_colors = dta()$colors,
             # Rows
             cluster_rows = check(input$row, "cluster"),
             show_rownames = check(input$row, "show_name"),
             # Columns
             cluster_cols = check(input$column, "cluster"),
             show_colnames = check(input$column, "show_name"),
             # Other preferences
             border_color = NA)
  })
  # 
  # Render output
  output$heatmap <- renderPlot({
    heatMap()
  })
  
  # Download plot
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0(input$file_label, '.', input$file_format) },
    content = function(file) {
      ggsave(file, plot = heatMap(), device = input$file_format)
    }
  )
  
  # Download csv of processed dataset
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$file_label, ".csv")
    },
    
    content = function(file) {
      df <- cbind(dta()$colAnn, t(df_mat_s()))
      write.csv(df, file, row.names = TRUE)
    }
  )
  
  # Link for second tab
  sc_Shiny_url <- a("Single-cell Shiny explorer", href="https://kazeera.shinyapps.io/single-cell-explorer/")
  output$sc_Shiny <- renderUI({
    tagList("Link: ", sc_Shiny_url)
  })
  
}

# Run app 
shinyApp(ui = ui, server = server)

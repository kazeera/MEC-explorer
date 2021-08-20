# Load libraries
library(shiny)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)

#' Acronyms:
#' hm = heatmap
#' ss = subset
#' ann = annotation

# Load helper functions
lapply(list.files(pattern = "helpers"), source)

# Import data
df_mat <- read.delim("data/human_MEC_proteome_combat_matrix.txt", row.names = 1)
df_row <- read.delim("data/human_MEC_proteome_row_info.txt", row.names = 1)
df_col <- read.delim("data/human_MEC_proteome_sample_info.txt", row.names = 1)  %>%
   .[,c("cellType", "patientAgeGroup", "hormoneStatus")]

# Color palettes
all_pals <- rownames(RColorBrewer::brewer.pal.info)

# Define UI ----
ui <- fluidPage(
  # App title ----
  titlePanel(title = "Heatmap Explorer"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # helpText("Query human mammary epithelial cell proteome first published in Mahendralingam et al., Nature Metabolism, 2021."),
    # Sidebar panel for inputs ----
    sidebarPanel( 
      # CSS
      id="tPanel", 
      style="overflow-y:scroll; max-height:600px; position=absolute",
      textInput("file_label", h4("Title for Plot:"),  value = ""),
      
      # Process data
      h4("Subset Rows"),
      textInput("ss_row_list",  h5("Input list here:"), value = "CD14 KRT18 KRT19 EGFR EPCAM VIM FOXA1 LTF GATA3 TP63 ITGA6 FOLR1"),
      textInput("ss_row_delim", h6("Delimited by:"), value = " "),
      
      h3("Perform data transformations:"),
      checkboxGroupInput("transf_checkbox", label = "", choices = list("z-score (row)" = "z_row"), selected = 0),
      
      # Heatmap visualization options
      h3("Heatmap Preferences"),
      # Scaling
      selectInput("hm_scale", label = "Scale by: ", choices = c("Row", "Column", "None"), selected = "Row"),
      # Clustering/showing name
      checkboxGroupInput("hm_row", label = "Rows", choices = list("cluster" = "cluster", "show_name" = "show_name"), selected = c("show_name", "cluster")),
      checkboxGroupInput("hm_column", label = "Columns", choices = list("cluster" = "cluster", "show_name" = "show_name"), selected = 0),
      # Cell size
      numericInput("hm_cell_h", h6("Cell height"), value = "", width = 100),
      numericInput("hm_cell_w", h6("Cell width"), value = "", width = 100),

      # Gradient color picker
      selectInput("hm_color_pal", label = "Pick heatmap color", choices = all_pals, selected = all_pals[3]),
      radioButtons("hm_color_rev", label = "Reverse palette?",  choices = c("No", "Yes"), selected = "No"),
      
      # Save Plot
      selectInput("file_format", label = "Save plot as", choices = c("png", "jpeg", "pdf", "svg", "tiff"), selected = "png"),
      downloadButton('downloadPlot', 'Download Plot'),
      
      # helpText("Data source: Human mammary epithelial cell proteome"),
      downloadButton('downloadData', 'Download Data')
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      # h3("Heatmap"),
      plotOutput(outputId = "hm")
    )
  ))


# Define server logic ----
server <- function(input, output) {
  # Process data
  df_mat_s <- reactive({
    # Return full data frame if input list is empty
    if(is.null(input$ss_row_list) | input$ss_row_list == "") 
         return(df_mat)
    
    # Subset based on input list
    df_s <- strsplit(input$ss_row_list, input$ss_row_delim) %>% 
      unlist %>% 
      get_to_keep(rownames(df_mat), list_to_match = .) %>% 
      subset_table(df_mat, rows=.)
    
    # Transform data
    if(is.null(input$transf_checkbox))
      return(df_s)
    # print(input$transf_checkbox)
    transform_data(df_s, input$transf_checkbox)
  })
  
  # Color of heatmap
  hm_color <- reactive({
      get_col_pal(brew_pal = input$hm_color_pal, 
                  rev =  ifelse(input$hm_color_rev == "No", F, T)) %>%
        get_col_grad(n=50) %>% 
      return
  })
  
  # Annotation colors
  ann_colors <- list(cellType=c("BC"="red", "LC"="darkblue", "LP"="skyblue"), 
                     patientAgeGroup=get_element_colors(v =  df_col[,"patientAgeGroup"] %>% unique, colRamp = get_col_pal("Greys"), rearr = F),
                     hormoneStatus = c("Luteal"="darkgoldenrod1", "Post.Menopausal"="darkgoldenrod4", "Follicular"="darkgoldenrod"))
  
  # Make heatmap
  heatMap <- reactive({
    pheatmap(mat = df_mat_s(), 
             scale = tolower(input$hm_scale),
             main = input$file_label,
             # Colors
             color = hm_color(),
             # Cell size
             cellwidth = ifelse(input$hm_cell_w == "", NA, input$hm_cell_w),
             cellheight = ifelse(input$hm_cell_h == "", NA, input$hm_cell_h),
             # Annotations
             annotation_row = NA, #TODO
             annotation_col = df_col,
             annotation_colors = ann_colors,
             # Rows
             cluster_rows = check(input$hm_row, "cluster"),
             show_rownames = check(input$hm_row, "show_name"),
             # Columns
             cluster_cols = check(input$hm_column, "cluster"),
             show_colnames = check(input$hm_column, "show_name"),
             # Other preferences
             border_color = NA)
  })
  # Render output
  output$hm <- renderPlot({
    heatMap()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0(input$file_label, '.', input$file_format) },
    content = function(file) {
      ggsave(file, plot = heatMap(), device = input$file_format)
    }
  )
  
    
  # Downloadable csv of processed dataset ----
  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste0("MEC_data", ".csv")
    },
    
    content = function(file) {
      df <- cbind(df_col, t(df_mat_s()))
      write.csv(df, file, row.names = TRUE)
    }
  )
}

# Run app 
shinyApp(ui = ui, server = server)

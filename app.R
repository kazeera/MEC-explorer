# Load libraries
library(shiny)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(Seurat)
library(viridis) #color palettes
library(cowplot) #plot_grid
#' Acronyms:
#' hm = heatmap
#' ss = subset
#' ann = annotation

# Load helper functions
lapply(list.files(pattern = "helpers"), source)

# Import proteome data
hprot_data <- readRDS("data/human_proteome.rds")

# Import single cell data
cell_colors <- readRDS("data/CellType_colors.rds")
# sc_data <- readRDS("data/Seurat_MEC_object.rds")

# Accessible colors for single cell plots
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
sc_colors  <- list(Gray="#EDEDED", Black="#000000", Orange="#E69F00", SkyBlue="#56B4E9", Green="#009E73", 
                   Yellow="#F0E442", Blue="#0072sc_data", Vermillion="#D55E00", Pink="#CC79A7")
# scales::show_col(unlist(color_blind_pal))
hm_colors <- rownames(RColorBrewer::brewer.pal.info)

# Color palettes
ui <- navbarPage("Khokha Lab MEC Explorer",
                 
                 # Tab 1
                 tabPanel("Human Proteome",  
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            # Sidebar panel for inputs ----
                            sidebarPanel( 
                              # CSS
                              id="scPanel", 
                              style="overflow-y:scroll; max-height:800px; position=absolute",
                              
                              helpText("Query human mammary epithelial cell (MEC) proteome (first published in Mahendralingam et al., Nature Metabolism, 2021)"),
                              textInput("file_label", h4("Title for Plot:"),  value = ""),
                              
                              # Make drop down menu with all gene names in single cell object
                              selectizeInput("hm_geneIDs", multiple = T, label = h4('Select proteins of interest'), choices = NULL, options = list(create = FALSE)),
                              
                              h4("Data Transformations"),
                              checkboxGroupInput("transf_checkbox", label = "", choices = list("z-score (row)" = "z_row", "z-score (column)" = "z_column"), selected = 0),
                              
                              # Heatmap visualization options
                              h4("Heatmap Preferences"),
                              # Scaling
                              selectInput("hm_scale", label = "Scale by: ", choices = c("Row", "Column", "None"), selected = "Row"),
                              # Clustering/showing name
                              checkboxGroupInput("hm_row", label = "Rows", choices = list("cluster" = "cluster", "show_name" = "show_name"), selected = c("show_name", "cluster")),
                              checkboxGroupInput("hm_column", label = "Columns", choices = list("cluster" = "cluster", "show_name" = "show_name"), selected = 0),
                              # Cell size
                              numericInput("hm_cell_h", h6("Cell height"), value = "", width = 100),
                              numericInput("hm_cell_w", h6("Cell width"), value = "", width = 100),
                              
                              # Gradient color picker
                              selectInput("hm_color_pal", label = "Pick heatmap color", choices = hm_colors, selected = hm_colors[3]),
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
                              plotOutput(outputId = "hm") # hm = heat map
                            ))),
                 # Tab 2
                 tabPanel("Human scRNA-seq",
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            # Sidebar panel for inputs ----
                            sidebarPanel( 
                              # CSS
                              id="hmPanel", 
                              style="height:800px;",
                              helpText("Query human mammary epithelial cell (MEC) single cell RNA seq data (first published in Mahendralingam et al., Nature Metabolism, 2021)"),
                              
                              # Make drop down menu with discrete meta data variables
                              selectInput("sc_group.by", label = h5("Select Annotation"), selected = "CellType", choices = c("CellType","CellStates","Seurat_clusters","Seurat_cellCyclePhase","Source","Patient")),# colnames(sc_data@meta.data)),
                              
                              # Make drop down menu with all gene names in single cell object
                              selectizeInput("sc_geneIDs", multiple = T, label = h5('Select genes'), choices = NULL, options = list(create = FALSE)), # if TRUE, allows newly created inputs),
                              
                              # Point size
                              numericInput("sc_pt.size", h6("Point Size"), value = "0.1", min = 0.0001, max=10, step = 0.1, width = 100),
                              
                              # # Gradient color picker
                              # selectInput("sc_color_low", label = "Select color (low)", choices = names(sc_colors), selected = "Gray"),
                              # selectInput("sc_color_high", label = "Select color (high)", choices = names(sc_colors), selected = "Vermillion"),
                              
                              # Save Plot
                              selectInput("file_format_sc", label = "Save plot as", choices = c("png", "jpeg", "pdf", "svg", "tiff"), selected = "png"),
                              downloadButton('downloadPlot_sc', 'Download Plots'),
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              plotOutput(outputId = "sc") # sc = single cell plot
                            )))
)


# Define server logic ----
server <- function(input, output, session) {
  ## Proteomics ----------------------------------------
  # For drop down menu with gene names, update and filter server-side to speed up program
  updateSelectizeInput(session, 'hm_geneIDs', 
                       choices = rownames(hprot_data$mat), 
                       selected = c("CD14","KRT18","KRT19","EGFR","EPCAM","VIM","FOXA1","LTF","GATA3","TP63","ITGA6","FOLR1"),
                       server = TRUE)
  
  # Process data
  df_mat_s <- reactive({
    df_s <-
      get_to_keep(rownames(hprot_data$mat), list_to_match = input$hm_geneIDs) %>%
      subset_table(hprot_data$mat, rows=.)
    # Transform data
    if(is.null(input$transf_checkbox))
      return(df_s)
    # print(input$transf_checkbox)
    transform_data(df_s, input$transf_checkbox)
  })
  
  # Color of heatmap
  hm_color <- reactive({
    get_col_pal(brew_pal = input$hm_color_pal, rev =  ifelse(input$hm_color_rev == "No", F, T)) %>%
      get_col_grad(n=50) %>% 
      return
  })
  
  # Annotation colors
  ann_colors <- list(cellType=c("BC"="red", "ML"="darkblue", "LP"="skyblue"), 
                     patientAgeGroup=get_element_colors(v = hprot_data$col[,"patientAgeGroup"] %>% unique, colRamp = get_col_pal("Greys"), rearr = F),
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
             annotation_col = hprot_data$col,
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
      paste0("human_MEC_proteome", ".csv")
    },
    
    content = function(file) {
      df <- cbind(hprot_data$col, t(df_mat_s()))
      write.csv(df, file, row.names = TRUE)
    }
  )
  
  ## scRNA ----------------------------------------
  sc_data <- readRDS("data/sc_Seurat_object.rds")
  
  # For drop down menu with gene names, update and filter server-side to speed up program
  updateSelectizeInput(session, 'sc_geneIDs', choices = rownames(sc_data), selected = c("KRT18","VIM"), server = TRUE)
  
  # Umap for annotations
  umap1 <- reactive({
    DimPlot(sc_data, pt.size = input$sc_pt.size, cols= turbo(length(unique(unlist(sc_data[[input$sc_group.by]])))), group.by = input$sc_group.by)
  })
  
  # Umap for features
  umap2 <- reactive({
    FeaturePlot(sc_data, features = input$sc_geneIDs, pt.size = input$sc_pt.size, reduction = "umapharmony", cols = viridis(3)) #c(sc_colors[[input$sc_color_low]], sc_colors[[input$sc_color_high]]),
  })
  
  # Stacked violin plot for features
  vln <- reactive({
    VlnPlot(sc_data, features = input$sc_geneIDs, stack = T, group.by = input$sc_group.by)
  })
  
  # Heatmap for features
  sc_hm <- reactive({
    DoHeatmap(sc_data, features = input$sc_geneIDs, group.colors = turbo(length(unique(unlist(sc_data[[input$sc_group.by]])))), group.by = input$sc_group.by) + 
      viridis::scale_fill_viridis() + 
      theme(axis.text.y = element_text(face="bold",size=15, color = "black"))# + guides(fill="none") # removes color bar
  })
  
  # Render output
  output$sc <- renderPlot({
    cowplot::plot_grid(umap1(), vln(), umap2(), sc_hm())
  })
  
  # Download plot
  output$downloadPlot_sc <- downloadHandler(
    filename = function() { paste0('MEC_single_cell.', input$file_format_sc)},
    content = function(filename) {
      ggsave(filename, plot = plot_grid(umap1(), vln(), umap2(), sc_hm()), device = input$file_format_sc)#,  width = units(8.5, "in"), height = units(5, "in"), device = input$file_format_sc)
    }
  )
}

# Run app 
shinyApp(ui = ui, server = server)

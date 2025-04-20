# 1) Packages
library(SeuratData)
library(shiny)
library(shinydashboard)
library(Seurat)
library(patchwork)
library(ggpubr)
library(viridis)
library(DT)
library(shinycssloaders)

# 2) Load & preprocess PBMC3K once on app start
InstallData("pbmc3k")
pbmc3k <- LoadData("pbmc3k")
pbmc3k <- subset(pbmc3k, subset = seurat_annotations != "NA")
pbmc3k <- SCTransform(pbmc3k)
pbmc3k <- RunPCA(pbmc3k, features = VariableFeatures(pbmc3k))
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:30)
pbmc3k <- FindClusters(pbmc3k)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:30)
Idents(pbmc3k) <- pbmc3k$seurat_annotations
markers <- FindAllMarkers(pbmc3k, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# UI
ui <- dashboardPage(skin = "purple",
  dashboardHeader(title = "PBMC3K Explorer"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Explore", tabName = "explore", icon = icon("chart-bar")),
      menuItem("Markers", tabName = "markers", icon = icon("table"))
    ),
    br(),
    textInput(
      "features",
      "Gene(s) or feature(s):",
      placeholder = "e.g. CD3D, MS4A1"
    ),
    actionButton(
      "makePlot",
      "Generate plots",
      icon = icon("play"),
      class = "btn-primary"
    )
  ),
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css?family=Roboto"),
      tags$style(HTML("
        /* Global font & bg */
        body { font-family: 'Roboto', sans-serif; }
        .content-wrapper { background-color: #f7f7f7; }

        /* Header & sidebar */
        .main-header .logo, .main-header .navbar { background-color: #6f42c1 !important; }
        .sidebar { background-color: #3f0f3f; }
        .sidebar .sidebar-menu > li.active > a { background-color: #2c0d2c !important; }

        /* Boxes */
        .box {
          box-shadow: 0 4px 10px rgba(0,0,0,0.15);
          border-radius: 8px;
          border: none;
        }

        /* Primary action button */
        .btn-primary {
          background-color: #6f42c1;
          border-color: #6f42c1;
          border-radius: 4px;
        }
        .btn-primary:hover {
          background-color: #59328f;
          border-color: #59328f;
        }
      "))
    ),
    tabItems(
      tabItem(
        tabName = "explore",
        fluidRow(
          box(
            title = "FeaturePlot",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            withSpinner(plotOutput("featurePlot", height = "auto"), type = 6)
          )
        ),
        fluidRow(
          box(
            title = "ViolinPlot",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            withSpinner(plotOutput("violinPlot", height = "auto"), type = 6)
          )
        ),
        fluidRow(
          box(
            title = "DotPlot",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            withSpinner(plotOutput("dotPlot", height = "auto"), type = 6)
          )
        )
      ),
      tabItem(
        tabName = "markers",
        fluidRow(
          box(
            title = "Marker Genes",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            DTOutput("markersTable")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Publication-ready theme
  theme_set(theme_classic(base_size = 14))
  update_geom_defaults("point", list(size = 2))

  # Markers table
  output$markersTable <- renderDT({
    datatable(
      markers,
      filter = "top",
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # Parse features
  feats <- eventReactive(input$makePlot, {
    req(input$features)
    trimws(unlist(strsplit(input$features, ",")))
  })

  # FeaturePlot
  output$featurePlot <- renderPlot({
    req(feats())
    plots <- FeaturePlot(
      pbmc3k,
      features = feats(),
      cols = c("lightgrey", "blue", "red"),
      pt.size = 1,
      combine = FALSE
    )
    wrap_plots(plots, ncol = min(length(plots), 3)) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
  }, height = function() {
    ceiling(length(feats()) / 3) * 300
  })

  # ViolinPlot
  output$violinPlot <- renderPlot({
    req(feats())
    VlnPlot(
      pbmc3k,
      features = feats(),
      pt.size = 0
    ) +
      stat_compare_means(aes(label = ..p.signif..)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }, height = function() {
    200 + 50 * length(feats())
  })

  # DotPlot
  output$dotPlot <- renderPlot({
    req(feats())
    DotPlot(
      pbmc3k,
      features = feats(),
      dot.scale = 8
    ) +
      RotatedAxis() +
      scale_color_viridis(option = "D") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }, height = function() {
    200 + 50 * length(feats())
  })
}

# Run app
shinyApp(ui, server)
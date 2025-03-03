library(shiny)
library(DT)
library(ggplot2)
setwd("/home/lera/Documents/2025-DAAD/Shiny/")

data_path0<-"0171619.1.1/0171619/1.1/data/0-data/"
headerphoto<-"NCEI 0171619.jpg"
data_path<-"0171619.1.1/0171619/1.1/data/1-data/"
datasets<-c("","attendance v1-0-DATA.csv","interval v1-0-DATA.csv","repro v1-1-DATA.csv","validation v1-1-DATA.csv")
metadatas<-c("","attendance v1-0-READ ME.txt","interval v1-0- READ ME.txt","repro v1-1-READ ME.txt","validation v1-1-READ ME.txt")
names(metadatas) <- datasets
page.title<-"Estimating nest-level phenology and reproductive success of colonial seabirds using time-lapse cameras"
paper.link<-"https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13015"

ui <- fluidPage(
  # Set theme
  theme = bslib::bs_theme(bootswatch = "darkly"),
  # Set title
  titlePanel(tags$a(page.title,href=paper.link)),
  # Set text input box
  textInput("wd", "What is your working directory?"),
  # Divide page into tabs
  tabsetPanel(
    # Set 1st tab
    tabPanel("Data",
      # Set row of panels
      fluidRow(
        # Set subrow column panels
        column(3,
               # Set drop down list
               selectInput("dataset", "Dataset", choices = datasets),
               # Set check box
               checkboxInput("show_metadata", "Show Metadata", value = FALSE)
          ),
        # Set subrow column panels
        column(9,
               # Set hidden panel
               conditionalPanel(
                 condition = "input.show_metadata == true",
                 verbatimTextOutput("metadata")
              ) 
        )
      ),
      # Set table
      DT::DTOutput("table"),
    ),
    
    # Set 2ed tab
    tabPanel("Plots",
      # Set file input
      fileInput("upload", "Load species_locality.csv", accept = c(".csv")),
      # Set page row
      fluidRow(
        # Set subrow column panel
        column(6,
               # Set plot
               plotOutput("species_locality_plot", width = "600px",height = "700px", click = "plot_click")
        ),
        column(6,
               verbatimTextOutput("locality_info")
        )
      )
    ),
    
    # Set 3ed tab
    tabPanel("Galery",
      # Divide panel in two
      sidebarLayout(
          sidebarPanel(
            # Set button
            actionButton("show_img", "Show me penguins!", class = "btn-lg btn-success")
          ),
          # Set conditional image
          mainPanel(
            uiOutput("image_ui")
          )
      )
    )
  )
)



server <- function(input, output, session) {
  
  # Save input working directory
  wd<-reactive({
    req(input$wd)
  })
  
  # Open and process selected file
  dataset <- reactive({
    # Ensure input is not NULL
    req(input$dataset)  
    # Load data
    df <- read.csv(file.path(wd(),data_path, input$dataset))
    # Drop columns with no values
    df <- df[,colSums(is.na(df))<nrow(df)]
    #Fill NA
    df[is.na(df)] <- "NA"
    # Convert only character columns to factors
    char_cols <- sapply(df, is.character)
    df[char_cols] <- lapply(df[char_cols], as.factor)
    # Convert string to date type
    if ("DATE" %in% names(df)){
      df$DATE<-as.Date(as.character(df$DATE), format = "%d/%m/%Y")
    }
    # Give back the data
    return(df)
  })
  
  # --------- Display data in table ----------
  # Display table
  output$table <- renderDT({
    datatable(dataset(), filter = "top", options = list(pageLength = 10))
  })
  
  # Get metadata
  output$metadata <- renderText({  
    # Ensure input is not NULL
    req(input$dataset)  
    # Assemble file path
    fileName <- file.path(wd(),data_path, metadatas[input$dataset])  
    # Read entire file as text
    if (file.exists(fileName)) {
      return(paste(readLines(fileName), collapse = "\n")) 
    } else {
      return("Metadata file not found.")
    }
  })
  
  # ---------- Display interactive plot ----------
  # Open file to plot
  species_locality_data <- reactive({
    # Ensure input is not NULL
    req(input$upload)
    # Get file extention
    ext <- tools::file_ext(input$upload$name)
    # Validate .csv file extention
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ","),
           validate("Invalid file; Please upload a .csv file")
    )
    # Load the data
    df <- read.csv(paste0(wd(),input$upload$name))
    # Transform string to numerical 
    df$Latitude <- as.numeric(gsub("−", "-", df$Latitude))
    df$Longitude <- as.numeric(gsub("−", "-", df$Longitude))
    # Give back the data
    return(df)
  })
  
  # Plot 
  output$species_locality_plot <- renderPlot({
    ggplot(species_locality_data(), aes(x = Longitude, y = Latitude, shape = Species, color = Species)) +
      geom_point(size = 4, position = "jitter") +  # Adjust point size
      scale_shape_manual(values = c("Adélie" = 16, "Gentoo" = 17, "Chinstrap" = 15)) +
      scale_color_manual(values = c("Adélie" = "black", "Gentoo" = "red", "Chinstrap" = "blue")) +
      theme_minimal() +
      labs(title = "Penguin Localities", x = "Longitude", y = "Latitude") +
      theme(legend.title = element_blank())
  })
  
  # Render text when clicked on the plot
  output$locality_info <- renderText({
    # Ensure that the click input is available
    req(input$plot_click)  
    # Get the closest point to the clicked coordinates
    click_data <- nearPoints(species_locality_data(), input$plot_click, threshold = 5, maxpoints = 1)
    # Give back data in text box
    if (nrow(click_data) > 0) {
      return(paste("Locality: ", click_data$Locality))
    } else {
      return("Click on a point to see the locality.")
    }
  })
  

  # ---------- Display local photo on demand ----------
  # Set default value
  image_displayed <- reactiveVal(FALSE)
  
  observeEvent(input$show_img, {
    # Toggle the image display state on button click
    image_displayed(!image_displayed())
  })
  
  # Render the image dynamically when the button is clicked
  output$image_ui <- renderUI({
    req(image_displayed())  # Make sure the image display state is TRUE
    
    # If image_displayed is TRUE, show the image
    if (image_displayed()) {
      # Provide the path to your local image
      img_path <- normalizePath(paste0(wd(),data_path0,headerphoto))
      
      # Use renderImage to show the image
      renderImage({
        list(src = img_path, width = 692*0.4,
              height = 922*0.4, alt = "penguins")
        }, deleteFile = FALSE)
    }
  })
}

shinyApp(ui, server)

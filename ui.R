# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com

# version : 0.1.3
# since: DEC.03.2015
# last fix : JAN.11.2016

# version 0.1.2
# connect server.R with ui.R and read tab deletion file form local storage
# it can be draw basic heatmap

# version 0.1.3
# back to bersion 0.1.2


###############################################################################
library(shiny)

## version 0.1.3

# Define UI for dataset viewer application
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Draw Heatmap by the Shiny"),
  
  # Sidebar with controls to select a dataset and specify the
  # number of observations to view
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose CSV File',
                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
      tags$hr(),
      
      numericInput("obs", "Number of observations to view:", 20),
      textInput("maintitle", label = h6("Enter main title"), value = ""),
      textInput("ylab", label = h6("Enter ylab's title"), value = ""),
      checkboxInput('fixed', 'Fix column', FALSE),
      checkboxInput('header', 'Header', TRUE),
      checkboxInput('rowname', 'Row name', TRUE),
      radioButtons('sep', 'Separator',
                   c("Comma (if you load the data with CSV check here)" =',',
                     Semicolon=';',
                     Tab='\t'),
                   'Tab'),
      radioButtons('quote', 'Quote',
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   'Double Quote'),
      radioButtons('color', 'heatmap Color',
        c( 'Spectral' = 1,
           'YlOrRd' = 2,
           'Greys'= 3,
         'RdYlBu' = 4),
           'RdYlBu')
      
      #downloadButton('downloadData', 'Download')
    ),
    
    # Show a summary of the dataset and an HTML table with the
    # requested number of observations
    mainPanel(
      # with grey box "summary" information is here
      #verbatimTextOutput("summary"),
      
      # view "view" data in here
      #tableOutput("view")
      
      # make tab panel in the main panel
      tabsetPanel(
        tabPanel("Contants", tableOutput('contents')),
        tabPanel("Heatmap", plotOutput("heatmap", height = 2400, width = 1600)),
        tabPanel("Help", includeHTML("help.html"))
        #        tabPanel("View", tableOutput("view"))
      )
      
    )
  )
))
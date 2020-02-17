# FERMI PROCESS MODEL for PUBLIC GOODS GAME
# << UI >>
# by Aydin Mohseni

# Load the shiny GUI library
library(shiny)
library(ggplot2)
library(ggthemes)

# Set encoding for special characters
Sys.setlocale("LC_ALL", "fr_FR.UTF-8")

# Define UI for application
shinyUI(fluidPage(
  
  # CSS for visual
  includeCSS("www/style.css"),
  
  # Title
  titlePanel("Tragedies of the Commons"),
  
  # Load MathJax 
  withMathJax(),
  fluidRow(
    style = "background-color:#F2F2F2; margin-top: 30px; margin-bottom: 30px; padding: 10px", 
    column(
      width = 4,
      # Introduction text:
      p(
        tags$b("Description"),
        "."
      )
    ),
    column(
      width = 4,
      # Introduction text:
      p(
        "Lipsum."
      )
    ),
    column(
      width = 4,
      p(
        tags$b("Result:"),
        "Dolor."
      )
    )
  ),
  
  # Sidebar for Parameter Input
  sidebarLayout(
    
    sidebarPanel(
      
      # Use MathJax
      withMathJax(),
      
      # Group size
      sliderInput("groupSize",
                  "Group size:",
                  min = 1,
                  max = 100,
                  value = 50),
      
      # Perceived risk of disaster
      sliderInput("perceivedRiskOfDisaster",
                  "Perceived risk of disaster:",
                  min = 0,
                  max = 1,
                  value = 0.25),
      
      # Perceived cost of disaster
      sliderInput("perceivedCostOfDisaster",
                  "Perceived cost of disaster:",
                  min = 0,
                  max = 1,
                  value = 0.75),
      
      # Cost of cooperation
      sliderInput("costOfCooperation",
                  "Cost of cooperation:",
                  min = 0,
                  max = 1,
                  value = 0.2),
      
      # Critical fraction of cooperators required to avoid disaster risk
      sliderInput("criticalFraction",
                  "Critical fraction of cooperators required to avoid disaster risk:",
                  min = 0,
                  max = 1,
                  value = 0.8),
      
      # Selection coefficient
      sliderInput("selection",
                  "Selection coefficient:",
                  min = 0,
                  max = 10,
                  value = 1),
      
      # Mutation rate
      sliderInput("mutation",
                  "Mutation rate:",
                  min = 0,
                  max = 0.5,
                  value = 0.05)
      
    ),
    
    # Main Panel with Stationary Distribution + Simulation & Stats Panels
    mainPanel(
      style = "padding: 10px; margin-bottom: 10px",
      plotOutput("payoffsPlot", height = "350px"),
      tags$br(),
      plotOutput("stationaryDistributionPlot", height = "350px"),
      tags$br(),
      plotOutput("simulationPlot", height = "350px")
    )
  )
))
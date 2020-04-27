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
  titlePanel("AI Tragedy of the Commons"),
  
  # Load MathJax 
  withMathJax(),
  fluidRow(
    style = "background-color:#F2F2F2; margin-top: 30px; margin-bottom: 30px; padding: 10px", 
    column(
      width = 6,
      # Introduction text:
      p(
        tags$b("Description"),
        "A population of size \\(Z\\) is organized into groups of size \\(N\\). Each individual can be understood as an AI researcher able to realize or ignore proposed norms of research and practice. Each group can be thought of as an organizational unit, such as a laboratory, university department, or collaborative, able to sustain shared norms.",
        tags$br(), tags$br(),
        "Each individual has an initial endowment \\(b\\), viewed as the asset value at stake. In the simplest scenario, participants choose either to cooperate (\\(C\\)) by adhering to the proposed norm, and so to contribute a fraction \\(cb\\) of their endowment to realizing adherence, or to defect (\\(D\\)), and contribute nothing. Contributing produces a cost in extra effort (e.g., taking additional precautions to ensure their research and business practices adhere to the proposed norm), or in restricted opportunities (e.g., by refusing to take on lucrative projects that would violate the proposed norm)."
      )
    ),
    column(
      width = 6,
      # Introduction text:
      p(
        "A successful agreement is reached if the overall number of contributors to the public group exceeds a certain threshold \\(n^*cb\\). In that case, all participants will keep whatever they have. Otherwise, with probability \\(r\\) corresponding to the risk of disaster if agreement is not attained, everyone in the group will lose a portion \\(d \\in [0,1]\\) of the payoff she has accrued. That is, if \\(d=1\\), individuals risk losing all of their accrued payoffs if they fail to come to an agreement, while \\(d=0\\) represents the case where failure to reach an agreement poses no risk of negative consequences.",
        tags$br(), tags$br(),
        "The payoffs to each action in a group of size \\(N\\) when \\(n_C\\) individuals are cooperating are given by
        $$ \\pi_C(n_C) = b \\Theta (n_C-n^*) + (1-r) b (1-d \\Theta(n_C-n^*)) - cb,$$
        $$ \\pi_D(n_C) = \\pi_D(n_C) + cb $$
        where \\(\\Theta\\) is the heaviside step function."
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
                  "Group size (\\(N\\)):",
                  min = 1,
                  max = 20,
                  value = 6),
      
      # Perceived risk of disaster
      sliderInput("perceivedRiskOfDisaster",
                  "Risk of disaster (\\(r\\)):",
                  min = 0,
                  max = 1,
                  value = 0.8),
      
      # Perceived cost of disaster
      sliderInput("perceivedCostOfDisaster",
                  "Cost of disaster (\\(d\\)):",
                  min = 0,
                  max = 1,
                  value = 0.8),
      
      # Cost of cooperation
      sliderInput("costOfCooperation",
                  "Cost of cooperation (\\(c\\)):",
                  min = 0,
                  max = 1,
                  value = 0.1),
      
      # Critical fraction of cooperators required to avoid disaster risk
      sliderInput("criticalFraction",
                  "Critical fraction of cooperators required to avoid disaster risk (\\(n^*\\)):",
                  min = 0,
                  max = 1,
                  value = 0.6),
      
      # Selection coefficient
      sliderInput("selection",
                  "Selection coefficient (\\(\\lambda\\)):",
                  min = 0,
                  max = 10,
                  value = 5),
      
      # Mutation rate
      sliderInput("mutation",
                  "Mutation rate (\\(\\mu\\)):",
                  min = 0,
                  max = 0.5,
                  value = 0.1)
      
    ),
    
    # Main Panel with Stationary Distribution + Simulation & Stats Panels
    mainPanel(
      style = "padding: 10px; margin-bottom: 0px",
      plotOutput("selectionGradientPlot", height = "300px"),
      plotOutput("payoffsPlot", height = "300px"),
      plotOutput("stationaryDistributionPlot", height = "300px")
    )
  )
))
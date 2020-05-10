# FERMI PROCESS MODEL for PUBLIC GOODS GAME
# << SERVER >>
# by Aydin Mohseni

# Load packages
library(shiny)
library(ggplot2)
library(expm)

options(shiny.sanitize.errors = FALSE)

Z <- 100 # The size of the population

# Define server logic
shinyServer(function(input, output, session) {

  # PROCESS: Compute the desired quantities:
  # (1) Payoffs, 
  # (2) Fermi Process Transition Matrix,
  # (3) Stationary Distribution.
  computeDynamics <- reactive({
    
    # GLOBAL VARIABLES
    # Import the required parameters fro the UI
    N <- as.numeric(input$groupSize) # The size of groups 
    r <- as.numeric(input$perceivedRiskOfDisaster) # The risk of disaster occuring in the absence of cooperation succeeding
    d <- as.numeric(input$perceivedCostOfDisaster) # The cost of disaster (in terms of te proportion of endownments lost) if cooperation fails and disaster occurs
    b <- 1 # The initial endowment of each individual
    c <- as.numeric(input$costOfCooperation) # The proportion of the endowment contributed by cooperators
    p_crit <- as.numeric(input$criticalFraction) # The critical fraction of cooperators at which cooperation succeeds
    n_crit <- N * p_crit # The critical number of cooperators at which cooperation succeeds
    R <- as.numeric(input$selection) # The rationlaity coefficient (0:= random selection; 1:= replicator dynamics; ∞ := best reponse dynamics)
    M <- as.numeric(input$mutation) # Mutation/error/noise parameter (portion of the time a random strategy is chosen)
    
    # PAYOFF FUNCTIONS    
    # Define the payoffs to each cooperators and defectors
    # as a function of the number of cooperators n_C in their group
    payoff_cooperator <-
      function(n_C) {
        b * (n_C - n_crit >= 0) + b * (1 - r *d) * (n_C - n_crit < 0) - (c * b)
      }
    payoff_defector <-
      function(n_C) {
        payoff_cooperator(n_C) + (c * b)
      }
    # Now, define the _average_ payoff to each type
    # as a function of the group sizes N
    # and the fraction of each type in the broader population n_C_pop
    mean_payoff_cooperator <-
      function(n_C_pop) {
        # Note: If there are no cooperators, just return 0
        if (n_C_pop == 0) {
          z <- 0
        } else {
          # Find the fraction of cooperators in the group
          x <- n_C_pop / Z
          # Compute the vector of (binomially distributed) probabilities for
          # each possible group composed of k cooperators and (N - k) defectors
          p <- dbinom(x = 0:N, size = N, prob = x)
          # Compute the vector of conditional probabilities (k/N) of being a cooperator in each combination
          q <- 0:N / N
          # Compute the payoffs to cooperators in each group
          y <- sapply(c(0:N), payoff_cooperator)
          # Find the expected payoff for cooperators across all groups
          z <- sum(y * q * p) / sum(p * q)
        }
        return(z)
      }
    mean_payoff_defector <-
      function(n_C_pop) {
        # Note: If there are no defectors, just return 0
        if (n_C_pop == Z) {
          z <- 0
        } else {
          # Find the fraction of cooperators in the group
          x <- n_C_pop / Z
          # Compute the vector of (binomially distributed) probabilities for
          # each possible group composed of k cooperators and (N - k) defectors
          p <- dbinom(x = 0:N, size = N, prob = x)
          # Compute the vector of conditional probabilities ((N-k)/N) of being a defector in each combination
          q <- N:0 / N
          # Compute the payoffs to defectors in each group
          y <- sapply(c(0:N), payoff_defector)
          # Find the expected payoff for defectors across all groups
          z <- sum(y * q * p) / sum(p * q)
        }
        return(z)
      }
    
    # Create payoff data frame
    payoff_cooperator_vec <- sapply(seq(from = 0, to = Z, by = 1), mean_payoff_cooperator)
    payoff_defector_vec <- sapply(seq(from = 0, to = Z, by = 1), mean_payoff_defector)
    PayoffsDF <- data.frame(
      N = rep(seq(
        from = 0, to = 1, by = 1 / Z
      ), times = 2),
      Payoff = c(payoff_cooperator_vec, payoff_defector_vec),
      Strategy = rep(c("Cooperate", "Defect"), each = Z + 1)
    )
    
    # Specify the transition probabilities of the population changing
    # from a state with n_C_pop cooperators in the population 
    # to one with n_C_pop + 1, n_C_pop - 1, or n_C_pop cooperators.
    
    # Pr(n_C_pop -> n_C_pop + 1)
    Prob_n_C_Increase <-
      function(k) {
        if (k < Z) {
          z <- (k / Z) * (Z - k) / (Z - 1) * (1 - M) * (1 + exp(R * (
            mean_payoff_defector(k) - mean_payoff_cooperator(k)
          ))) ^ -1 + M / 2
        } else {
          z <- 0
        }
        return(z)
      }
    # Pr(n_C_pop -> n_C_pop - 1)
    Prob_n_C_Decrease <- function(k) {
      if (k > 0) {
        z <- (k / Z) * (Z - k) / (Z - 1)  * (1 + exp(R * (
          mean_payoff_cooperator(k) - mean_payoff_defector(k)
        ))) ^ -1 + M / 2
      } else {
        z <- 0
      }
      return(z)
    }
    # Pr(n_C_pop -> n_C_pop)
    Prob_n_C_Stay <- function(k) {
      1 - Prob_n_C_Increase(k) - Prob_n_C_Decrease(k)
    }
    
    # TRANSITION MATRIX
    # We compute the transition matrix by recursive applications of the transition probabilities.
    # Create an empty transition matrix.
    FP <- matrix(
      data = NA,
      nrow = Z + 1,
      ncol = Z + 1,
      byrow = TRUE
    ) 
    # Fill the transition matrix
    for (i in 1:(Z + 1)) { # For each row 
      for (j in 1:(Z + 1)) {# For each column
        if (j == (i - 1)) { # If the state has one fewer Cooperator
          FP[i, j] <- Prob_n_C_Decrease(i - 1)
        }
        if (j == i) { # If the state is the same
          FP[i, j] <- Prob_n_C_Stay(i - 1)
        }
        if (j == (i + 1)) { # If the state has one more Cooperator
          FP[i, j] <- Prob_n_C_Increase(i - 1)
        }
        if ((j != (i - 1)) & (j != i) & (j != (i + 1))) {
          FP[i, j] <- 0
        }
      }
    }
    
    # SELECTION GADIENT
    # Now, use the transition matrix to compute the gradient of selection 
    # defined as G(i) = Prob_n_C_Increase(i) - Prob_n_C_Decrease(i).
    Prob_n_C_Increase_vec <- rep(0, times = (Z - 1))
    Prob_n_C_Decrease_vec <- rep(0, times = (Z - 1))
    for (i in 2:Z) {
      Prob_n_C_Increase_vec[i - 1] <- FP[i, i + 1]
    }
    for (i in 2:Z) {
      Prob_n_C_Decrease_vec[i - 1] <- FP[i, i - 1]
    }
    selectionGradient <- Prob_n_C_Increase_vec - Prob_n_C_Decrease_vec
    
    # STATIONARY DISTRIBUTION  
    # Note: We know that the stationary distribution µ exists because the addition of mutation makes the Markov process Ergodic.
    # That is: Where the Fermi Process is already finite, and aperiodic,
    # it is now irreducible (has only one recursive class) as every state is reachable from any other,
    # and every state will be visited an infinite number of times in the limit.
    # In the absence of mutation, there are were two (singleton) recursive classes
    # corresponding to the two absorbing states where the population is all cooperators or all defectors.
    # It follows from being ergodic that the limit distribution π is independent of any initial distribution.
    
    # We compute the stationary distribution of the process
    # First, we calculate µ_0 = ( Σ^{N}_{k=0} Π^{k}_{i=1} PforwardM_{i-1} / PbackM_{i} ) ^ -1
    Mu0Vector <- c() # Empty vector, in which to store each product
    # Empty matrix, in which to store each element of the products
    Mu0Matrix <- matrix(data = 1,
                        nrow = Z + 1,
                        ncol = Z + 1)
    for (k in 2:(Z + 1)) {
      for (i in 2:k) {
        Mu0Matrix[k, i - 1] <- FP[i - 1, i] / FP[i, i - 1]
      }
    }
    # Take the product of the rows of Mu0Matrix
    for (i in 1:(Z + 1)) {
      Mu0Vector[i] <- prod(Mu0Matrix[i,])
    }
    # Compute µ_0
    Mu0 <- sum(Mu0Vector) ^ -1
    
    # Now, we calculate µ_k = µ_0 * Π^{k}_{i=1} PforwardM_{i-1} / PbackM_{i} )
    MuVector <- c()
    # Empty matrix, in which to store each element of the products
    MuMatrix <- matrix(data = 1,
                       nrow = Z + 1,
                       ncol = Z + 1)
    for (k in 2:(Z + 1)) {
      for (i in 2:k) {
        MuMatrix[k, i - 1] <- FP[i - 1, i] / FP[i, i - 1]
      }
    }
    for (i in 1:(Z + 1)) {
      MuVector[i] <- Mu0 * prod(MuMatrix[i,])
    }
    
    # OUTPUT the results of our computations 
    # to be accessed by other reactive contexts
    return(
      list(
        PayoffsDF,
        FP,
        MuVector,
        selectionGradient
      )
    )
    
  }) # END COMPUTATIONS 
  
  # OUTPUT: Payoff functions plot
  output$payoffsPlot <- renderPlot({
    
    # Import relevant variables
    N <- as.numeric(input$groupSize)
    PayoffsDF <- computeDynamics()[[1]]
    
    # Trim {0, 1} endopoints of C & D payoffs (for aesthetics)
    PayoffsDF[which(PayoffsDF$N == 0 & PayoffsDF$Strategy == "Cooperate"), 2] <- NA
    PayoffsDF[which(PayoffsDF$N == 1 & PayoffsDF$Strategy == "Defect"), 2] <- NA
    
    # Plot payoff functions
    plot_PayoffsDF <- ggplot(data = PayoffsDF, aes(x = N, y = Payoff, group = Strategy)) +
      geom_line(aes(color = Strategy), size = 2) + 
      scale_color_manual(values = c("#3576BD", "Black")) +
      ggtitle("Average Payoffs of Each Strategy") +
      labs(x = "Fraction of Cooperators", y = "Payoff") +
      theme_minimal() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          margin = margin(t = 30, b = 20, unit = "pt"),
          lineheight = 1.15
        ),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.spacing.x = unit(10, 'pt'),
        legend.text = element_text(size = 14),
        axis.title.x =  element_text(margin = margin(t = 15, unit = "pt")),
        axis.title.y =  element_text(margin = margin(r = 15, unit = "pt")),
        text = element_text(size = 16)
      )
    print(plot_PayoffsDF)
    
  }) # END PLOT OUTPUT
  
  # OUTPUT: Stationary distribution plot
  output$stationaryDistributionPlot <- renderPlot({
    
    # Import relevant variables
    MuVector <- computeDynamics()[[3]]
    
    # Plot the stationary distribution
    # Print the stationary distribution µ
    MuDF <-
      data.frame(N = seq(from = 0, to = 1, by = 1 / Z), Probability = MuVector)
    plot_MuDF <- ggplot(data = MuDF, aes(x = N, y = Probability)) +
      geom_density(
        stat = "identity",
        fill = "#3576BD",
        colour = "#3576BD",
        size = 1
      ) +
      ggtitle("Stationary Distribution") +
      labs(x = "Fraction of Cooperators", y = bquote('Probability')) +
      theme_minimal() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          margin = margin(t = 30, b = 20, unit = "pt"),
          lineheight = 1.15
        ),
        axis.title.x =  element_text(margin = margin(t = 15, unit = "pt")),
        axis.title.y =  element_text(margin = margin(r = 15, unit = "pt")),
        text = element_text(size = 16)
      ) +
      geom_vline(
        aes(
        colour = "Critical Threshold"
        ),
        xintercept = as.numeric(input$criticalFraction),
        colour = "Black",
        linetype="dotted",
        size = 2
      )
    print(plot_MuDF)
  }) # END PLOT OUTPUT
  
  # OUTPUT: Selection gradient plot
  output$selectionGradientPlot <- renderPlot({
    
    # Import relevant variables
    selectionGradient <- computeDynamics()[[4]]
    
    # Find all of the (stable & unstable) fixed points using the selection gradient.
    # Creat the empty vectors in which to store the fixed points once we have determined their location
    stableFixedPoint <- c()
    unstableFixedPoint <- c()
    # Determine the sign of the selection gradient for each state
    selectionGradientSign <- (selectionGradient > 0)
    # Find where the sign of the selection gradient flips.
    # These will correspond to the interior fixed points.
    for (i in 1:(Z - 2)) {
      # Find the (interior) stable states by locating where the gradient flips from positive to negative.
      if (selectionGradientSign[i] == TRUE &
          selectionGradientSign[i + 1] == FALSE) {
        stableFixedPoint <-
          append(stableFixedPoint, i + 1 / 2, after = length(stableFixedPoint))
      }
      # Find the (interior) unstable states by locating where the gradient flips from negative to positive
      if (selectionGradientSign[i] == FALSE &
          selectionGradientSign[i + 1] == TRUE) {
        unstableFixedPoint <-
          append(unstableFixedPoint, i + 1 / 2, after = length(unstableFixedPoint))
      }
    }
    # Determine the exterior (stable & unstable) fixed points using the selection gradient
    # First, we consider the two simple cases where one strategy dominates the other:
    # 1. The dynamics is such that Cooperate dominates Defect,
    # then the fixed pointe structure simplifies to 0-unstable & 1-stable.
    if (all(selectionGradient >= 0)) {
      stableFixedPoint <- Z
      unstableFixedPoint <- 0
    }
    # 2. The dynamics is such that Defect dominates Cooperate,
    # then the fixed pointe structure is 1-unstable & 0-stable.
    if (all(selectionGradient <= 0)) {
      stableFixedPoint <- 0
      unstableFixedPoint <- Z
    }
    # If neither strategy dominates the other, 
    # then we find whether the first and last inteior fixed point are stable or unstable,
    # and deduce that their correponding exterior fixed point must exhibit the opposite stability.
    if (!all(selectionGradient >= 0) &
        !all(selectionGradient <= 0)) {
      # Check if there is only an interior stable unstable state,
      # then both exterior fixed points are stable.
      if (is.null(stableFixedPoint)) {
        stableFixedPoint <- c(0, Z)
      } else {
        # Find the fixed point closest to 0, and infer that the fixed point at 0 has the opposite stability
        if (stableFixedPoint[1] < unstableFixedPoint[1]) {
          unstableFixedPoint <- append(unstableFixedPoint, 0, after = 0)
        } else {
          stableFixedPoint <- append(stableFixedPoint, 0, after = 0)
        }
        # Find the fixed point closest to Z, and infer that the fixed point at Z has the opposite stability
        if (stableFixedPoint[length(stableFixedPoint)] > unstableFixedPoint[length(unstableFixedPoint)]) {
          unstableFixedPoint <-
            append(unstableFixedPoint, Z, after = length(unstableFixedPoint))
        } else {
          stableFixedPoint <-
            append(stableFixedPoint, Z, after = length(stableFixedPoint))
        }
        # Alternatively, if there is only an interior stable unstable state,
        # then both exterior fixed points are stable.
        if (is.null(stableFixedPoint)) {
          stableFixedPoint <- c(0, Z)
        }
      }
    }
    
    # Plot the stationary distribution µ
    GradientDF <-
      data.frame(N = seq(from = 1 / Z, to = 1 - 1 / Z, by = 1 / Z), Gradient = selectionGradient)
    plot_GradientDF <- ggplot(data = GradientDF, aes(x = N, y = Gradient)) +
      geom_line(size = 2, color = "#3576BD") +
      scale_color_manual(values = c("#3576BD")) +
      ggtitle("Gradient of Selection") +
      labs(x = "Fraction of Cooperators", y = bquote('Selection for Cooperation')) +
      theme_minimal() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          margin = margin(t = 30, b = 20, unit = "pt"),
          lineheight = 1.15
        ),
        axis.title.x =  element_text(margin = margin(t = 15, unit = "pt")),
        axis.title.y =  element_text(margin = margin(r = 15, unit = "pt")),
        text = element_text(size = 16)
      ) + geom_point(
        # Plot stable fixed points
        data = data.frame(x = stableFixedPoint / Z, y = rep(0, times = length(stableFixedPoint))),
        aes(x, y),
        shape = 21,
        colour = "#3576BD",
        fill = "#3576BD",
        size = 5,
        stroke = 2
      ) + geom_point(
        # Plot unstable fixed points
        data = data.frame(x = unstableFixedPoint / Z, y = rep(0, times = length(unstableFixedPoint))),
        aes(x, y),
        shape = 21,
        colour = "#3576BD",
        fill = "white",
        size = 5,
        stroke = 2
      )
    print(plot_GradientDF)
  }) # END PLOT OUTPUT
  
}) ## EOD

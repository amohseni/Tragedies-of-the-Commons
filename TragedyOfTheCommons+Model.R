# FERMI PROCESS MODEL for PUBLIC GOODS GAME
# by Aydin Mohseni

# Load packages
library(shiny)
library(ggplot2)

# Define global variables
Z <- 100 # The size of the population
N <- 10 # The size of groups 
# n_C <- 10 # The number of cooperators
# n_D <- N - n_C # The number of defectors
r <- .5 # The risk of disaster occuring in the absence of cooperation succeeding
d <- .5 # The cost of disaster (in terms of te proportion of endownments lost) if cooperation fails and disaster occurs
b <- 1 # The initial endowment of each individual
c <- .1 # The proportion of the endowment contributed by cooperators
p_crit <- 0.8 # The critical number of cooperators at which cooperation succeeds
n_crit <- N * p_crit # The critical number of cooperators at which cooperation succeeds
R <- 1 # The rationlaity coefficient (0:= random selection; 1:= replicator dynamics; ∞ := best reponse dynamics)
M <- 0.1 # Mutation/error/noise parameter (portion of the time a random strategy is chosen)

# Define the payoffs to each cooperators and defectors as follows,
# as a function of the number of cooperators n_c, their 
payoff_cooperator <-
  function(n_C) {
    b * (n_C - n_crit >= 0) + b * (1 - r *d) * (n_C - n_crit < 0) - (c * b)
  }
payoff_defector <-
  function(n_C) {
    payoff_cooperator(n_C) + (c * b)
  }
# Now, define the _average_ payoff to each type
# as a function of the number of group sizes 
# and the fraction of each type in the broader population
mean_payoff_cooperator <-
  function(n_C_pop) {
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
    return(z)
  }
mean_payoff_defector <-
  function(n_C_pop) {
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
    return(z)
  }

# Create payoff data frame
payoff_cooperator_vec <- sapply(c(seq(from = 0, to = Z, by = 1)), mean_payoff_cooperator)
payoff_defector_vec <- sapply(seq(from = 0, to = Z, by = 1), mean_payoff_defector)
PayoffsDF <- data.frame(
  N = rep(seq(
    from = 0, to = 1, by = 1 / Z
  ), times = 2),
  Payoff = c(payoff_cooperator_vec, payoff_defector_vec),
  Strategy = rep(c("Cooperate", "Defect"), each = Z + 1)
)


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

# Specify the transition probabilities of the population changing
# from a state with n_C_pop cooperators in the population 
# to one with n_C_pop + 1, n_C_pop - 1, or n_C_pop cooperators.

# Pr(n_C_pop -> n_C_pop + 1)
Prob_n_C_Increase <-
  function(k) {
    (Z - k) / Z * (1 - M) * (1 + exp(R * (
      mean_payoff_defector(k) - mean_payoff_cooperator(k)
    ))) ^ -1 + M * ((Z - k) / (Z - 1))
  }
# Pr(n_C_pop -> n_C_pop - 1)
Prob_n_C_Decrease <- function(k) {
  k / Z * (1 - M)  * (1 + exp(R * (
    mean_payoff_cooperator(k) - mean_payoff_defector(k)
  ))) ^ -1 + M * (k / (Z - 1))
}
# Pr(n_C_pop -> n_C_pop)
Prob_n_C_Stay <- function(k) {
  1 - Prob_n_C_Increase(k) - Prob_n_C_Decrease(k)
}

# TRANSITION MATRIX
# We compute the transition matrix by recursive applications of the transition probabilities
MPM <- outer(
  0:Z,
  0:Z,
  FUN = function(row, column)
    ifelse(column == row + 1, Prob_n_C_Increase(row),
           ifelse(
             column == row, Prob_n_C_Stay(row),
             ifelse(column == row - 1, Prob_n_C_Decrease(row),
                    0)
           ))
)

# STATIONARY DISTRIBUTION  

# Note: We know that the stationary distribution µ exists because the addition of mutation makes the Markov process Ergodic
# That is: Where the Fermi Process is already finite, and aperiodic,
# it is now irreducible (has only one recursive class) as every state communicates with every other (is in a single class),
# and every state will be visited an infinite number of times in the limit.
# In the absence of mutation, there are were two (singleton) recursive classes
# corresponding to the two absorbing states where the population is all A-type or B-type.
# It follows from being ergodic that the limit distribution π is independent of any initial distribution π(0).

# We compute the stationary distribution of the process
# First, we calculate µ_0 = ( Σ^{N}_{k=0} Π^{k}_{i=1} PforwardM_{i-1} / PbackM_{i} ) ^ -1
Mu0Vector <- c() # Empty vector, in which to store each product
# Empty matrix, in which to store each element of the products
Mu0Matrix <- matrix(data = 1,
                    nrow = Z + 1,
                    ncol = Z + 1)
for (k in 2:(Z + 1)) {
  for (i in 2:k) {
    Mu0Matrix[k, i - 1] <- MPM[i - 1, i] / MPM[i, i - 1]
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
    MuMatrix[k, i - 1] <- MPM[i - 1, i] / MPM[i, i - 1]
  }
}
for (i in 1:(Z + 1)) {
  MuVector[i] <- Mu0 * prod(MuMatrix[i,])
}

# Print the stationary distribution µ
MuDF <-
  data.frame(N = seq(from = 0, to = 1, by = 1 / Z), Probability = MuVector)
plot_MuDF <- ggplot(data = MuDF, aes(x = N, y = Probability)) +
  geom_area(
    stat = "identity",
    #width = 1,
    fill = "#3576BD",
    colour = "#3576BD",
    size = 0.1
  ) +
  ggtitle("Stationary Distribution") +
  labs(x = "Fraction of Cooperators", y = bquote('Probability')) +
  #xlim(0, 1) +
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
  )
print(plot_MuDF)

# SELECTION GADIENT
# Now, use the transition matrix to compute the gradient of selection 
# defined as G(i) = Prob_n_C_Increase(i) - Prob_n_C_Decrease(i).
Prob_n_C_Increase_vec <- sapply(c(seq(from = 1, to = Z - 1, by = 1)), Prob_n_C_Increase)
Prob_n_C_Decrease_vec <- sapply(seq(from = 1, to = Z - 1, by = 1), Prob_n_C_Decrease)
selectionGradient <- Prob_n_C_Increase_vec - Prob_n_C_Decrease_vec


# We find all (stable & unstable) fixed points using the selection gradient
# find the sign of the 
selectionGradientSign <- (selectionGradient > 0)
stableFixedPoint <- c()
unstableFixedPoint <- c()
for (i in 1:(Z - 2)) {
  if (selectionGradientSign[i] == TRUE &
      selectionGradientSign[i + 1] == FALSE) {
    stableFixedPoint <-
      append(stableFixedPoint, i + 1 / 2, after = length(stableFixedPoint))
  }
  if (selectionGradientSign[i] == FALSE &
      selectionGradientSign[i + 1] == TRUE) {
    unstableFixedPoint <-
      append(unstableFixedPoint, i + 1 / 2, after = length(unstableFixedPoint))
  }
}

# We find all (stable & unstable) fixed points using the selection gradient
# find the sign of the 
selectionGradientSign <- (selectionGradient > 0)
stableFixedPoint <- c()
unstableFixedPoint <- c()
for (i in 1:(Z - 1)) {
  if (selectionGradientSign[i] == TRUE &
      selectionGradientSign[i + 1] == FALSE) {
    stableFixedPoint <-
      append(stableFixedPoint, i + 1 / 2, after = length(stableFixedPoint))
  }
  if (selectionGradientSign[i] == FALSE &
      selectionGradientSign[i + 1] == TRUE) {
    unstableFixedPoint <-
      append(unstableFixedPoint, i + 1 / 2, after = length(unstableFixedPoint))
  }
}

# Plot the stationary distribution
# Print the stationary distribution µ
GradientDF <-
  data.frame(N = seq(from = 1 / Z, to = 1 - 1 / Z, by = 1 / Z), Gradient = selectionGradient)
plot_GradientDF <- ggplot(data = GradientDF, aes(x = N, y = Gradient)) +
  geom_line(size = 2, color = "#3576BD") +
  scale_color_manual(values = c("#3576BD")) +
  ggtitle("Gradient of Selection") +
  labs(x = "Fraction of Cooperators", y = bquote('Selection for Cooperation')) +
  #ylim(-0.5, 0.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      margin = margin(t = 30, b = 20, unit = "pt"),
      lineheight = 1.15
    ),
    axis.title.x =  element_text(margin = margin(t = 15, unit = "pt")),
    axis.title.y =  element_text(margin = margin(r = 15, unit = "pt")),
    legend.title = element_blank(),
    legend.position="bottom",
    legend.spacing.x = unit(10, 'pt'),
    legend.text = element_text(size = 14),
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
  ) + scale_shape_manual(c("Stable", "Unstable"), aesthetics = "shape", breaks = waiver())
print(plot_GradientDF)


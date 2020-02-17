# FERMI PROCESS MODEL for PUBLIC GOODS GAME
# by Aydin Mohseni

# Load packages
library(shiny)
library(ggplot2)

# Define global variables
Z <- 10 # The size of the population
N <- Z # The size of groups 
n_C <- 5 # The number of cooperators
n_D <- N - n_c # The number of defectors
b <- 1 # The initial endowment of each individual
c <- 0.5 # The proportion of the endowment contributed by cooperators
n_crit <- N / 2 # The critical number of cooperators at which cooperation succeeds
r <- 0.5 # The risk of disaster occuring in the absence of cooperation succeeding
d <- 1 # The cost of disaster (in terms of te proportion of endownments lost) if cooperation fails and disaster occurs
M <- 0.1 # Mutation/error/noise parameter (portion of the time a random strategy is chosen)
R <- 1 # The rationlaity coefficient (0:= random selection; 1:= replicator dynamics; ∞ := best reponse dynamics)

# Define the payoffs to each cooperators and defectors as follows,
# as a function of the number of cooperators n_c, their 
payoff_cooperator <-
  function(n_C) {
    b * (n_C - n_crit >= 0) + (1 - r) * b * (1 - d * (n_C - n_crit >= 0)) - (c * b)
  }
payoff_defector <-
  function(n_C) {
    payoff_cooperator(n_C) + (c * b)
  }

# Create payoff data frame
payoff_cooperator_vec <- sapply(0:N, payoff_cooperator)
payoff_defector_vec <- sapply(0:N, payoff_defector)
PayoffsDF <- data.frame(N = rep(0:N, times = 2), Payoff = c(payoff_cooperator_vec, payoff_defector_vec), Strategy = rep(c("Cooperate", "Defect"), each = N + 1))

# Plot payoff functions
p <- ggplot(data = PayoffsDF, aes(x = N, y = Payoff, group = Strategy)) +
  geom_line(aes(color = Strategy), size = 2) + 
  scale_color_manual(values = c("Black", "Red")) +
  ggtitle("Strategy Payoffs") +
  labs(x = "Number of Cooperators", y = "Payoff") +
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
print(p)

# Specify the transition probabilities of an individual 
# changing from one strategy i to another j (where i,j in {C, D})
# given the current count of individuals playing each strategy (n_C, n_D).
# Pr(n_C -> n_C + 1) 
Prob_n_C_Increase <-
  function(n_C)
    ((N - n_C) / N) * ((n_C / (N - 1)) * (1 - M) / (1 + exp(
      R * (payoff_defector(n_C) - payoff_cooperator(n_C))
    )) + M * ((N - n_C) / (N - 1)))
# Pr(n_D -> n_D + 1)
Prob_n_C_Decrease <- function(n_C)
  ((n_C) / N) * (((N - n_C) / (N - 1)) * (1 - M) / (1 + exp(
    R * (payoff_cooperator(n_C) - payoff_defector(n_C))
  )) + M * (n_C / (N - 1)))
# Pr(j -> j)
Prob_n_C_Stay <- function(n_C)
  1 - Prob_n_C_Increase(n_C) - Prob_n_D_Derease(n_C)

# Pr(j -> k)

# TRANSITION MATRIX
# We compute the transition matrix by recursive applications of the transition probabilities
MPM <- outer(
  0:N,
  0:N,
  FUN = function(r, c)
    ifelse(c == r + 1, Prob_n_C_Increase(r),
           ifelse(
             c == r, Prob_n_C_Stay(r),
             ifelse(c == r - 1, Prob_n_C_Decrease(r),
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
                    nrow = N + 1,
                    ncol = N + 1)
for (k in 2:(N + 1)) {
  for (i in 2:k) {
    Mu0Matrix[k, i - 1] <- MPM[i - 1, i] / MPM[i, i - 1]
  }
}
# Take the product of the rows of Mu0Matrix
for (i in 1:(N + 1)) {
  Mu0Vector[i] <- prod(Mu0Matrix[i,])
}
# Compute µ_0
Mu0 <- sum(Mu0Vector) ^ -1

# Now, we calculate µ_k = µ_0 * Π^{k}_{i=1} PforwardM_{i-1} / PbackM_{i} )
MuVector <- c()
# Empty matrix, in which to store each element of the products
MuMatrix <- matrix(data = 1,
                   nrow = N + 1,
                   ncol = N + 1)
for (k in 2:(N + 1)) {
  for (i in 2:k) {
    MuMatrix[k, i - 1] <- MPM[i - 1, i] / MPM[i, i - 1]
  }
}
for (i in 1:(N + 1)) {
  MuVector[i] <- Mu0 * prod(MuMatrix[i,])
}

# Print the stationary distribution µ
MuDF <- data.frame(N = c(0:N), Probability = MuVector)
q <- ggplot(data = MuDF, aes(x = N, y = Probability)) +
  geom_bar(
    stat = "identity",
    width = 1,
    fill = "gray",
    colour = "black",
    size = 0.1
  ) +
  ggtitle("Stationary Distribution") +
  labs(x = "Number of Cooperators", y = bquote('Probability')) +
  # ylim(c(0, 1)) +
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
    xintercept = which.max(MuVector) - 1,
    colour = "red",
    size = 2
  )
print(q)


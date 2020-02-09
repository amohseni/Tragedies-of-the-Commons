# General Tragedy of the Commons Model


# Define the key variables
N <- c(10) # The size of the group 
n_C <- c(5) # The number of cooperators
n_D <- N - n_c # The number of defectors
b <- c(1) # The initial endowment of each individual
c <- c(0.5) # The proportion of the endowment contributed by cooperators
n_crit <- N / 2 # The critical number of cooperators at which cooperation succeeds
d <- c(1) # The proportion of endownments lost if both cooperation fails and disaster occurs
r <- c(0.5) # The risk of disaster occuring in the absence of cooperation succeeding

}
# Define the payoffs to each cooperators and defectors as follows,
# as a function of the number of cooperators n_c, their 
payoff_cooperator <-
  function(N, n_C, c, d, r) {
    z <- b * (n_C - n_crit >= 0) + (1 - r) * b * (1 - d * (n_C - n_crit >= 0)) - (c * b)
    print(paste("Payoff to cooperator", z, sep = " "))
    return(z)
  }
payoff_defector <-
  function(N, n_C, c, d, r) {
    z <- payoff_cooperator(N, n_C, c, d, r) + (c * b)
    return(z)
  }

# Specify the matrix of transition probability in terms of the Fermi dynamics Master-Equation



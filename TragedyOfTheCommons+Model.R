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
# Pr(j -> j+1)
PforwardM <-
  function(i)
    (1 -  η) * ((i * fA(i) / P(i)) * ((N - i) / N)) +   η  * (((N - i) * fB(i) / P(i)) * ((N -
                                                                                             i) / N))
# Pr(j -> j-1)
PbackM <-
  function(i)
    (1 -  η) * (((N - i) * fB(i) / P(i)) * (i / N)) +   η  * ((i * fA(i) / P(i)) * (i / N))
# Pr(j -> j)
PstayM <- function(i)
  1 - (PforwardM(i) + PbackM(i))

# The Fermi Process with Mutation (MPM) Transition Matrix 
MPM <- matrix(nrow = N + 1,
              ncol = N + 1,
              byrow = TRUE)

# We relabel the row and columns to correspond to the number of A-types in the N-size population
rownames(MPM) = c(0:N)
colnames(MPM) = c(0:N)

# TRANSITION MATRIX
# We compute the transition matrix by recursive applications of the transition probabilities
MPM <- outer(
  0:N,
  0:N,
  FUN = function(r, c)
    ifelse(c == r - 1, PbackM(r),
           ifelse(
             c == r, PstayM(r),
             ifelse(c == r + 1, PforwardM(r),
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

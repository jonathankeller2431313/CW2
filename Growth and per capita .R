setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Studium/Pop Ecology/CW2/R scripts/ competition version ")
library(deSolve)

# Model function (as before)
model_stage <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Bird species 1: eggs and adults
    dE1 <- b1 * A1 * (1 - (A1 + alpha12 * A2) / K1) - sigma1 * E1 - beta1 * E1 * R
    dA1 <- sigma1 * s * E1 - m1 * A1
    
    # Bird species 2: eggs and adults
    dE2 <- b2 * A2 * (1 - (A2 + alpha21 * A1) / K2) - sigma2 * E2 - beta2 * E2 * R
    dA2 <- sigma2 * s * E2 - m2 * A2
    
    # Rats (predator) feeding on eggs only
    dR <- delta1 * beta1 * E1 * R + delta2 * beta2 * E2 * R - gamma * R
    
    list(c(dE1, dA1, dE2, dA2, dR))
  })
}

# Define parameters
parameters <- c(
  b1      = 2.0,    # max egg production rate per adult for species 1
  K1      = 180,    # carrying capacity for adults of species 1
  alpha12 = 0.2,    # competition coefficient: impact of species 2 adults on species 1
  sigma1  = 0.3,    # maturation rate from egg to adult for species 1
  m1      = 0.05,   # adult mortality rate for species 1
  beta1   = 0.1,    # predation rate on eggs of species 1
  
  b2      = 1.5,    # max egg production rate per adult for species 2
  K2      = 200,    # carrying capacity for adults of species 2
  alpha21 = 0.5,    # competition coefficient: impact of species 1 adults on species 2
  sigma2  = 0.3,    # maturation rate from egg to adult for species 2
  m2      = 0.05,   # adult mortality rate for species 2
  beta2   = 0.1,    # predation rate on eggs of species 2
  
  delta1  = 0.3,    # conversion efficiency from species 1 eggs to rat growth
  delta2  = 0.3,    # conversion efficiency from species 2 eggs to rat growth
  gamma   = 0.2,    # natural mortality rate of rats
  
  s       = 0.6     # survival from eggs to adult
)

# Define initial state
state <- c(
  E1 = 20,    # initial number of eggs for species 1
  A1 = 100,   # initial number of adults for species 1
  E2 = 20,    # initial number of eggs for species 2
  A2 = 100,   # initial number of adults for species 2
  R  = 150    # initial number of rats
)

# Time vector for the simulation
times <- seq(0, 100, by = 1)

# Solve the ODE system
out <- ode(y = state, times = times, func = model_stage, parms = parameters)
out_df <- as.data.frame(out)

# ------------------------------
# Plot 1: Population Growth Plot
# ------------------------------
# Calculate discrete population growth for species 1 adults (A1), species 2 adults (A2), and Rats (R)
# Use diff() to compute changes over time; add NA for the first time point.
out_df$A1_growth <- c(NA, diff(out_df$A1))
out_df$A2_growth <- c(NA, diff(out_df$A2))
out_df$R_growth  <- c(NA, diff(out_df$R))

# Plot the population growth rates over time
matplot(out_df$time, out_df[, c("A1_growth", "A2_growth", "R_growth")],
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "Population Growth (Change per Time Unit)",
        col = c("blue", "red", "green"),
        main = "Population Growth over Time")
legend("topright", legend = c("Species 1 Adults", "Species 2 Adults", "Rats"),
       col = c("blue", "red", "green"), lty = 1, lwd = 2)

# ------------------------------
# Plot 2: Per Capita Growth Plot
# ------------------------------
# Calculate per capita growth: (growth/population)
# We avoid division by zero by using NA when population is 0.
out_df$A1_pc_growth <- with(out_df, ifelse(A1 > 0, A1_growth / A1, NA))
out_df$A2_pc_growth <- with(out_df, ifelse(A2 > 0, A2_growth / A2, NA))
out_df$R_pc_growth  <- with(out_df, ifelse(R > 0, R_growth / R, NA))

# Remove the first row (NA) for plotting purposes
pc_growth_df <- out_df[-1, ]

matplot(pc_growth_df$time, pc_growth_df[, c("A1_pc_growth", "A2_pc_growth", "R_pc_growth")],
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "Per Capita Growth (Change / Population)",
        col = c("blue", "red", "green"),
        main = "Per Capita Growth over Time")
legend("topright", legend = c("Species 1 Adults", "Species 2 Adults", "Rats"),
       col = c("blue", "red", "green"), lty = 1, lwd = 2)

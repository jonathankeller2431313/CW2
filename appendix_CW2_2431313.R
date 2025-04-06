setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Studium/Pop Ecology/CW2/R scripts/ competition version ")

#install these packages if not done before
#install.packages(ggplot2)
#install.packages(deSolve)

require(ggplot2)
require(deSolve) 

################## defining model via the equations from the methods section (9) (10) (11) #########################

model_stage <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dE1 <- b1 * A1 * (1 - (A1 + alpha12 * A2) / K1) - sigma1 * E1 - beta1 * E1 * R  #equation 9 for egg stage of bird species 1
    dA1 <- sigma1 * s * E1 - m1* A1   #equation 10 for adult stage of bird species 1
    
    dE2 <- b2 * A2 * (1 - (A2 + alpha21 * A1) / K2) - sigma2 * E2 - beta2 * E2 * R #equation 9 for egg stage of bird species 2
    dA2 <- sigma2 * s * E2 - m2*A2    ##equation 10 for adult stage of bird species 2
    
    dR <- delta1 * beta1 * E1 * R + delta2 * beta2 * E2 * R - gamma * R #equation 11 for rat population R
    
    list(c(dE1, dA1, dE2, dA2, dR))
  })
}

#setting parameters for the model
parameters <- c(
  #Bird species 1
  b1      = 2.0,    #egg production rate per
  K1      = 180,    #carrying capacity 
  alpha12 = 0.2,   #competition coefficient: impact of species 2 adults on species 1
  sigma1  = 0.3,    #maturation rate 
  m1      = 0.05,    #mortality rate 
  beta1   = 0.1,  #predation rate
  
  #Bird species 2
  b2      = 1.5,    #egg production rate 
  K2      = 200,     #carrying capacity
  alpha21 = 0.5,   #competition coefficient: impact of species 1 adults on species 2
  sigma2  = 0.3,    #maturation rate 
  m2      = 0.05,    #adult mortality
  beta2   = 0.1,  #predation rate
  
  # Rat parameters
  delta1  = 0.3,    #conversion efficiency from species 1 
  delta2  = 0.3,    #conversion efficiency from species 2
  gamma   = 0.2,    #natural mortality rate of rats
  
  s       = 0.6     #survival from eggs to adult
)

#set initial numbers
state <- c(
  E1 = 20,   #initial number of eggs for species 1
  A1 = 100,   #initial number of adults for species 1
  E2 = 20,    #initial number of eggs for species 2
  A2 = 100,   #initial number of adults for species 2
  R  = 150     #initial number of rats
)


times <- seq(0, 100, by = 1)  #simulate over 100 time units


out <- ode(y = state, times = times, func = model_stage, parms = parameters) #solve differential equation with ode function





#######################plot figure 1 (time series plot with default parameters)#################


out_df <- as.data.frame(out) #convert to data frame for plotting

#look at results 
print(out_df)

#Plot the results of the model defined in lines 9-22: Eggs and adults for both bird species and rats
matplot(out_df$time, out_df[,2:6],
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "Population size",
        col = c("blue", "darkblue", "red", "darkred", "darkgreen"))
legend("topright",
       legend = c("Species 1 Eggs", "Species 1 Adults",
                  "Species 2 Eggs", "Species 2 Adults",
                  "Rats"),
       col = c("blue", "darkblue", "red", "darkred", "darkgreen"),
       lty = 1, lwd = 2, cex = 0.38)



########################plot figure 2 (population growth rate)#####################



#Calculate discrete population growth for species 1 adults (A1), species 2 adults (A2), and Rats (R)
#based on model defined above in lines 9-22
out_df$A1_growth <- c(NA, diff(out_df$A1))
out_df$A2_growth <- c(NA, diff(out_df$A2))
out_df$R_growth  <- c(NA, diff(out_df$R))

#look at resulting values 
print(out_df)

#Plot the population growth rates over time
matplot(out_df$time, out_df[, c("A1_growth", "A2_growth", "R_growth")],
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "Population Growth (Change per Time Unit)",
        col = c("blue", "red", "green"),
        main = "Population Growth over Time")
legend("topright", legend = c("Species 1 Adults", "Species 2 Adults", "Rats"),
       col = c("blue", "red", "green"), lty = 1, lwd = 2)




#########################plot figure 3 (per-capita grwoth rate)####################



#Calculate per capita growth: (growth/population). Based on model defined above in lines 9-22
out_df$A1_pc_growth <- with(out_df, ifelse(A1 > 0, A1_growth / A1, NA)) #avoid division by zero by using NA when population is 0.
out_df$A2_pc_growth <- with(out_df, ifelse(A2 > 0, A2_growth / A2, NA))
out_df$R_pc_growth  <- with(out_df, ifelse(R > 0, R_growth / R, NA))

#look at resulting values 
print(out_df)

#Remove the first row (NA) for plotting 
pc_growth_df <- out_df[-1, ]

matplot(pc_growth_df$time, pc_growth_df[, c("A1_pc_growth", "A2_pc_growth", "R_pc_growth")],
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "Per Capita Growth (Change / Population)",
        col = c("blue", "red", "green"),
        main = "Per Capita Growth over Time")
legend("topright", legend = c("Species 1 Adults", "Species 2 Adults", "Rats"),
       col = c("blue", "red", "green"), lty = 1, lwd = 2)




###################figure 4 (Bifurcation plot for initial rat population)#####################

#Range of initial rat populations R
rat_range <- seq(0, 150, by = 5)

#Store results
results <- data.frame(R_init = numeric(), A1_final = numeric(), A2_final = numeric(), R_final = numeric())

#for loop running model which has already been defined in lines 9-22 with each value of R within the defined range
for (R_init in rat_range) {
  state <- c(E1 = 20, A1 = 100, E2 = 20, A2 = 100, R = R_init)
  out <- ode(y = state, times = times, func = model_stage, parms = parameters)
  out_df <- as.data.frame(out)
  
  #store final population sizes
  results <- rbind(results, data.frame(
    R_init = R_init,
    A1_final = tail(out_df$A1, 1),
    A2_final = tail(out_df$A2, 1),
    R_final = tail(out_df$R, 1)
  ))
}

#Plot bifurcation 
plot(results$R_init, results$A1_final, type = "b", col = "blue", pch = 16, ylim = c(0, max(results$A1_final, results$A2_final, results$R_final)),
     xlab = "Initial Rat Population", ylab = "Final Population Size", main = "Bifurcation Plot for Initial Rat Population")
points(results$R_init, results$A2_final, type = "b", col = "red", pch = 16)
points(results$R_init, results$R_final, type = "b", col = "green", pch = 16)
legend("topright", legend = c("Species 1 Adults", "Species 2 Adults", "Rats"), col = c("blue", "red", "green"), pch = 16, lty = 1)



################################figure 5 (time series with initial R=5)##################################


# change initial state so that R=5
state <- c(
  E1 = 20,   # initial number of eggs for species 1
  A1 = 100,   # initial number of adults for species 1
  E2 = 20,    # initial number of eggs for species 2
  A2 = 100,   # initial number of adults for species 2
  R  = 5     # initial number of rats
)


#Solve the differential equation with ODE funtion 
out <- ode(y = state, times = times, func = model_stage, parms = parameters)

#Convert output to a data frame plotting
out_df <- as.data.frame(out)

#look at resulting values 
print(out_df)

#Plot the results
matplot(out_df$time, out_df[,2:6],
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "Population size",
        col = c("blue", "darkblue", "red", "darkred", "darkgreen"))
legend("topright",
       legend = c("Species 1 Eggs", "Species 1 Adults",
                  "Species 2 Eggs", "Species 2 Adults",
                  "Rats"),
       col = c("blue", "darkblue", "red", "darkred", "darkgreen"),
       lty = 1, lwd = 2, cex = 0.38)



#################################figure 6 (time series with initial R=0)################################


# change initial state so that R=0
state <- c(
  E1 = 20,   # initial number of eggs for species 1
  A1 = 100,   # initial number of adults for species 1
  E2 = 20,    # initial number of eggs for species 2
  A2 = 100,   # initial number of adults for species 2
  R  = 0     # initial number of rats
)



# Solve differential equation with ODE function 
out <- ode(y = state, times = times, func = model_stage, parms = parameters)

# Convert output to data frame for plotting
out_df <- as.data.frame(out)

#look at resulting values 
print(out_df)



# Plot the results: Eggs and adults for both bird species and rats
matplot(out_df$time, out_df[,2:6],
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "Population size",
        col = c("blue", "darkblue", "red", "darkred", "darkgreen"))
legend("topright",
       legend = c("Species 1 Eggs", "Species 1 Adults",
                  "Species 2 Eggs", "Species 2 Adults",
                  "Rats"),
       col = c("blue", "darkblue", "red", "darkred", "darkgreen"),
       lty = 1, lwd = 2, cex = 0.38)





###############################figure 7 (time series with rat mortality 0.5)########################


parameters <- c(
  #Bird species 1
  b1      = 2.0,    #egg production rate per
  K1      = 180,    #carrying capacity 
  alpha12 = 0.2,   #competition coefficient: impact of species 2 adults on species 1
  sigma1  = 0.3,    #maturation rate 
  m1      = 0.05,    #mortality rate 
  beta1   = 0.1,  #predation rate
  
  #Bird species 2
  b2      = 1.5,    #egg production rate 
  K2      = 200,     #carrying capacity
  alpha21 = 0.5,   #competition coefficient: impact of species 1 adults on species 2
  sigma2  = 0.3,    #maturation rate 
  m2      = 0.05,    #adult mortality
  beta2   = 0.1,  #predation rate
  
  # Rat parameters
  delta1  = 0.3,    #conversion efficiency from species 1 
  delta2  = 0.3,    #conversion efficiency from species 2
  gamma   = 0.5,    #natural mortality rate of rats
  
  s       = 0.6     #survival from eggs to adult
)

#set initial numbers back to R=150
state <- c(
  E1 = 20,   
  A1 = 100,   
  E2 = 20,    
  A2 = 100,
  R  = 150
)


# Solve differential equation with ODE function 
out <- ode(y = state, times = times, func = model_stage, parms = parameters)

# Convert output to data frame for plotting
out_df <- as.data.frame(out)

#look at resulting values 
print(out_df)



# Plot the results: Eggs and adults for both bird species and rats
matplot(out_df$time, out_df[,2:6],
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "Population size",
        col = c("blue", "darkblue", "red", "darkred", "darkgreen"))
legend("topright",
       legend = c("Species 1 Eggs", "Species 1 Adults",
                  "Species 2 Eggs", "Species 2 Adults",
                  "Rats"),
       col = c("blue", "darkblue", "red", "darkred", "darkgreen"),
       lty = 1, lwd = 2, cex = 0.38)


##########################figure 8 (bifurcation plot for varying rat mortality)#####################


#Define a range for rat mortality (gamma) 
gamma_values <- seq(0.1, 1.0, by = 0.05)

#empty data frame to store results in 
results <- data.frame(
  gamma = gamma_values,
  E1_final = NA,
  A1_final = NA,
  E2_final = NA,
  A2_final = NA,
  R_final  = NA
)

# Loop over the gamma values, run the simulation based on model defined in 
#lines 9-22, and record the final state variables
for (i in seq_along(gamma_values)) {
  params_temp <- parameters
  params_temp["gamma"] <- gamma_values[i]
  
  out <- ode(y = state, times = times, func = model_stage, parms = params_temp)
  out_df <- as.data.frame(out)
  
  # Record final values (last row of the output)
  final_row <- tail(out_df, 1)
  results$E1_final[i] <- final_row$E1
  results$A1_final[i] <- final_row$A1
  results$E2_final[i] <- final_row$E2
  results$A2_final[i] <- final_row$A2
  results$R_final[i]  <- final_row$R
}

# Create a long-format data frame for plotting
# For each gamma, repeat the gamma value for each state variable.
results_long <- data.frame(
  gamma = rep(results$gamma, times = 5),
  StateVariable = rep(c("E1", "A1", "E2", "A2", "R"), each = nrow(results)),
  FinalValue = c(results$E1_final, results$A1_final, results$E2_final, results$A2_final, results$R_final)
)

#Create the bifurcation plot 
ggplot(results_long, aes(x = gamma, y = FinalValue, color = StateVariable)) +
  geom_line() +
  geom_point() +
  labs(x = "Rat Mortality Rate (γ)",
       y = "Final Population Size",
       title = "Bifurcation Plot: Final Populations vs. Rat Mortality") +
  theme_bw()+
  theme(
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black") 
  )


###################################Table 1 (testing statistical significance of linear fitted model)############



# Fit a linear model for A1 vs. gamma
lm_A1 <- lm(A1_final ~ gamma, data = results)
summary(lm_A1) #give out result

# Fit a linear model for A2 vs. gamma
lm_A2 <- lm(A2_final ~ gamma, data =results)
summary(lm_A2)

# Fit a linear model for E1 vs. gamma
lm_E1 <- lm(E1_final ~ gamma, data =results)
summary(lm_E1)

# Fit a linear model for E2 vs. gamma
lm_E2 <- lm(E2_final ~ gamma, data =results)
summary(lm_E2)

# Fit a linear model for R vs. gamma
lm_R <- lm(R_final ~ gamma, data =results)
summary(lm_R)



#################################figure 9 (isocline plot)################################


#Define parameters
K1 <- 180         # carrying capacity for species 1 
K2 <- 200         # carrying capacity for species 2 
alpha12 <- 0.2    # effect of species 2 on species 1
alpha21 <- 0.5    # effect of species 1 on species 2

#Create a sequence for A2 for the blue isocline: from 0 to K1/alpha12 = 900
A2_vals <- seq(0, K1/alpha12, length.out = 100)
# Blue isocline: A1 = K1 - alpha12 * A2
A1_iso1 <- K1 - alpha12 * A2_vals

#For the red isocline, valid for A2 in [0, K2] (0 to 200)
A2_vals2 <- seq(0, K2, length.out = 100)
# Red isocline: A1 = (K2 - A2) / alpha21
A1_iso2 <- (K2 - A2_vals2) / alpha21

# Plot the isoclines,
plot(A2_vals, A1_iso1, type = "l", col = "blue", lwd = 2,
     xlab = "Species 2 Adult Population (A2)",
     ylab = "Species 1 Adult Population (A1)",
     xlim = c(0, 900),   # extend x-axis to show blue isocline intersection at A2 = 900
     ylim = c(0, 450),   # extend y-axis to show red isocline intersection at A1 = 400
     main = "Isocline Plot for Bird Species Competition")
lines(A2_vals2, A1_iso2, col = "red", lwd = 2)
abline(h = 0, v = 0, col = "gray", lty = 2)  # add reference lines at 0
legend("topright", legend = c(expression(A[1] + alpha[12]*A[2] == K[1]),
                              expression(A[2] + alpha[21]*A[1] == K[2])),
       col = c("blue", "red"), lwd = 2)





##########################figure 10 (bifurcation plot for varying alpha12)##############################

#set parameters back to default 
parameters <- c(
  #Bird species 1
  b1      = 2.0,    #egg production rate per
  K1      = 180,    #carrying capacity 
  alpha12 = 0.2,   #competition coefficient: impact of species 2 adults on species 1
  sigma1  = 0.3,    #maturation rate 
  m1      = 0.05,    #mortality rate 
  beta1   = 0.1,  #predation rate
  
  #Bird species 2
  b2      = 1.5,    #egg production rate 
  K2      = 200,     #carrying capacity
  alpha21 = 0.5,   #competition coefficient: impact of species 1 adults on species 2
  sigma2  = 0.3,    #maturation rate 
  m2      = 0.05,    #adult mortality
  beta2   = 0.1,  #predation rate
  
  # Rat parameters
  delta1  = 0.3,    #conversion efficiency from species 1 
  delta2  = 0.3,    #conversion efficiency from species 2
  gamma   = 0.2,    #natural mortality rate of rats
  
  s       = 0.6     #survival from eggs to adult
)


# Define range for varying alpha12 (competition from species 2 to species 1)
alpha12_values <- seq(0, 1, by = 0.1)

#empty data frame to store results
results <- data.frame(alpha12 = alpha12_values, E1_final = NA, A1_final = NA, E2_final = NA, A2_final = NA, R_final = NA)

#for loop over varying alpha12 running simulationss of the model defined in lines 9-22
for (i in 1:length(alpha12_values)) {
  params_temp <- parameters
  params_temp["alpha12"] <- alpha12_values[i]  # Update alpha12 while keeping alpha21 constant
  
  out <- ode(y = state, times = times, func = model_stage, parms = params_temp)
  out_df <- as.data.frame(out)
  
  #store final population values
  final_row <- tail(out_df, 1)
  results$E1_final[i] <- final_row$E1
  results$A1_final[i] <- final_row$A1
  results$E2_final[i] <- final_row$E2
  results$A2_final[i] <- final_row$A2
  results$R_final[i]  <- final_row$R
}

#Convert results to long format for plotting
results_long <- data.frame(
  alpha12 = rep(results$alpha12, times = 5),
  StateVariable = rep(c("E1", "A1", "E2", "A2", "R"), each = length(alpha12_values)),
  FinalValue = c(results$E1_final, results$A1_final, results$E2_final, results$A2_final, results$R_final)
)

#Plot: Final population vs competition coefficient alpha12
ggplot(results_long, aes(x = alpha12, y = FinalValue, color = StateVariable)) +
  geom_line() +
  geom_point() +
  labs(x = "Competition Coefficient α12",
       y = "Final Population Size",
       title = "Effect of α12 (Competition from Species 2 on Species 1)") +
  theme_bw() +  # White background with black border
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black") 
  )



######################figure 11 (bifurcation plot for varying alpha21)##########################


# Define range for varying alpha12 (competition from species 2 to species 1)
alpha21_values <- seq(0, 1, by = 0.1)

#emtpty  data frame to store results
results <- data.frame(alpha21 = alpha21_values, E1_final = NA, A1_final = NA, E2_final = NA, A2_final = NA, R_final = NA)

#for Loop over varying alpha12 simulating the model defined in lines 9-22
for (i in 1:length(alpha21_values)) {
  params_temp <- parameters
  params_temp["alpha21"] <- alpha21_values[i]  # Update alpha12 while keeping alpha21 constant
  
  out <- ode(y = state, times = times, func = model_stage, parms = params_temp)
  out_df <- as.data.frame(out)
  
  #Store final population values
  final_row <- tail(out_df, 1)
  results$E1_final[i] <- final_row$E1
  results$A1_final[i] <- final_row$A1
  results$E2_final[i] <- final_row$E2
  results$A2_final[i] <- final_row$A2
  results$R_final[i]  <- final_row$R
}

#convert results to long format for plotting
results_long <- data.frame(
  alpha21 = rep(results$alpha21, times = 5),
  StateVariable = rep(c("E1", "A1", "E2", "A2", "R"), each = length(alpha21_values)),
  FinalValue = c(results$E1_final, results$A1_final, results$E2_final, results$A2_final, results$R_final)
)

# Plot: Final population vs competition coefficient alpha21
ggplot(results_long, aes(x = alpha21, y = FinalValue, color = StateVariable)) +
  geom_line() +
  geom_point() +
  labs(x = "Competition Coefficient α21",
       y = "Final Population Size",
       title = "Effect of α21 ") +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black") 
  )




###########################figure 12 (interplay between gamma and alpha12)############################


#define ranges for rat mortality (gamma) and competition coefficient (alpha21)
gamma_values   <- seq(0, 1.0, length.out = 20)  # x-axis values
alpha21_values <- seq(0, 1.0, length.out = 20)    # y-axis values

# Create a grid of parameter combinations
results <- expand.grid(gamma = gamma_values, alpha21 = alpha21_values)
results$A2_final <- NA  # final population of species 2 adults

#for loop over each combination and simulate the model defined in lines 9-22
for (i in 1:nrow(results)) {
  params_temp <- parameters
  params_temp["gamma"]   <- results$gamma[i]
  params_temp["alpha21"] <- results$alpha21[i]
  
  out <- ode(y = state, times = times, func = model_stage, parms = params_temp)
  out_df <- as.data.frame(out)
  
  # Record the final species 2 adult population (A2) at time = 100
  results$A2_final[i] <- tail(out_df$A2, 1)
}

#put in long format data frame for plotting
results_long <- data.frame(
  gamma = rep(results$gamma, times = 1),
  alpha21 = results$alpha21,
  A2_final = results$A2_final
)

#plot the heatmap
ggplot(results_long, aes(x = gamma, y = alpha21, fill = A2_final)) +
  geom_tile() +
  scale_fill_gradient(name = "Final A2", low = "white", high = "red") +
  labs(x = expression("Rat Mortality (" * gamma * ")"),
       y = expression("Competition Coefficient (" * alpha[21] * ")"),
       title = "Heatmap: Final Species 2 Adults vs. Rat Mortality and " ~ alpha[21]) +
  theme_bw()




###########################figure 13 (interplay between gamma and alpha21)############################


#Define the range for rat mortality (gamma) and competition coefficient (alpha12)
gamma_values   <- seq(0.1, 1.0, length.out = 20)   # rat mortality on x-axis
alpha12_values <- seq(0.1, 0.5, length.out = 20)       # alpha12 on y-axis

#Create a grid of parameter combinations
results <- expand.grid(gamma = gamma_values, alpha12 = alpha12_values)
results$A1_final <- NA  # final population of Species 1 adults

#for Loop over each combination and simulate the model as defined in lines 9-11
for (i in 1:nrow(results)) {
  params_temp <-parameters
  params_temp["gamma"]   <- results$gamma[i]
  params_temp["alpha12"] <- results$alpha12[i]
  
  out <- ode(y = state, times = times, func = model_stage, parms = params_temp)
  out_df <- as.data.frame(out)
  
  # Record final Species 1 adult population at time = 100
  results$A1_final[i] <- tail(out_df$A1, 1)
}

#put in long format data frame for plotting
results_long <- data.frame(
  gamma = results$gamma,
  alpha12 = results$alpha12,
  A1_final = results$A1_final
)

#plot the heatmap
ggplot(results_long, aes(x = gamma, y = alpha12, fill = A1_final)) +
  geom_tile() +
  scale_fill_gradient(name = "Final A1", low = "white", high = "red") +
  labs(x = expression("Rat Mortality (" * gamma * ")"),
       y = expression("Competition Coefficient (" * alpha[12] * ")"),
       title = "Heatmap: Final Species 1 Adults vs. Rat Mortality and " ~ alpha[12]) +
  theme_bw()




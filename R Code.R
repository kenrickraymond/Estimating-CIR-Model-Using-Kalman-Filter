SimulateCIR = function(beta, alpha, sigma, lambda, dt, ratestart, months, tau) {
  
  # This function simulates one-factor CIR short rate and returns zero-coupon rates.
  # Parameters:
  #   beta: Long-term mean level of the short rate
  #   alpha: Speed of reversion to the mean
  #   sigma: Volatility of the short rate
  #   lambda: Risk premium parameter
  #   dt: Time step in years
  #   ratestart: Initial short rate
  #   months: Number of months to simulate
  #   tau: Vector of maturities (time to maturity)
  
  # CIR parameters (default values)
  # beta = 0.10
  # alpha = 0.05
  # sigma = 0.075
  # lambda = -0.4
  # dt = 1/12
  # months = 120
  # ratestart = 0.10
  # tau = c(3/12, 6/12, 2, 5)  
  
  # Initialize vectors/matrices
  srt <- vector(mode="numeric", length=480)  # short rate vector sampled weekly
  Rt = matrix(, 120, length(tau))  # Zero-coupon rate matrix
  
  # Simulation of the short rate following the CIR model
  srt[1] = ratestart
  for (i in 1:480) {
    srt[i + 1] = srt[i] + alpha * (beta - srt[i]) * dt / 4 + 
      sqrt(srt[i]) * sigma * sqrt(dt / 4) * rnorm(1)
  }
  
  # Calculate zero-coupon rates based on term structure dynamics
  for (i in 1:months) {
    Rttemp = srt[i * 4 - 3]  # short rate at the start of the month
    for (j in 1:length(tau)) {
      
      # Term structure parameters
      Gamma = sqrt((alpha + lambda)^2 + 2 * sigma^2) 
      CIRB = 2 * (exp(Gamma * tau[j]) - 1) / 
        ((Gamma + alpha + lambda) * (exp(Gamma * tau[j]) - 1) + 2 * Gamma)
      CIRA = 2 * alpha * beta / (sigma^2) * 
        log(2 * Gamma * exp((Gamma + alpha + lambda) * tau[j] / 2) / 
              ((Gamma + alpha + lambda) * (exp(Gamma * tau[j]) - 1) + 2 * Gamma))
      
      # Measurement equation for Kalman filter: Rt = A + B * srt
      A = -CIRA / tau[j]
      B = CIRB / tau[j]
      Rt[i, j] = A + B * Rttemp
    }
  }
  
  # Return zero-coupon rates matrix
  return(Rt)
}

LogLikelihoodCIR = function(para, Y, tau, nrow, ncol) {
  # This function computes the log-likelihood value of CIR parameters given the zero-coupon rates Y
  # Parameters:
  #   para: Vector of parameters [beta, alpha, sigma, lambda, measurement noise]
  #   Y: Matrix of observed zero-coupon rates
  #   tau: Vector of maturities (time to maturity)
  #   nrow: Number of rows in the observed rates matrix Y
  #   ncol: Number of columns in the observed rates matrix Y
  
  # Extract parameters
  beta = para[1]
  alpha = para[2]
  sigma = para[3]
  lambda = para[4]
  sigmai = para[5:length(para)]  # Measurement noise
  
  # Initialize measurement noise covariance matrix
  R = diag(ncol)
  for (i in 1:ncol) {
    R[i, i] = sigmai[i]^2
  }
  
  dt = 1/12
  
  # System matrices initialization
  C = beta * (1 - exp(-alpha * dt))  # Transition coefficient
  F = exp(-alpha * dt)                # Transition coefficient
  A = as.matrix(t(rep(0, ncol)))      # Measurement coefficient
  H = as.matrix(A)
  
  # Create A and B for term structure dynamics
  for (j in 1:ncol) {
    Gamma = sqrt((alpha + lambda)^2 + 2 * sigma^2) 
    CIRB = 2 * (exp(Gamma * tau[j]) - 1) / 
      ((Gamma + alpha + lambda) * (exp(Gamma * tau[j]) - 1) + 2 * Gamma)
    CIRA = 2 * alpha * beta / (sigma^2) * 
      log(2 * Gamma * exp((Gamma + alpha + lambda) * tau[j] / 2) / 
            ((Gamma + alpha + lambda) * (exp(Gamma * tau[j]) - 1) + 2 * Gamma))
    A = -CIRA / tau[j]
    B = CIRB / tau[j]
  }
  
  # Kalman Filter
  # Initialize state vector
  AdjR = as.matrix(beta)  # Unconditional mean of state variable, r_{t-1}
  VarR = as.matrix(sigma^2 * beta / (2 * alpha))  # Unconditional variance
  LL = t(rep(0, nrow))  # Initialize likelihood vector
  
  for (i in 1:nrow) {
    
    # Update state variable and its variance
    PredR = C + F * AdjR  # E[r_{t_i}|F_{t_{i-1}}]
    Q = beta * sigma^2 * (1 - exp(-alpha * dt))^2 / (2 * alpha) + 
      sigma^2 / alpha * (exp(-alpha * dt) - exp(-2 * alpha * dt)) * AdjR 
    VarR = F * VarR * F + Q # var[r_{t_{i+1}|F_{t_i}}]
    
    # Forecast measurement and its variance
    PredZ = A + H * c(PredR) # E[z_{t_i}| F_{t_{i-1}}]
    VarZ = (t(H) * c(VarR)) %*% H + R # var[z_t_i|F_{t_{i-1}}]
    
    # Update state inference
    PredError = as.matrix(Y[i, ] - PredZ) # Zeta
    KalmanGain = (c(VarR) * H) %*% solve(VarZ) 
    AdjR = PredR + KalmanGain %*% t(PredError)  
    VarR = c(VarR) %*% (c(1) - KalmanGain %*% c(t(H)))  
    
    # Construct the likelihood function
    DetY = det(VarZ)
    LL[i] = -(ncol / 2) * log(2 * pi) - 0.5 * log(DetY) - 
      0.5 * PredError %*% solve(VarZ) %*% t(PredError)  
  }
  
  sumll = -sum(LL)
  return(sumll)
}

MinimizeLL = function(n_sim=100) {
  # This function minimizes the log-likelihood function to estimate CIR model parameters
  
  CIRsim = matrix(NA, n_sim, 8)  # Create a matrix to store parameter estimates
  
  tau = c(3/12, 6/12, 2, 5)  # Maturities
  
  for (i in 1:n_sim) {
    # Simulate zero-coupon rates
    Ya = SimulateCIR(0.10, 0.05, 0.075, -0.4, 1/12, 0.10, 120, tau)
    Y = Ya
    nrow = dim(Y)[1]
    ncol = dim(Y)[2]
    
    # Initial parameter guesses
    para = c(0.10, 0.05, 0.075, -0.4, 0.1 * rnorm(ncol))
    
    # Minimize negative log-likelihood
    x = nlm(LogLikelihoodCIR, para, Ya, tau, nrow, ncol)
    CIRsim[i, ] = x$estimate
  }
  
  # Calculate mean and standard deviation of parameter estimates
  mean_estimates = colMeans(CIRsim, na.rm = TRUE)
  sd_estimates = apply(CIRsim, 2, sd, na.rm = TRUE)
  
  return(list(mean = mean_estimates, sd = sd_estimates))  # Return mean and standard deviation
}

# Run the function and store the result
data = MinimizeLL(n_sim = 200)

# Print the result
print(data)


# # Plotting Zero Coupon Rates over Time for Different Maturities
library(ggplot2)

# Simulate zero-coupon rates
tau = c(3/12, 6/12, 2, 5)  
zero_coupon_rates = SimulateCIR(beta = 0.10, alpha = 0.05, sigma = 0.075, 
                                lambda = -0.40, dt = 1/12, ratestart = 0.10, 
                                months = 120, tau = tau)

# Prepare data for plotting
time_points = seq(1, 120)
data_list = list()
for (i in 1:length(tau)) {
  data_list[[i]] = data.frame(Time = time_points, 
                              Rate = zero_coupon_rates[, i], 
                              Maturity = paste(tau[i], "years", sep = " "))
}

# Combine all data frames into one
plot_data = do.call(rbind, data_list)

# Plot using ggplot2
ggplot(plot_data, aes(x = Time, y = Rate, color = Maturity)) +
  geom_line() +
  labs(title = "Zero-Coupon Rates Over Time for Different Maturities",
       x = "Time (Months)",
       y = "Zero-Coupon Rate",
       color = "Maturity") +
  theme_minimal()


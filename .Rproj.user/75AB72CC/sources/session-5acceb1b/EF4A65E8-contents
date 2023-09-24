# This script uses TMB to estimate length-at-age of a simulated population using a three parameter Von Bertalanffy Growth Function fit using maximum likelihood.

require(TMB)
# NOTE: You need RTools to make TMB function!!!!
# Follow the online instructions here:
# https://github.com/kaskr/adcomp


###### Simulate data
N_obs = 1000 # Number of observations
Min_age = 0 # Minimum age in the population
Max_age = 10 # Maximum age in the population


# Parameters
Obs_sd = 20 # Standard deviation in 
Linf = 350 # Asymptoptic length/mean length at maximum age (mm)
k = .3 # Growth rate (yr^-1)
t0 = -1.5 # Estimated age at length 0 (yr)


# Simulation
Age_i = runif(N_obs, Min_age, Max_age)

Obs_i = rnorm(N_obs, # 1000 observations
              Linf * (1 - exp(-k * (Age_i - t0))), # Each with a length
              Obs_sd) # standard deviation


# Plot length-at-age
plot(Age_i,Obs_i , xlab = "Age (yr)", ylab = "Length (mm)")


# Load model template
compile("tmb/vbgf.cpp") # Compile a C++ templated into a shared object file
dyn.load(dynlib("tmb/vbgf"))


# Run the model
obj = MakeADFun(
  data = list(Length = Obs_i, Age = Age_i), 
  parameters = list(Linf = 400, Kappa = .3, T_zero = 0, LogSigma = 0), DLL = "vbgf") # Prepare model object, data, and initial parameter estimates for fitting


# Fit model using nlminb optimiser
opt1 = nlminb(obj$par, 
              obj$fn, 
              obj$gr) 

summ = summary(sdreport(obj)) # Get summary of parameter estimates
summ


# Plot the fit
curve(summ[1,1] * (1 - exp(-summ[2,1] * (x - summ[3,1]))), 
      from = 0, 
      to = Max_age, 
      col = 2, 
      add = T)

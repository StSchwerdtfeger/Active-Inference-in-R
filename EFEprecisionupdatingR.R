###########################################################
#     Example code for simulated expected free energy     #
#             precision (beta/gamma) updates              #
# (associated with dopamine in the neural process theory) #
###########################################################

# R conversion of: Supplementary Code for: A Step-by-Step Tutorial
# on Active Inference Modelling and its Application to Empirical Data
# By: Ryan Smith, Karl J. Friston, Christopher J. Whyte

# This script will reproduce the simulation results in Figure 9
# Following the example from p. 71ff.

softmax <- function(probs){
  if(is.null(dim(probs))) probs <- matrix(probs,ncol= length(probs))
  exp(probs)/apply(probs,1, function(x) sum(exp(x)))
}

# Here you can set the number of policies and the distributions that
# contribute to prior and posterior policy precision:

# Set a fixed-form prior distribution over policies (habits)
E = c(1, 1, 1, 1, 1)                            

# Set an example expected free energy distribution over policies:
G = c(12.505, 9.51, 12.5034, 12.505, 12.505)

# Set an example variational free energy distribution over policies
# after a new observation:
VFE = c(17.0207, 1.7321, 1.7321, 17.0387, 17.0387) # F = not possible in R!

gamma_0 = 1            # Starting expected free energy precision value
gamma = gamma_0        # Initial expected free energy precision to be updated
beta_prior = 1/gamma   # Initial prior on expected free energy precision
beta_posterior = beta_prior # Initial posterior on expected free energy precision

psi = 2                # Step size parameter (promotes stable convergence) 

varup = 16             # number of variational updates (16)

# Set up lists to store values:
gamma_dopamineLIST <- list()  
gamma_dopamineLIST

policies_neuralLIST <- list()
policies_neuralLIST

for (ni in 1:varup){  
  
  # calculate prior and posterior over policies 
  # (see main text for explanation of equations) 
  
  # Prior over policies
  pi_0 = exp(log(E) - (gamma%*%G))/(sum(exp(log(E) - (gamma%*%G)))) 
  
  # Posterior over policies
  pi_posterior = softmax(log(E)-(gamma%*%G)-VFE) 
  
  # Calculate expeted free energy precision:
  G_error = (pi_posterior - pi_0)%*%(-G) # somewhat different to Matlab concerening t()
  
  # change in beta: gradient of VFE with respect to gamma
  # (recall gamm = 1/beta!)
  beta_update = beta_posterior - beta_prior + G_error

  # update posterior precision estimate (with step size of psi = 2, 
  # which reduces the magnitude of each update and can promote
  # stable convergence) 
  beta_posterior = beta_posterior-beta_update/psi

  # update expected free energy precision
  gamma = 1/beta_posterior

  # simulate dopamine responses:
  
  # simulated neural encoding of precision (beta_posterior^-1) at 
  # each iteration of variational updating    
  gamma_dopamineLIST[[ni]] <- gamma 
  
  # neural encoding of posterior over policies at 
  # each iteration of variational updating
  policies_neuralLIST[[ni]] <- pi_posterior   
  
}

# Final policy prior:
pi_0

  #            [,1]      [,2]       [,3]       [,4]       [,5]
  # [1,] 0.02214869 0.9113612 0.02219271 0.02214869 0.02214869

# Final policy posterior:
pi_posterior

  #              [,1]      [,2]       [,3]         [,4]         [,5]
  # [1,] 5.438184e-09 0.9762277 0.02377229 5.341173e-09 5.341173e-09

# Final difference vector:
diffvec =  pi_posterior - pi_0
diffvec

#             [,1]       [,2]        [,3]        [,4]        [,5]
# [1,] -0.02214868 0.06486647 0.001579574 -0.02214868 -0.02214868


# Negativ expected free energy:
NegEFE = -G
NegEFE

# [1] -12.5050  -9.5100 -12.5034 -12.5050 -12.5050

# Prior G Precision (Prior Gamma):
gamma_0

# [1]  1

# Posterior G Precision (Gamma):
gamma

#          [,1]
# [1,] 1.241122


########################
# Gamma Dopamine Plot: #
########################

# Add prior value: 
# Create matrix of added gamma_0 
gammaADDprior = matrix(c(1,1,1), nrow = 1, ncol=3, byrow=TRUE)
gammaADDprior

# Use rbind() to combine and eventually add the priors.
gamma_dopamine_preplot= rbind(t(gamma_dopamineLIST))
gamma_dopamine_preplot
gamma_dopamine_plot = rbind(t(gammaADDprior), t(gamma_dopamine_preplot))
gamma_dopamine_plot

# Caculate gradients for plot:

library(pracma)
rbind(gamma_dopamine_plot[,1])

gamma_dopa_grad=gradient(unlist(gamma_dopamine_plot))
gamma_dopa_grad

# Plot Rate of Change in Precision (Phasic Dopamine):

plot(gamma_dopa_grad, type = "l", main = "Rate of Change in Precision (Phasic Dopamine)",
     xlab = "Updates", ylab = "Gamma Gradient")


########################################################
# Plot Expected Free Energy Precision (Tonic Dopamine) #
########################################################

plot(gamma_dopamine_plot, type = "l", main = "Plot Expected Free Energy Precision (Tonic Dopamine)",
     xlab = "Updates", ylab = "Gamma")


###################################################
# Firing rates encoding beliefs about each policy #
###################################################

# plot firing rates encoding beliefs about each
# policy (columns = policies, rows = updates over time)

policies_neural_plot = rbind(t(policies_neuralLIST))
policies_neural_plot


# Converting into plotable form:
pn_plot_pre=rbind(policies_neural_plot)
pn_plot <- data.frame(matrix(unlist(pn_plot_pre), nrow=5, ncol=16,byrow=FALSE),stringsAsFactors = FALSE)
x=as.matrix(pn_plot)

# Not very elegant code, but I will optimize this in the future
plot(x[,1], type = "l", main = "Firing rates encoding beliefs about each policy",
     xlim=c(1,5))
lines(x[,2]);lines(x[,3]);lines(x[,4]);lines(x[,5]);lines(x[,6]);
lines(x[,7]);lines(x[,8]);lines(x[,9]);lines(x[,10]);lines(x[,11]);
lines(x[,12]);lines(x[,13]);lines(x[,14]);lines(x[,15]);lines(x[,16])



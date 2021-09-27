########################################################################
#----------------------------------------------------------------------#
# Example code for simulating state and outcome prediction errors in R #
#----------------------------------------------------------------------#
########################################################################

# R conversion of the Supplementary Code Matlab for: 
# "A Step-by-Step Tutorial on Active Inference Modelling and its 
# Application to Empirical Data", by Ryan Smith, Karl J. Friston, 
# Christopher J. Whyte

# Following page 44.

# Softmax function for matrices (!) from https://gist.github.com/bnicenboim/6fdf8f64a5840b74ae3b279e2a56dd07
# (Different function than in the "Pencil and Paper" script)
# The matlab script uses exp(x)/sum(exp(x)) and spm_softmax(x,k).

softmax <- function(probs){
  if(is.null(dim(probs))) probs <- matrix(probs,ncol= length(probs))
  exp(probs)/apply(probs,1, function(x) sum(exp(x)))
}

# Non-Zero-Log
nonZ = exp(-16)


#### Set up model to calculate state prediction errors:
##

# A = likelihood = matrix:

A = matrix( c(.8, .4,
              .2, .6),
            nrow = 2,
            ncol = 2,
            byrow = TRUE)
A


# B transition matrix (Note that the previous transition matrix
# functions as a prior):
# transition prior from previous timestep

Bt1 = matrix( c(.9, .2, 
                .1, .8),
            nrow = 2,
            ncol = 2,
            byrow = TRUE)
Bt1


# B transition matrix:
# transition prior from current timestep

Bt2 = matrix( c(.2, .3, 
                .8, .7),
            nrow = 2,
            ncol = 2,
            byrow = TRUE)
Bt2

# Observation
o = as.matrix(c(1,0))

# Prior distrubution over states (note that we will use the same value
# for sPItauminus1/plus1/tau)
sPItau = c(.5, .5)
sPItaumin1 = c(.5, .5)
sPItauplus1 = c(.5, .5)

# Depolarisation term (initial value)
vO = log(sPItau)
  
# Normalize Bt2 vis softmax function. Does not need to be transposed 
# different to the matlab code.  
Bt2cross = softmax(Bt2)
Bt2cross

#           [,1]      [,2]
# [1,] 0.4750208 0.5249792
# [2,] 0.5249792 0.4750208


#### Calculate state prediction error (single iteration):
##

stateError = (.5*((log(Bt1%*%sPItaumin1) + log(Bt2cross%*%sPItauplus1))) 
+ log(t(A)%*%o) - log(sPItau))
stateError

#            [,1]
# [1,] -0.1754885
# [2,] -0.9689710


# Depolarisation
v = vO+stateError
v

#            [,1]
# [1,] -0.8686356
# [2,] -1.6621182


# Updated Distribution over states:
# Via softmax
s = (exp(v)/sum(exp(v)))
s

#           [,1]
# [1,] 0.6885786
# [2,] 0.3114214



#### Set up model to calculate outcome prediction errors
## This minimizes expected free energy (maximizes reward and 
## information-gain)

### Calculate risk (reward-seeking) term under two policies

A = matrix(c(.9, .1,
             .1, .9),
           ncol = 2,
           nrow = 2,
           byrow = TRUE)

# State under policy 1
S1 = c(.9, .1)

# State under polcy 2
S2 = c(.5, .5)

# Preferred outcome
C = c(1, 0)

# Predicted poutcome under policy 1 (here dot prduct, different to Matlab)
o1 = A%*%S1
o1

# Predicted poutcome under policy 2 ( - " - )
o2 = A%*%S2
o2

# Add small number to avoid log(0)
z = exp(-16)

# Dotproduct:     https://www.statology.org/dot-product-in-r/
# Given a vector a = [a1, a2, a3] and a vector b = [b1, b2, b3]
# the dot product =  a · b  =  a1 * b1  +  a2 * b2  +  a3 * b3

# Risk under policy 1
risk1 = sum(o1*(log(o1)-log(C+z)))
risk1

# [1] 2.408606

# Risk under policy 1
risk2 = sum(o2*(log(o2)-log(C+z)))
risk2

# [1] 7.306853


### Calculate ambiguity (information-seeking) term under two policies:

A = matrix(c(.4, .2,
             .6, .8),
           ncol = 2,
           nrow = 2,
           byrow = TRUE)

# State under policy 1
S1 = c(.9, .1)

# State under policy 2
S2 = c(.1, .9)

# Ambiguity 1
x = diag(t(A)%*%log(A))
amb1 = sum((x)*S1)
amb1      # [1] -0.6557507

# Ambiguity 2
x = diag(t(A)%*%log(A))
amb2 = sum((x)*S2)
amb2      # [1] -0.5176633

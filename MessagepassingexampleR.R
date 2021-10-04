##################################
#--------------------------------#
#  Message Passing Examples in R #
#--------------------------------#
##################################

# CONVERSION of MATlab code Message_passing_example.m, including commentaries.

# Supplementary Code for: A Tutorial on Active Inference 
# Modelling and its Application to Empirical Data by Ryan Smith 
# and Christopher J. Whyte. We also acknowledge Samuel Taylor for contributing 
# to this example code

# This script provides two examples of (marginal) message passing, based on
# the steps described in the main text. Each of the two examples (sections)
# need to be run separately. The first example fixes all observed
# variables immediately and does not include variables associated with the
# neural process theory. The second example provides observations
# sequentially and also adds in the neural process theory variables. To
# remind the reader, the message passing steps in the main text are:
  
# 	1. Initialize the values of the approximate posteriors q(s_(?,?) ) 
#      for all hidden variables (i.e., all edges) in the graph. 

# 	2. Fix the value of observed variables (here, o_?).

# 	3. Choose an edge (V) corresponding to the hidden variable you want to 
#      infer (here, s_(?,?)).

# 	4. Calculate the messages, ?(s_(?,?)), which take on values sent by 
#      each factor node connected to V.

# 	5. Pass a message from each connected factor node N to V (often written 
#      as ?_(N?V)). 

# 	6. Update the approximate posterior represented by V according to the 
#      following rule: q(s_(?,?) )? ? ?(s_(?,?))? ?(s_(?,?)). The arrow 
#      notation here indicates messages from two different factors arriving 
#      at the same edge. 
#        6A. Normalize the product of these messages so that q(s_(?,?) ) 
#            corresponds to a proper probability distribution. 
#        6B. Use this new q(s_(?,?) ) to update the messages sent by 
#            connected factors (i.e., for the next round of message passing).

# 	7. Repeat steps 4-6 sequentially for each edge.

# 	8. Steps 3-7 are then repeated until the difference between updates 
#      converges to some acceptably low value (i.e., resulting in stable 
#      posterior beliefs for all edges). 

library(pracma) # for grad(), numerical gradiant Matlab-Style


##########################################################
# Example 1: Fixed observation and message passing steps #
##########################################################

# This section carries out marginal message passing on a graph with beliefs
# about states at two time points. In this first example, both observations 
# are fixed from the start (i.e., there are no ts as in full active inference
# models with sequentially presented observations) to provide the simplest
# example possible. We also highlight where each of the message passing
# steps described in the main text are carried out.

# Note that some steps (7 and 8) appear out of order when they involve loops that
# repeat earlier steps

library(dplyr)

# Softmax function
softmax <- function(par){             # THANKS TO: https://rpubs.com/FJRubio/softmax
  n.par <- length(par)                # equivalent: exp(result)/sum(exp(result))
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

# Specify generative model and initialize variables

#  priors

D = c(.5, .5)

# likelihood mapping
A = matrix(c(.9, .1,
             .1, .9),
           nrow = 2,
           ncol = 2,
           byrow= TRUE) 
 
# transitions
B = matrix(c(1, 0,
             0, 1),
            nrow = 2,
            ncol = 2,
            byrow = TRUE)

# number of timesteps
Timepoints = 2

# number of iterations of message passing
NumIterations = 16

# at exp(-16) to log inputs as log(0) is not defined
z = exp(-16)

# Initialize posterior STEP 1
qs = c(.5,.5)

# List to store results
qsLIST <- vector("list", 16*2)
dim(qsLIST) <- matrix(c(16,2))


for (Ni in 1:NumIterations){
  # fix observations (STEP 2) [using a list]
  o<- vector("list", 1*2)
  dim(o) = matrix(c(1,2))
  # Observations
  o[[1]] = matrix(c(1, 0))
  o[[2]] = matrix(c(1, 0))

    # For each edge (hidden state) - STEP 7 
    for (tau in 1:Timepoints){
      # Choose an edge
      q = log(qs[tau]+z)   
      
      # compute MEssage sent by D and B (STEP 4) using the posterior
      # computed in 6B
      if (tau == 1){  # first time point
        lnD = log(D+z)  # Message 1
        lnBs= log(t(B)%*%qs[]+z)
      }
      else if (tau == Timepoints){
        lnBs= log(B%*%qs[]+z)
      }
      lnAo = log((t(A)%*%as.vector(o[[tau]]))+z)  # as.vector(o...) 
      
      # STEPS 5-6 (Pass messages and update the posterior)
      # Since all terms are in log space, this is addition instead of
      # multiplication. This corresponds to equation 16 in the main
      # text (within the softmax)
      if (tau == 1){
        q = .5*lnD+.5*lnBs + as.vector(lnAo)
      }
      else if (tau == Timepoints){
        q = .5*(lnBs) + as.vector(lnAo)
      }
      # normalize using a softmax function to find posterior (STEP 6A)
      qs = softmax (q)
      qsLIST[[Ni,tau]]<- qs # store value for each iteration
    } 
  #  Repeat for remaining edges (STEP 7)
print(qs) # Repeat until convergence/for fixed number of iterations (STEP 8)
#print(qsResults)
}

# Final posterior beliefs over states
qs

# qsLIST of each iteration, for each tau. 
qsLIST


###################################
# Plotting firing rates (traces): #
###################################

# Adding prior to starting values:
addD = vector("list", 1*2)
dim(addD) = matrix(c(1,2))
addD[[1,1]] <- c(.5,.5)
addD[[1,2]] <- c(.5,.5)
addD
# Use rbind to place D up front which results in a new list of dim 17*2
qsLISTplot = rbind(addD, qsLIST)

# Transform upper list into a dataframe. Not that the list neeeds
# to be transposed for this via t(a). This results in a dataframe
# with a dim of 17*4
qsPLOT <- data.frame(matrix(unlist(t(qsLISTplot)), nrow=17, ncol=4,byrow=TRUE),stringsAsFactors = TRUE)
qsPLOT

# Plot: Approximate posteriors (1 per edge per time point)
plot(qsPLOT$X1, type = "l", col= "green", 
     ylim = c(0,1), # Set dimension y axis 0 to 1
    xlab = "Message passing iterations", ylab = "qs_tau") # Name labels
lines(qsPLOT$X2, col="blue")     # add the other columns of qsPLOT
lines(qsPLOT$X3, col="red")
lines(qsPLOT$X4, col="yellow")


#####################################################
# Example 2: Sequential observations and simulation #
#             of firing rates and ERPs              #
#####################################################

# This script performs state estimation using the message passing 
# algorithm introduced in Parr, Markovic, Kiebel, & Friston (2019).
# This script can be thought of as the full message passing solution to 
# problem 2 in the pencil and paper exercises. It also generates
# simulated firing rates and ERPs in the same manner as those shown in
# figs. 8, 10, 11, 14, 15, and 16. Unlike example 1, observations are
# presented sequentially (i.e., two ts and two taus).

# Specify generative model and initialise variables

#  priors

D = c(.5, .5)

# likelihood mapping
A = matrix(c(.9, .1,
             .1, .9),
           nrow = 2,
           ncol = 2,
           byrow= TRUE) 

# transitions
B = matrix(c(1, 0,
             0, 1),
           nrow = 2,
           ncol = 2,
           byrow = TRUE)

# number of timesteps
Timesteps = 2

# number of iterations of message passing
NumIterations = 16

# at exp(-16) to log inputs as log(0) is not defined
z = exp(-16)

# Initialize posterior STEP 1
qs2 = c(.5,.5)

# List to store results
xn <- vector("list", 16*2*2)
dim(xn) <- matrix(c(16,2,2))


# Message passing:

for (t in 1:Timesteps){
  
      for (Ni in 1:NumIterations){   # STEP 8 loop of VMP
          
          # Likelihood
          #lnAo = log((t(A)%*%as.vector(o[[t,tau]]))+z)  # as.vector(o...) + multiple places in the loop necessary!!
  
          for (tau in 1:Timesteps){ # STEP 7 loop of VMP
            # fix observations sequentially (STEP 2) [using a list of lists]
            # Each tau for Timestep1
            o1<- vector("list", 1*2)
            dim(o1) = matrix(c(1, 2))
            o1[[1]] = matrix(c(1, 0))
            o1[[2]] = matrix(c(0, 0))
            
            # Each tau for Timestep1
            o2<- vector("list", 1*2)
            dim(o2) = matrix(c(1,2))
            o2[[1]] = matrix(c(1, 0))       
            o2[[2]] = matrix(c(1, 0))               
            
            # list of lists
            o <- vector("list", 2*2)
            dim(o) = matrix(c(2,2))
            o[1,] <- o1 # assign to list
            o[2,] <- o2 #      -''-

          # initialise depolarization variable: v = ln(s)
          # choose an edge (Step 3 of VMP)
          v = log(qs2[])
            
            # get correct D and B for each time point (Steps 4-5 of VMP)
            # using using the posterior computed in Step 6B  
            if (tau == 1){  # first time point
            lnD = log(D+z)  # pas(Message 1)
            
            lnBs= log(t(B)%*%qs2[]+z) # future (Message 2)
           }
            else if (tau == Timesteps){ # last time point
            lnBs= log(B%*%qs2[]+z)      # no contribution from future (only
                                        # Message 1)
           } 
             # Likelihood Message 3
            lnAo = log((t(A)%*%as.vector(o[[t,tau]]))+z)  # as.vector(o...) + multiple places in the loop as below!!
    
             # calculate state prediction error: equation 24
            if (tau == 1){
            epsilon = .5*lnD+.5*lnBs + as.vector(lnAo) - v
           }
            else if (tau == Timesteps){
            epsilon = .5*(lnBs) + as.vector(lnAo) - v
      }
            
            # (Step 6 of VMP)
            # update depolarization variable: equation 25
            v = v + epsilon[]

            # normalize using a softmax function to find posterior (STEP 6A of VMP)
            # equation 26 (Step 6A of VMP)
            qs2 = softmax (v)
            xn[[Ni,tau,t]]<- qs2 # store qs for firing rate plots
      } 
           #  Repeat for remaining edges (STEP 7 of VMP)
    }  # Repeat until convergence/for number of iterations (Step 8 of VMP)
  print(qs2)
}

# Final posterior beliefs over states
qs2

#            [,1]
# [1,] 0.98780474
# [2,] 0.01219526

# xn of each iteration, for each t and each tau. 
xn


###################################
# Plotting firing rates (traces): #
###################################


# Note: the Matlab code uses functions to transform xn into
# "firing_rate". A different approach to matlab is used here.

# Adding prior to starting values:
addD2 = vector("list", 1*2*1)
dim(addD2) = matrix(c(1,2,1))
addD2[[1,1,1]] <- c(.5,.5)
addD2[[1,2,1]] <- c(.5,.5)
addD2

# Use rbind to place D up front of [,,1] which results in a new list [x,y]
xnplotpre = rbind(addD2[,,1], xn[,,1])
xnplotpre[1,2]
# Use rbind to combine lists
xnplotpre2 = rbind(xnplotpre, xn[,,2])
xnplotpre2

# Transform upper list into a dataframe. Not that the list neeeds
# to be transposed for this via t(a). This results in a dataframe
# with a dim of 17*4
xnplotpredf <- data.frame(matrix(unlist(t(xnplotpre2)), nrow=4, ncol=33,byrow=FALSE),stringsAsFactors = FALSE)
xnplot = t(xnplotpredf)
xnplot

# Plot firing rate (traces)
plot(xnplot[,1], type = "l", col = "green",
     ylim = c(0,1), # Set dimension y axis 0 to 1
     xlab = "Message passing iterations", ylab = "Firing rates") # Name labes
lines(xnplot[,2], col="blue")     # add the other columns of qsPLOT
lines(xnplot[,3], col="red")
lines(xnplot[,4], col="yellow")


library(pracma) # for gradient(), numerical gradiant Matlab-Style

num_states = 2
num_epochs = 2
time_tau = matrix(c(1, 2, 1, 2,
                    1, 1, 2, 2),
                  ncol =4,
                  nrow=2,
                  byrow=TRUE) 

# Sort of, but not the right function yet...
for (t_tau in 1:ncol(time_tau)){
    for (epoch in 1:num_epochs){
    # firing rate 
     firing_rate[] = unlist(xn)
      ERP = gradient(firing_rate)
    }
}


# Plotting ERP - still not working right, but got a little further:
grad1=gradient(xnplot[,1])
grad1
grad1=rbind(c(0),as.matrix(grad1)) # add 0 to starting value:

grad2=gradient(xnplot[,2])
grad2
grad2=rbind(c(0),as.matrix(grad2)) #        --   "   --


grad3=gradient(xnplot[,3])
grad3
grad3=rbind(c(0),as.matrix(grad3))

grad4=gradient(xnplot[,4])
grad4
grad4=rbind(c(0),as.matrix(grad4))


plot(grad4, ylim=c(-0.5,0.5), type = "l")
lines(grad3)
lines(grad2)
lines(grad1)

# Firing rates:
image((64*as.matrix(xnplot)), col = grey.colors(256))





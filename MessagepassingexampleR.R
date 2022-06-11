##################################
#--------------------------------#
#  Message Passing Examples in R #
#--------------------------------#
##################################

# Supplementary Code for: A Tutorial on Active Inference 
# Modelling and its Application to Empirical Data by Ryan Smith 
# and Christopher J. Whyte. We also acknowledge Samuel Taylor for contributing 
# to this example code.

# CONVERSION of Matlab code Message_passing_example.m. The commentary
# is included as well and extendet mostly by those that refer to R code
# exclusively. Conversion by Steffen Schwerdtfeger https://github.com/StSchwerdtfeger

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

library(pracma) # for gradient(), numerical gradiant Matlab-Style


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

# Add exp(-16) to log inputs as log(0) is not defined.
# A function is usually defined for this purpose. 
# However, we will just add the term z to every log-input.
z = exp(-16)

for (Ni in 1:NumIterations){ # START of loop
  if (Ni == 1){ # set up qs and corresponding list
    # Initialize posterior STEP 1
    qs = c(.5,.5)
    
    # List to store results
    qsLIST <- vector("list", 16*2)
    dim(qsLIST) <- matrix(c(16,2))
    
    # fix observations (STEP 2) [using a list]
    o<- vector("list", 1*2)
    dim(o) = matrix(c(1,2))
    # Observations
    o[[1]] = matrix(c(1, 0)) # 
    o[[2]] = matrix(c(1, 0)) # same observation for each tau!
    
  }
  
  # For each edge (hidden state) - STEP 7 
  for (tau in 1:Timepoints){
    # Choose an edge
    q = log(qs[tau]+z)   
    
    # compute Message sent by D and B (STEP 4) using the posterior
    # computed in 6B
    if (tau == 1){  # first time point
      lnD = log(D+z)  # Message 1
      lnBs= log(t(B)%*%qs[]+z)
    }
    else if (tau == Timepoints){
      lnBs= log(B%*%qs[]+z)
    }
      lnAo = log((t(A)%*%as.vector(o[[tau]]))+z)  # Likelihood (as.vector(o...)) 
      
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
  # Repeat until convergence/for fixed number of iterations (STEP 8)
  print(qs) 
} # END of loop!

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
# (((The Matlab script uses: qs_plot = [D' D';qs(:,:,1) qs(:,:,2)])))
qsLISTplot = rbind(addD, qsLIST)

# Transform upper list into a dataframe to use $. Note that the elememts of
# the list have to be unlisted and to be transposed for this fist via t(a). 
# This results in a dataframe with a dim of 17*4. 
qsPLOT <- data.frame(matrix(unlist(t(qsLISTplot)), nrow=17, ncol=4,byrow=TRUE),stringsAsFactors = TRUE)
qsPLOT

# Plot: Approximate posteriors (1 per edge per time point)
plot(qsPLOT$X1, type = "l", col= "green", 
     ylim = c(0,1), # Set dimension y axis 0 to 1
    xlab = "Message passing iterations", ylab = "qs_tau") # Name labels
lines(qsPLOT$X2, col="blue")     # add the other columns of qsPLOT
lines(qsPLOT$X3, col="red")
lines(qsPLOT$X4, col="yellow")

# GREEN + BLUE: tau = 1
# RED + YELLOW: tau = 2

# see p. 29, Fig. 8 "Firing rates (upper-right) correspond to the 
# magnitude of posteriors over each state[...]".
# RECALL that the observation o was the same over all t.


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

for (t in 1:Timesteps){   # START of loop
  if(t == 1){
    # Initialize posterior STEP 1
    qs2 = c(.5,.5)
    
    # List to store results
    xn <- vector("list", 16*2*2)
    dim(xn) <- matrix(c(16,2,2))
    
    # fix observations sequentially (STEP 2) [using a list of lists]
    # Each tau for Timestep 1
    o1<- vector("list", 1*Timesteps)
    dim(o1) = matrix(c(1, 2))
    o1[[1]] = matrix(c(1, 0)) # o at t = 1, tau = 1
    o1[[2]] = matrix(c(0, 0)) # o at t = 1, tau = 2
    
    # Each tau for Timestep 2
    o2<- vector("list", 1*Timesteps)
    dim(o2) = matrix(c(1,2))
    o2[[1]] = matrix(c(1, 0)) # o at t = 1, tau = 1       
    o2[[2]] = matrix(c(1, 0)) # o at t = 2, tau = 2              
    
    # list of lists for all o
    o <- vector("list", Timesteps*Timesteps)
    dim(o) = matrix(c(2,2))
    o[1,] <- o1 # assign to list
    o[2,] <- o2 #      -''-
  } # END if t == 1 for setup lists etc.
  
  for (Ni in 1:NumIterations){             # STEP 8 loop of VMP
      
    for (tau in 1:Timesteps){              # STEP 7 loop of VMP
      # initialise depolarization variable: v = ln(s)
      # choose an edge (Step 3 of VMP)
      v = log(qs2) # different to Matlab qs2 is always the current t 
                   # and not obtained from a list such as e.g. qs2[[t]]
            
      # get correct D and B for each time point (Steps 4-5 of VMP)
      # using using the posterior computed in Step 6B  
      if (tau == 1){                  # first time point
        lnD = log(D+z)                # past (Message 1)
        lnBs= log(t(B)%*%qs2+z)       # future (Message 2)
      } 
      else if (tau == Timesteps){     # last time point
        lnBs= log(B%*%qs2+z)          # no contribution from future (only
                                      # Message 1)
      } 
      
      lnAo = log((t(A)%*%as.vector(o[[t,tau]]))+z)  # Likelihood/present 
                                                    #    => Message 3
    
      # Calculate state prediction error: equation 24
      if (tau == 1){
        epsilon = .5*lnD + .5*lnBs + as.vector(lnAo) - v
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
      } # END LOOP OVER tau
      #  Repeat for remaining edges (STEP 7 of VMP)
    } # END LOOP OVER Ni 
    # Repeat until convergence/for number of iterations (Step 8 of VMP)
  print(qs2)
} # END of loop (over t)

# Final posterior beliefs over states
qs2

#            [,1]
# [1,] 0.98780474
# [2,] 0.01219526

xn


##################################
# ERP - Event Related Potentials #
##################################


# Getting Data in Form. In the Matlab code this is, e.g., done
# via implemented functions like spm_cat... # https://github.com/neurodebian/spm12/blob/master/spm_cat.m

library(pracma) # library for gradient(), calculates numerical gradient 
                # in Matlab-Style. gradient() for central difference!


# Example:
x= 1:10
# 1 2 3 4 5 6 7 8 9 10
gradient(x)
# 1 1 1 1 1 1 1 1 1 1 


# Set up List for ERP plot:
ERPplot=vector("list", 2*4)
dim(ERPplot) = matrix(c(2,4))
ERPplot

# Grad 1
grad1pre = rbind(xn[,,1])
grad1predf <- data.frame(matrix(unlist(t(grad1pre)), nrow=4, ncol=16,byrow=FALSE),stringsAsFactors = FALSE)
grad1preMX=as.matrix(grad1predf)
grad1preMX

grad1grad1=gradient(grad1preMX[1,])
grad1grad1
ERPplot[[1,1]]= grad1grad1 # Add to ERP list

grad1grad2=gradient(grad1preMX[2,])
grad1grad2
ERPplot[[1,2]]= grad1grad2

grad1grad3=gradient(grad1preMX[3,])
grad1grad3
ERPplot[[1,3]]= grad1grad3

grad1grad4=gradient(grad1preMX[4,])
grad1grad4
ERPplot[[1,4]]= grad1grad4

# Grad 2
grad2pre = rbind(xn[,,2])
grad2predf <- data.frame(matrix(unlist(t(grad2pre)), nrow=4, ncol=16,byrow=FALSE),stringsAsFactors = FALSE)
grad2preMX=as.matrix(grad2predf)
grad2preMX

grad2grad1=gradient(grad2preMX[1,])
grad2grad1
ERPplot[[2,1]] = grad2grad1 # add to list
  
grad2grad2=gradient(grad2preMX[2,])
grad2grad2
ERPplot[[2,2]] = grad2grad2 

grad2grad3=gradient(grad2preMX[3,])
grad2grad3
ERPplot[[2,3]] = grad2grad3 

grad2grad4=gradient(grad2preMX[4,])
grad2grad4
ERPplot[[2,4]] = grad2grad4 

# ERPplot as data.frame, in order to ADD ZERO value!
ERP =  data.frame(matrix(unlist(ERPplot), 
                                    nrow=32, 
                                    ncol=4,
                                    byrow=FALSE),
                             stringsAsFactors = FALSE)

ERPplusZERO = rbind(c(0,0,0,0), ERP) # Add Zero!
ERPplusZERO  

# FINALLY plot ERP:

plot(ERPplusZERO[,1], type = "l", col = "green",
     ylim = c(-.05,.05), # Set dimension y axis 0 to 1
     xlab = "Message passing iterations", ylab = "Response", # Name labels
             main = "Event related potentials")
lines(ERPplusZERO[,2], col="blue")     # add the other columns of xnplot
lines(ERPplusZERO[,3], col="red")
lines(ERPplusZERO[,4], col="orange")



###################################
# Plotting firing rates (traces): #
###################################


# Note: the Matlab code uses functions to transform xn into
# "firing_rate". A different approach as to Matlab is used here.

# Adding prior to starting values:
addD2 = vector("list", 1*2*1)
dim(addD2) = matrix(c(1,2,1))
addD2[[1,1,1]] <- c(.5,.5)
addD2[[1,2,1]] <- c(.5,.5)
addD2

# Use rbind to place D up front of [,,1] which results in a new list [x,y]
xnplotpre = rbind(addD2[,,1], xn[,,1]) 
xnplotpre

# Use rbind to combine lists
xnplotpre2 = rbind(xnplotpre, xn[,,2])
xnplotpre2

# Transform upper list into a dataframe. Note that the list needs
# to be transposed for this via t(a). This results in a dataframe
# with a dim of 17*4
xnplotpredf <- data.frame(matrix(unlist(t(xnplotpre2)), nrow=4, ncol=33,byrow=FALSE),stringsAsFactors = FALSE)
xnplot = t(xnplotpredf)
xnplot

# Plot firing rate (traces)
plot(xnplot[,1], type = "l", col = "green",
     ylim = c(0,1), # Set dimension y axis 0 to 1
     xlab = "Message passing iterations", ylab = "Firing rates") # Name labes
lines(xnplot[,2], col="blue")     # add the other columns of xnplot
lines(xnplot[,3], col="red")
lines(xnplot[,4], col="yellow")



#########################
# Plotting Firing rates #
#########################


library(matlab) # for imagesc() 

firepre= rbind(xn[,,1], xn[,,2])
firepre
# Add prior D:
prior = vector("list", 1*2)
dim(prior) = matrix(c(1,2))
priorD = matrix(c(.5,.5), ncol=1, nrow=2, byrow=TRUE)
prior[[1,1]] = priorD
prior[[1,2]] = priorD
prior 

firepre2=rbind(prior, firepre)
firepre2

firepre3=rbind(firepre2)
firepre3

# Convert to df to make data plotable:
fireplot <- data.frame(matrix(unlist(t(firepre3)), nrow=4, ncol=33,byrow=FALSE),stringsAsFactors = FALSE)
fireplot

# FINALLY: Plot fire rate!

imagesc(1-as.matrix(fireplot), col = grey.colors(256))


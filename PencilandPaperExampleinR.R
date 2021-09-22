###################################################
##-----------------------------------------------##
## Code/solutions for pencil and paper exercises ##
##-----------------------------------------------##
###################################################


# Conversion of the Supplementary Code for: 
# "A Step-by-Step Tutorial on Active Inference Modelling and its Application 
# to Empirical Data, by Ryan Smith, Karl J. Friston, Christopher J. Whyte"



#####################
#-------------------#
# Static perception #
#-------------------#
#####################


# We will work with this equation
# s = σ(log(prior) + log(t(likelihood)%*%obs))

# As log(0) is not defined, 0.01 was added to all log-input (usually exp(-16)).
# Note: The example 1 in the paper on p. 134 did not add .01 to logs.
# The matlab code may will (keep in mind, if the answers in the paper
# may differ from here, or from the matlab skript).

# Observations
obs        = c(1, 0)         # VECTOR

# Our Prior:                 # MATRIX/VECTOR
prior      = matrix(
  c(.75,.25),
  nrow = 2,
  ncol = 1,
  byrow = TRUE)


likelihood = matrix(         # MATRIX
  c(.8, .2, 
    .2, .8),           # the data elements 
  nrow=2,              # number of rows 
  ncol=2,              # number of columns 
  byrow = TRUE)        # fill matrix by rows 


# express generative model in terms of update equations
logs = log(prior+.01) + log(t(likelihood +.01)%*%obs)
logs

# Normalized via softmax function:
s = exp(logs+.01)/(sum(exp(logs+.01)))
s

# Alternative: via this function:
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

# Use formula directly:
s = softmax( log(prior+.01) + log(t(likelihood)%*%obs+.01) )
s



########################
#----------------------#
# DYNAMICAL PERCEPTION #
#----------------------#
########################



# As log(0) is not defined, 0.01 was added to all log-input.
# (usually exp(-16))

# For the Pencil and paper example we use the 0.01!
nzlog = .01

# Alternative
nonzerolog = exp(-16)

# Our softmax function
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


# GENERATIVE MODEL (FOR EXAMPEL 2 Dyn. perception)

# D = prior = matrix:

D = matrix( c(.75, .25),
            nrow = 2,
            ncol = 1,
            byrow = TRUE)
D


# A = likelihood = matrix:

A = matrix( c(.8, .2,.2,.8),
            nrow = 2,
            ncol = 2,
            byrow = TRUE)
A


# B transition matrix

B = matrix( c(0, 1, 1, 0),
            nrow = 2,
            ncol = 2,
            byrow = TRUE)
B



##################
## Without loop ##
##################


# Only first and last equation applied, as there are only two timesteps given
# in the example.

# 𝑠𝜏=1=𝜎(12(ln𝐃+ln𝐁𝜏†*𝑠𝜏+1)+ln𝐀T𝑜𝜏) 
# 𝑠1<𝜏<𝑇=𝜎(12(ln𝐁𝜏−1s𝜏−1+ln𝐁𝜏†𝑠𝜏+1)+ln𝐀T𝑜𝜏) 
# 𝑠𝜏=𝑇=𝜎(12(ln𝐁𝜏−1s𝜏−1)+ln𝐀T𝑜𝜏)

# True observation at t-1 for tau 1 and 2.
ot11 = c(1, 0)
ot12 = c(0, 0)

# True observation at t-2 for ...
ot21 = c(1, 0)       
ot22= c(0, 1)               

# Initialize Approx. Posterior:
sINt = c(.5, .5)


# Equations:
# p. 135, or 31 (Bayes network...):

################################
# Model inversion time step 1! #
################################

# S t=1
st1 = softmax( ((.5*(log(D+nzlog)) )+ (.5*(log(t((B+nzlog))%*%sINt))) + ((log(t(A)+nzlog)%*%ot11))) )
st1  

# > st1
#           [,1]
# [1,] 0.8683268    
# [2,] 0.1316732


# For S t= 2

stall = softmax ( (.5*(log((B%*%st1)+nzlog))) +  (log((t(A)%*%ot12)+nzlog))    )
stall   
#           [,1]
# [1,] 0.2865401
# [2,] 0.7134599



###############################
# Model inversion time step 2 #
###############################


st21 = softmax( ((.5*(log(D+nzlog)) )+ (.5*(log(t((B+nzlog))%*%stall))) + ((log(t(A)+nzlog)%*%ot21))) )
st21  

#            [,1]
# [1,] 0.91150706   
# [2,] 0.08849294



st2all = softmax( (.5*(log((B%*%st21)+nzlog))) +  ((log(t(A)+nzlog)%*%ot22))    )
st2all

#            [,1]
# [1,] 0.07813653
# [2,] 0.92186347 # POSTERIOR OVER STATES


#############
#############
# With loop #
#############
#############

# This time for excercise 2 p. 136

# You may clear your environment for this, or make sure to address the
# right generative model:

# Prior

D = matrix( c(.5, .5),
            nrow = 2,
            ncol = 1,
            byrow = TRUE)
D


# A = likelihood = matrix:

A = matrix( c(.9, .1,
              .1, .9),
            nrow = 2,
            ncol = 2,
            byrow = TRUE)
A


# B transition matrix

B = matrix( c(1, 0, 
              0, 1),
            nrow = 2,
            ncol = 2,
            byrow = TRUE)
B




# Timesteps 
n = 2

# Initial approx posterior
qs= matrix(c(.5,.5))

# Result = all qs for all timesteps and each tau!
for(t in 1:n){
  for (tau in 1:n){
    o1<- vector("list", 1*2)
    # True observation at t 1
    dim(o1) = matrix(c(1, 2))
    o1[[1]] = matrix(c(1, 0))
    o1[[2]] = matrix(c(0, 0))
    
    o2<- vector("list", 1*2)
    dim(o2) = matrix(c(1,2))
    # True observation at t2
    o2[[1]] = matrix(c(1, 0))       
    o2[[2]] = matrix(c(1, 0))               
    
    # list of lists
    o <- vector("list", 2*2)
    dim(o) = matrix(c(2,2))
    
    o[1,] <- o1 # assign to list
    o[2,] <- o2 #      -''-
    logAo = log(t(A)%*%as.vector(o[[t,tau]])+.01) # if likelihood UP HERE DOESNT CHANGE RESULT?
    
    if (tau == 1){
      
      logD  <- log(D +.01)          # past
      logBs <- log(t(B)%*%(qs[,])+.01) # future
    }
    else if(tau == 2){
      logBs <- log(t(B)%*%(qs[,])+.01) # no contribution in future
    }
    logAo = log(t(A)%*%as.vector(o[[t,tau]])+.01) #likelihood
    if (tau == 1){
      logs = .5*logD + .5*logBs + as.vector(logAo)
    }
    else if (tau == 2) {
      logs = .5*logBs + as.vector(logAo)
    }
    qs = softmax(logs)   # Eqivalent: qs =  (exp(logs)/sum(exp(logs)))
    print(qs)}
}

#  Result for exercise 2
#          [,1]
# [1,] 0.7178485
# [2,] 0.2821515
#          [,1]
# [1,] 0.6121623
# [2,] 0.3878377
#          [,1]
# [1,] 0.9118586
# [2,] 0.0881414
#          [,1]
# [1,] 0.96205579
# [2,] 0.03794421


# Results for generative model of example 2 (not exercise 2) 
#           [,1]
# [1,] 0.8683268
# [2,] 0.1316732
#           [,1]
# [1,] 0.2865401
# [2,] 0.7134599
#           [,1] 
# [1,] 0.91150706
# [2,] 0.08849294
#           [,1]
# [1,] 0.07813653
# [2,] 0.92186347

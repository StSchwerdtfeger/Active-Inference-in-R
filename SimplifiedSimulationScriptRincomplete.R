#####################################
#-----------------------------------#
# Simplified Simulation Script in R #
#-----------------------------------#
#####################################

# Supplementary Code for: A Step-by-Step Tutorial on Active Inference 
# Modelling and its Application to Empirical Data

# By: Ryan Smith, Karl J. Friston, Christopher J. Whyte

############## Converted by Steffen Schwerdtfeger (incomplete so far!)

library(pracma)  # needed?
library(rlist)   # for list.append() function NEEDED?
library(imager)  # for squeeze()

# This code simulates a single trial of the explore-exploit task introduced 
# in the active inference tutorial using a stripped down version of the 
# modelinversion scheme implemented in the spm_MDP_VB_X.m script. 

# Note that this implementation uses the marginal message passing scheme
# described in (Parr et al., 2019), and will return very slightly 
# (negligably) different values than the spm_MDP_VB_X.m script in 
# simulation results.

# Parr, T., Markovic, D., Kiebel, S., & Friston, K. J. (2019). Neuronal 
# message passing using Mean-field, Bethe, and Marginal approximations. 
# Scientific Reports, 9, 1889.


## Functions:

# Reg. Softmax function
softmax <- function(par){             # THANKS TO: https://rpubs.com/FJRubio/softmax
  n.par <- length(par)                # EQUIVALENT: exp(result)/sum(exp(result))
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

# Alternative Softmax function, columnwise:
softmaxCOL <- function(probs){
  if(is.null(dim(probs))) probs <- matrix(probs,ncol= length(probs))
  exp(probs)/apply(probs,1, function(x) sum(exp(x)))
}

# natural log that replaces zero values with very small values for 
# numerical reasons (as log(0) is not defined).
nat_log = function (x) {
x = log(x+exp(-16))
}

# Isfield function, Matlab-Style
isfield  <- function(structure, field){
  any(field %in% names(structure) )
} 

# Normalizing vectors. So far this seems to be enough in R
# Compare with Matlab code for col_norm(B). 
# So far I used the same function for B_norm() in the Matlab script. 
col_norm = function(x){t(t(x)/colSums(x))}

# B_norm (1 and 2 necesserry????)
B_norm = function(x){
  bb = x
  z  = sum(bb)        # create normalizing constant from sum of columns
  bb = bb/z           # divide columns by constant
  bb[is.nan(bb)] = 0  # replace NaN with zero
  b = bb
  return(b)
}

B_norm2 = function(y){
  x=col_norm(t(y))
  x[is.nan(x)] = 0
  x
}

# B_norm3 (does both B_norm and B_norm2)
B_norm3 = function(x){
  if (ncol(x) == 1){
    bb = x
    z  = sum(bb)        # create normalizing constant from sum of columns
    bb = bb/z           # divide columns by constant
    bb[is.nan(bb)] = 0  # replace NaN with zero
    b = bb
    return(b)
  }
  else {
    x = B_norm2(x)
    z=t(x)
    return(z)
  }
}

# dot product along dimension f (maybe %*%?)
md_dot = function (A,s,f){ 
  if (f == 1){
    B = t(A)%*%s 
  }
  else if (f == 2){
    B = A%*%s
  }
}

# Select last element:
last <- function(x) { return( x[length(x)] ) }


# epistemic value term (Bayesian surprise) in expected free energy 
G_epistemic_value = function(X1,X2) {
  
  # probability distribution over the hidden causes: i.e., Q(s)
  
  qx = spm_cross1(X2) # this is the outer product of the posterior over states
  # calculated with respect to itself
  
  # accumulate expectation of entropy: i.e., E[lnP(o|s)]
  G     = 0
  qo    = 0
  find = which(qx > exp(-16))
  
  for (i in find){
    # probability over outcomes for this combination of causes
    po   = matrix(c(1))
    for (g in 1:numel(X1)){
      po = as.matrix(spm_cross2(po,X1[[g]][,i]))
    }  
    qo = qo + t(qx[[i]]%*%t(po))
    G  = G  + qx[[i]]*po*nat_log(po)
  }
  
  # subtract entropy of expectations: i.e., E[lnQ(o)]
  G  = G - qo*nat_log(qo)
  G[!(apply(G, 1, function(y) any(y == 0))),] # Delete rows with zero values (THIS MIGHT JUST BE WRONG TEST WITH OTHER EXAMPLES!!!)
  G = min(G)                                  # Otherwise one gets a vector 4x1, with 3 rows with "0"... 
} 


### SPM Functions
##

# For vectors:
spm_wnorm = function(A){
  # This uses an equiv function to Matlab bsxfun to subtract the 
  # inverse of each column entry from the inverse of the sum of the 
  # columns and then divide by 2. Pracma bsxfun does not work here,
  # therefore spmwnorm2 function for matrices. 
  
  A   = A + exp(-16)
  A   = arrayfun("-",1/sum(A),1/A)/2
} 

# For matrices
spm_wnorm2 = function(A){
  # This uses the bsxfun/arrayfun(mapply) function to subtract the 
  # inverse of each column entry from the inverse of the sum of the 
  # columns and then divide by 2.
  
  A   = A + exp(-16)
  A1 = 1/colSums(A)
  A2 = (1/A)
  A = mapply("-",A1,A2)/2
} 

spm_cross1 = function (INPUT) {
  
  X = as.matrix(INPUT[[1]])
  x = as.matrix(INPUT[[2]])
  
  # outer product of first pair of arguments
  #--------------------------------------------------------------------------
  A <- X
  dim(A[]) <- c(numel(X), ones(1,ndims(x)))
  B <- x
  dim(B) <- c(ones(1,ndims(X)),numel(x))
  Y = squeeze(A%*%B) # Squeeze actually not necessary in this example
}

spm_cross2 = function (X,x) {
  X = as.matrix(X)
  x = as.matrix(x)
  # outer product of first pair of arguments
  #--------------------------------------------------------------------------
  A2 <- X
  #A2[] <- reshape(A2,c(ones(1,ndims(X)),numel(x)))
  #A2 = reshape(as.matrix(A2),c(dim(X),ones(1,ndims(x))))
  dim(A2) <- c(dim(X), ones(1,ndims(X)))
  
  B2 <- x
  #B2 <- reshape(B2,c(ones(1,ndims(X)),dim(x)))
  dim(B2) <- c(ones(1,ndims(X)),numel(x))
  
  Y = as.matrix(A2)%*%t(as.matrix(B2)) # Squeeze actually not necessary in this example
  #Y = as.matrix(Y)
  #dim(Y) = matrix(c(numel(Y),1))
  Findim = as.matrix(Y)
  dim(Findim) = matrix(c(numel(Findim),1))
  Y = Findim
}


## Simulation Settings

# To simulate the task when prior beliefs (d) are separated from the 
# generative process, set the 'Gen_model' variable directly
# below to 1. To do so for priors (d), likelihoods (a), and habits (e), 
# set the 'Gen_model' variable to 2:
  
  Gen_model = 1 # as in the main tutorial code, many parameters 
                # can be adjusted in the model setup, within the 
                # explore_exploit_model function [[starting on line 810]](MAtlab code).
                # This includes, among others (similar to in the main 
                # tutorial script):
    
## Set up POMDP model structure
#=====================================================================
  
# Please note that the main tutorial script ('Step_by_Step_AI_Guide.m') 
# has more thorough descriptions of how to specify this generative 
# model and the other parameters that might be included. Below we 
# only describe the elements used to specify this specific model. 
# Also, unlike the main tutorial script which focuses on learning 
# initial state priors (d), this version also enables habits (priors
# over policies; e) and separation of the generative process from the 
# generative model for the likelihood function (a).
  

## Specify Generative Model
  
#MDP = explore_exploit_model(Gen_model)
  
# Model specification is reproduced at the bottom of this script 
# (starting on line 810), but see main tutorial script for more 
# complete walk-through


# Number of time points or 'epochs' within a trial: T
# =========================================================================
  
# Here, we specify 3 time points (T), in which the agent 1) starts 
# in a 'Start' state, 2) first moves to either a 'Hint' state or a 
# 'Choose Left' or 'Choose Right' slot machine state, and 3) either 
# moves from the Hint state to one of the choice states or moves from
# one of the choice states back to the Start state.

Time = 3

# Priors about initial states: D and d
# =========================================================================
  
#--------------------------------------------------------------------------
# Specify prior probabilities about initial states in the generative 
# process (D)
# Note: By default, these will also be the priors for the generative 
# model
#--------------------------------------------------------------------------

# Setup a list for D  
D = vector("list", 2*1)
dim(D) =  matrix(c(2,1))
D
 
# For the 'context' state factor, we can specify that the 'left better' 
# context (i.e., where the left slot machine is more likely to win)
# is the true context:
  
D[[1,1]] = t(matrix(c(1, 0)))  # c('left better','right better')

# For the 'behavior' state factor, we can specify that the agent 
# always begins a trial in the 'start' state (i.e., before choosing
# to either pick a slot machine or first ask for a hint:

D[[2,1]] = t(matrix(c(1, 0, 0, 0))) # c('start','hint','choose-left','choose-right')


#--------------------------------------------------------------------------
# Specify prior beliefs about initial states in the generative model 
# (d) Note: This is optional, and will simulate learning priors over 
# states if specified.
#--------------------------------------------------------------------------
  
# Note that these are technically what are called 'Dirichlet 
# concentration paramaters', which need not take on values between
# 0 and 1. These values are added to after each trial, based on 
# posterior beliefs about initial states. For example, if the agent
# believed at the end of trial 1 that it was in the 'left better' 
# context, then d{1} on trial 2 would be d{1} = [1.5 0.5]' (although
# how large the increase in value is after each trial depends on a 
# learning rate). In general, higher values indicate more confidence
# in one's beliefs about initial states, and entail that beliefs will
# change more slowly (e.g., the shape of the distribution encoded by 
# d{1} = [25 25]' will change much more slowly than the shape of the 
# distribution encoded by d{1} = [.5 0.5]' with each new observation).

# Setup a list for d
d = vector("list", 2*1)
dim(d) = matrix(c(2,1))
d

# For context beliefs, we can specify that the agent starts out believing 
# that both contexts are equally likely, but with somewhat low confidence 
# in these beliefs:

d[[1,1]] = t(matrix(c(.25,.25)))# c('left better','right better')

# For behavior beliefs, we can specify that the agent expects with 
# certainty that it will begin a trial in the 'start' state:

d[[2,1]] = t(matrix(c(1, 0, 0, 0))) # c('start','hint','choose-left','choose-right')


# State-outcome mappings and beliefs: A and a
#=========================================================================
  
#--------------------------------------------------------------------------
# Specify the probabilities of outcomes given each state in the 
# generative process (A)
# This includes one matrix per outcome modality
# Note: By default, these will also be the beliefs in the generative 
# model
#--------------------------------------------------------------------------
  
# First we specify the mapping from states to observed hints (outcome
# modality 1). Here, the rows correspond to observations, the columns
# correspond to the first state factor (context), and the third 
# dimension corresponds to behavior. Each column is a probability 
# distribution that must sum to 1.

# Setup a list:
A = vector("list", 3*4)
dim(A) = matrix(c(3,4))
A


# We start by specifying that both contexts generate the 'No Hint'
# observation across all behavior states:
  
## number of states in each state factor (2 and 4)
Ns = t(matrix(c(ncol(D[[1]]), ncol(D[[2]])))) 
Ns


### A[1,1:4]

for (i in 1:Ns[,2]){
  # [[i,,]] for all behavior states
  A[[1,i]] <- matrix(c(1,1,     # No Hint
                       0,0,     # Machine-Left Hint
                       0,0),    # Machine-Right Hint
                     ncol = 2, nrow = 3, byrow = TRUE)
}


### A[1,2]

# Then we specify that the 'Get Hint' behavior state generates a hint that
# either the left or right slot machine is better, depending on the context
# state. In this case, the hints are accurate with a probability of pHA. 

pHA = 1 # By default we set this to 1, but try changing its value to 
# see how it affects model behavior

A[[1,2]] = matrix(c(0,     0,    # No Hint
                    pHA, (1-pHA),   # Machine-Left Hint
                    (1-pHA), pHA),  # Machine-Right Hint
                  nrow = 3, ncol = 2, byrow = TRUE)

# Next we specify the mapping between states and wins/losses. The first 
# two behavior states ('Start' and 'Get Hint') do not generate either 
# win or loss observations in either context:

for (i in 1:2){
  
  A[[2,i]] = matrix(c( 1, 1,   # Null
                       0, 0,   # Loss
                       0, 0),  # Win
                    ncol = 2, nrow = 3, byrow = TRUE) 
}

# Choosing the left machine (behavior state 3) generates wins with
# probability pWin, which differs depending on the context state (columns):

pWin = .8  # By default we set this to .8, but try changing its value to 
# see how it affects model behavior

A[[2,3]] = matrix(c(  0,      0,        # Null        
                      (1-pWin),   pWin,     # Loss
                      pWin,    (1-pWin)), # Win
                  ncol = 2, nrow = 3, byrow=TRUE) 


# Choosing the right machine (behavior state 4) generates wins with
# probability pWin, with the reverse mapping to context states from 
# choosing the left machine:

A[[2,4]] = matrix(c(  0,      0,       # Null
                      pWin, (1-pWin),   # Loss
                      (1-pWin),  pWin),    # Win
                  ncol = 2, nrow = 3, byrow = TRUE)

# Finally, we specify an identity mapping between behavior states and
# observed behaviors, to ensure the agent knows that behaviors were carried
# out as planned. Here, each row corresponds to each behavior state.

for (i in 1:Ns[,2]){ 
  
  # x will be added to a 2*4 Matrix  
  x = matrix(c(1,1), ncol = 2, nrow = 1, byrow=TRUE)                                
  A3 = matrix(c(0,0, 0,0 ,0,0, 0,0), # c('start','hint','choose-left','choose-right')
              ncol = 2, nrow = 4, byrow=TRUE)  
  #y = A3[col = i]
  A3[i,] = 1
  A[[3,i]] <- A3
  
}

# See to understand the consequence of the loop
A[3,]
A 


#--------------------------------------------------------------------------
# Specify prior beliefs about state-outcome mappings in the generative model 
# (a)
# Note: This is optional, and will simulate learning state-outcome mappings 
# if specified.
#--------------------------------------------------------------------------
  
# Similar to learning priors over initial states, this simply
# requires specifying a matrix (a) with the same structure as the
# generative process (A), but with Dirichlet concentration parameters
# that can encode beliefs (and confidence in those beliefs) that need
# not match the generative process. Learning then corresponds to
# adding to the values of matrix entries, based on what outcomes were 
# observed when the agent believed it was in a particular state. For
# example, if the agent observed a win while believing it was in the 
# 'left better' context and the 'choose left machine' behavior state,
# the corresponding probability value would increase for that 
# location in the state outcome-mapping (i.e., a{2}(3,1,3) might 
# change from .8 to 1.8).

# One simple way to set up this matrix is by:

#  1. initially identifying it with the generative process 
#  2. multiplying the values by a large number to prevent learning all
#     aspects of the matrix (so the shape of the distribution changes 
#     very slowly)
#  3. adjusting the elements you want to differ from the generative 
#     process.

# Setup a list of lists:
a = vector("list", 3*4)
dim(a) = matrix(c(3,4))
a

# As another example, to simulate learning the hint accuracy one
# might specify:

for (i in 1:4){
  apre = as.matrix(A)
  a[[1,i]] = apre[[1,i]]*200
  a[[2,i]] = apre[[2,i]]*200
  a[[3,i]] = apre[[3,i]]*200
}

a[[1,2]] =  matrix(c(0,     0,    # No Hint
                   .25,   .25,    # Machine-Left Hint
                   .25,   .25),   # Machine-Right Hint
                       nrow = 3, ncol = 2, byrow = TRUE)
a
# Controlled transitions and transition beliefs : B{:,:,u} and b(:,:,u)
#==========================================================================
  
#--------------------------------------------------------------------------
# Next, we have to specify the probabilistic transitions between 
# hidden states under each action (sometimes called 'control states'). 
# Note: By default, these will also be the transitions beliefs 
# for the generative model
#--------------------------------------------------------------------------
  
# Columns are states at time t. Rows are states at t+1.

# Setup list:
B = vector("list", 2*4)
dim(B) = matrix(c(2,4))
B


# The agent cannot control the context state, so there is only 1 'action',
#  indicating that contexts remain stable within a trial:

B[[1,1]] = matrix(c(1, 0,   # 'Left Better' Context
                    0, 1),  # 'Right Better' Context
                  nrow = 2, ncol = 2, byrow = TRUE)

# The agent can control the behavior state, and we include 4 possible 
# actions:

# Move to the Start state from any other state
B[[2,1]] = matrix(c(1, 1, 1, 1,  # Start State
                    0, 0, 0, 0,  # Hint
                    0, 0, 0, 0,  # Choose Left Machine
                    0, 0, 0, 0), # Choose Right Machine
                  nrow = 4, ncol = 4, byrow = TRUE)

# Move to the Hint state from any other state
B[[2,2]] = matrix(c(0, 0, 0, 0,  # Start State
                    1, 1, 1, 1,  # Hint
                    0, 0, 0, 0,  # Choose Left Machine
                    0, 0, 0, 0), # Choose Right Machine
                  nrow = 4, ncol = 4, byrow = TRUE)

# Move to the Choose Left state from any other state
B[[2,3]] = matrix(c(0, 0, 0, 0,  # Start State
                    0, 0, 0, 0,  # Hint
                    1, 1, 1, 1,  # Choose Left Machine
                    0, 0, 0, 0), # Choose Right Machine
                  nrow = 4, ncol = 4, byrow = TRUE)

# Move to the Choose Right state from any other state
B[[2,4]] = matrix(c(0, 0, 0, 0,  # Start State
                    0, 0, 0, 0,  # Hint
                    0, 0, 0, 0,  # Choose Left Machine
                    1, 1, 1, 1), # Choose Right Machine        
                  nrow = 4, ncol = 4, byrow = TRUE)


#--------------------------------------------------------------------------
# Specify prior beliefs about state transitions in the generative model
# (b). This is a set of matrices with the same structure as B.
# Note: This is optional, and will simulate learning state transitions 
# if specified.
#--------------------------------------------------------------------------
  
# For this example, we will not simulate learning transition beliefs. 
# But, similar to learning d and a, this just involves accumulating
# Dirichlet concentration parameters. Here, transition beliefs are 
# updated after each trial when the agent believes it was in a given 
# state at time t and and another state at t+1.

# Preferred outcomes: C and c
#==========================================================================
  
#--------------------------------------------------------------------------
# Next, we have to specify the 'prior preferences', encoded here as log
# probabilities. 
#--------------------------------------------------------------------------
# One matrix per outcome modality. Each row is an observation, and each
# columns is a time point. Negative values indicate lower preference,
# positive values indicate a high preference. Stronger preferences promote
# risky choices and reduced information-seeking.

# We can start by setting a 0 preference for all outcomes:

# No: number of outcomes in each outcome modality
No = matrix(c(nrow(A[[1,1]]), nrow(A[[2,1]]), nrow(A[[3,1]]) ))
No # Unused here so far, as a different appoach than in the Matlabcode
# was used at this point.  
# Though this tells us the number of rows needed for C

# Setup list for C
C = vector("list", 3*1)
dim(C) = matrix(c(3,1))
C

# Hints
C[[1,1]]      = matrix(c(0, 0, 0,   # Not hint    
                         0, 0, 0,   # Machine-Left Hint
                         0, 0, 0),  # Machine-Left Hint
                     nrow = 3, ncol = 3, byrow = TRUE)

# Wins/Losses
C[[2,1]]      = matrix(c(0, 0, 0,   # Null      
                         0, 0, 0,   # Loss
                         0, 0, 0),  # Win
                     nrow = 3, ncol = 3, byrow = TRUE)

# Observed Behaviors

C[[3,1]]      = matrix(c(0, 0, 0, 0,   # Start State     
                         0, 0, 0, 0,   # Hint
                         0, 0, 0, 0,   # Choose Left Machine
                         0, 0, 0, 0),  # Choose Right Machine
                     nrow = 4, ncol = 4, byrow = TRUE)

# Then we can specify a 'loss aversion' magnitude (la) at time points 2 
# and 3, and a 'reward seeking' (or 'risk-seeking') magnitude (rs). 
# Here, rs is divided by 2 at the third time point to encode a smaller 
# win ($2 instead of $4) if taking the hint before choosing a slot 
# machine.

la = 1   # By default we set this to 1, but try changing its value to 
# see how it affects model behavior

rs = 4  # By default we set this to 4, but try changing its value to 
        # see how it affects model behavior

C[[2,1]] =  matrix(c(0,  0,   0,      # Null
                     0, -la, -la,     # Loss
                     0,  rs,  rs/2),  # win
                 nrow = 3, ncol = 3, byrow = TRUE)


#--------------------------------------------------------------------------
# One can also optionally choose to simulate preference learning by
# specifying a Dirichlet distribution over preferences (c). 
#--------------------------------------------------------------------------
  
# This will not be simulated here. However, this works by increasing 
# the reference magnitude for an outcome each time that outcome is 
# observed. The assumption here is that preferences naturally increase
# for entering situations that are more familiar.

# Allowable policies: U or V. 
#==========================================================================
  
#--------------------------------------------------------------------------
# Each policy is a sequence of actions over time that the agent can 
# consider. 
#--------------------------------------------------------------------------
  
# For our simulations, we will specify V, where rows correspond to 
# time points and should be length T-1 (here, 2 transitions, from 
# time point 1 to time point 2, and time point 2 to time point 3):
  
NumPolicies = 5 # Number of policies
NumFactors = 2 # Number of state factors

# V         = ones(Time-1,NumPolicies,NumFactors)
# V        # actually works, but used lists for R
V = vector("list", 1*2)
dim(V) = matrix(c(1,2))


# Deep policies v[[]]

V[[1]]         = matrix(c(1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1), # Context state is not controllable
                        nrow = 2, ncol = 5, byrow = TRUE)


V[[2]]         = matrix(c(1, 2, 2, 3, 4,
                          1, 3, 4, 1, 1), 
                        nrow = 2, ncol = 5, byrow = TRUE)

# For V[[2]], columns left to right indicate policies allowing: 
# 1. staying in the start state 
# 2. taking the hint then choosing the left machine
# 3. taking the hint then choosing the right machine
# 4. choosing the left machine right away (then returning to start 
#    state)
# 5. choosing the right machine right away (then returning to start 
#    state)


# Habits: E and e. 
#==========================================================================
  
#--------------------------------------------------------------------------
# Optional: a columns vector with one entry per policy, indicating 
# the prior probability of choosing that policy (i.e., independent 
# of other beliefs). 
#--------------------------------------------------------------------------
  
# We will not equip our agent with habits with any starting habits 
# (flat distribution over policies):
  
E = matrix(c(.1, .1, .1, .1, .1),
            nrow = 1, ncol = 5, byrow = TRUE)

# To incorporate habit learning, where policies become more likely 
# after each time they are chosen, we can also specify concentration
# parameters by specifying e:

e = matrix(c(1, 1, 1, 1, 1), 
             nrow = 1, ncol = 5, byrow = TRUE)


# Additional optional parameters. 
#==========================================================================
  
# Eta: learning rate (0-1) controlling the magnitude of concentration 
# parameter updates after each trial (if learning is enabled).

eta = 1 # Default (maximum) learning rate

# Omega: forgetting rate (0-1) controlling the magnitude of 
# reduction in concentration parameter values after each trial 
# (if learning is enabled).

omega = 1 # Default value indicating there is no forgetting 
          # (values < 1 indicate forgetting)

# Beta: Expected precision of expected free energy (G) over 
# policies (a positive value, with higher values indicating lower 
# expected precision). Lower values increase the influence of 
# habits (E) and otherwise make policy selection less deteriministic.

beta = 1 # By default this is set to 1, but try increasing its value 
         # to lower precision and see how it affects model behavior

# Alpha: An 'inverse temperature' or 'action precision' parameter that 
# controls how much randomness there is when selecting actions 
# (e.g., how  often the agent might choose not to take the hint, 
# even if the model assigned the highest probability to that action. 
# This is a positive number, where higher values indicate less 
# randomness. Here we set this to a fairly high value:

alpha = 32 # fairly low randomness in action selection


## Define POMDP Structure # Maybe not necessary in R
#==========================================================================
  
mdp.T = Time                 # Number of time steps
mdp.V = V                    # allowable (deep) policies

mdp.A = A                    # state-outcome mapping
mdp.B = B                    # transition probabilities
mdp.C = C                    # preferred states
mdp.D = D                    # priors over initial states
mdp.d = d                    # enable learning priors over initial states


if (Gen_model == 1){    # NECESSARY IN R???#################################################################################################################################################
mdp.E = E }               # prior over policies

if (Gen_model == 2){
mdp.a = a                # enable learning state-outcome mappings
mdp.e = e                # enable learning of prior over policies
}

mdp.E
mdp.eta = eta                 # learning rate
mdp.omega = omega             # forgetting rate
mdp.alpha = alpha             # action precision
mdp.beta = beta               # expected free energy precision


# respecify for use in inversion script (specific to this tutorial example)
mdp.NumPolicies = NumPolicies  # Number of policies
mdp.NumFactors = NumFactors    # Number of state factors

# MDP list/structure to work with
# MDP = list(mdp)

MDP = list(a,A,B,C,D,d,E,V,Time,eta,alpha,beta,omega,NumPolicies,NumFactors)     

# (Re)name items of list
names(MDP) = c("a","A","B","C","D","d","E","V","Time","eta","alpha","beta","omega","NumPolicies","NumFactors")
MDP$a

# Store initial paramater values of generative model for free energy 
# calculations after learning
#--------------------------------------------------------------------------

## 'complexity' of d vector concentration paramaters

if (isfield(MDP, "d")){
  x = numel(MDP$d)
  # Setup lists
  d_prior = vector("list", 1*x)  
  dim(d_prior) = matrix(c(1,x))
  
  d_complexity = vector("list", 1*x)  
  dim(d_complexity) = matrix(c(1,x))
  
  for (factor in 1:x){
    # store d vector values before learning
    y=MDP$d
    d_prior[1,factor] = y[factor]
  }
  for (factor in 1:x){
    # compute "complexity" - lower concentration paramaters have
    # smaller values creating a lower expected free energy thereby
    # encouraging 'novel' behaviour 
    z=as.vector(d_prior[[1,factor]])
    d_complexity[[factor]] = spm_wnorm(z)
  } 
  for (i in 1:x){
    dim(d_complexity[[1,i]]) = c(1,ncol(d_prior[[1,i]]))
  }
} 
# Check results
d_prior[1,2]
d_complexity[1,2] # Both correct!


## Similar for a:

if (isfield(MDP,"a")){
  # complexity of a maxtrix concentration parameters
  x=numel(MDP$a)
  
  #Setup list
  a_prior = vector("list", 1*x)
  dim(a_prior) = matrix(c(1,x))
  
  for (modality in 1:x){
    y = MDP$a
    a_prior[modality] = y[modality]
  }
  dim(a_prior) = matrix(c(nrow(MDP$a), ncol(MDP$a)))
  
  # Setup list for a_cpomplexity
  a_complexity = vector("list", nrow(a_prior)*ncol(a_prior))
  dim(a_complexity) = matrix(c(nrow(a_prior),ncol(a_prior)))
  
  for (modality in 1:ncol(a_prior)){
    for (rowdim in 1:nrow(a_prior)){
      z=as.matrix(a_prior[[rowdim,modality]])
      a_complexity[[rowdim,modality]] = spm_wnorm2(z)*(z>0)
    }
  }
}  

a_prior[2,4]         # Correct
a_complexity[1,2]    # Correct!


# Normalise matrices before model inversion/inference
#--------------------------------------------------------------------------

## normalize a matrix

if(isfield(MDP,"a")){
  a[] = lapply(MDP$a,col_norm)      
} else{
  a[] = lapply(MDP$A,col_norm)
}                                   # Works!

## Col_norm. of A 

A[] = lapply(MDP$A,col_norm)


## normalize b matrix

if(isfield(MDP,"b")){     
  b = lapply(MDP$b, \(x) if (!is.null(x)) col_norm(x) else NULL)
  }else{
  b = lapply(MDP$B, \(x) if (!is.null(x)) col_norm(x) else NULL)
}                                   # WORKS!
dim(b) = matrix(c(nrow(MDP$B), ncol(MDP$B))) # Dim B
b[2,]

# Normalize B matrix (in matlab above together with A and D)
B[] = lapply(MDP$B, \(x) if (!is.null(x)) col_norm(x) else NULL)
B

# normalize C and transform into log probability
for (ii in 1:numel(MDP$C)){              
  C[[ii]] <- MDP$C[[ii]] + 1/32 
}
for(i in 1:Time){
  C[[i,1]] <- t(log(softmaxCOL(t(MDP$C[[i,1]][,]))))
  C[[i,1]] = MDP$C[[i,1]][,-(Time+1)] # Crop matrix, so ncol is same 
}                                     # with all elements of MDP$C
 

## normalize d matrix
# Transpose first
d1 = vector("list", numel(d)*1)
dim(d1) = matrix(c(numel(d),1))
D1 = vector("list", numel(D)*1)
dim(D1) = matrix(c(numel(D),1))

for (i in 1:numel(d)){
d1[[i,1]]=t(d[[i,1]])
D1[[i,1]]=t(D[[i,1]])
}

if(isfield(MDP,"d")){
  d[] = lapply(d1,col_norm)
} else{
  d[] = lapply(D1,col_norm)
}                                   

## Col. norm of D
D[] = lapply(D1,col_norm)


e
# normalize E vector
if (isfield(MDP,"e")){
E = MDP$e
E = E./sum(E)
}else if (isfield(MDP,"E")){
E = MDP$E
E = E/sum(E)
}else{
E = col_norm(ones(NumPolicies,1))
E = E/sum(E)
}

# Initialize variables
#--------------------------------------------------------------------------
  
# numbers of transitions, policies and states
NumModalities = nrow(a)       # number of outcome factors
NumFactors = numel(d)          # number of hidden state factors
NumPolicies = ncol(V[[1,1]])        # number of allowable policies
NumPolicies
# Setup list:
NumStates = vector("list", NumFactors*1)
dim(NumStates) = matrix(c(NumFactors,1))
NumControllable_transitions = vector("list", NumFactors*1)
dim(NumControllable_transitions) = matrix(c(NumFactors,1))

for (factor in 1:NumFactors){
  # number of hidden states
  NumStates[[factor,1]] = ncol(b[[factor]])   
  # number of hidden controllable hidden states for each 
  # factor (number of B matrices)
  # [Not very elegant solution below. Problem: NULL's in b] 
  NumControllable_transitions[[factor]] = prod(dim(b[[factor]]))
  NumControllable_transitions[[factor]] = as.matrix(NumControllable_transitions[[factor]])/prod(dim(b[[1,1]]))
    if (factor == NumFactors){
      dim(NumControllable_transitions) = matrix(c(2,1))
    }
}

# initialize the approximate posterior over states conditioned on 
# policies for each factor as a flat distribution over states at 
# each time point

# Setup list for state_posterior:

state_posterior = vector("list", NumFactors*NumPolicies)
dim(state_posterior) = matrix(c(NumFactors,NumPolicies))

for (policy in 1:NumPolicies){
  for (factor in 1:NumFactors){
  NumStates[[factor,1]] = length(D[[factor,1]]) # number of states in each hidden state factor
  state_posterior[[factor,policy]] = matrix(1,NumStates[[factor]],Time,policy)/NumStates[[factor,1]]
  }  
} 
vv=state_posterior   # Still necessary? 
vvv=state_posterior  # vv and vvv to keep track of changes in v_depolarization...
state_posterior1=state_posterior # Still necessary?

# initialize the approximate posterior over policies as a flat 
# distribution over policies at each time point
policy_posteriors = matrix(1,NumPolicies,Time)/NumPolicies

# initialize posterior over actions
# Setup list:
chosen_action = vector("list", nrow(B)*(Time-1))
dim(chosen_action) = matrix(c(nrow(B), (Time-1)))
chosen_action = zeros(nrow(B),Time-1)
chosen_action

# if there is only one policy
for (factors in 1:NumFactors){ 
  if (NumControllable_transitions[[factors]] == 1){
  chosen_action[factors,] = ones(1)
  }
}
chosen_action

# Add chosen_action to MDP list [NOT very elegant, but works so far:]
MDP = list.append(MDP, chosen_action)
MDP
names(MDP)[numel(MDP)] <- "chosen_action"
MDP

# Setup list for Gamma[[ ]]
gamma = list()

# initialize expected free energy precision (beta)
posterior_beta = 1
gamma[[1]] = (1/posterior_beta) # expected free energy precision

# initialize expected free energy precision (beta)
posterior_beta = 1
gamma[[2]] = 1/posterior_beta   # expected free energy precision


# message passing variables
TimeConst = 4        # time constant for gradient descent
NumIterations  = 16  # number of message passing iterations

# Setup lists:
true_states = vector("list", NumFactors*Time)
dim(true_states) = matrix(c(NumFactors,Time))

outcomes = vector("list", NumModalities*Time)
dim(outcomes) = matrix(c(NumModalities, Time))

prob_state = vector("list", nrow(D[[2]])*ncol(D[[2]]))
dim(prob_state) = matrix(c(nrow(D[[2]]),ncol(D[[2]])))

O = vector("list", NumModalities*Time)
dim(O) = matrix(c(NumModalities, Time))

# Re-list a as a cell array in order to be able to permute.
a1 <- c(rep(list(array(0, c(3, 2, 4))), 2), list(array(0, c(4, 2, 4))))
for (i in 1:ncol(A)){
  for (ii in 1:nrow(A)){
    a1[[ii]][,,i] = a[[ii,i]] 
  }
}

# This is for the function G_epsitemic_value. Not very elegent 
# solution here...
a2=list()
a2
a2[[1]] = cbind(a[[1,1]][,1],a[[1,1]][,2],a[[1,2]][,1],a[[1,2]][,2],a[[1,3]][,1],a[[1,3]][,2],a[[1,4]][,1],a[[1,4]][,2])
a2[[2]] = cbind(a[[2,1]][,1],a[[2,1]][,2],a[[2,2]][,1],a[[2,2]][,2],a[[2,3]][,1],a[[2,3]][,2],a[[2,4]][,1],a[[2,4]][,2])
a2[[3]] = cbind(a[[3,1]][,1],a[[3,1]][,2],a[[3,2]][,1],a[[3,2]][,2],a[[3,3]][,1],a[[3,3]][,2],a[[3,4]][,1],a[[3,4]][,2])

# Change NULL to NA, as R is more comfortable with this
b[b == "NULL"] = NA
# Re-list a as a cell array in order to be able to permute.
# I have found this option way to late, as it seems closer
# to the cell array structure for which I used matrix lists.
b1 <- c(rep(list(array(0, c(2, 2, 1))), 1), list(array(0, c(4, 4, 4))))
b1[[1]][,,1] = b[[1,1]]
b1[[2]][,,1] = b[[2,1]]
b1[[2]][,,2] = b[[2,2]]
b1[[2]][,,3] = b[[2,3]]
b1[[2]][,,4] = b[[2,4]]
b1[[2]][,,1]

normalized_firing_rates = vector("list", NumFactors*Time*Time*NumPolicies*NumIterations)
dim(normalized_firing_rates) = matrix(c(NumFactors,Time,Time,NumPolicies,NumIterations))

prediction_error = vector("list", NumFactors*Time*Time*NumPolicies*NumIterations)
dim(prediction_error) = matrix(c(NumFactors,Time,Time,NumPolicies,NumIterations))

Ft = vector("list", Time*NumIterations*Time*NumFactors)
dim(Ft) = matrix(c(Time,NumIterations,Time,NumFactors))

Free = vector("list", NumPolicies*Time)
dim(Free) = matrix(c(NumPolicies,Time))

Expected_states = vector("list", 1*NumFactors)
dim(Expected_states) = matrix(c(1,NumFactors))

Gintermediate = vector("list", NumPolicies*1)
dim(Gintermediate) = matrix(c(NumPolicies,1))

policy_posteriors = matrix(1,NumPolicies,Time)/NumPolicies
policy_posteriors

policy_priors = vector("list", NumPolicies*Time)
dim(policy_priors) = matrix(c(NumPolicies,Time))

G = list()


### Lets go! Message passing and policy selection 
#--------------------------------------------------------------------------

# Execute the following to get the loop running in its current state
chosen_actionTEST = matrix(c(1,1,1,1), ncol = 2, nrow = 2, byrow = TRUE)

# Execute the short nested loop to reset state_posterior
state_posterior = vector("list", NumFactors*NumPolicies)
dim(state_posterior) = matrix(c(NumFactors,NumPolicies))
#Reset state posterior:
for (policy in 1:NumPolicies){
  for (factor in 1:NumFactors){
    NumStates[[factor,1]] = length(D[[factor,1]]) # number of states in each hidden state factor
    state_posterior[[factor,policy]] = matrix(1,NumStates[[factor]],Time,policy)/NumStates[[factor,1]]
  }  
} 


for (t in 1:Time){  # loop over time points
  
  # sample generative process
  #------------------------------------------------------------------------
  for (factor in 1:NumFactors){ # number of hidden state factors
    # Here we sample from the prior distribution over states to obtain the
    # state at each time point. At T = 1 we sample from the D vector, and at
    # time T > 1 we sample from the B matrix. To do this we make a vector 
    # containing the cumulative sum of the columns (which we know sum to one), 
    # generate a random number (0-1),and then use the find function to take 
    # the first number in the cumulative sum vector that is >= the random number. 
    # For example if our D vector is [.5 .5] 50% of the time the element of the 
    # vector corresponding to the state one will be >= to the random number. 
    
    # sample states 
    if (t == 1){
      prob_state <- D[[factor]] # sample initial state T = 1 # WORKS! 1 0 0 0 Transp.
    }
    else if (t > 1){
      prob_state <- as.matrix(B[[factor, true_states[[factor, t-1]]]][,chosen_actionTEST[[factor, t-1]]])
    } 
    true_states[[factor,t]] = which(cumsum(as.matrix(prob_state))>= runif(1))[1]
  } 
  
  # sample observations
  for (modality in 1:NumModalities){ # loop over number of outcome modalities
    outcomes[[modality,t]] = which(cumsum(a[[modality,true_states[[2,t]]]][,true_states[[1,t]]])>=runif(1))[1]
  } 
  
  # express observations as a structure containing a 1 x observations 
  # vector for each modality with a 1 in the position corresponding to
  # the observation recieved on that trial
  for (modality in 1:NumModalities){
    vec = zeros(1,nrow(a[[modality]]))
    index = outcomes[[modality,t]]
    vec[[1,index]] = 1
    O[[modality,t]] = vec       
  }
  
  # marginal message passing (minimize F and infer posterior over states)
  #----------------------------------------------------------------------
  
  for (policy in 1:NumPolicies){
    for (Ni in 1:NumIterations){ # number of iterations of message passing  
      for (factor in 1:NumFactors){
        lnAopre = as.matrix(dim(state_posterior[[factor]]))
        lnAo <- zeros(lnAopre[[1,1]],lnAopre[[2,1]]) # initialise matrix containing the log likelihood of observations
        for (tau in 1:3){ # loop over tau
          v_depolarization = nat_log(as.matrix(state_posterior[[factor, policy]][,tau])) # convert approximate posteriors into depolarisation variable v 
          vvv[[factor,policy]][,tau] = v_depolarization
          if (tau < t+1){ # Collect an observation from the generative process when tau <= t
            for (modal in 1:NumModalities){ # loop over observation modalities
              # this line uses the observation at each tau to index
              # into the A matrix to grab the likelihood of each hidden state
              #lnA = aperm(nat_log(a[[modal]][outcomes[[modal,tau]],],c(2, 3, 4, 5, 6, 1))) in Matlab                          
              lnA <- nat_log(a1[[modal]][outcomes[[modal, tau]],,])
              for (fj in 1:NumFactors){
                # dot product with state vector from other hidden state factors 
                # (this is what allows hidden states to interact in the likleihood mapping)    
                if (fj != factor){        
                  lnAs = md_dot((lnA),as.matrix(state_posterior[[fj]][,tau]),fj)
                  lnA = lnAs
                }
              }
              lnAo[,tau] = lnAo[,tau] + as.vector(lnA)
              }
            }

          # 'forwards' and 'backwards' messages at each tau
          if (tau == 1){ # first tau
            lnD = nat_log(d[[factor]]) # forward message
            lnBs = nat_log(t(B_norm3(as.matrix(b[[factor,V[[factor]][tau,policy]]])))%*%as.matrix(state_posterior[[factor,policy]][,tau+1])) # backward message
            }
          else if (tau == Time){ # last tau                    
            lnD  = nat_log(as.matrix(b[[factor,V[[factor]][tau-1,policy]]])%*%as.matrix(state_posterior[[factor,policy]][,tau-1])) # forward message 
            lnBspre = as.matrix(dim(d[[factor]])) 
            lnBs = zeros(lnBspre[1],lnBspre[2]) # backward message
          }
          else { # 1 > tau > Time
            lnD  = nat_log(b[[factor,V[[factor]][tau-1,policy]]]%*%as.matrix(state_posterior[[factor,policy]][,tau-1])) # forward message
            lnBs = nat_log(t(B_norm3(as.matrix(b[[factor,V[[factor]][tau,policy]]])))%*%as.vector(state_posterior[[factor,policy]][,tau+1])) # backward message
          }
          
          # variational free energy at each time point
          Ft[[tau,Ni,t,factor]] = t(state_posterior[[factor,policy]][,tau])%*%(.5*as.matrix(lnD) + .5*as.matrix(lnBs) + as.matrix(lnAo[,tau]) - nat_log(state_posterior[[factor,policy]][,tau]))
          
          # here we both combine the messages and perform a gradient
          # descent on the posterior. 
          v_depolarization = v_depolarization + (.5*(lnD) + .5*(lnBs) + as.matrix(lnAo[,tau]) - v_depolarization)/TimeConst
          # update posterior by running v through a softmax 
          state_posterior[[factor,policy]][,tau] = softmax(v_depolarization) # t(exp(v_depolarization)/sum(exp(v_depolarization))) # 
          # store v (non-normalized log posterior or 'membrane potential') 
          # from each epoch of gradient descent for each tau
          prediction_error[[factor,tau,t,policy,Ni]] = t(v_depolarization)
          # store state_posterior (normalised firing rate) from each epoch of
          # gradient descent for each tau
          normalized_firing_rates[[factor,t,tau,policy,Ni]] = t(state_posterior[[factor,policy]][,tau])    
          
        }
        
      }
    }
    
    # variational free energy for each policy (F)
    Ftarr <- array(as.numeric(unlist(Ft)), dim = c(3, 16, 3, 2))
    Fintermediate = rowSums(Ftarr,dims = 3) # sum over hidden state factors (Fintermediate is an intermediate F value)
    Fintermediate = colSums(Fintermediate) # sum over tau and squeeze into 16x3 matrix 
                                                   # in R squeeze() not necessary (at least for standard example)
    # store variational free energy at last iteration of message passing
    Free[[policy,t]] = Fintermediate[NumIterations,t]
  }
 
  # expected free energy (G) under each policy
  #----------------------------------------------------------------------
    
  # initialize intermediate expected free energy variable (Gintermediate) for each policy
  Gintermediate = zeros(NumPolicies,1)  
  # policy horizon for 'counterfactual rollout' for deep policies (described below)
  horizon = Time
  
    # loop over policies
    for (policy1 in 1:NumPolicies){
  
      # Bayesian surprise about 'd'
      if (isfield(MDP, "d")){
        for (factor1 in 1:NumFactors){
        Gintermediate[[policy1]] = Gintermediate[[policy1]] - d_complexity[[factor1]]%*%state_posterior[[factor1,policy1]][,1]
        }       
      }
    # This calculates the expected free energy from time t to the
    # policy horizon which, for deep policies, is the end of the trial T.
    # We can think about this in terms of a 'counterfactual rollout'
    # that asks, "what policy will best resolve uncertainty about the 
    # mapping between hidden states and observations (maximize
    # epistemic value) and bring about preferred outcomes"?
        
      for (timestep in t:horizon){
         # grab expected states for each policy and time
         for (factor2 in 1:NumFactors){
           Expected_states[[factor2]] = as.matrix(state_posterior[[factor2, policy1]][,timestep])
         }
        
         # calculate epistemic value term (Bayesian Surprise) and add to
         # expected free energy
         Gintermediate[[policy1]] = Gintermediate[[policy1]] + G_epistemic_value(a2,Expected_states) 
         
         }
      }
  # store expected free energy for each time point and clear intermediate
  # variable
  G[[t]] = Gintermediate 
  
}  

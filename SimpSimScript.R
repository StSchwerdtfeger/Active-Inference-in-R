#####################################
#-----------------------------------#
# Simplified Simulation Script in R #
#-----------------------------------#
#####################################

# Supplementary Code for: A Step-by-Step Tutorial on Active Inference 
# Modelling and its Application to Empirical Data

# By: Ryan Smith, Karl J. Friston, Christopher J. Whyte

# CONVERSION by Steffen Schwerdtfeger

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

#############
# Libraries #
#############

library(pracma)


#############
# Functions # 
#############

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
  log(x+exp(-16))
}

# Replication of the gradient function Matlab-Style
# for calcualting the 1D numerical gradient.  
# Alternative use function with same name via library("pracma")
gradient = function (x) {
  step = length(x)
  upLim  = step-1 # needed for differences on interior points
  
  # List for result
  result = vector("list",1*length(x))
  dim(result) = c(1,length(x))
  for (i in 1:step){
    if (i == 1){
      result[1] = (x[1+1]-x[1])/1  # Take forward differences on left... 
    }
    if (i == step){                # ....... and right edges.
      result[i] = (x[step]-x[step-1])/1
    }
    for (i in 2:upLim) { # needed for differences on interior points
      result[i] = (x[i+1]-x[i-1])/2
    }
  }
  # Print
  print(result)
} # END of function

# Normalizing vectors. So far this seems to be enough in R
# Compare with Matlab code for col_norm(B). 
# So far I used the same function for B_norm() in the Matlab script. 
col_norm = function(x){
  lapply(x,lapply, function(x)(t(t(x)/colSums(x))))
}
testoar = c(1,1)

# B_norm 
B_norm = function(x){
  bb = x
  z  = colSums(t(bb))     # create normalizing constant from sum of columns
  bb = t(bb/z)            # divide columns by constant
  # z is a vector [1] 1 1; using "/" z
  bb[is.nan(bb)] = 0   # replace NaN with zero (probably not necessary in R)
  b = bb
  return(b)
}

# used to only keep those values in a list that are greater than zero 
# refers to line 104 in the Matlab script:
GreaterZero = function(x){
  TrueGreatZero = function(x){x>0}
  checkLogic = lapply(x,TrueGreatZero) # Works both!
  trueISone = function(x){ x = x*1}
  checkNum = lapply(checkLogic,trueISone)
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

# Select last element of a vector, or 1D matrix:
last <- function(x) { return( x[[length(x)]] ) }


#################
# SPM Functions #
#################

spm_wnorm = function(x) { # Start of Function:
  
  # This is a replication of the bsxfun function to subtract the 
  # inverse of each column entry from the inverse of the sum of the 
  # columns and then divide by 2.
  if(ncol(as.matrix(x[[1]])) == 1){
    
    X = x
    
    X1   = lapply(X,"+", exp(-16))
    X2   = lapply(X1, colSums)
    
    X3 = lapply(X2, function(x)(1/x))
    X4 = lapply(X1, function(x)(1/x))
    
    X5 = lapply(X4,"-",as.numeric(X3[[1]]))
    x = lapply(X5, "/",2)
    
  } # End if
  else{
    x   = lapply(x,"+", exp(-16))
    for(i in 1:length(x)){
      A = 1/colSums(x[[i]]) # Correct up to here! Matlab:  1./sum(A,1)
      B = 1/x[[i]]
      x2 =  (A-B)/2 #lapply(X2,"-",X3)
      x[[i]] = x2
    } # End of loop
    x=x # Necessary to assign the values obtained from the loop!
  } # End else if
} # End of Function


#######################
# Simulation Settings #
#######################

# To simulate the task when prior beliefs (d) are separated from the 
# generative process, set the 'Gen_model' variable directly
# below to 1. To do so for priors (d), likelihoods (a), and habits (e), 
# set the 'Gen_model' variable to 2:

Gen_model = 2 # as in the main tutorial code, many parameters 
# can be adjusted in the model setup, within the 
# explore_exploit_model function [[starting on line 810]]
# (Matlab code). This includes, among others (similar to 
# in the main tutorial script):

# prior beliefs about context (d): alter line 876
# beliefs about hint accuracy in the likelihood (a): alter lines 996-998
# to adjust habits (e), alter line 1155

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

############################
# Specify Generative Model #
############################

# THE following refers to the Matlab script, as it has a slightyl different 
# structure, using a function toi switch between Gen_model 1 and 2.
# This is not needed here.
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

# Setup vector list for D  
D = c(list())

# Predim list, otherwise assigning via e.e D[[1]][[1]] is not possible, but
# necessary for the loop later on:
D[[1]] = c(rep(list(zeros(2,1)),1)) # can be random.
D[[2]] = c(rep(list(zeros(2,1)),1)) 

# For the 'context' state factor, we can specify that the 'left better' 
# context (i.e., where the left slot machine is more likely to win)
# is the true context:

# matrix(c('left better','right better'))
D[[1]][[1]] = matrix(c(1, 0))  

# For the 'behavior' state factor, we can specify that the agent 
# always begins a trial in the 'start' state (i.e., before choosing
# to either pick a slot machine or first ask for a hint:

# matrix(c('start','hint','choose-left','choose-right'))
D[[2]][[1]] = matrix(c(1, 0, 0, 0))

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

# Setup vector list for d
d = c(list())
d[[1]] = c(rep(list(zeros(2,1)),1)) # can be random.
d[[2]] = c(rep(list(zeros(2,1)),1)) # can be random.


# For context beliefs, we can specify that the agent starts out believing 
# that both contexts are equally likely, but with somewhat low confidence 
# in these beliefs:

# matrix(c('left better','right better'))
d[[1]][[1]] = matrix(c(.25,.25)) 

# For behavior beliefs, we can specify that the agent expects with 
# certainty that it will begin a trial in the 'start' state:

# matrix(c('start','hint','choose-left','choose-right'))
d[[2]][[1]] = matrix(c(1, 0, 0, 0))


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

# Setup a vector list:
A = c(list())

# We start by specifying that both contexts generate the 'No Hint'
# observation across all behavior states:

## number of states in each state factor (2 and 4)
Ns = matrix(c(length(D[[1]][[1]]), length(D[[2]][[1]])))      
Ns

### A_bhav / context (context  = left or right)
A_bhavCont = matrix(c(1,1,     # No Hint
                      0,0,     # Machine-Left Hint
                      0,0),    # Machine-Right Hint
                    nrow = 3, byrow = TRUE)

# Assign to vector list:
# Alternative to the loop in the Matlab script (line 905, 12.06.2022)
A[[1]] = c(rep(list(A_bhavCont),Ns[2]))


# TESTGROUND:
#A[[1]][[4]][1,1] # works!!! Add in double brakets to fixate an element
# and then 
#A[[1]][[4]]*c(1,2) # Vector and matrix multiplications work fine!!


# Then we specify that the 'Get Hint' behavior state generates a hint that
# either the left or right slot machine is better, depending on the context
# state. In this case, the hints are accurate with a probability of pHA. 

pHA = 1 # By default we set this to 1, but try changing its value to 
# see how it affects model behavior

A_hintAcc = matrix(c( 0   ,    0   ,   # No Hint
                      pHA  , (1-pHA),   # Machine-Left Hint
                      (1-pHA),   pHA) ,   # Machine-Right Hint
                   nrow = 3, byrow = TRUE)
# Assign to vector list:
A[[1]][[2]] = A_hintAcc 

# Next we specify the mapping between states and wins/losses. The first 
# two behavior states ('Start' and 'Get Hint') do not generate either 
# win or loss observations in either context:

A_contWL = matrix(c( 1, 1,   # Null
                     0, 0,   # Loss
                     0, 0),  # Win
                  nrow = 3, byrow = TRUE) 
# Assign to list:
# Alt to loop in Matlab script (line 927, 12.06.2022) 
A[[2]] = c(rep(list(A_contWL),2))

# Choosing the left machine (behavior state 3) generates wins with
# probability pWin, which differs depending on the context state (columns):

pWin = .8  # By default we set this to .8, but try changing its value to 
# see how it affects model behavior

# Does not need an extra name to be assigned: 
A[[2]][[3]] = matrix(c(      0   ,     0    ,  # Null        
                             (1-pWin),   pWin   ,  # Loss
                             pWin  , (1-pWin)),  # Win
                     nrow = 3, byrow = TRUE) 

# Choosing the right machine (behavior state 4) generates wins with
# probability pWin, with the reverse mapping to context states from 
# choosing the left machine:

# Assign to list:
A[[2]][[4]] = matrix(c(    0   ,      0   ,   # Null
                           pWin  ,  (1-pWin),   # Loss
                           (1-pWin),    pWin) ,   # Win
                     ncol = 2, nrow = 3, byrow = TRUE)

# Finally, we specify an identity mapping between behavior states and
# observed behaviors, to ensure the agent knows that behaviors were carried
# out as planned. Here, each row corresponds to each behavior state.

A_Id = matrix(c(0,0,  # Start
                0,0,  # Hint
                0,0,  # Choose-left
                0,0), # Choose-right
              nrow = 4, byrow = TRUE)  

# Assign to vector list
A[[3]] = c(rep(list(A_Id),nrow(A_Id)))

# Adds a c(1,1) vector to the respective line (Matlab line 956)
for (i in 1:Ns[2]){
  A[[3]][[i]][i,] = c(1,1)
}

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

# Setup vector list (here a modification of the list A)
a = A

# As another example, to simulate learning the hint accuracy one
# might specify:
a[[1]] = lapply(a[[1]],"*", 200) 
a[[2]] = lapply(a[[2]],"*", 200) 
a[[3]] = lapply(a[[3]],"*", 200) 


a[[1]][[2]] =  matrix(c(0,     0,    # No Hint
                        .25,   .25,    # Machine-Left Hint
                        .25,   .25),   # Machine-Right Hint
                      nrow = 3, byrow = TRUE)


# Controlled transitions and transition beliefs : B{:,:,u} and b(:,:,u)
#==========================================================================

#--------------------------------------------------------------------------
# Next, we have to specify the probabilistic transitions between 
# hidden states under each action (sometimes called 'control states'). 
# Note: By default, these will also be the transitions beliefs 
# for the generative model
#--------------------------------------------------------------------------

# Columns are states at time t. Rows are states at t+1.

# Setup vector list:
B = c(list())

# The agent cannot control the context state, so there is only 1 'action',
# indicating that contexts remain stable within a trial:

# Predim list, otherwise assigning via e.e B[[1]][[1]] is not possible, but
# necessary for the loop later on:
B[[1]] = c(rep(list(zeros(2,2)),4))

B[[1]][[1]] = matrix(c(1, 0,   # 'Left Better' Context
                       0, 1),  # 'Right Better' Context
                     nrow = 2, byrow = TRUE)

# For the loop, we have consider B[[1]][] having the same dim as B[[2]][].
# The reason is that B[[1]] refers to the context state, which is an identity
# matrix. In R we have to represent the identity, such that the the length
# of B[[1]] == B[[2]]. Otherwise the loop spits out funny values. Compare
# page 26 in the ATUT on B[[1]] being identity matrix by intention. 

for (i in 1:4){
  B[[1]][[i]] = B[[1]][[1]]
}

# The agent can control the behavior state, and we include 4 possible 
# actions:

# Pre dim again:
B[[2]] = c(rep(list(zeros(4,4)),4))

# Adds a c(1,1,1,1) vector to respective line, similar as before
for (i in 1:Ns[2]){
  B[[2]][[i]][i,] = ones(1,4)
}

# THE BELOW IS JUST FOR AN OVERVIEW, but otherwise not necessary so far
# in this script.

# Move to the Start state from any other state
B[[2]][[1]] = matrix(c(1, 1, 1, 1,  # Start State
                       0, 0, 0, 0,  # Hint
                       0, 0, 0, 0,  # Choose Left Machine
                       0, 0, 0, 0), # Choose Right Machine
                     nrow = 4, byrow = TRUE)

# Move to the Hint state from any other state
B[[2]][[2]] = matrix(c(0, 0, 0, 0,  # Start State
                       1, 1, 1, 1,  # Hint
                       0, 0, 0, 0,  # Choose Left Machine
                       0, 0, 0, 0), # Choose Right Machine
                     nrow = 4, byrow = TRUE)

# Move to the Choose Left state from any other state
B[[2]][[3]] = matrix(c(0, 0, 0, 0,  # Start State
                       0, 0, 0, 0,  # Hint
                       1, 1, 1, 1,  # Choose Left Machine
                       0, 0, 0, 0), # Choose Right Machine
                     nrow = 4, byrow = TRUE)

# Move to the Choose Right state from any other state
B[[2]][[4]] = matrix(c(0, 0, 0, 0,  # Start State
                       0, 0, 0, 0,  # Hint
                       0, 0, 0, 0,  # Choose Left Machine
                       1, 1, 1, 1), # Choose Right Machine        
                     nrow = 4, byrow = TRUE)


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
No = matrix(c(nrow(A[[1]][[1]]), nrow(A[[2]][[1]]), nrow(A[[3]][[1]]) ))
No 

# Setup vector list for C
C = c(list())

C[[1]] = c(rep(list(zeros(No[1],Time)),1))
C[[2]] = c(rep(list(zeros(No[2],Time)),1))
C[[3]] = c(rep(list(zeros(No[3],Time)),1))


# Alternative for an overview below, but otherwise not necessary:

# Hints
C[[1]][[1]]      = matrix(c(0, 0, 0,   # Not hint    
                            0, 0, 0,   # Machine-Left Hint
                            0, 0, 0),  # Machine-Left Hint
                          nrow = 3, byrow = TRUE)

# Wins/Losses
C[[2]][[1]]      = matrix(c(0, 0, 0,   # Null      
                            0, 0, 0,   # Loss
                            0, 0, 0),  # Win
                          nrow = 3, byrow = TRUE)

# Observed Behaviors

C[[3]][[1]]      = matrix(c(0, 0, 0,   # Start State     
                            0, 0, 0,   # Hint
                            0, 0, 0,   # Choose Left Machine
                            0, 0, 0),  # Choose Right Machine
                          nrow = 4, byrow = TRUE)

# Then we can specify a 'loss aversion' magnitude (la) at time points 2 
# and 3, and a 'reward seeking' (or 'risk-seeking') magnitude (rs). 
# Here, rs is divided by 2 at the third time point to encode a smaller 
# win ($2 instead of $4) if taking the hint before choosing a slot 
# machine.

la = 1   # By default we set this to 1, but try changing its value to 
# see how it affects model behavior

rs = 4  # By default we set this to 4, but try changing its value to 
# see how it affects model behavior

C[[2]][[1]] =  matrix(c(0,  0,   0,      # Null
                        0, -la, -la,     # Loss
                        0,  rs,  rs/2),  # win
                      nrow = 3, byrow = TRUE)

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

V = c(list())

# Pre dim
V[[1]] = c(rep(list(zeros(NumFactors,NumPolicies)),(Time-1)))

# Deep policies v[[]]

V[[1]][[1]]         = matrix(c(1, 1, 1, 1, 1,
                               1, 1, 1, 1, 1), # Context state is not controllable
                             nrow = 2, byrow = TRUE)


V[[1]][[2]]         = matrix(c(1, 2, 2, 3, 4,
                               1, 3, 4, 1, 1), 
                             nrow = 2, byrow = TRUE)

# For V[[1]][[2]], columns left to right indicate policies allowing: 
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

E = matrix(c(1, 1, 1, 1, 1),
           nrow = 5, ncol = 1, byrow = TRUE)

# To incorporate habit learning, where policies become more likely 
# after each time they are chosen, we can also specify concentration
# parameters by specifying e:

e = matrix(c(1, 1, 1, 1, 1), 
           nrow = 5, ncol = 1, byrow = TRUE)


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


## Define POMDP Structure
#==========================================================================

mdp.T = Time                 # Number of time steps
mdp.V = V                    # allowable (deep) policies

mdp.A = A                    # state-outcome mapping
mdp.B = B                    # transition probabilities
mdp.C = C                    # preferred states
mdp.D = D                    # priors over initial states
mdp.d = d                    # enable learning priors over initial states


if (Gen_model == 1){
  mdp.E = E }              # prior over policies

if (Gen_model == 2){
  mdp.a = a                # enable learning state-outcome mappings
  mdp.e = e                # enable learning of prior over policies
}

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

# Col_norm of A
A = col_norm(A)

# Col_norm of B
B = col_norm(B)

# Col_norm of D
D = col_norm(D)


# Store initial paramater values of generative model for free energy 
# calculations after learning
#--------------------------------------------------------------------------

## 'complexity' of d vector concentration paramaters

# Set vector list using d: as in the Matlab script, just 
# outside any loopin'
d_prior = MDP$d

# MDP$d set just for dimensional reasons
d_complexity = MDP$d  

for (factor in 1:length(d)){
  # compute "complexity" - lower concentration paramaters have
  # smaller values creating a lower expected free energy thereby
  # encouraging 'novel' behaviour 
  d_complexity[[factor]] = spm_wnorm(d_prior[[factor]])
} 

d_complexity 

# complexity of a maxtrix concentration parameters
# similar to d:
a_prior = MDP$a
a_complexity = c(list())

for (modality in 1:length(a)){
  a_complexity[[modality]] = spm_wnorm(a_prior[[modality]]) # *GreaterZero(a_prior[[modality]]) 
  for (i in length(a_complexity[[modality]])){
    a_complexity[[modality]][[i]]=a_complexity[[modality]][[i]]*GreaterZero(a_prior[[modality]])[[i]]
  }
} # Checked on all GreaterZero results and seems to be fine

# Normalise matrices before model inversion/inference
#--------------------------------------------------------------------------

## normalize a matrix

if(Gen_model==2){
  a = col_norm(a)      # For check: MDP$a = a
} 
if(Gen_model==1){
  a = col_norm(A)
}                                   

# No field for b defined in the SimpSimScript, so no coditions:
b = col_norm(B)      # For check: MDP$a = a

for(modality in 1:length(C)){
  C[[modality]] = lapply(C[[modality]],"+",1/32) # So far definitly right!
  for (t in 1:Time){
    C[[modality]][[1]][,t] = nat_log(softmaxCOL(C[[modality]][[1]][,t]))
  }
}

## normalize d matrix

if(Gen_model==2){
  d = col_norm(d)      # For check: MDP$a = a
} 
if(Gen_model==1){
  d = col_norm(D)
}                                   


# normalize E vector
if (Gen_model == 1){
  E = MDP$e
  E = E/sum(E) 
}
if (Gen_model==2){
  E = MDP$E
  E = E/sum(E)
}
#else{  # NOT NEEDED, but maybe for the actual Tut script.
#  E = col_norm(ones(NumPolicies,1))
#  E = E/sum(E)
#}

# Initialize variables
#--------------------------------------------------------------------------

# numbers of transitions, policies and states
NumModalities = length(a)         # number of outcome factors
NumFactors = length(d)            # number of hidden state factors
NumPolicies = ncol(V[[1]][[1]])   # number of allowable policies

# Setup matrix
NumStates = matrix(0,nrow=NumFactors)
NumControllable_transitions = matrix(0,nrow=NumFactors)

for(factor in 1:NumFactors){
  NumStates[[factor]] = nrow(b[[factor]][[1]])
  NumControllable_transitions[[factor]] = length(b[[factor]])
}

# Setup vector list for state_posterior:
state_posterior = c(list())
for(factor in 1: NumFactors){
  state_posterior[[factor]] = c(rep(list(),1))
}

# initialize the approximate posterior over states conditioned on 
# policies for each factor as a flat distribution over states at 
# each time point

for(policy in 1:NumPolicies){
  for(factor in 1: NumFactors){
    state_posterior[[factor]][[policy]] = matrix(1, 
                                                 nrow = NumStates[[factor]],
                                                 ncol = Time)
    state_posterior[[factor]][[policy]] = state_posterior[[factor]][[policy]]/NumStates[[factor]]
  }
}

# initialize the approximate posterior over policies as a flat 
# distribution over policies at each time point
policy_posteriors = matrix(1,NumPolicies,Time)/NumPolicies

# initialize posterior over actions:
chosen_action = matrix(0, nrow = length(B), ncol = Time-1)

# if there is only one policy
for(factor in 1:NumFactors){
  if(NumControllable_transitions[[factor]]==1){
    chosen_action[factor,] = matrix(1,nrow =1,ncol=Time-1)
  }
}

# Add to MDP list:   
# [[LOOKING FOR DIFFERENT WAY]] => Rlist package list.append() can be used
# but adds another package needed... 
MDP = list(a,A,B,C,D,d,E,V,Time,eta,alpha,beta,omega,NumPolicies,NumFactors,chosen_action)

# (Re)name items of list
names(MDP) = c("a","A","B","C","D","d","E","V","Time","eta","alpha","beta","omega","NumPolicies","NumFactors", "chosen_action")

# Set list for gamma:
gamma = list()

# initialize expected free energy precision (beta)
posterior_beta = 1
gamma[[1]] = (1/posterior_beta) # expected free energy precision

# message passing variables
TimeConst = 4        # time constant for gradient descent
NumIterations  = 16  # number of message passing iterations

### Lets go! Message passing and policy selection 
#--------------------------------------------------------------------------

# Pre-set of lists and matrices

# Set simple list for prob_state:    
prob_state = c(list())
# Set matrix for true_states:
true_states = matrix(1, nrow = NumFactors, ncol = Time)
# Set matrix fro outcomes:
outcomes = matrix(1, nrow = NumModalities, ncol = Time)
# Set up vector list, dimmed as matrix for O:
O = vector("list", NumModalities*Time)
dim(O) = matrix(c(NumModalities, Time))

# Re-list "a" as a cell array in order to be able to perform a
# work aroung of the permute() function in Matlab.
# Still looking for a consistent approach, when enlisting..
# Orginial matlab line: 209 i  question (within the loop) 
# lnA = permute(nat_log(a{modal}(outcomes(modal,tau),:,:,:,:,:)),[2 3 4 5 6 1])
aa <- c(rep(list(array(0, c(3, 2, 4))), 2), list(array(0, c(4, 2, 4))))
for (i in 1:length(A)){
  for (ii in 1:length(A[[1]])){
    aa[[i]][,,ii] = a[[i]][[ii]] 
  }
}

# Setup vector list for Ft
# Ft[[1]][tau,Ni,t,factor]
Ft = list()
Ft <- c(rep(list(array(0, c(Time,NumIterations, Time,NumFactors)))))

# List for VFE:
VFE <- c(rep(list(array(0, c(NumPolicies,Time)))))

# List for Fintermediate:
#Fintermediate <- c(rep(list(array(0, c(NumIterations, Time)))))


# List for prediction error and normalized_firign_rates
# prediction_error[[factor]][[1]][Ni,nrow(state_posterior[[factor]][[1]]),tau,t,policiy]
prediction_error = list()
for(factor in 1:NumFactors){
  prediction_error[[factor]] <- c(rep(list(array(0, c(NumIterations,nrow(state_posterior[[factor]][[1]]), Time, Time,NumPolicies)))))
}
normalized_firing_rates=prediction_error



vv= state_posterior
#state_posterior=vv
### Adjustments for testing:
###
chosen_action = matrix(c(1,1,1,1), ncol = 2, nrow = 2, byrow = TRUE)

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
      prob_state = D[[factor]] # sample initial state T = 1 # WORKS! 1 0 0 0 Transp.
    } #
    else if (t > 1){
      prob_state <- B[[factor]][[true_states[[factor, t-1]]]][,chosen_action[[factor, t-1]]]
    } #
    true_states[[factor,t]] = which(cumsum(prob_state[[1]])>= runif(1))[1]
  } # End loop factor
  
  # sample observations
  for (modality in 1:NumModalities){ # loop over number of outcome modalities
    outcomes[[modality,t]] = which(cumsum(a[[modality]][[true_states[[2,t]]]][,true_states[[1,t]]])>=runif(1))[1]
  }# End loop modality
  
  # express observations as a structure containing a 1 x observations 
  # vector for each modality with a 1 in the position corresponding to
  # the observation recieved on that trial
  for (modality in 1:NumModalities){
    vec = matrix(0,ncol=nrow(a[[modality]][[1]]))
    index = outcomes[[modality,t]]
    vec[[1,index]] = 1
    O[[modality,t]] = vec 
    
  } # End loop modality
  
  # marginal message passing (minimize F and infer posterior over states)
  #----------------------------------------------------------------------
  
  for (policy in 1:NumPolicies){
    for (Ni in 1:NumIterations){ # number of iterations of message passing  
      for (factor in 1:NumFactors){
        # initialise matrix containing the log likelihood of observations:
        lnAo = rep(list(matrix(0, nrow = nrow(state_posterior[[factor]][[1]]), 
                               ncol = ncol(state_posterior[[factor]][[1]]))), 
                   length(state_posterior[[factor]]))
        for (tau in 1:Time){ # loop over tau
          # convert approximate posteriors into depolarisation variable v
          v_depolarization = nat_log(state_posterior[[factor]][[policy]][,tau])
          if (tau < t+1){ # Collect an observation from the generative process when tau <= t
            for (modal in 1:NumModalities){ # loop over observation modalities
              # this line uses the observation at each tau to index
              # into the A matrix to grab the likelihood of each hidden state
              # NOTE: makes use of array list "aa" as alternative to "a",
              # to work around the permute() function in Matlab:
              lnA = nat_log(aa[[modal]][outcomes[[modal, tau]],,])
              for (fj in 1:NumFactors){
                # dot product with state vector from other hidden state factors 
                # (this is what allows hidden states to interact in the likleihood mapping)    
                if (fj != factor){        
                  lnAs = md_dot((lnA),state_posterior[[fj]][[1]][,tau],fj)
                  lnA = lnAs # Matlab uses "clear" here as well - needed? 
                } # End if fj != factor
              } # End loop fj
              lnAo[[1]][,tau] = mapply("+",lnAo[[1]][,tau], lnA) # very different to Matlab, but results in the same
              # Matlab: size(lnAo(:,:)) % dim is 4 15 and it is just all collumns next to each other!
              # size of lnAo itself though is 4 3 5; considering tau=3, the change only concerns the cols of lnAo[[1]]
            } # End loop modal 
          } # End if tau < t+1
          # 'forwards' and 'backwards' messages at each tau
          if (tau == 1){ # first tau
            # forward message
            lnD = nat_log(d[[factor]][[1]]) 
            # backward message
            lnBs = nat_log(t(B_norm(b[[factor]][[V[[1]][[tau]][factor,policy]]]))%*%as.matrix(state_posterior[[factor]][[policy]][,tau+1]))
          } # End if tau == 1
          else if (tau == Time){ # last tau  
            # foward message
            lnD  = nat_log(b[[factor]][[V[[1]][[tau-1]][factor,policy]]]%*%as.matrix(state_posterior[[factor]][[policy]][,tau-1]))  
            # Backward message:
            lnBs = matrix(0, nrow = nrow(d[[factor]][[1]]), ncol = ncol(d[[factor]][[1]]), byrow = TRUE)
          }
          else { # 1 < tau < Time
            # foward message
            lnD  = nat_log(b[[factor]][[V[[1]][[tau-1]][factor,policy]]]%*%as.matrix(state_posterior[[factor]][[policy]][,tau-1]))  
            lnBs = nat_log(t(B_norm(b[[factor]][[V[[1]][[tau]][factor, policy]]]))%*%as.matrix(state_posterior[[factor]][[policy]][,tau+1]))
          }
          
        } # End loop tau
        # here we both combine the messages and perform a gradient
        # descent on the posterior.
        v_depolarization = v_depolarization + ((.5*(lnD) + .5*(lnBs) + as.matrix(lnAo[[1]][,tau])) - v_depolarization)/TimeConst
        # variational free energy at each time point
        Ft[[1]][tau,Ni,t,factor] = t(state_posterior[[factor]][[policy]][,tau])%*%as.matrix(as.matrix(.5*(lnD)) + as.matrix(.5*(lnBs)) + as.matrix(lnAo[[1]][,tau])-as.matrix(nat_log(state_posterior[[factor]][[policy]][,tau])))
        # update state_posterior running v through a softmax:
        state_posterior[[factor]][[policy]][,tau] = softmax(v_depolarization) 
        # store state_posterior (normalised firing rate) from each epoch of
        # gradient descent for each tau
        normalized_firing_rates[[factor]][[1]][Ni,,tau,t,policy] = state_posterior[[factor]][[policy]][,tau]
        # store v (non-normalized log posterior or 'membrane potential') 
        # from each epoch of gradient descent for each tau
        prediction_error[[factor]][[1]][Ni,,tau,t,policy] = as.vector(v_depolarization)
      } # End loop factor
    } # End loop Ni
    # variational free energy for each policy (F)
    Fintermediate = rowSums(Ft[[1]],dims=3)            # sum over hidden state factors (Fintermediate is an intermediate F value)
    Fintermediate1 = colSums(Fintermediate,dims=1)     # sum over tau and squeeze into 16x3 matrix (squeeze here in R not needed)
    # store variational free energy at last iteration of message passing
    VFE[[1]][policy,t] = last(Fintermediate1[,t])
  } # End loop policies
} # End loop Time

                
# Now works up to VFE, i.e., F(policy,t) in Matlab.
                

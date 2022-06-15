
% Replication attempt 06.2022

rng('shuffle')
close all
clear

Gen_model = 2

MDP = explore_exploit_model(Gen_model);

% generative process
A = MDP.A;         % Likelihood matrices
B = MDP.B;         % Transition matrices
C = MDP.C;         % Preferences over outcomes
D = MDP.D;         % Priors over initial states    
T = MDP.T;         % Time points per trial
V = MDP.V;         % Policies
beta = MDP.beta;   % Expected free energy precision
alpha = MDP.alpha; % Action precision
eta = MDP.eta;     % Learning rate
omega = MDP.omega; % Forgetting rate


A = col_norm(A);
B = col_norm(B);
D = col_norm(D);

% generative model (lowercase matrices/vectors are beliefs about capitalized matrices/vectors)

NumPolicies = MDP.NumPolicies; % Number of policies
NumFactors = MDP.NumFactors;   % Number of state factors

% Store initial paramater values of generative model for free energy 
% calculations after learning
%--------------------------------------------------------------------------

% 'complexity' of d vector concentration paramaters

if isfield(MDP,'d')
    for factor = 1:numel(MDP.d)
        % store d vector values before learning
        d_prior{factor} = MDP.d{factor};
        % compute "complexity" - lower concentration paramaters have
        % smaller values creating a lower expected free energy thereby
        % encouraging 'novel' behaviour 
        d_complexity{factor} = spm_wnorm(d_prior{factor});
    end 
end 


if isfield(MDP,'a')
    % complexity of a maxtrix concentration parameters
    for modality = 1:numel(MDP.a)
        a_prior{modality} = MDP.a{modality};
        a_complexity{modality} = spm_wnorm(a_prior{modality}).*(a_prior{modality} > 0);
    end
end 

% normalize A matrix
if isfield(MDP,'a')
    a = col_norm(MDP.a);
else 
    a = col_norm(MDP.A);
end 

% normalize B matrix
if isfield(MDP,'b')
    b = col_norm(MDP.b);
else 
    b = col_norm(MDP.B);
end 

% normalize C and transform into log probability
for ii = 1:numel(C)
    C{ii} = MDP.C{ii} + 1/32;
    for t = 1:T
        C{ii}(:,t) = nat_log(exp(C{ii}(:,t))/sum(exp(C{ii}(:,t))));
    end 
end 



% normalize D vector
if isfield(MDP,'d')
    d = col_norm(MDP.d);
else 
    d = col_norm(MDP.D);
end 

% normalize E vector
if isfield(MDP,'e')
    E = MDP.e;
    E = E./sum(E);
elseif isfield(MDP,'E')
    E = MDP.E;
    E = E./sum(E);
else
    E = col_norm(ones(NumPolicies,1));
    E = E./sum(E);
end

% Initialize variables
%--------------------------------------------------------------------------

% numbers of transitions, policies and states
NumModalities = numel(a);                    % number of outcome factors
NumFactors = numel(d);                       % number of hidden state factors
NumPolicies = size(V,2);                     % number of allowable policies

test= size(b{1},3)
for factor = 1:NumFactors
    NumStates(factor) = size(b{factor},1);   % number of hidden states
    NumControllable_transitions(factor) = size(b{factor},3); % number of hidden controllable hidden states for each factor (number of B matrices)
end

% initialize the approximate posterior over states conditioned on policies
% for each factor as a flat distribution over states at each time point
for policy = 1:NumPolicies
    for factor = 1:NumFactors
        NumStates(factor) = length(D{factor}); % number of states in each hidden state factor
        state_posterior{factor} = ones(NumStates(factor),T,policy)/NumStates(factor); 
    end  
end 



% initialize the approximate posterior over policies as a flat distribution 
% over policies at each time point
policy_posteriors = ones(NumPolicies,T)/NumPolicies; 

% initialize posterior over actions
chosen_action = zeros(ndims(B),T-1);
test=ones(1,T-1)
 
% if there is only one policy
for factors = 1:NumFactors 
    if NumControllable_transitions(factors) == 1
        chosen_action(factors,:) = ones(1,T-1);
    end
end
%MDP.chosen_action = chosen_action;
chosen_action = [1 1 ; 1 1] 
MDP.chosen_action = chosen_action;
% initialize expected free energy precision (beta)
posterior_beta = 1;
gamma(1) = 1/posterior_beta; % expected free energy precision
    
% message passing variables
TimeConst = 4; % time constant for gradient descent
NumIterations  = 16; % number of message passing iterations



% Lets go! Message passing and policy selection 
%--------------------------------------------------------------------------

for t = 1:T % loop over time points  
    
    % sample generative process
    %----------------------------------------------------------------------
    
    for factor = 1:NumFactors % number of hidden state factors
        % Here we sample from the prior distribution over states to obtain the
        % state at each time point. At T = 1 we sample from the D vector, and at
        % time T > 1 we sample from the B matrix. To do this we make a vector 
        % containing the cumulative sum of the columns (which we know sum to one), 
        % generate a random number (0-1),and then use the find function to take 
        % the first number in the cumulative sum vector that is >= the random number. 
        % For example if our D vector is [.5 .5] 50% of the time the element of the 
        % vector corresponding to the state one will be >= to the random number. 

        % sample states 
        if t == 1
            prob_state = D{factor}; % sample initial state T = 1
        elseif t>1
            prob_state = B{factor}(:,true_states(factor,t-1),MDP.chosen_action(factor,t-1));
        end 
            true_states(factor,t) = find(cumsum(prob_state)>= rand,1);
    end 
        % sample observations
    for modality = 1:NumModalities % loop over number of outcome modalities
        outcomes(modality,t) = find(cumsum(a{modality }(:,true_states(1,t),true_states(2,t)))>=rand,1);
    end
    % express observations as a structure containing a 1 x observations 
    % vector for each modality with a 1 in the position corresponding to
    % the observation recieved on that trial
    for modality = 1:NumModalities
        vec = zeros(1,size(a{modality},1));
        index = outcomes(modality,t);
        vec(1,index) = 1;
        O{modality,t} = vec;
        clear vec
    end 
    
    % marginal message passing (minimize F and infer posterior over states)
    %----------------------------------------------------------------------
    
    for policy = 1:NumPolicies
        for Ni = 1:NumIterations % number of iterations of message passing  
            for factor = 1:NumFactors
            lnAo = zeros(size(state_posterior{factor})); % initialise matrix containing the log likelihood of observations
                for tau = 1:T % loop over tau
                    v_depolarization = nat_log(state_posterior{factor}(:,tau,policy)); % convert approximate posteriors into depolarisation variable v
                    if tau<t+1 % Collect an observation from the generative process when tau <= t
                        for modal = 1:NumModalities % loop over observation modalities
                            % this line uses the observation at each tau to index
                            % into the A matrix to grab the likelihood of each hidden state
                            lnA = permute(nat_log(a{modal}(outcomes(modal,tau),:,:,:,:,:)),[2 3 4 5 6 1]);
                            for fj = 1:NumFactors
                                % dot product with state vector from other hidden state factors 
                                % (this is what allows hidden states to interact in the likleihood mapping)    
                                if fj ~= factor        
                                    lnAs = md_dot((lnA),state_posterior{fj}(:,tau),fj);
                                    clear lnA
                                    lnA = lnAs; 
                                    clear lnAs
                                end
                            end
                            lnAo(:,tau) = lnAo(:,tau) + lnA;
                        end
                    end
                end
            end
        end
    end
end
test = lnAo(:,2)
test= size(lnAo(:,:)) % dim is 4 15 and it is just all collumns next to each other!
%lnAotest(:,1) = lnAo(:,1) + lnA;

%%%%%%%%%%%%%%%% Functions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% normalise vector columns
function b = col_norm(B)
numfactors = numel(B);
for f = 1:numfactors
    bb{f} = B{f}; 
    z = sum(bb{f},1); %create normalizing constant from sum of columns
    bb{f} = bb{f}./z; %divide columns by constant
end 
b = bb;
end 


function A  = spm_wnorm(A)
% This uses the bsxfun function to subtract the inverse of each column
% entry from the inverse of the sum of the columns and then divide by 2.
% 
A   = A + exp(-16);
A   = bsxfun(@minus,1./sum(A,1),1./A)/2;
end 

% natural log that replaces zero values with very small values for numerical reasons.
function y = nat_log(x)
y = log(x+exp(-16));
end 

% dot product along dimension f
function B = md_dot(A,s,f)
if f == 1
    B = A'*s;
elseif f == 2
    B = A*s;
end 
end


function MDP = explore_exploit_model(Gen_model)

% Number of time points or 'epochs' within a trial: T
% =========================================================================

% Here, we specify 3 time points (T), in which the agent 1) starts in a 'Start'
% state, 2) first moves to either a 'Hint' state or a 'Choose Left' or 'Choose
% Right' slot machine state, and 3) either moves from the Hint state to one
% of the choice states or moves from one of the choice states back to the
% Start state.

T = 3; % R Script it is Time

% Priors about initial states: D and d
% =========================================================================

%--------------------------------------------------------------------------
% Specify prior probabilities about initial states in the generative 
% process (D)
% Note: By default, these will also be the priors for the generative model
%--------------------------------------------------------------------------

% For the 'context' state factor, we can specify that the 'left better' context 
% (i.e., where the left slot machine is more likely to win) is the true context:

D{1} = [1 0]';  % {'left better','right better'}

% For the 'behavior' state factor, we can specify that the agent always
% begins a trial in the 'start' state (i.e., before choosing to either pick
% a slot machine or first ask for a hint:

D{2} = [1 0 0 0]'; % {'start','hint','choose-left','choose-right'}

%--------------------------------------------------------------------------
% Specify prior beliefs about initial states in the generative model (d)
% Note: This is optional, and will simulate learning priors over states 
% if specified.
%--------------------------------------------------------------------------

% Note that these are technically what are called 'Dirichlet concentration
% paramaters', which need not take on values between 0 and 1. These values
% are added to after each trial, based on posterior beliefs about initial
% states. For example, if the agent believed at the end of trial 1 that it 
% was in the 'left better' context, then d{1} on trial 2 would be 
% d{1} = [1.5 0.5]' (although how large the increase in value is after 
% each trial depends on a learning rate). In general, higher values 
% indicate more confidence in one's beliefs about initial states, and 
% entail that beliefs will change more slowly (e.g., the shape of the 
% distribution encoded by d{1} = [25 25]' will change much more slowly 
% than the shape of the distribution encoded by d{1} = [.5 0.5]' with each 
% new observation).

% For context beliefs, we can specify that the agent starts out believing 
% that both contexts are equally likely, but with somewhat low confidence in 
% these beliefs:

d{1} = [.25 .25]';  % {'left better','right better'}

% For behavior beliefs, we can specify that the agent expects with 
% certainty that it will begin a trial in the 'start' state:

d{2} = [1 0 0 0]'; % {'start','hint','choose-left','choose-right'}


% State-outcome mappings and beliefs: A and a
% =========================================================================

%--------------------------------------------------------------------------
% Specify the probabilities of outcomes given each state in the generative 
% process (A)
% This includes one matrix per outcome modality
% Note: By default, these will also be the beliefs in the generative model
%--------------------------------------------------------------------------

% First we specify the mapping from states to observed hints (outcome
% modality 1). Here, the rows correspond to observations, the columns
% correspond to the first state factor (context), and the third dimension
% corresponds to behavior. Each column is a probability distribution
% that must sum to 1.

% We start by specifying that both contexts generate the 'No Hint'
% observation across all behavior states:

Ns = [length(D{1}) length(D{2})]; % number of states in each state factor (2 and 4)

for i = 1:Ns(2) 

    A{1}(:,:,i) = [1 1; % No Hint
                   0 0; % Machine-Left Hint
                   0 0];% Machine-Right Hint
end

% Then we specify that the 'Get Hint' behavior state generates a hint that
% either the left or right slot machine is better, depending on the context
% state. In this case, the hints are accurate with a probability of pHA. 

pHA = 1; % By default we set this to 1, but try changing its value to 
          % see how it affects model behavior

A{1}(:,:,2) = [0     0;      % No Hint
               pHA 1-pHA;    % Machine-Left Hint
               1-pHA pHA];   % Machine-Right Hint

% Next we specify the mapping between states and wins/losses. The first two
% behavior states ('Start' and 'Get Hint') do not generate either win or
% loss observations in either context:

for i = 1:2

    A{2}(:,:,i) = [1 1;  % Null
                   0 0;  % Loss
                   0 0]; % Win
end
           
% Choosing the left machine (behavior state 3) generates wins with
% probability pWin, which differs depending on the context state (columns):

pWin = .8; % By default we set this to 1, but try changing its value to 
          % see how it affects model behavior
           
A{2}(:,:,3) = [0      0;     % Null        
               1-pWin pWin;  % Loss
               pWin 1-pWin]; % Win

% Choosing the right machine (behavior state 4) generates wins with
% probability pWin, with the reverse mapping to context states from 
% choosing the left machine:
           
A{2}(:,:,4) = [0      0;     % Null
               pWin 1-pWin;  % Loss
               1-pWin pWin]; % Win
           
% Finally, we specify an identity mapping between behavior states and
% observed behaviors, to ensure the agent knows that behaviors were carried
% out as planned. Here, each row corresponds to each behavior state.
           
for i = 1:Ns(2) 

    A{3}(i,:,i) = [1 1];

end

%--------------------------------------------------------------------------
% Specify prior beliefs about state-outcome mappings in the generative model 
% (a)
% Note: This is optional, and will simulate learning state-outcome mappings 
% if specified.
%--------------------------------------------------------------------------
           
% Similar to learning priors over initial states, this simply
% requires specifying a matrix (a) with the same structure as the
% generative process (A), but with Dirichlet concentration parameters that
% can encode beliefs (and confidence in those beliefs) that need not
% match the generative process. Learning then corresponds to
% adding to the values of matrix entries, based on what outcomes were 
% observed when the agent believed it was in a particular state. For
% example, if the agent observed a win while believing it was in the 
% 'left better' context and the 'choose left machine' behavior state,
% the corresponding probability value would increase for that location in
% the state outcome-mapping (i.e., a{2}(3,1,3) might change from .8 to
% 1.8).

% One simple way to set up this matrix is by:
 
% 1. initially identifying it with the generative process 
% 2. multiplying the values by a large number to prevent learning all
%    aspects of the matrix (so the shape of the distribution changes very slowly)
% 3. adjusting the elements you want to differ from the generative process.

% To simulate learning the hint accuracy we
% can specify:

a{1} = A{1}*200;
a{2} = A{2}*200;
a{3} = A{3}*200;

a{1}(:,:,2) =  [0     0;     % No Hint
               .25   .25;    % Machine-Left Hint
               .25   .25];   % Machine-Right Hint
    

% Controlled transitions and transition beliefs : B{:,:,u} and b(:,:,u)
%==========================================================================

%--------------------------------------------------------------------------
% Next, we have to specify the probabilistic transitions between hidden states
% under each action (sometimes called 'control states'). 
% Note: By default, these will also be the transitions beliefs 
% for the generative model
%--------------------------------------------------------------------------

% Columns are states at time t. Rows are states at t+1.

% The agent cannot control the context state, so there is only 1 'action',
% indicating that contexts remain stable within a trial:

B{1}(:,:,1) = [1 0;  % 'Left Better' Context
               0 1]; % 'Right Better' Context
           
% The agent can control the behavior state, and we include 4 possible 
% actions:

% Move to the Start state from any other state
B{2}(:,:,1) = [1 1 1 1;  % Start State
               0 0 0 0;  % Hint
               0 0 0 0;  % Choose Left Machine
               0 0 0 0]; % Choose Right Machine
           
% Move to the Hint state from any other state
B{2}(:,:,2) = [0 0 0 0;  % Start State
               1 1 1 1;  % Hint
               0 0 0 0;  % Choose Left Machine
               0 0 0 0]; % Choose Right Machine

% Move to the Choose Left state from any other state
B{2}(:,:,3) = [0 0 0 0;  % Start State
               0 0 0 0;  % Hint
               1 1 1 1;  % Choose Left Machine
               0 0 0 0]; % Choose Right Machine

% Move to the Choose Right state from any other state
B{2}(:,:,4) = [0 0 0 0;  % Start State
               0 0 0 0;  % Hint
               0 0 0 0;  % Choose Left Machine
               1 1 1 1]; % Choose Right Machine        
           
%--------------------------------------------------------------------------
% Specify prior beliefs about state transitions in the generative model
% (b). This is a set of matrices with the same structure as B.
% Note: This is optional, and will simulate learning state transitions if 
% specified.
%--------------------------------------------------------------------------
          
% For this example, we will not simulate learning transition beliefs. 
% But, similar to learning d and a, this just involves accumulating
% Dirichlet concentration parameters. Here, transition beliefs are updated
% after each trial when the agent believes it was in a given state at time
% t and and another state at t+1.

% Preferred outcomes: C and c
%==========================================================================

%--------------------------------------------------------------------------
% Next, we have to specify the 'prior preferences', encoded here as log
% probabilities. 
%--------------------------------------------------------------------------

% One matrix per outcome modality. Each row is an observation, and each
% columns is a time point. Negative values indicate lower preference,
% positive values indicate a high preference. Stronger preferences promote
% risky choices and reduced information-seeking.

% We can start by setting a 0 preference for all outcomes:

No = [size(A{1},1) size(A{2},1) size(A{3},1)]; % number of outcomes in 
                                               % each outcome modality

C{1}      = zeros(No(1),T); % Hints
C{2}      = zeros(No(2),T); % Wins/Losses
C{3}      = zeros(No(3),T); % Observed Behaviors

% Then we can specify a 'loss aversion' magnitude (la) at time points 2 
% and 3, and a 'reward seeking' (or 'risk-seeking') magnitude (rs). Here,
% rs is divided by 2 at the third time point to encode a smaller win ($2
% instead of $4) if taking the hint before choosing a slot machine.

la = 1; % By default we set this to 1, but try changing its value to 
        % see how it affects model behavior

rs = 4; % By default we set this to 4, but try changing its value to 
        % see how it affects model behavior

C{2}(:,:) =    [0  0   0   ;  % Null
                0 -la -la  ;  % Loss
                0  rs  rs/2]; % win
            
%--------------------------------------------------------------------------
% One can also optionally choose to simulate preference learning by
% specifying a Dirichlet distribution over preferences (c). 
%--------------------------------------------------------------------------

% This will not be simulated here. However, this works by increasing the
% preference magnitude for an outcome each time that outcome is observed.
% The assumption here is that preferences naturally increase for entering
% situations that are more familiar.

% Allowable policies: U or V. 
%==========================================================================

%--------------------------------------------------------------------------
% Each policy is a sequence of actions over time that the agent can 
% consider. 
%--------------------------------------------------------------------------

% For our simulations, we will specify V, where rows correspond to time 
% points and should be length T-1 (here, 2 transitions, from time point 1
% to time point 2, and time point 2 to time point 3):

NumPolicies = 5; % Number of policies
NumFactors = 2; % Number of state factors

V         = ones(T-1,NumPolicies,NumFactors);

V(:,:,1) = [1 1 1 1 1;
            1 1 1 1 1]; % Context state is not controllable

V(:,:,2) = [1 2 2 3 4;
            1 3 4 1 1];
        
% For V(:,:,2), columns left to right indicate policies allowing: 
% 1. staying in the start state 
% 2. taking the hint then choosing the left machine
% 3. taking the hint then choosing the right machine
% 4. choosing the left machine right away (then returning to start state)
% 5. choosing the right machine right away (then returning to start state)


% Habits: E and e. 
%==========================================================================

%--------------------------------------------------------------------------
% Optional: a columns vector with one entry per policy, indicating the 
% prior probability of choosing that policy (i.e., independent of other 
% beliefs). 
%--------------------------------------------------------------------------

% We will not equip our agent with habits with any starting habits 
% (flat distribution over policies):

E = [1 1 1 1 1]';

% To incorporate habit learning, where policies become more likely after 
% each time they are chosen, we can also specify concentration parameters
% by specifying e:

 e = [1 1 1 1 1]';

% Additional optional parameters. 
%==========================================================================

% Eta: learning rate (0-1) controlling the magnitude of concentration parameter
% updates after each trial (if learning is enabled).

     eta = 1; % Default (maximum) learning rate
     
% Omega: forgetting rate (0-1) controlling the magnitude of reduction in concentration
% parameter values after each trial (if learning is enabled).

     omega = 1; % Default value indicating there is no forgetting (values < 1 indicate forgetting)

% Beta: Expected precision of expected free energy (G) over policies (a 
% positive value, with higher values indicating lower expected precision).
% Lower values increase the influence of habits (E) and otherwise make
% policy selection less deteriministic.

     beta = 1; % By default this is set to 1, but try increasing its value 
               % to lower precision and see how it affects model behavior

% Alpha: An 'inverse temperature' or 'action precision' parameter that 
% controls how much randomness there is when selecting actions (e.g., how 
% often the agent might choose not to take the hint, even if the model 
% assigned the highest probability to that action. This is a positive 
% number, where higher values indicate less randomness. Here we set this to 
% a fairly high value:

    alpha = 32; % fairly low randomness in action selectionend
    
    %% Define POMDP Structure
%==========================================================================

mdp.T = T;                    % Number of time steps
mdp.V = V;                    % allowable (deep) policies

mdp.A = A;                    % state-outcome mapping
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred states
mdp.D = D;                    % priors over initial states
mdp.d = d;                    % enable learning priors over initial states

if Gen_model == 1
    mdp.E = E;                % prior over policies
elseif Gen_model == 2
    mdp.a = a;                % enable learning state-outcome mappings
    mdp.e = e;                % enable learning of prior over policies
end 

mdp.eta = eta;                % learning rate
mdp.omega = omega;            % forgetting rate
mdp.alpha = alpha;            % action precision
mdp.beta = beta;              % expected free energy precision

%respecify for use in inversion script (specific to this tutorial example)
mdp.NumPolicies = NumPolicies; % Number of policies
mdp.NumFactors = NumFactors; % Number of state factors
    
   
% We can add labels to states, outcomes, and actions for subsequent plotting:

label.factor{1}   = 'contexts';   label.name{1}    = {'left-better','right-better'};
label.factor{2}   = 'choice states';     label.name{2}    = {'start','hint','choose left','choose right'};
label.modality{1} = 'hint';    label.outcome{1} = {'null','left hint','right hint'};
label.modality{2} = 'win/lose';  label.outcome{2} = {'null','lose','win'};
label.modality{3} = 'observed action';  label.outcome{3} = {'start','hint','choose left','choose right'};
label.action{2} = {'start','hint','left','right'};
mdp.label = label;

MDP = mdp;
end

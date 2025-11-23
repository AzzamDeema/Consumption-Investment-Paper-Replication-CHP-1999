%% Replication of "Machine Replacement and the Business Cycle: Lumps and Bumps" (1999)
% Authors: Russell Cooper, John Haltiwanger and Laura Power 

%% 1) Calibration of Parameters
beta=0.9;               % discount rate
delta_k=0.1;            % depreciation of capital
lambda=0.75;            % replacement opportunity cost
                        % A is aggregate shock
A_L=0.75;               % A low state
A_H=1.25;               % A high state    
P=[0.9,0.1;0.1,0.9];    % markov matrix procss % rows sum to 1
fc=0.2;                 % fixed cost
% Grid Pre-allocation
k_n=12;                 % capital (1 + the cut-off level; 11 for low state, 6 for high state, use larger of the two states so 11+1=12)
a_n=2;                  % aggregate shock
productivity_n=20;      % idiosyncratic poductivity 
productivity=linspace(0.4,1.6,productivity_n); %  Set an equally spaced grid of 20 gridpoints for ε productivity. pg. 926 row vector 1x20
% Inicialization of grid
comp_delta_k=1-delta_k;
for ik=0:k_n-1
    k(ik+1)=comp_delta_k^ik;       % Array indices must be positive 
                                   % integers or logical values so +1 in k(ik+1)
                                   % (1-delta_k)^ ik
                                   
end

% Revenue 3-D
revenue=productivity'*A_L*k;  % pg 924, e transpose for correct dimensions
                   % productivity'--> 20x1 column vector, A_l--> 1x1 scalar, k 1x8 row
                   % vector--> 20x1*1x1*1*8, revenue 3-D 20*8*2 double
revenue(:,:,a_n)=productivity'*A_H*k;
% More pre-allocations, 3-D
V_replace=zeros(productivity_n,k_n,a_n); % same dimenions as revenue
V_noreplace=zeros(productivity_n,k_n,a_n);
V_inicial=max(V_replace,V_noreplace);
z_function=zeros(productivity_n,k_n,a_n); % The choice variable z is an indicator variable 
                          % taking the value 1 if capital replacement occurs
                          % and 0 if the firm does not replace its capital 
                          % in the period under consideration.

%% 2) Write the value function iteration loop (inlcluding policy function).
toler=1e-9;         % torlerance. 
max_iter=10000; % maximum number of iterations.
current_iter=0; % Inicilaize iterations.
norm_V=inf;     % Inicialize Norm ie difference.

tic
while (norm_V>toler) && (current_iter<max_iter)
    z_function=zeros(productivity_n,k_n,a_n);    % Pre-allocation to jumpstart VFI, probability
    for ia=1:a_n
        V_replace(:,:,ia)=revenue(:,:,ia)*lambda-fc+...
            beta*P(ia,1)*repmat(V_inicial(:,1,1),1,k_n)+...
            beta*P(ia,a_n)*repmat(V_inicial(:,1,a_n),1,k_n);
        for jk=1:k_n-1
            V_noreplace(:,jk,ia)=revenue(:,jk+1,ia)+...
                beta*P(ia,1)*V_inicial(:,jk+1,1)+...
                beta*P(ia,a_n)*V_inicial(:,jk+1,a_n);
        end
    end
    V_1=max(V_replace,V_noreplace);
    z_function(V_1==V_replace)=1; % Probability = 1
    norm_V=max(max(max(abs(V_1-V_inicial)))); % triple max due to 3 dimensions
    current_iter=current_iter+1;
    V_inicial=V_1;
end
toc

%% 3) Examine the policy function. There should be a cut-off level of capital 
% below which capital is always replaced. Reduce the size of the k dimension
% to 1 plus the cut-off level and repeat questions 1. and 2.

%% 4) Plot the policy function using a spy plot.

figure(799)
% Low State
subplot(2,1,1); spy(z_function(:,:,1)')
title('Policy Function: Low State'); xlabel('Lapse of Time Since Replacement'); ylabel('Idiosyncratic Productivity');
% High State
subplot(2,1,2); spy(z_function(:,:,2)')
title('Policy Function: High State'); xlabel('Lapse of Time Since Replacement'); ylabel('Idiosyncratic Productivity');

%% Question: Comment on the important features of the policy function.

% Answer:
% In high state of idiosyncratic productivity, firms will engage in more
% capital replacement than in the low state. Specifically, in the low state,
% the cut-off level of capital is 11, under which firms always replace 
% capital for every idiosyncratic productivity level. For high 
% aggregate states, this is 6. Replacement is more probable in high states.
% With higher levels of idiosyncratic productivity, firms face lower costs
% of replacement, thus replace more often. Additionally, as lapse of time
% rises, replacement of capital increases weakly with it.

% --> Capital k_n=12
% 1 + the cut-off level; 11 for low state, 6 for high state, use larger of
% the two states so 11+1=12.


%% 5) Replication of Figure 1: THEORETICAL HAZARD FOR MACHINE REPLACEMENT
% Plot the hazard function of capital replacement for the two values of A. 
% Essentially, replicate Figure 1 in Cooper, Haltiwanger, & Power (1999)("the paper").

% (pg's 926-927) Figure 1 in CHP 1999 depicts "...the hazard functions for 
% this economy conditional on the aggregate state. The horizontal axis 
% measures the time since last replacement. Thus, from Proposition 2, 
% the hazard is increasing in the time since last replace ment, given A. 
% In the discussion that follows, we term this an "increasing hazard" to 
% focus on the relationship between the probability of replacement and 
% the time since last replacement given A.17 Further, fixed costs and 
% persistent shocks combine to imply a procyclical hazard: the probability 
% of replacement increases during periods of high profitability."
% PROPOSITION 2: H(k, A) is decreasing in k. pg 925

hazard=zeros(k_n, a_n);
for ik=1:k_n
    for ja=1:a_n
        hazard(ik,ja)=sum(z_function(:,ik,ja))/productivity_n; % sum --> probability of replacement over time
    end
end

figure (800);
plot((1:k_n),hazard(:,a_n),"Color",'k');                % plot high state
hold on
plot((1:k_n),hazard(:,1),'LineWidth',2,'Color','k');    % plot low state
title('THEORETICAL HAZARD FOR MACHINE REPLACEMENT');
xlabel('Lapse of Time Since Replacement');
ylabel('Probability of Replacement');
legend('Low state','High State');

%% Question: Comment on the important features of the hazard function.

% Answer: 
% Hazard for both levels of A increases with the lapse of time since
% replacement. Moreover, Hazard for low state is above the high state until
% period five, when they both become equal. This implies that replacement is
% more likely to occur with low states. Once capital becomes "old
% enough", firms replace capital quickly. This makes sense once accounting
% for increasing revenues as well. Once both hazard functions meet, this
% implies that at that age of capital, firms will certainly replace capital
% no matter the level of productivity state.

%% 6) Simulate a time series for the evolution of one firm in the model using 
% the policy function and the stochastic processes for the exogenous 
% variables. Plot sample paths for something like 40 periods with the 
% output of the firm and the capital stock of the firm.
T = 40;

% Setup
A_vals = [A_L, A_H];  % numeric values for aggregate states
A_idx  = zeros(T,1); 
A_idx(1)=1;           % start in low state
eps_idx = zeros(T,1);
k_idx   = zeros(T,1); 
k_idx(1)=1;           % start with new capital

z_sim  = zeros(T,1);
k_path = zeros(T,1);
y_path = zeros(T,1);

for t = 1:T
    % draw idiosyncratic productivity from the 20 point uniform grid
    eps_idx(t) = randi(productivity_n); % stochastic random process for simulation of productivity over t time periods.
    je = eps_idx(t);
    ia = A_idx(t);
    ik = k_idx(t);

    % Policy decision
    z_sim(t) = z_function(je, ik, ia);

    % realized output this period (lambda if replacing, otherwise not apply)
    k_path(t)=k(ik);
    if z_sim(t) == 1
        y_path(t)=A_vals(ia)*productivity(je)*lambda*k_path(t);
    else
        y_path(t)=A_vals(ia)*productivity(je)*k_path(t);
    end

    % Transitions
    if t < T
        % Capital transition
        if z_sim(t)==1
            k_idx(t+1)=1;
        else
            k_idx(t+1)=min(ik + 1, k_n);
        end
        % Aggregate state transition (markov matrix P)
        u=rand;
        if u<=P(ia,1)
            A_idx(t+1)=1;
        else
            A_idx(t+1)=2;
        end
    end
end

% Plot sample paths --> output and capital stock
figure(801)
subplot(2,1,1)
plot(1:T, y_path, 'k', 'LineWidth', 1.5)
title('Output over Time (One Firm)')
xlabel('Time')
ylabel('Output')

subplot(2,1,2)
plot(1:T, k_path, 'k', 'LineWidth', 1.5)
title('Capital Stock over Time (One Firm)')
xlabel('Time')
ylabel('Capital')

%% Question: Comment on the behavior of firms in the model.

% Answer:
% Overall, there is a zigzag for capital, most likely due to depreciation 
% and aging of capital. There are long periods of low output, or inaction, 
% for firms, before a large spike, probably due to large investment. 
% Here, investment is lumpy. Output variation is partly due to shocks of 
% capital stock variation. Moreover, capital replacement varies wherein a 
% high state pushes replacement time farther. A low state draw induces 
% firms to replace them earlier.

%% 7) Replication of Figure 3:  CONVERGENCE WITHOUT AGGREGATE SHOCKS-BASELINE PARAMETERS
% Assume that A is fixed. Replicate Figure 3 in the paper.

% (pg's 927-928) Figure 3 in CHP 1999 depicts "...the behavior of aggregate
% investment starting from a uniform distribution of producers across the 
% capital state space in an economy without aggregate uncertainty. 
% Initially there is a burst of investment since many producers start with 
% relatively old capital. Over time the distribution evolves until a 
% stationary distribution is obtained and the level of aggregate investnent
% is constant. Each producer undertakes machine replacement at stochastic 
% intervals of time: the deterministic machine replacement system is made 
% stochastic by the producer-specific shocks. Yet, due to the assumed large
%  number of plants, the economy has a stationary distribution with 
% constant aggregate investnent." 
% "...the initial burst of investment shown in period 1 of Figure 3 creates
% a relatively large fraction of "young plants" in periods 2 and 3. Since, 
% from the hazards shown in Figure 1, young plants have low spike 
% probabilities, the initial burst of investment is followed by a couple of
% periods of below-average investment."

A=1;                            % Fixed A, no uncertainties
periods=50;                     % Period specification in CHP 1999.
investment=zeros(periods,1);    % Investment rate
initial=zeros(k_n,periods+1);   % Note: initial(j_k, i_t)-->fraction of 
                                % firms that are at age j_k in period i_t 
                                % i.e. their capital is at grid point k(j_k).
initial(:,1)=ones(k_n,1);              
for i_t=1:periods

    for j_k=1:k_n
    investment(i_t)=investment(i_t)+hazard(j_k,A)*initial(j_k,i_t);
    end

    initial(1,i_t+1)=investment(i_t);

    for j_a=a_n:k_n
        initial(j_a,i_t+1)=(1-hazard(j_a-1,A))*initial(j_a-1,i_t);
    end

end

rate_investment=investment/k_n;
figure(802);
inv_periods=(1:periods);
plot(inv_periods,rate_investment,'Color','k');
title('CONVERGENCE WITHOUT AGGREGATE SHOCKS-BASELINE PARAMETERS');
xlabel('Period');
ylabel('Investment Rate');

%% 8) Now let A follow the Markov process specified earlier. (Some uncertainties)
% Replicate Figure 4 in the paper.

markov_t=160; % period specification in CHP 1999.

markov_investment=zeros(markov_t-1,1);
markov_initial=zeros(k_n, markov_t);
markov_initial(:,1)=ones(k_n,1);

A_idx_markov=zeros(markov_t,1);
A_idx_markov(1)=1; % start in low state
for t=1:markov_t-1
    i_a=A_idx_markov(t);
    u=rand;
    if u<=P(i_a,1)
        A_idx_markov(t+1)=1;
    else
        A_idx_markov(t+1)=2;
    end
end

for i_t=1:markov_t-1

    % aggregate investment this period
    for j_k=1:k_n
        markov_investment(i_t)=markov_investment(i_t)+hazard(j_k, A_idx_markov(i_t))*markov_initial(j_k,i_t);
    end

    % replacers move to age 1 next period
    markov_initial(1,i_t+1)=markov_investment(i_t);

    % keep survivors at last age
    for j_k=2:k_n
        markov_initial(j_k,i_t+1)=markov_initial(j_k,i_t+1)+(1-hazard(j_k-1, A_idx_markov(i_t)))*markov_initial(j_k-1,i_t);
    end
    markov_initial(k_n,i_t+1)=markov_initial(k_n,i_t+1)+(1-hazard(k_n, A_idx_markov(i_t)))*markov_initial(k_n,i_t);

end

markov_rate_investment=markov_investment/k_n;

figure(803);
t_minus_1=1:markov_t-1;
just_t=1:markov_t;
[AX,line1,line2]=plotyy(t_minus_1,markov_rate_investment,just_t,A_idx_markov);
ylabel(AX(1),'Investment Rate','Color','k')
ylabel(AX(2),'Aggregate State','Color','k')
set(line1,'Color','k')
set(line2,'Color','k')
line2.Marker='x';
line2.LineStyle='--';
xlabel('period');
title('AGGREGATE INVESTMENT FLUCTUATIONS-BASELINE SIMULATION');

%% 9) Question: Comment on the implications for firm investment behavior that you 
% see as important from your replications of Figures 3 and 4.

% Answer:
% In figure 3, with A fixed, firms high no uncertainty therefore the
% investment rate begins with a large spike where most firms immediately
% replace capital. The replacement rate then drops due to most machines
% having little hazard to be replaced. This pattern then moves in a
% declining waves as capital becomes older until settling on a stationary
% level. 
% Meanwhile, in figure 4, firm face A following a Markov process, resetting 
% the investment rate depending on the endowed state. Within each
% aggregate state, some firms will choose to invest, before displaying the
% same oscillating wages to a near stationary level, before resetting with a
% new aggregate state switch.
% Therefore, without aggregate uncertainty shocks, capital investment lumpiness
% creates an initial spike in investment before oscillating back to a
% stationary level. That said, when incorporating aggregate shocks in the
% form of a Markov process for A, hazard changes inducing new phases for
% investment rate moves with reoccurring and resetting spikes.
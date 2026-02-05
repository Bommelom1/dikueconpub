% This solves the bellman equation for a given set of prices and parameters

addpath('matlabinclude');
addpath('autotrade');

clear
clc

% read default parameters
mp0=trmodel.default_params;

% Parameter adjustments (relative to defaults)
mp0.ntypes = 1;
mp0.ncartypes=2;
mp0.abar=10;

% populate mp with remaining fields and update model parameter dependencies
mp = trmodel.setparams(mp0);

% model indexing
s = trmodel.index(mp);

% initialize prices where to solve bellman
p = mp.pnew.*0.85.^(1:mp.abar-1)'; 

% pre-compute age transition matrices
F = trmodel.age_transition(mp, s);

ap=dpsolver.setup;
ap.printfxp=1;  % show output

t=1;  fprintf('Solving model for sonsumer type, tau=%d\n ', t);

% pre-compute utility
util = trmodel.utility(mp, s, t, p);   % update the current period, post-trade payoffs stored in the global util_tau

% step 1: solve the DP problem for the price vector p for each consumer type
bellman=@(ev) trmodel.bellman(mp, s, util, F, ev);

ev0=zeros(s.ns, 1); % initial guess on expected value function
ev_sa= dpsolver.sa(bellman, ev0, ap);
ev_poly= dpsolver.poly(bellman, ev0(:,t), ap, mp.bet);




% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline

% This solves the bellman equation for a given set of prices and parameters

addpath('matlabinclude');
addpath('autotrade');

clear
clc

% set parameters
mp = setparams.default(); % default parameters used for illustration

% set number of types
mp.ntypes = 2000;
mp.ncartypes=3;

% update model parameter dependencies
mp = trmodel.update_mp(mp);

% model indexing
s = trmodel.index(mp);

% initialize prices where to solve bellman
p = equilibrium.price_init(mp,s);

ev0=zeros(s.ns, mp.ntypes); % initial guess on expected value function

% initialize arrays
util=nan(s.ns,s.nd, mp.ntypes); 	% utility function
ev_sa=nan(s.ns,mp.ntypes);	% expected value function solved by successive approximations (dpsolver.sa)
ev_poly=nan(s.ns, mp.ntypes); 	% expected value function solved by poly algorithm (dpsolver.poly)

% pre-compute age transition matrices
F = trmodel.age_transition(mp, s);

ap=dpsolver.setup;
ap.printfxp=2;  % show output

for t=1:mp.ntypes;
	fprintf('tau=%d\n ', t)

	% pre-compute utility
	util(:,:,t) = trmodel.utility(mp, s, t, p);   % update the current period, post-trade payoffs stored in the global util_tau

    % step 1: solve the DP problem for the price vector p for each consumer type
    bellman=@(ev) trmodel.bellman(mp, s, util(:,:,t), F, ev);

    ev_sa(:,t)= dpsolver.sa(bellman, ev0(:,t), ap);

    ev_poly(:,t)= dpsolver.poly(bellman, ev0(:,t), ap, mp.bet);

end


% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline
